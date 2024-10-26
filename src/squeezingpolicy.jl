struct StateChoice{Z,A<:AbstractVector{ItemState}}# <: AbstractVector{ItemState}
	lb::Z
	ub::Z
	x::A
end

Base.size(sc::StateChoice) = size(sc.x)
Base.@propagate_inbounds Base.getindex(sc::StateChoice, i::Int) = sc.x[i]

struct Policy{Z,A} <: AbstractVector{StateChoice{Z,A}}
	cutoffs::Vector{Z}
	xs::Vector{A}
    ub::Z
end

Base.size(p::Policy) = size(p.xs)
Base.@propagate_inbounds Base.getindex(p::Policy, i::Int) =
	StateChoice(p.cutoffs[i], i+1>length(p.xs) ? p.ub : p.cutoffs[i+1], p.xs[i])

struct DiffObj{Obj}
    obj1::Obj
    obj2::Obj
end

function (d::DiffObj)(z, fcall)
    # x should have been set
    f1, obj1 = value(d.obj1, z)
    f2, obj2 = value(d.obj2, z)
    fcall[] += 2
    return f2 - f1
end

struct Equal_Obj{M,KW<:NamedTuple}
    m::M
    kwargs::KW
    fcall::RefValue{Int}
end

Equal_Obj(m, kwargs::NamedTuple=NamedTuple()) = Equal_Obj(m, kwargs, Ref(0))

struct Default_Zero_Margin{F, O}
    equal_obj::F
    obj2::O
end

_copyx(obj::Objective{<:Any, A}, x::A) where A<:SVector =
    Objective(obj.f, x, obj.fcall)

# Fallback method assumes A is a mutable array
_copyx(obj::Objective{<:Any, A}, x::A) where A =
    Objective(obj.f, copyto!(obj.x, x), obj.fcall)

function (zm::Default_Zero_Margin)(obj::Objective, i, lb, ub)
    # The order of true and false matters in case there is no interior solution
    obj = _setx(obj, true, i)
    obj2 = _setx(_copyx(zm.obj2, obj.x), false, i)
    return zm.equal_obj(obj, obj2, lb, ub)
end

struct SqueezingPolicyTrace{Z}
	i::Int
    z::Z
    lb::Z
    ub::Z
	s::ItemState
end

struct SqueezingPolicy{Z,A,AO,F1,F2,O,TR} <: CDCPSolver
	scdca::Bool
    allu::A
    pool::Vector{StateChoice{Z,A}}
    seen::Set{StateChoice{Z,A}}
    branching::Vector{Int}
    lookup_zero_margin::Dict{Tuple{Int,AO},Z}
    zero_margin::F1
    equal_obj::F2
    obj2::O
    zero_margin_call::RefValue{Int}
    equal_obj_call::RefValue{Int}
    trace::TR
    nobranching::Bool
end

# A special value for indicating no new cutoff
const _nextunknown = (Inf, Inf, 0, undetermined)

function _init(::Type{<:SqueezingPolicy}, obj, scdca::Bool, equal_obj, zbounds::Tuple{Z,Z};
		zero_margin=nothing, x0=nothing, trace::Bool=false,
        nobranching::Bool=false, kwargs...) where Z
	S = length(obj.x)
    allu = obj.x isa SVector ? _fillstate(SVector{S,ItemState}, undetermined) :
        fill(undetermined, S)
	if x0 === nothing
		if obj.x isa SVector
        	x = Policy([zbounds[1]], [allu], zbounds[2])
    	else
            x = Policy([zbounds[1]], [allu], zbounds[2])
		end
	else
		x = x0
    end
    tr = trace ? [SqueezingPolicyTrace{Z}[]] : nothing
    obj2 = deepcopy(obj)
    if zero_margin === nothing
        zero_margin = Default_Zero_Margin(equal_obj, obj2)
    end
    pool = [x[i] for i in eachindex(x.xs)]
	return SqueezingPolicy(scdca, allu, pool, Set(pool), collect(1:length(x.xs)),
        Dict{Tuple{Int,typeof(obj.x)},Z}(),
        zero_margin, equal_obj, obj2, Ref(0), Ref(0), tr, nobranching), x
end

_squeeze(x::StateChoice, s::ItemState, i::Int) =
    StateChoice(x.lb, x.ub, _squeeze(x.x, s, i))

_setitemstate(x::StateChoice, s::ItemState, i::Int) =
    StateChoice(x.lb, x.ub, _setitemstate(x.x, s, i))

function squeeze!(p::CDCP{<:SqueezingPolicy}, x::StateChoice, i::Int)
	obj, scdca, tr = p.obj, p.solver.scdca, p.solver.trace
    lookup = p.solver.lookup_zero_margin
    obj = _setchoice(obj, scdca ? setsup(x.x) : setsub(x.x))
    key = (i, obj.x)
    z = get(lookup, key, nothing)
    if z === nothing
        z, obj = p.solver.zero_margin(obj, i, x.lb, x.ub)
        lookup[key] = z
        p.solver.zero_margin_call[] += 1
        p.obj = obj
    end
    if z <= x.lb # Include the whole interval
        x = _squeeze(x, included, i)
        tr === nothing || push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, included))
        return x, _nextunknown, included
    elseif z >= x.ub # Cannot include any part of the interval
        tr === nothing ||
            push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, undetermined))
        return squeeze_exclude!(p, x, i)
    else
        tr === nothing ||
            push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, undetermined))
        ub0 = x.ub # ub may be replaced by squeeze_exclude!
        x, next, state = squeeze_exclude!(p, x, i)
        next = next[1] < z ? (next[1], z, next[3:4]...) : (z, ub0, i, included)
        # next[1] can be smaller than x.ub if cannot exclude any part by squeeze_exclude!
        x = StateChoice(x.lb, min(x.ub, next[1]), x.x)
        return x, next, state
    end
end

function squeeze_exclude!(p::CDCP{<:SqueezingPolicy}, x::StateChoice, i::Int)
    obj, scdca, tr = p.obj, p.solver.scdca, p.solver.trace
    lookup = p.solver.lookup_zero_margin
    obj = _setchoice(obj, scdca ? setsub(x.x) : setsup(x.x))
    key = (i, obj.x)
    z = get(lookup, key, nothing)
    if z === nothing
        z, obj = p.solver.zero_margin(obj, i, x.lb, x.ub)
        lookup[key] = z
        p.solver.zero_margin_call[] += 1
        p.obj = obj
    end
    if z >= x.ub # Exclude the whole interval
        x = _squeeze(x, excluded, i)
        tr === nothing || push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, excluded))
        return x, _nextunknown, excluded
    elseif z <= x.lb # Cannot exclude any part of the interval
        x = _setitemstate(x, aux, i) # Don't use _squeeze
        tr === nothing || push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, aux))
        return x, _nextunknown, aux
    else
        tr === nothing || push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, excluded))
        ub0 = x.ub
        x = StateChoice(x.lb, z, _squeeze(x.x, excluded, i))
        return x, (z, ub0, i, aux), excluded
    end
end

function branch(x::StateChoice, i::Int)
    xin = StateChoice(x.lb, x.ub, _squeeze(x.x, included, i))
	xex = StateChoice(x.lb, x.ub, _squeeze(x.x, excluded, i))
	return xin, xex
end

function _squeezenext(::SVector{S}, s::ItemState, k::Int) where S
    if @generated
        ex = :(())
        for i in 1:S
            push!(ex.args, :(ifelse($i==k, s, undetermined)))
        end
        return :(SVector{S,ItemState}($ex))
    else
        SVector{S,ItemState}(ntuple(i->ifelse(i==k, s, undetermined), S))
    end
end

function squeeze!(p::CDCP{<:SqueezingPolicy})
	S = length(p.x.xs[1])
    pool, seen, branching = p.solver.pool, p.solver.seen, p.solver.branching
    while !isempty(branching)
        k = pop!(branching)
        x = pool[k]
        ub0 = x.ub # Save the original one before it's replaced
        next = _nextunknown
        lastaux = nothing
        i = findfirst(==(undetermined), x.x)
        while i !== nothing
            p.obj.fcall < p.maxfcall || return maxfcall_reached
	    	x, next1, itemstate = squeeze!(p, x, i)
            next = ifelse(next1[1] < next[1], next1, next)
	    	lastaux = ifelse(itemstate == aux, i, nothing)
	    	if i < S
	    		i = findnext(==(undetermined), x.x, i+1)
	    		i === nothing && (i = findfirst(==(undetermined), x.x))
	    	else
	    		i = findfirst(==(undetermined), x.x)
	    	end
	    end

        # Avoid accumulating duplicates
        N = length(seen)
        push!(seen, x)
        length(seen) > N || continue

        if lastaux === nothing || p.solver.nobranching
            pool[k] = x
        else
            xin, xout = branch(x, lastaux)
            pool[k] = xin
            push!(pool, xout)
            push!(branching, k, length(pool))
        end
        if next != _nextunknown
            lb, ub, i, s = next
            sc = StateChoice(lb, ub, p.solver.allu)#_squeezenext(x.x, s, i))
            push!(pool, sc)
            push!(branching, length(pool))
            # Must make sure no interval is left out
            if ub < ub0
                sc = StateChoice(ub, ub0, p.solver.allu)
                push!(pool, sc)
                push!(branching, length(pool))
            end
        end
    end
    return inprogress
end

_lb(x::StateChoice) = x.lb

# By default, assume any type is valid unless it contains Inf/-Inf
checktype(::Objective, z) = !(Inf in z || -Inf in z)

function combine_branch!(p::CDCP{<:SqueezingPolicy})
    obj, pool = p.obj, p.solver.pool
    sort!(pool, by=_lb)
    # Overwrite p.x
    cutoffs = empty!(p.x.cutoffs)
    xs = empty!(p.x.xs)
    z = pool[1].lb
    # Handle the case where lb is something like -Inf
    if !checktype(obj, z)
        push!(cutoffs, z)
        push!(xs, _fillstate(eltype(xs), excluded))
        for x in pool
            z1 = x.lb
            if z1 > z
                z = z1
                break
            end
        end
    end
    znext = z
    xmax = nothing
    while true
        fx = -Inf
        # Search through overlapping intervals and determine the next z
        for x in pool
            if x.lb > z
                znext = x.lb
                break
            elseif z < x.ub
                obj.fcall < p.maxfcall || return maxfcall_reached
                obj = _setchoice(obj, setsub(x.x))
                fx1, obj = value(obj, z)
                if fx1 > fx
                    fx = fx1
                    xmax = x.x
                end
            end
        end
        if isempty(xs)
            push!(cutoffs, z)
            push!(xs, xmax)
        elseif xmax != xs[end]
            obj = _setchoice(obj, setsub(xs[end]))
            obj2 = _setchoice(p.solver.obj2, setsub(xmax))
            lastz = cutoffs[end]
            znew, obj = p.solver.equal_obj(obj, obj2, lastz, z)
            p.solver.equal_obj_call[] += 1
            push!(cutoffs, ifelse(lastz < znew < z, znew, z))
            push!(xs, xmax)
        end
        if znext > z
            z = znext
        else
            break
        end
    end
    obj = p.obj
    return success
end

function solve!(p::CDCP{<:SqueezingPolicy}; restart::Bool=false)
	# restart && (p = _reinit!(p))
	p.state = squeeze!(p)
    if p.state == maxfcall_reached
        @warn "maxfcall is reached before convergence"
        return p
    end
    p.solver.nobranching || (p.state = combine_branch!(p))
	return p
end

# Compatibility with old interface
# Some keyword arguments may simply be dummies
# ! Rewrite for the new interface is strongly recommended
function policy!((cutoffs, policies), containers, C::Integer; restart::Bool=true, kwargs...)
    containers isa CDCP{<:SqueezingPolicy} ||
        throw(ArgumentError("containers of type $(typeof(containers)) is no longer supported; please switch to the new interface"))
    p = containers
    solve!(p; restart=restart)
    resize!(cutoffs, length(p.x.cutoffs)+1)
    copyto!(cutoffs, p.x.cutoffs)
    cutoffs[end] = p.x.ub
    resize!(policies, length(p.x.xs))
    for i in eachindex(policies)
        policies[i] = BitVector(setsub(p.x.xs[i]))
    end
    return cutoffs, policies
end

function policy(C::Integer;
        scdca::Bool, obj, equalise_obj,
        D_j_obj=nothing, zero_D_j_obj=nothing, show_time::Bool=false,
        emptyset=falses(C), restart::Bool=true, kwargs...)
    p = solve(SqueezingPolicy, obj, C, scdca, equalise_obj, (-Inf, Inf);
        zero_margin=zero_D_j_obj, kwargs...)
    cutoffs = copy(p.x.cutoffs)
    push!(cutoffs, p.x.ub)
    policies = [BitVector(setsub(x)) for x in p.x.xs]
    return cutoffs, policies
end
