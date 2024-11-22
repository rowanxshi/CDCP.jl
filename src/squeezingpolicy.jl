struct IntervalChoice{Z,A<:AbstractVector{ItemState}}# <: AbstractVector{ItemState}
	lb::Z
	ub::Z
	x::A
end

Base.size(ic::IntervalChoice) = size(ic.x)
Base.@propagate_inbounds Base.getindex(ic::IntervalChoice, i::Int) = ic.x[i]

struct Policy{Z,A} <: AbstractVector{IntervalChoice{Z,A}}
	cutoffs::Vector{Z}
	xs::Vector{A}
    ub::Z
end

Base.size(p::Policy) = size(p.xs)
Base.@propagate_inbounds Base.getindex(p::Policy, i::Int) =
	IntervalChoice(p.cutoffs[i], i+1>length(p.xs) ? p.ub : p.cutoffs[i+1], p.xs[i])

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

struct SqueezingPolicy{Z,A,AO,F1,F2,O,S,TR} <: CDCPSolver
	scdca::Bool
    pool::Vector{IntervalChoice{Z,A}}
    squeezing::Vector{Int}
    branching::Vector{Int}
    lookup_zero_margin::Dict{Tuple{Int,AO},Z}
    zero_margin::F1
    matcheds::Vector{Vector{IntervalChoice{Z,A}}}
    singlesolvers::Vector{S}
    equal_obj::F2
    obj2::O
    zero_margin_call::RefValue{Int}
    equal_obj_call::RefValue{Int}
    trace::TR
    nobranching::Bool
end

function _init(::Type{<:SqueezingPolicy}, obj, scdca::Bool, equal_obj, zbounds::Tuple{Z,Z};
		zero_margin=nothing, x0=nothing, ntasks=1, trace::Bool=false,
        nobranching::Bool=false, singlekw=NamedTuple(), kwargs...) where Z
	S = length(obj.x)
	if x0 === nothing
		if obj.x isa SVector
            allu = _fillstate(SVector{S,ItemState}, undetermined)
        	x = Policy([zbounds[1]], [allu], zbounds[2])
    	else
            allu = fill(undetermined, S)
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
    A = eltype(x.xs)
    matcheds = [IntervalChoice{Z,A}[] for _ in 1:ntasks]
    ss = [init(Squeezing, obj, S, scdca; z=zero(Z), singlekw...) for _ in 1:ntasks]
	return SqueezingPolicy(scdca, pool, collect(1:length(x.xs)), Int[],
        Dict{Tuple{Int,typeof(obj.x)},Z}(), zero_margin, matcheds, ss,
        equal_obj, obj2, Ref(0), Ref(0), tr, nobranching), x
end

_squeeze(x::IntervalChoice, s::ItemState, i::Int) =
    IntervalChoice(x.lb, x.ub, _squeeze(x.x, s, i))

_setitemstate(x::IntervalChoice, s::ItemState, i::Int) =
    IntervalChoice(x.lb, x.ub, _setitemstate(x.x, s, i))

function squeeze!(p::CDCProblem{<:SqueezingPolicy}, x::IntervalChoice, i::Int)
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
        return (x,)
    elseif z >= x.ub # Cannot include any part of the interval
        tr === nothing ||
            push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, undetermined))
        return squeeze_exclude!(p, x, i)
    else
        tr === nothing ||
            push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, undetermined))
        x1 = IntervalChoice(x.lb, z, x.x)
        xs = squeeze_exclude!(p, x1, i)
        x2 = IntervalChoice(z, x.ub, x.x)
        return xs..., x2
    end
end

function squeeze_exclude!(p::CDCProblem{<:SqueezingPolicy}, x::IntervalChoice, i::Int)
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
        return (x,)
    elseif z <= x.lb # Cannot exclude any part of the interval
        x = _setitemstate(x, aux, i) # Don't use _squeeze
        tr === nothing || push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, aux))
        return (x,)
    else
        tr === nothing || push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, excluded))
        x1 = IntervalChoice(x.lb, z, _squeeze(x.x, excluded, i))
        x2 = IntervalChoice(z, x.ub, _setitemstate(x.x, aux, i))
        return x1, x2
    end
end

function squeeze!(p::CDCProblem{<:SqueezingPolicy})
    pool, squeezing, branching = p.solver.pool, p.solver.squeezing, p.solver.branching
    while !isempty(squeezing)
        k = pop!(squeezing)
        x = pool[k]
        i = findfirst(==(undetermined), x.x)
        if i === nothing
            findfirst(==(aux), x.x) === nothing || push!(branching, k)
        else
            p.obj.fcall < p.maxfcall || return maxfcall_reached
	    	xs = squeeze!(p, x, i)
            pool[k] = xs[1]
            push!(squeezing, k)
            for j in 2:length(xs)
                push!(pool, xs[j])
                push!(squeezing, length(pool))
            end
        end
    end
    return inprogress
end

function _reset!(p::CDCProblem{<:Squeezing}, z, x)
    p.state = inprogress
    p.solver = Squeezing(p.solver.scdca, empty!(p.solver.branching), z, p.solver.trace)
    p.x = x
    p.fx = -Inf
    return p
end

function setx0scdcb(x0::SVector{S,ItemState}, xl::SVector{S,ItemState}) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(ifelse(xl[$i]==included, included, x0[$i])))
		end
		return :(SVector{S,ItemState}($ex))
	else
		return SVector{S,ItemState}(
            ntuple(i->ifelse(x[i]==included, included, x0[i]), S))
	end
end

function search!(p::CDCProblem{<:Squeezing}, matched, zl0, xl0, zr0, xr0, x0, obj2, equal_obj)
    zl, xl, zr, xr = zl0, xl0, zr0, xr0
    while true
        # xl and xr should always be different here
        obj1 = _setchoice(p.obj, setsub(xl))
        obj2 = _setchoice(obj2, setsub(xr))
        z, obj = equal_obj(obj1, obj2, zl, zr)
        # p.solver.equal_obj_call[] += 1
        if zl < z < zr # Additional cutoff points in between
            x0next = p.solver.scdca ? x0 : setx0scdcb(x0, xl)
            xnew = _solvesingle!(p, z, x0next)
            if xnew == xl || xnew == xr
                push!(matched, IntervalChoice(zl, z, xl), IntervalChoice(z, zr, xr))
                if zr != zr0 # Move to the interval on the right
                    zl, xl = zr, xr
                    zr, xr = zr0, xr0
                else
                    return
                end
            else # More cutoff points on the left
                zr, xr = z, xnew
            end
        else # No split of interval
            if z == zl
                push!(matched, IntervalChoice(zl, zr, xr))
            elseif z == zr
                push!(matched, IntervalChoice(zl, zr, xl))
            else
                error("z is not in between")
            end
            if zr != zr0 # Move to the interval on the right
                zl, xl = zr, xr
                zr, xr = zr0, xr0
            else
                return
            end
        end
    end
end

function _solvesingle!(p::CDCProblem{<:Squeezing}, z, x0)
    _reset!(p, z, x0)
    solve!(p)
    p.state == success || error("single-agenet solver fails with z = ", p.solver.z)
    return p.x
end

function branching!(p::CDCProblem{<:SqueezingPolicy}, k::Int, itask::Int)
    pool, matched = p.solver.pool, p.solver.matcheds[itask]
    x = pool[k]
    sp = p.solver.singlesolvers[itask]
    xl = _solvesingle!(sp, x.lb, x.x)
    xr = _solvesingle!(sp, x.ub, x.x)
    if xl == xr
        pool[k] = IntervalChoice(x.lb, x.ub, xl)
    else
        search!(sp, matched, x.lb, xl, x.ub, xr, x.x, p.solver.obj2, p.solver.equal_obj)
        # Overwrite the old interval with aux
        pool[k] = pop!(matched)
    end
end

function branching!(p::CDCProblem{<:SqueezingPolicy})
    pool, branching, matcheds = p.solver.pool, p.solver.branching, p.solver.matcheds
    ntasks = length(matcheds)
    for m in matcheds
        empty!(m) # Just to be safe
    end
    if ntasks == 1 # No multithreading
        for k in branching
            branching!(p, k, 1)
        end
    else
        @sync for itask in 1:ntasks
            Threads.@spawn for ik in itask:ntasks:length(branching)
                branching!(p, branching[ik], itask)
            end
        end
    end
    append!(pool, matcheds...)
    return inprogress
end

_lb(x::IntervalChoice) = x.lb

function concat!(p::CDCProblem{<:SqueezingPolicy})
    pool = p.solver.pool
    sort!(pool, by=_lb)
    cutoffs = resize!(p.x.cutoffs, 1)
    xs = resize!(p.x.xs, 1)
    cutoffs[1] = pool[1].lb
    xs[1] = xlast = pool[1].x
    for x in pool
        # Filter out potential singletons
        if x.lb < x.ub && x.x != xlast
            push!(cutoffs, x.lb)
            push!(xs, x.x)
            xlast = x.x
        end
    end
    return success
end

function solve!(p::CDCProblem{<:SqueezingPolicy}; restart::Bool=false)
	# restart && (p = _reinit!(p))
	p.state = squeeze!(p)
    if p.state == maxfcall_reached
        @warn "maxfcall is reached before convergence"
        return p
    end
    p.solver.nobranching || (p.state = branching!(p))
    concat!(p)
	return p
end

