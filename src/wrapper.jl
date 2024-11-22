function _containers(C::Integer)
	sub = falses(C)
	sup = trues(C)
	aux = falses(C)
	return (sub, sup, aux)
end

struct Interval{T <: AbstractVector{Bool}, N <: Real}
	sub::T
	sup::T
	aux::T
	l::N
	r::N
end
Interval(J::AbstractVector{Bool}, l::Real, r::Real) =
    Interval(J, deepcopy(J), deepcopy(J), l, r)
Interval(Vs, l::Real, r::Real) = Interval(Vs..., l, r)

function _checklength(::SVector{S,ItemState}, sub, sup, aux) where S
    length(sub) == length(sup) == length(aux) == S ||
        throw(ArgumentError("length of sub, sup, aux are inconsistent with C"))
end

function _parse_state(sub::Bool, sup::Bool, isaux::Bool)
    if isaux
        return aux
    elseif sub != sup
        return undetermined
    elseif sub
        return included
    else
        return excluded
    end
end

function _parse_triplet(x::SVector{S,ItemState}, sub, sup, aux) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(_parse_state(
                @inbounds(sub[$i]), @inbounds(sup[$i]), @inbounds(aux[$i]))))
		end
		return :(_checklength(x, sub, sup, aux); SVector{S,ItemState}($ex))
	else
        _checklength(x, sub, sup, aux)
		return SVector{S,ItemState}(ntuple(i->_parse_state(sub[i], sup[i], aux[i]), S))
	end
end

function _invert_state!(sub::AbstractVector{Bool}, sup::AbstractVector{Bool},
        isaux::AbstractVector{Bool}, x::SVector{S,ItemState}) where S
    @inbounds for i in 1:S
        xi = x[i]
        if xi == aux
            isaux[i] = true
        else
            isaux[i] = false
            if xi == undetermined
                sub[i] = false
                sup[i] = true
            elseif xi == included
                sub[i] = true
                sup[i] = true
            else
                sub[i] = false
                sup[i] = false
            end
        end
    end
end

function D_j(obj)
	D_j_obj(j::Integer, J::AbstractArray{Bool}, z::Real = 1) = let obj = obj
		bool_j = J[j]
		J = setindex!(J, true, j)
		marg = obj(J, z)
		J = setindex!(J, false, j)
		marg -= obj(J, z)
		J = setindex!(J, bool_j, j)
		return marg
	end
end

function zero_D_j(equalise_obj, holder::AbstractVector{Bool})
	zero_D_j_obj(j::Integer, J::AbstractVector{Bool}, extras...) = let equalise_obj = equalise_obj, holder = holder
		J2 = copyto!(holder, J)
		J2 = setindex!(J2, !J[j], j)
		equalise_obj((J, J2), extras...)
	end
end

struct Wrapped_Equalise_Obj{F}
    f::F
end

function (w::Wrapped_Equalise_Obj)(obj1, obj2, l, r)
    obj1 = addfcall(obj1, 2)
    return w.f((obj1.x, obj2.x), l, r), obj1
end

struct Wrapped_Zero_D_j_Obj{F}
    f::F
    J::Vector{Bool}
end

function (w::Wrapped_Zero_D_j_Obj)(obj, i, lb, ub)
    copyto!(w.J, obj.x)
    z = w.f(i, w.J, lb, ub)
    obj = addfcall(obj, 1)
    return z, obj
end

function restart!((working, converged, done), C)
	empty!.((working, converged, done))
	int = Interval(_containers(C), -Inf, Inf)
	push!(working, int)
end

function naive!(J::AbstractVector{Bool}; obj, z=nothing)
    C = length(J)
    wobj = Objective(obj, copy(J))
    p = solve(BruteForce, wobj, C; z=z)
    copyto!(J, p.x)
end

naive(C::Integer; obj, z=nothing) = naive!(falses(C); obj, z)

# Wrappers for singleagent method
function solve!((sub, sup, aux); scdca::Bool, obj,
        containers = nothing, restart::Bool = true, z=nothing, kwargs...)
    C = length(sub)
    wobj = Objective(obj, SVector{C,Bool}(trues(C)))
    p = init(Squeezing, wobj, C, scdca; z=z, restart=restart)
    # Allow setting initial choice by (sub, sup, aux)
    restart || (p.x = _parse_triplet(p.x, sub, sup, aux))
    solve!(p)
    # Translate result
    _invert_state!(sub, sup, aux, p.x)
end

function solve(C::Integer; obj, kwargs...)
	wobj = Objective(obj, SVector{C,Bool}(trues(C)))
    p = solve(Squeezing, wobj, scdca; z=z, restart=restart)
    _invert_state!(sub, sup, aux, p.x)
end

# D_j_obj is not used and hence ignored
function policy!(cutoffspolicies, containers, C::Integer;
        scdca::Bool, obj, equalise_obj, zero_D_j_obj = zero_D_j(equalise_obj, falses(C)),
        show_time::Bool = false, emptyset = falses(C), restart::Bool = true,
        ntasks=1, trace::Bool=false,
        nobranching::Bool=false, singlekw=NamedTuple(), kwargs...)

    restart && restart!(containers, C)

    weq_obj = Wrapped_Equalise_Obj(equalise_obj)
    wobj = Objective(obj, SVector{C,Bool}(trues(C)))
    wzero_dj = Wrapped_Zero_D_j_Obj(zero_D_j_obj, fill(false, C))

    p = init(SqueezingPolicy, wobj, C, scdca, weq_obj, (-Inf, Inf);
        zero_margin=wzero_dj, ntasks=ntasks, trace=trace, nobranching=nobranching,
        singlekw=singlekw, kwargs...)

    pool, squeezing = p.solver.pool, p.solver.squeezing

    # Handle initial choices passed via (cutoffs, policies)
    _initialized = false
    if containers !== nothing
        working, converged, done = containers
        if !isempty(working)
            resize!(pool, length(working))
            resize!(squeezing, length(working))
            for (i, v) in enumerate(working)
                pool[i] = IntervalChoice(v.l, v.r,
                    _parse_triplet(pool[1].x, v.sub, v.sup, v.aux))
                squeezing[i] = i
            end
            _initialized = true
        end
    end
    if cutoffspolicies !== nothing
        if _initialized
            @warn "Both working and policies are non-empty; contents in policies are ignored"
        else
            cutoffs, policies = cutoffspolicies
            resize!(pool, length(policies))
            resize!(squeezing, length(policies))
            for i in 1:length(policies)
                pool[i] = IntervalChoice(cutoffs[i], cutoffs[i+1], policies[i])
                squeezing[i] = i
            end
        end
    else
        policies = Union{Nothing, BitVector}[]
	    cutoffs = Float64[-Inf, Inf]
    end

    solve!(p)

    # Copy results back
    resize!(cutoffs, length(p.x.cutoffs)+1)
    copyto!(cutoffs, p.x.cutoffs)
    cutoffs[end] = p.x.ub
    resize!(policies, length(p.x.xs))
    for i in eachindex(policies)
        policies[i] = BitVector(setsub(p.x.xs[i]))
    end
    return cutoffs, policies
end

function policy(C::Integer; kwargs...)
    working = Vector{Interval{BitVector, Float64}}(undef, 0)
	converged = similar(working)
	done = similar(working)
    policy!(nothing, (working, converged, done), C; kwargs...)
end
