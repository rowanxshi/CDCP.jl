"""
    naive(C::Integer; obj)

Solve a combinatorial discrete choice problem over `C` choices with exhaustion. (Generally used for testing or time-trial exercises.) The solver an objective function `π(J)`.

See also: [`naive!`](@ref), [`solve!`](@ref), [`policy`](@ref)

!!! warning

	This method exists only for the sake of backward compatibility.
	Future use should prefer the interface based on `solve(Naive, ...)`.
"""
function naive(C::Integer; obj, z=nothing)
	naive!(falses(C); obj, z)
end

function naive!(J::AbstractVector{Bool}; obj, z=nothing)
	Base.depwarn("Consider the new interface for solution with exhaustion using the problem with `Naive`", :naive!)
	C = length(J)
	wobj = Objective(obj, copy(J))
	cdcp = solve(Naive, wobj, C; z=z)
	copyto!(J, cdcp.x)
end

"""
    solve(C::Integer; scdca::Bool, obj, [; containers])

Solve a combinatorial discrete choice problem over `C` choices with SCD-C from above if `scdca` is `true` (otherwise, from below). The solver uses the objective function `obj(J)` which must accept as argument a Boolean vector with length corresponding to the number of items in the problem.

See also: [`solve!`](@ref), [`policy`](@ref)

!!! warning

	This method exists only for the sake of backward compatibility.
	Future use should prefer the interface based on `solve(Squeezing, ...)`.
"""
function solve(C::Integer; obj, kwargs...)
	wobj = Objective(obj, SVector{C,Bool}(trues(C)))
	return solve(Squeezing, wobj, scdca; z=z, restart=restart)
end

function solve!((sub, sup, aux); scdca::Bool, obj, D_j_obj = nothing, containers = nothing, restart::Bool = true, z=nothing, kwargs...)
	Base.depwarn("Consider the new interface for solving the single-agent problem with `Squeezing`", :solve!)
	D_j_obj===nothing || @warn "D_j_obj doesn't need to be specified explicitly"
	C = length(sub)
	wobj = Objective(obj, SVector{C,Bool}(trues(C)))
	cdcp = init(Squeezing, wobj, C, scdca; z=z, restart=restart)
	# Allow setting initial choice by (sub, sup, aux)
	restart || (cdcp.x = _parse_triplet(cdcp.x, sub, sup, aux))
	solve!(cdcp)
	# Translate result
	_invert_state!(sub, sup, aux, cdcp.x)
	return sup
end

"""
    policy(C::Integer; scdca::Bool, obj, equalise_obj, [zero_D_j_obj])

Find the policy function for a combinatorial discrete choice problem over `C` choices and one dimension of heterogeneity with SCD-C from above if `scdca` is `true` (otherwise, from below) and SCD-T. The solver uses the objective function `obj(J, z)`, which must accept as argument a Boolean vector with length `C` and the heterogeneous component `z` and the `equalise_obj((J1, J2), l, r)` function, which identifies the `z` where the agent is indifferent between the pair `(J1, J2)`. The routine provides the interval `[l, r]` along which this search is performed; if the marginal type is not in the interval, it is sufficient to return `nothing`.

The solver can optionally take `zero_D_j_obj(j, J, l, r)`, a user-supplied function that identifies the `z` where the marginal value of item `j` to set `J` is zero. The solver provides the interval `[l, r]` within which this marginal type is located. Please return `NaN` instead of `nothing` if the marginal type is not within the interval. If not provided, the solver automatically constructs these using the `equalise_obj` function.

See also: [`solve!`](@ref), [`solve`](@ref)

!!! warning

	This method exists only for the sake of backward compatibility.
	Future use should prefer the interface based on `solve(SqueezingPolicy, ...)`.
"""
function policy(C::Integer; kwargs...)
	working = Vector{Interval{BitVector, Float64}}(undef, 0)
	converged = similar(working)
	done = similar(working)
	policy!(nothing, (working, converged, done), C; kwargs...)
end

# D_j_obj is not used and hence ignored
function policy!(cutoffspolicies, containers, C::Integer; scdca::Bool, obj, equalise_obj, D_j_obj=nothing, zero_D_j_obj = zero_D_j(equalise_obj, falses(C)), show_time::Bool = false, emptyset = falses(C), restart::Bool = true, ntasks=1, nobranching::Bool=false, singlekw=NamedTuple(), kwargs...)
	Base.depwarn("Consider the new interface for solving the policy problem with `SqueezingPolicy`", :policy!)

	D_j_obj===nothing || @warn "D_j_obj doesn't need to be specified explicitly"

	restart && restart!(containers, C)

	weq_obj = Wrapped_Equalise_Obj(equalise_obj)
	wobj = Objective(obj, SVector{C,Bool}(trues(C)))
	wzero_dj = Wrapped_Zero_D_j_Obj(zero_D_j_obj, fill(false, C))

	cdcp = init(SqueezingPolicy, wobj, C, scdca, weq_obj, (-Inf, Inf); zero_margin=wzero_dj, ntasks=ntasks, nobranching=nobranching, singlekw=singlekw, kwargs...)

	intervalchoices, squeezing = cdcp.solver.intervalchoices, cdcp.solver.squeezing_indices

	# Handle initial choices passed via (cutoffs, policies)
	_initialized = false
	if containers !== nothing
		working, converged, done = containers
		if !isempty(working)
			resize!(intervalchoices, length(working))
			resize!(squeezing, length(working))
			for (i, v) in enumerate(working)
				intervalchoices[i] = IntervalChoice(v.l, v.r, _parse_triplet(intervalchoices[1].itemstates, v.sub, v.sup, v.aux))
				squeezing[i] = i
			end
			_initialized = true
		end
	end
	policies = Union{Nothing, BitVector}[]
	cutoffs = Float64[-Inf, Inf]
	if cutoffspolicies !== nothing
		if !_initialized
			cutoffs, policies = cutoffspolicies
			resize!(intervalchoices, length(policies))
			resize!(squeezing, length(policies))
			for i in 1:length(policies)
				intervalchoices[i] = IntervalChoice(cutoffs[i], cutoffs[i+1], policies[i])
				squeezing[i] = i
			end
		end  
	end

	solve!(cdcp)

	# Copy results back
	resize!(cutoffs, length(cdcp.x.cutoffs)+1)
	copyto!(cutoffs, cdcp.x.cutoffs)
	cutoffs[end] = cdcp.x.zright
	resize!(policies, length(cdcp.x.itemstates_s))
	for i in eachindex(policies)
		policies[i] = BitVector(to_sub(cdcp.x.itemstates_s[i]))
	end
	return cutoffs, policies
end
