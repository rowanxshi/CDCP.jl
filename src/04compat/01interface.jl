"""
    naive(C::Integer; obj)

Solve a combinatorial discrete choice problem with objective function `π(J)` over `C` choices with exhaustion.

!!! warning

    This method exists only for the sake of backward compatibility. Future use should prefer the interface based on `solve(Naive, ...)`.
"""
function naive(C::Integer; obj, z=nothing)
	naive!(falses(C); obj, z)
end
"""
    naive!(J::AbstractVector{Bool}; obj)

Similar to [`naive`](@ref), but accepts a pre-allocated vector `J`.

!!! warning

    This method exists only for the sake of backward compatibility. Future use should prefer the interface based on `solve(Naive, ...)`.
"""
function naive!(J::AbstractVector{Bool}; obj)
	Base.depwarn("Consider the new interface for solution with exhaustion using the problem with `Naive`", :naive!)
	C = length(J)
	wrapped_obj = Objective(obj, copy(J))
	cdcp = solve(Naive, wrapped_obj, C)
	copyto!(J, cdcp.x)
end

"""
    solve(C::Integer; scdca::Bool, obj)

Solve a combinatorial discrete choice problem over `C` choices with SCD-C from above if `scdca` is `true` (otherwise, from below).

The solver uses the objective function `obj(J)` which must accept as argument a Boolean vector of length `C`.

!!! warning

    This method exists only for the sake of backward compatibility. Future use should prefer the interface based on `solve(Squeezing, ...)`.
"""
function solve(C::Integer; scdca::Bool, obj)
	wrapped_obj = Objective(obj, SVector{C,Bool}(trues(C)))
	return solve(Squeezing, wrapped_obj, scdca)
end
"""
    solve!((sub, sup, aux); scdca::Bool, obj, restart::Bool=true)

Solve in-place a combinatorial discrete choice problem with SCD-C from above if `scdca` is true (otherwise, from below). The solver uses the objective function `obj(J)` and initiates using the Boolean vectors `(sub, sup, aux)` if `restart=false`. The objective function `obj` must accept as argument a Boolean vector with length corresponding to the number of items in the problem.

!!! warning

    This method exists only for the sake of backward compatibility. Future use should prefer the interface based on `solve(Squeezing, ...)`.
"""
function solve!((sub, sup, aux); scdca::Bool, obj, restart::Bool=true, kwargs...)
	Base.depwarn("Consider the new interface for solving the single-agent problem with `Squeezing`", :solve!)
	C = length(sub)
	wrapped_obj = Objective(obj, SVector{C,Bool}(trues(C)))
	cdcp = init(Squeezing, wrapped_obj, C, scdca)
	# allow setting initial choice by (sub, sup, aux)
	restart || (cdcp.x = invert_triplet(cdcp.x, sub, sup, aux))
	solve!(cdcp)
	# translate result
	invert_state!(sub, sup, aux, cdcp.x)
	return sup
end

"""
    policy(C::Integer; scdca::Bool, obj, equalise_obj, [zero_D_j_obj])

Find the policy function for a combinatorial discrete choice problem over `C` choices and one dimension of heterogeneity with SCD-C from above if `scdca` is `true` (otherwise, from below) and SCD-T. The solver uses the objective function `obj(J, z)`, which must accept as argument a Boolean vector with length `C` and the heterogeneous component `z` and the `equalise_obj((J1, J2), l, r)` function, which identifies the `z` where the agent is indifferent between the pair `(J1, J2)`. The routine provides the interval `[l, r]` along which this search is performed; if the marginal type is not in the interval, it is sufficient to return `nothing`.

The solver can optionally take `zero_D_j_obj(j, J, l, r)`, a user-supplied function that identifies the `z` where the marginal value of item `j` to set `J` is zero. The solver provides the interval `[l, r]` within which this marginal type is located. Please return `NaN` instead of `nothing` if the marginal type is not within the interval. If not provided, the solver automatically constructs these using the `equalise_obj` function.

!!! warning

    This method exists only for the sake of backward compatibility. Future use should prefer the interface based on `solve(SqueezingPolicy, ...)`.
"""
function policy(C::Integer; kwargs...)
	working = Vector{Interval{BitVector, Float64}}(undef, 0)
	converged = similar(working)
	done = similar(working)
	cutoffs = Float64[-Inf, Inf]
	policies = Union{Nothing, BitVector}[]
	policy!((cutoffs, policies), (working, converged, done), C; kwargs...)
end
"""
    policy!((cutoffs, policies), containers0, C::Integer; scdca::Bool, obj, equalise_obj, zero_D_j_obj=zero_D_j(equalise_obj, falses(C)), restart::Bool=true)

The similar to [`policy`](@ref), but the solver initiates with `working0` in `containers0 = (working0, _, _)` if `restart=false`.

!!! warning

    This method exists only for the sake of backward compatibility. Future use should prefer the interface based on `solve(SqueezingPolicy, ...)`.
"""
function policy!((cutoffs, policies), containers, C::Integer; scdca::Bool, obj, equalise_obj, zero_D_j_obj = zero_D_j(equalise_obj, falses(C)), restart::Bool=true, singlekw=NamedTuple(), kwargs...)
	Base.depwarn("Consider the new interface for solving the policy problem with `SqueezingPolicy`", :policy!)

	restart && restart!(containers, C)
	wrapped_eq_obj = Wrapped_Equalise_Obj(equalise_obj)
	wrapped_obj = Objective(obj, SVector{C,Bool}(trues(C)))
	wrapped_zero_dj = Wrapped_Zero_D_j_Obj(zero_D_j_obj, fill(false, C))
	cdcp = init(SqueezingPolicy, wrapped_obj, C, scdca, wrapped_eq_obj; zero_margin=wrapped_zero_dj, singlekw, kwargs...)
	invert_cutoffspolicies!(cdcp, containers) # handle initial choices passed via containers

	solve!(cdcp)
	resize!(cutoffs, length(cdcp.x.cutoffs)+1)
	copyto!(cutoffs, cdcp.x.cutoffs)
	cutoffs[end] = cdcp.x.zright
	resize!(policies, length(cdcp.x.itemstates_s))
	for i in eachindex(policies)
		policies[i] = BitVector(to_sub(cdcp.x.itemstates_s[i]))
	end
	return cutoffs, policies
end
