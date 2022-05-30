"""
	policy(C::Integer; scdca::Bool, obj, equalise_obj, [D_j_obj, zero_D_j_obj])

Find the policy function for a combinatorial discrete choice problem over `C` choices and one dimension of heterogeneity with SCD-C from above if `scdca` is `true` (otherwise, from below) and SCD-T. The solver uses the objective function `obj(J, z)`, which must accept as argument a Boolean vector with length `C` and the heterogeneous component `z` and the `equalise_obj((J1, J2), l, r)` function, which identifies the `z` where the agent is indifferent between the pair `(J1, J2)`. The routine provides the interval `[l, r]` along which this search is performed; if the marginal type is not in the interval, it is sufficient to return `nothing`.

The solver can optionally take `D_j_obj(j, J, z)` and `zero_D_j_obj(j, J, l, r)`, user-supplied functions. The first is the marginal value function, while the second identifies the `z` where the marginal value of item `j` to set `J` is zero. For the second, the solver provides the interval `[l, r]` within which this marginal type is located. It is sufficient to return `nothing` is the marginal type is not within the interval. If not provided, the solver automatically constructs these using the `equalise_obj` function.

See also: [`solve!`](@ref), [`solve`](@ref)
"""
function policy(C::Integer; e...)
	working = Vector{interval{Float64, Float64}}(undef, 0)
	converged = similar(working)
	done = similar(working)
	
	policies = Vector{Union{Nothing, BitVector}}(nothing, 1)
	cutoffs = Float64[-Inf, Inf]
	
	policy!((cutoffs, policies), (working, converged, done), C; e...)
end

function policy!((cutoffs, policies), (working, converged, done), C::Integer; scdca::Bool, obj, equalise_obj, D_j_obj = D_j(obj), zero_D_j_obj = zero_D_j(equalise_obj, falses(C)), show_time::Bool = false, emptyset = falses(C))
	empty!.((working, converged, done))
	int = interval(_containers(C), -Inf, Inf)
	push!(working, int)
	
	cdcp = (; scdca, obj, D_j_obj, zero_D_j_obj, equalise_obj, emptyset)
	
	# initial converge
	t1 = @elapsed converge!(working, converged; cdcp...)
	show_time && @info "initial converge: $t1"
	
	# branch if necessary: local optimal stored in done
	t2 = @elapsed converge_branches!(working, converged, done; cdcp...)
	show_time && @info "branching: $t2"
	
	# brute force among local optima; final policy function stored in converged
	policy = (cutoffs, policies)
	
	t3 = @timed brute!(policy, working, converged, done; cdcp...)
	show_time && @info "brute forcing: $t3"
	
	any(isnothing, policies) && @warn("Some intervals do not have associated policies.")
	policy
end