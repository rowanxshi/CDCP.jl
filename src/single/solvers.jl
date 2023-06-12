"""
	solve!((sub, sup, aux); scdca::Bool, obj, [D_j_obj; containers])

Solve in-place a combinatorial discrete choice problem with SCD-C from above if `scdca` is `true` (otherwise, from below). The solver uses the preallocated Boolean vectors `(sub, sup, aux)` as well as the objective function `obj(J)`. The objective function `obj` must accept as argument a Boolean vector with length corresponding to the number of items in the problem.

The solver can optionally take `D_j_obj(j, J)`, a user-supplied marginal value function; otherwise it will construct one automatically given `π`. It may also optionally take preallocated `containers = (working, converged)`, where `working` and `converged` are both `Vectors` with element type matching `(sub, sup, aux)`. These are used for the branching step and will be automatically allocated if not supplied.

See also: [`solve`](@ref), [`policy`](@ref)
"""
function solve!(Vs; scdca::Bool, obj, D_j_obj = D_j(obj), containers = _containers(Vs), restart::Bool = true)
	sub, sup, aux = Vs
	if restart
		sub = fill!(sub, false)
		sup = fill!(sup, true)
		aux = fill!(aux, false)
	end

	# squeeze
	(sub, sup, aux) = converge!((sub, sup, aux); D_j_obj, scdca)
	isequal(sub, sup) && return sub

	# if no convergence from squeezing, branch
	(working, converged) = containers
	empty!(working)
	append!(working, collect(branch((sub, sup, aux))))
	empty!(converged)
	converge_branches!((working, converged); D_j_obj, scdca)

	# among results in converged, choose best
	i_argmax = 0; max_prof = -Inf
	for (i, option) in enumerate(converged)
		prof = obj(first(option))
		if prof > max_prof
			max_prof = prof
			i_argmax = i
		end
	end
	argmax = first(converged[i_argmax])
	sub = copyto!(sub, argmax)
	sup = copyto!(sup, argmax)
end

"""
	solve(C::Integer; scdca::Bool, obj, [D_j_obj; containers])

Solve a combinatorial discrete choice problem over `C` choices with SCD-C from above if `scdca` is `true` (otherwise, from below). The solver uses the objective function `obj(J)` which must accept as argument a Boolean vector with length corresponding to the number of items in the problem.

The solver can optionally take `D_j_obj(j, J)`, a user-supplied marginal value function; otherwise it will construct one automatically given `obj`. It may also optionally take preallocated `containers = (working, converged)`, where `working` and `converged` are both `Vectors` with element type matching `(sub, sup, aux)`. These are used for the branching step and will be automatically allocated if not supplied.

See also: [`solve!`](@ref), [`policy`](@ref)
"""
function solve(C::Integer; obj, D_j_obj = D_j(obj), containers = _containers(_containers(C))) 
	Vs = isempty(first(containers)) ? _containers(C) : first(first(containers))
	solve!(Vs; obj, D_j_obj, containers)
end

"""
	naive!(J::AbstractVector{Bool}; obj)

Solve in-place a combinatorial discrete choice problem with simple brute force. (Generally used for testing or time-trial exercises.) The solver expectes a pre-allocated Boolean vector `J` and objective function `obj(J)`.

See also: [`naive`](@ref), [`solve!`](@ref), [`policy`](@ref)
"""
function naive!(J::AbstractVector{Bool}; obj)
	i_max = -1; max = -Inf
	
	for n in 0:(2^length(J)-1)
		J = digits!(J, n, base = 2)
		obj(J) > max && begin
			max = obj(J)
			i_max = n
		end
	end
	
	J = digits!(J, i_max, base = 2)
end

"""
	naive(C::Integer; obj)

Solve a combinatorial discrete choice problem over `C` choices with simple brute force. (Generally used for testing or time-trial exercises.) The solver an objective function `π(J)`.

See also: [`naive!`](@ref), [`solve!`](@ref), [`policy`](@ref)
"""
naive(C::Integer; obj) = naive!(falses(C); obj)
