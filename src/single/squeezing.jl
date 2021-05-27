function update!((sub, sup, aux)::NTuple{3, Any}, j::Integer, D_j_π::Function, scdca::Bool)
	# for excluding: look at best case
	exclude = (scdca && D_j_π(j, sub) ≤ 0.0) || (!scdca && D_j_π(j, sup) < 0.0)
	if exclude
		sup[j] = false
		aux .= false
		return true
	end

	# for including: look at worst case
	include = (scdca && D_j_π(j, sup) > 0.0) || (!scdca && D_j_π(j, sub) ≥ 0.0)
	if include
		sub[j] = true
		aux .= false
		return true
	end
	
	aux[j] = true
	return false
end

function converge!((sub, sup, aux)::NTuple{3, Any}, D_j_π::Function, scdca::Bool)
	converged = false
	
	@inbounds while !converged
		j = next_undetermined((sub, sup, aux))
		iszero(j) && break
		update!((sub, sup, aux), j, D_j_π, scdca)
	end
	
	(sub, sup, aux)
end
