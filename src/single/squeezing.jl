function update!((sub, sup, aux), j::Integer; D_j_obj, scdca::Bool)
	# for excluding: look at best case
	exclude = (scdca && D_j_obj(j, sub) ≤ 0) || (!scdca && D_j_obj(j, sup) < 0)
	if exclude
		sup = setindex!(sup, false, j)
		aux = fill!(aux, false)
		return (sub, sup, aux)
	end

	# for including: look at worst case
	include = (scdca && D_j_obj(j, sup) > 0) || (!scdca && D_j_obj(j, sub) ≥ 0)
	if include
		sub = setindex!(sub, true, j)
		aux = fill!(aux, false)
		return (sub, sup, aux)
	end
	
	aux = setindex!(aux, true, j)
	return (sub, sup, aux)
end

function converge!(Vs; D_j_obj, scdca::Bool)
	@inbounds while true
		j = next_undetermined(Vs)
		iszero(j) && break
		Vs = update!(Vs, j; D_j_obj, scdca)
	end
	
	Vs
end
