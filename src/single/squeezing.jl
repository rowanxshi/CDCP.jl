function update!((sub, sup, aux), j::Integer; D_j_obj, scdca::Bool)
	# for excluding: look at best case
	exclude = if scdca
		Djbest = D_j_obj(j, sub)
		signbit(Djbest) || iszero(Djbest)
	else
		signbit(D_j_obj(j, sup))
	end
	if exclude
		sup = setindex!(sup, false, j)
		aux = fill!(aux, false)
		return (sub, sup, aux)
	end

	# for including: look at worst case
	include = if scdca
		Djworst = D_j_obj(j, sup)
		!signbit(Djworst) && !iszero(Djwosrt)
	else
		!signbit(D_j_obj(j, sub))
	end
	if include
		sub = setindex!(sub, true, j)
		aux = fill!(aux, false)
		return (sub, sup, aux)
	end
	
	aux = setindex!(aux, true, j)
	return (sub, sup, aux)
end

function converge!(Vs; kw...)
	@inbounds while true
		j = next_undetermined(Vs)
		iszero(j) && break
		Vs = update!(Vs, j; kw...)
	end
	
	Vs
end
