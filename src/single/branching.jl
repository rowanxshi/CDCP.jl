function branch((sub, sup, aux))
	aux_in = fill!(copy(aux), false)
	aux = fill!(aux, false)
	
	j = next_undetermined((sub, sup, aux_in))
	
	sub_in = setindex!(copy(sub), true, j)
	sup_out = setindex!(copy(sub), false, j)
	
	return (sub_in, sup, aux_in), (sub, sup_out, aux)
end

function converge_branches!((working, converged); cdcp...)
	while !isempty(working)
		(sub, sup, aux) = pop!(working)
		if isequal(sub, sup)
			push!(converged, (sub, sup, aux))
			continue
		end
		(sub, sup, aux) = converge!((sub, sup, aux); cdcp...)
		isequal(sub, sup) ? push!(converged, (sub, sup, aux)) : append!(working, collect(branch((sub, sup, aux))))
	end
	return converged
end
