function branch((sub, sup, aux)::NTuple{3, Any})
	aux_in = similar(sub)
	fill!(aux_in, false)
	fill!(aux, false)
	
	j = next_undetermined((sub, sup, aux_in))
	
	sub_in = similar(sub)
	copyto!(sub_in, sub)
	sub_in[j] = true
	
	sup_out = similar(sup)
	copyto!(sup_out, sup)
	sup_out[j] = false
	
	return (sub_in, sup, aux_in), (sub, sup_out, aux)
end

function converge_branches!((working, converged), D_j_π::F, scdca::Bool) where F <: Function
	while !isempty(working)
		(sub, sup, aux) = pop!(working)
		if isequal(sub, sup)
			push!(converged, (sub, sup, aux))
			continue
		end
		converge!((sub, sup, aux), D_j_π, scdca)
		isequal(sub, sup) ? push!(converged, (sub, sup, aux)) : append!(working, collect(branch((sub, sup, aux))))
	end
	return converged
end
