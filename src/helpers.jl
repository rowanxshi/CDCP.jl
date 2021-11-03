function next_undetermined((sub, sup, aux)::Tuple{<: AbstractVector{Bool}, <: AbstractVector{Bool}, <: AbstractVector{Bool}})
	j = zero(first(eachindex(sub)))
	for j in eachindex(sub)
		(aux[j] || isequal(sub[j], sup[j])) && continue
		return j
	end
	return j
end

function containers(C::Int)
	sub = falses(C)
	sup = trues(C)
	aux = falses(C)
	
	(sub, sup, aux)
end