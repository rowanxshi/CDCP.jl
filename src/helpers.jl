function next_undetermined((sub, sup, aux))
	j = zero(first(eachindex(sub)))
	for j in eachindex(sub)
		(aux[j] || isequal(sub[j], sup[j])) && continue
		return j
	end
	return j
end

function _containers(C::Integer)
	sub = falses(C)
	sup = trues(C)
	aux = falses(C)
	
	(sub, sup, aux)
end