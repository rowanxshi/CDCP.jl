function next_undetermined((sub, sup, aux)::Tuple{<: AbstractVector{Bool}, <: AbstractVector{Bool}, <: AbstractVector{Bool}})
	j = zero(first(eachindex(sub)))
	for j in eachindex(sub)
		(aux[j] || isequal(sub[j], sup[j])) && continue
		return j
	end
	return j
end

function D_j(π::Function)
	function D_j_π(j::Integer, J::A) where A <: AbstractArray{Bool}
		bool_j = J[j]
		J[j] = true
		marg = π(J)
		J[j] = false
		marg -= π(J)
		J[j] = bool_j
		return marg
	end
end

