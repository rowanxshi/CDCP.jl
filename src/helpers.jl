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

struct interval{N1 <: Real, N2 <: Real}
	sub::BitVector
	sup::BitVector
	aux::BitVector
	l::N1
	r::N2
end
interval(J::BitVector, l::N1, r::N2) where {N1 <: Real, N2 <: Real} = interval(J, J, J, l, r)
function int_isless(int1::interval, int2::interval)
	int1.l < int2.l
end
function next_undetermined(int::interval)
	next_undetermined((int.sub, int.sup, int.aux))
end
