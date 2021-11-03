struct interval{N1 <: Real, N2 <: Real}
	sub::BitVector
	sup::BitVector
	aux::BitVector
	l::N1
	r::N2
end
interval(J::BitVector, l::N1, r::N2) where {N1 <: Real, N2 <: Real} = interval(J, J, J, l, r)
interval(triad::NTuple{3, BitVector}, l::N1, r::N2) where {N1 <: Real, N2 <: Real} = interval(triad[1], triad[2], triad[3], l, r)
int_isless(int1::interval, int2::interval) =int1.l < int2.l
next_undetermined(int::interval) = next_undetermined((int.sub, int.sup, int.aux))

function zero_D_j(C::Int, equalise_π::F) where F <: Function
	holder = falses(C)
	function zero_D_j_π(j::Integer, J::AbstractVector{Bool})
		holder .= J
		bool_j = J[j]
		J[j] = true
		holder[j] = false
		z = equalise_π((J, holder))
		J[j] = bool_j
		return z
	end
	
	zero_D_j_π
end

function branch(int::interval)
	((a,b,c), (d,e,f)) = branch((int.sub, int.sup, int.aux))
	return interval(a, b, c, int.l, int.r), interval(d, e, f, int.l, int.r)
end