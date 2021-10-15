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
function branch(int::interval)
	((a,b,c), (d,e,f)) = branch((int.sub, int.sup, int.aux))
	return interval(a, b, c, int.l, int.r), interval(d, e, f, int.l, int.r)
end