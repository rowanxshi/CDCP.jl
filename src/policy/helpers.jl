struct interval{N1 <: Real, N2 <: Real}
	sub::BitVector
	sup::BitVector
	aux::BitVector
	l::N1
	r::N2
end
interval(J::AbstractVector{Bool}, l::Real, r::Real) = interval(J, J, J, l, r)
interval(Vs, l::Real, r::Real) = interval(Vs..., l, r)
int_isless(int1::interval, int2::interval) = int1.l < int2.l
next_undetermined(int::interval) = next_undetermined((int.sub, int.sup, int.aux))

function zero_D_j(equalise_obj, holder::AbstractVector{Bool})
	zero_D_j_obj(j::Integer, J::AbstractVector{Bool}, extras...) = let equalise_obj = equalise_obj, holder = holder
		J2 = copyto!(holder, J)
		J2 = setindex!(J2, !J[j], j)
		equalise_obj((J, J2), extras...)
	end
end

function branch(int::interval)
	((a,b,c), (d,e,f)) = branch((int.sub, int.sup, int.aux))
	return interval(a, b, c, int.l, int.r), interval(d, e, f, int.l, int.r)
end