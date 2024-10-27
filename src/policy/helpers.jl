struct Interval{T <: AbstractVector{Bool}, N <: Real}
	sub::T
	sup::T
	aux::T
	l::N
	r::N
end
Interval(J::AbstractVector{Bool}, l::Real, r::Real) = Interval(J, deepcopy(J), deepcopy(J), l, r)
Interval(Vs, l::Real, r::Real) = Interval(Vs..., l, r)
int_isless(int1::Interval, int2::Interval) = int1.l < int2.l
next_undetermined(int::Interval) = next_undetermined((int.sub, int.sup, int.aux))

function zero_D_j(equalise_obj, holder::AbstractVector{Bool})
	zero_D_j_obj(j::Integer, J::AbstractVector{Bool}, extras...) = let equalise_obj = equalise_obj, holder = holder
		J2 = copyto!(holder, J)
		J2 = setindex!(J2, !J[j], j)
		equalise_obj((J, J2), extras...)
	end
end

function branch(int::Interval)
	((a,b,c), (d,e,f)) = branch((int.sub, int.sup, int.aux))
	return Interval(a, b, c, int.l, int.r), Interval(d, e, f, int.l, int.r)
end

middle(int::Interval) = middle(int.l, int.r)
function middle(l, r)
	all(isinf, (l, r)) && return zero(l)
	all(isfinite, (l, r)) && return (l + r)/2
	isfinite(l) && return 2l
	isfinite(r) && return r/2
end
