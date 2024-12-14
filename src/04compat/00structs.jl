struct Wrapped_Equalise_Obj{F}
	f::F
end

function (w::Wrapped_Equalise_Obj)(obj1, obj2, l, r)
	obj1 = addfcall(obj1, 2)
	return w.f((obj1.x, obj2.x), l, r), obj1
end

struct Wrapped_Zero_D_j_Obj{F}
	f::F
	J::Vector{Bool}
end

function (w::Wrapped_Zero_D_j_Obj)(obj, i, lb, ub)
	copyto!(w.J, obj.x)
	z = w.f(i, w.J, lb, ub)
	obj = addfcall(obj, 1)
	return z, obj
end

struct Interval{T <: AbstractVector{Bool}, N <: Real}
	sub::T
	sup::T
	aux::T
	l::N
	r::N
end
Interval(J::AbstractVector{Bool}, l::Real, r::Real) =
	Interval(J, deepcopy(J), deepcopy(J), l, r)
Interval(Vs, l::Real, r::Real) = Interval(Vs..., l, r)
