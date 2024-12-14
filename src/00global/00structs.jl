"""
    CDCPSolver

Abstract type for all solution algorithms for a [`CDCProblem`](@ref).
"""
abstract type CDCPSolver end

"""
    Objective{F,A}

A wrapped objective function for solving a [`CDCProblem`](@ref). This facilitates maintaining an internal interface for dealing with objective functions that is independent from user interface.

Users are *not* required to construct `Objective` unless there is a need for fine-grained control.

# Constructor
	Objective(f, x, [z=0])

Construct an instance of `Objective` with objective function `f` and an input vector `x`.
`f` must always accept `x` as the first argument
and may additionally accept an optional argument `z` for parameter
(e.g., productivity).
"""
struct Objective{F,A}
	f::F
	x::A
	fcall::Int
end

function Objective(f, x)
	Objective(f, x, 0)
end

function _setx(obj::Objective{<:Any,<:SVector}, v, i)
	Objective(obj.f, setindex(obj.x, v, i), obj.fcall)
end

# Fallback method assumes x is mutable
function _setx(obj::Objective, v, i)
	Objective(obj.f, setindex!(obj.x, v, i), obj.fcall)
end

"""
    addfcall(obj::Objective, n=1)

Add `n` to the counter for function call for `obj`.
"""
function addfcall(obj::Objective, n=1)
	Objective(obj.f, obj.x, obj.fcall+n)
end
function _clearfcall(obj::Objective)
	Objective(obj.f, obj.x, 0)
end

"""
    value(obj::Objective, z)

Evaluate `obj` with parameter `z`.
The parameter is ignored if `z` is `nothing`.
"""
function value(obj::Objective, z)
	obj.f(obj.x, z), addfcall(obj)
end
function value(obj::Objective, ::Nothing)
	obj.f(obj.x), addfcall(obj)
end

"""
    margin(obj::Objective, i::Int, z)

Evalutate the change in `obj` with optional parameter `z`
when the `i`th item is included or not.
This corresponds to the `D_j` function in earlier implementation.
"""
function margin(obj::Objective, i::Int, z)
	obj = _setx(obj, true, i)
	f1, obj = value(obj, z)
	obj = _setx(obj, false, i)
	f0, obj = value(obj, z)
	return f1, f0, obj
end

@enum SolverState::Int8 begin
	inprogress
	success
	maxfcall_reached
end

@enum ItemState::Int8 begin
	undetermined
	included
	excluded
	aux
end

"""
    CDCProblem{M<:CDCPSolver, O<:Objective, A, TF<:AbstractFloat}

Results from solving a combinatorial discrete choice problem.
When a solution is attained, it can be retrived from the field `x`.
"""
mutable struct CDCProblem{M<:CDCPSolver, O<:Objective, A, TF<:AbstractFloat}
	solver::M
	obj::O
	x::A
	fx::TF
	maxfcall::Int
	state::SolverState
end
