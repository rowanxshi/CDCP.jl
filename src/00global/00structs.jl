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
	Objective(f, ℒ, [z=0])

Construct an instance of `Objective` with objective function `f` and an input vector `ℒ`. `f` must always accept `ℒ` as the first argument and may additionally accept an optional argument `z` for parameter (e.g., productivity).
"""
struct Objective{F,A}
	f::F
	ℒ::A
	fcall::Int
end

function Objective(f, ℒ)
	Objective(f, ℒ, 0)
end

function _setℒ(obj::Objective{<:Any,<:SVector}, v, i)
	Objective(obj.f, setindex(obj.ℒ, v, i), obj.fcall)
end

# Fallback method assumes ℒ is mutable
function _setℒ(obj::Objective, v, i)
	Objective(obj.f, setindex!(obj.ℒ, v, i), obj.fcall)
end

"""
    addfcall(obj::Objective, n=1)

Add `n` to the counter for function call for `obj`.
"""
function addfcall(obj::Objective, n=1)
	Objective(obj.f, obj.ℒ, obj.fcall+n)
end
function _clearfcall(obj::Objective)
	Objective(obj.f, obj.ℒ, 0)
end

function (obj::Objective)(z)
	obj.f(obj.ℒ, z), addfcall(obj)
end
function (obj::Objective)(::Nothing)
	obj.f(obj.ℒ), addfcall(obj)
end

"""
    margin(obj::Objective, i::Int, z)

Evalutate the change in `obj` with optional parameter `z`
when the `i`th item is included or not.
This corresponds to the `D_j` function in earlier implementation.
"""
function margin(obj::Objective, i::Int, z)
	obj = _setℒ(obj, true, i)
	f1, obj = obj(z)
	obj = _setℒ(obj, false, i)
	f0, obj = obj(z)
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
    CDCProblem{M<:CDCPSolver, O<:Objective, A, F<:AbstractFloat}

Results from solving a combinatorial discrete choice problem.
When a solution is attained, it can be retrived from the field `x`.
"""
mutable struct CDCProblem{M<:CDCPSolver, O<:Objective, A, F<:AbstractFloat}
	solver::M
	obj::O
	x::A
	value::F
	maxfcall::Int
	state::SolverState
end
