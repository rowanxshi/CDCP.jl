"""
    CDCPSolver

Abstract type for all solution algorithms for a [`CDCProblem`](@ref).
"""
abstract type CDCPSolver end

"""
    Objective{F,V <: AbstractVector}

A wrapped objective function for solving a [`CDCProblem`](@ref). This facilitates maintaining an internal interface for dealing with objective functions that is independent from user interface.

Users are *not* required to construct `Objective` unless there is a need for fine-grained control.

# Constructor
	Objective(f, ℒ, fcall = 0)

Construct an instance of `Objective` with objective function `f` and an input vector `ℒ`. `f` must always accept `ℒ` as the first argument and may additionally accept an optional argument `z` for parameter (e.g., productivity).
"""
struct Objective{F,V <: AbstractVector}
	f::F
	ℒ::V
	fcall::Int
end

function Objective(f, ℒ)
	Objective(f, ℒ, 0)
end

function setℒ(obj::Objective{<:Any,V}, ℒ::V) where V
	Objective(obj.f, ℒ, obj.fcall)
end

function StaticArrays.setindex(obj::Objective{<:Any,<:SVector}, value, index)
	Objective(obj.f, setindex(obj.ℒ, value, index), obj.fcall)
end
function StaticArrays.setindex(obj::Objective, value, index) # fallback method assumes ℒ is mutable
	Objective(obj.f, setindex!(obj.ℒ, value, index), obj.fcall)
end

function init_Objective(obj, S::Integer)
	S = Int(S)
	(S > 0) || throw(ArgumentError("the number of items S must be a positive integer"))
	if obj isa Objective
		obj = clearfcall(obj)
		(length(obj.ℒ) == S) || throw(ArgumentError("length of obj.ℒ is not $S"))
	else
		ℒ = (S < static_threshold()) ? SVector{S, Bool}(ntuple(i->false, S)) : Vector{Bool}(undef, S)
		obj = Objective(obj, ℒ)
	end
	obj
end

# default threshold for determining whether SVector is used for a choice
function static_threshold()
	256
end

"""
    addfcall(obj::Objective, n=1)

Add `n` to the counter for function call for `obj`.
"""
function addfcall(obj::Objective, n=1)
	Objective(obj.f, obj.ℒ, obj.fcall+n)
end
function clearfcall(obj::Objective)
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

Evalutate the change in `obj` with optional parameter `z` when the `i`th item is included or not. This corresponds to the `D_j` function in earlier implementation.
"""
function margin(obj::Objective, i::Int, z)
	obj = setindex(obj, true, i)
	value1, obj = obj(z)
	obj = setindex(obj, false, i)
	value0, obj = obj(z)
	return value1, value0, obj
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
    CDCProblem{M<:CDCPSolver, O<:Objective, T, F<:AbstractFloat}

Results from solving a combinatorial discrete choice problem. When a solution is attained, it can be retrieved from the field `x`.
"""
mutable struct CDCProblem{M<:CDCPSolver, O<:Objective, T, F<:AbstractFloat}
	solver::M
	obj::O
	x::T
	value::F
	state::SolverState
end
