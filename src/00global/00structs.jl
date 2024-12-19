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
	Objective(f, ℒ, fcall = 0)

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

function _setchoice(obj::Objective{<:Any,A}, ℒ::A) where A
	Objective(obj.f, ℒ, obj.fcall)
end

function _setℒ(obj::Objective{<:Any,<:SVector}, value, index)
	Objective(obj.f, setindex(obj.ℒ, value, index), obj.fcall)
end
function _setℒ(obj::Objective, value, index) # fallback method assumes ℒ is mutable
	Objective(obj.f, setindex!(obj.ℒ, value, index), obj.fcall)
end

function init_Objective(obj, S::Integer)
	S = Int(S)
	(S > 0) || throw(ArgumentError("the number of items S must be a positive integer"))
	if obj isa Objective
		obj = _clearfcall(obj)
		(length(obj.ℒ) == S) || throw(ArgumentError("length of obj.ℒ is not $S"))
	else
		ℒ = (S < _static_threshold()) ? SVector{S, Bool}(ntuple(i->false, S)) : Vector{Bool}(undef, S)
		obj = Objective(obj, ℒ)
	end
	obj
end

# default threshold for determining whether SVector is used for a choice
function _static_threshold()
	256
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

Evalutate the change in `obj` with optional parameter `z` when the `i`th item is included or not. This corresponds to the `D_j` function in earlier implementation.
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

function allundetermined!(itemstates)
	if itemstates isa SVector
		S = length(itemstates)
		itemstates = _fillstate(SVector{S,ItemState}, undetermined)
	else
		fill!(itemstates, undetermined)
	end
	itemstates
end
function allundetermined(obj::Objective)
	allundetermined(obj.ℒ)
end
function allundetermined(ℒ::AbstractVector)
	S = length(ℒ)
	if ℒ isa SVector
		itemstates = _fillstate(SVector{S,ItemState}, undetermined)
	else
		itemstates = fill(undetermined, S)
	end
end

"""
    CDCProblem{M<:CDCPSolver, O<:Objective, A, F<:AbstractFloat}

Results from solving a combinatorial discrete choice problem. When a solution is attained, it can be retrieved from the field `x`.
"""
mutable struct CDCProblem{M<:CDCPSolver, O<:Objective, A, F<:AbstractFloat}
	solver::M
	obj::O
	x::A
	value::F
	state::SolverState
end
