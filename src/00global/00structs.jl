"""
    CDCPSolver

Abstract type for all solution algorithms for a [`CDCProblem`](@ref).

See also [`Naive`](@ref), [`Squeezing`](@ref), [`SqueezingPolicy`](@ref).
"""
abstract type CDCPSolver end

"""
    Objective{F,V <: AbstractVector}

A wrapped objective function for solving a [`CDCProblem`](@ref). This facilitates maintaining an internal interface for dealing with objective functions that is independent from user interface.

Users are *not* required to construct `Objective` unless there is a need for fine-grained control.

Fields
===
* `f::F`: objective function
* `ℒ::V`: the decision vector at which to evaluate the objective function
* `fcall::Int = 0`: the number of times the objective has been called, used to track whether the maximum number of objective calls has been reached (see also [`CDCProblem`](@ref))
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

function init_Objective(obj::Objective, S::Integer)
	S = Int(S)
	(S > 0) || throw(ArgumentError("the number of items S must be a positive integer"))
	(length(obj.ℒ) == S) || throw(ArgumentError("length of obj.ℒ is not $S"))
	obj = clearfcall(obj)
end
function init_Objective(obj, S::Integer)
	S = Int(S)
	(S > 0) || throw(ArgumentError("the number of items S must be a positive integer"))
	ℒ = (S < static_threshold()) ? SVector{S, Bool}(ntuple(i->false, S)) : Vector{Bool}(undef, S)
	Objective(obj, ℒ)
end

# default threshold for determining whether SVector is used for a choice
function static_threshold()
	256
end

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
    margin(obj::Objective, ℓ::Int, z)

Evaluate the change in `obj` with optional parameter `z` when the `ℓ`th item is included or not. This corresponds to the `D_j` function in earlier implementation.
"""
function margin(obj::Objective, i::Int, z)
	obj = setindex(obj, true, i)
	value1, obj = obj(z)
	obj = setindex(obj, false, i)
	value0, obj = obj(z)
	return value1, value0, obj
end

"""
    SolverState::Int8

An `Enum` type with values:
* `inprogress`
* `success`
* `maxfcall_reached`

See also [`solve`](@ref), [`CDCProblem`](@ref).
"""
@enum SolverState::Int8 begin
	inprogress
	success
	maxfcall_reached
end

"""
    ItemState::Int8

An `Enum` type describing the state of item `ℓ` with values:
* `undetermined`: not yet checked
* `included`: from squeezing, definitely included
* `excluded`: from squeezing, definitely excluded
* `aux`: checked, but cannot definitively include or exclude based on squeezing

See also [`solve`](@ref), [`CDCProblem`](@ref).
"""
@enum ItemState::Int8 begin
	undetermined
	included
	excluded
	aux
end
function Base.convert(::Type{Bool}, itemstate::ItemState)
	if itemstate == included
		return true
	elseif itemstate == excluded
		return false
	else
		error("Attempted to convert $itemstate to boolean")
	end
end

"""
    CDCProblem{M<:CDCPSolver, O<:Objective, T, F<:AbstractFloat}

Results from solving a combinatorial discrete choice problem with [`solve`](@ref). When a solution is attained, it can be retrieved from the field `x`.

Fields
===
* `solver::`[`CDCPSolver`](@ref)
* `obj::`[`Objective`]: the objective function
* `x`: the solution, once [`solve`](@ref)d
* `value`: the maximum value found by the solver
* `state::`[`SolverState`](@ref)
"""
mutable struct CDCProblem{M<:CDCPSolver, O<:Objective, T, F<:AbstractFloat}
	solver::M
	obj::O
	x::T
	value::F
	state::SolverState
end
