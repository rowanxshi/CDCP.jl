"""
    Objective{F,A}

A wrapped objective function for solving a [`CDCProblem`](@ref).
This facilitates maintaining an internal interface for dealing with objective functions
that is independent from user interface.

Users are *not* required to construct `Objective`
unless there is a need for fine-grained control.

# Constructor
    Objective(f, x, [z=0])

Construct an instance of `Objective` with objective function `f`
and an input vector `x`.
`f` must always accept `x` as the first argument
and may additionally accept an optional argument `z` for parameter
(e.g., productivity).
"""
struct Objective{F,A}
    f::F
    x::A
    fcall::Int
end

Objective(f, x) = Objective(f, x, 0)

"""
    CDCPSolver

Abstract type for all solution algorithms for a [`CDCProblem`](@ref).
"""
abstract type CDCPSolver end

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
