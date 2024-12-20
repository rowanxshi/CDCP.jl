"""
    Squeezing{V <: AbstractVector,Z} <: CDCPSolver

A type for solving a [`CDCProblem`](@ref) with a single-agent method as in Arkolakis, Eckert and Shi (2023).

# Usage
	solve(Squeezing, obj, scdca::Bool; kwargs...)
	solve!(cdcp::CDCProblem{<:Squeezing}; restart::Bool=false)

Pass the type `Squeezing` as the first argument to `solve` indicates the use of the single-agent method for the problem. Users are required to specify the objective function `obj` that returns the value evaluated at a choice vector `ℒ` with an optional parameter `z` that is typically a number. `obj` must have a method of either `obj(ℒ)` or `obj(ℒ, z)` with `ℒ` being a Boolean choice vector. `obj` must not restrict the specific type of `ℒ` but only assume `ℒ` is a vector with element type being `Bool`. Specifically, `obj` must *not* try to modify the elements in `ℒ` when it is called. It should only read from `ℒ` with `getindex`. The problem should satisfy SCD-C from above if `scdca` is `true` and SCD-C from below if `scdca` is `false`.
"""
struct Squeezing{V <: AbstractVector,Z} <: CDCPSolver
	scdca::Bool
	branching::Vector{V}
	z::Z
end
