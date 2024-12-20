"""
    Squeezing{V <: AbstractVector,Z} <: CDCPSolver

A type for solving a [`CDCProblem`](@ref) with single-agent squeezing.

See also [`solve`](@ref).

Fields
===
* `scdca::Bool`: does the problem obey single crossing differences in choices from above? (if `false`, then from below is assumed)
* `branching::Vector{V}`: container holding branching outcomes
* `z`: the parameter (e.g. productivity) of the problem
"""
struct Squeezing{V <: AbstractVector,Z} <: CDCPSolver
	scdca::Bool
	branching::Vector{V}
	z::Z
end
