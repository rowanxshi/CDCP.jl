"""
    Naive{Z,I} <: CDCPSolver

A type for solving a [`CDCProblem`](@ref) with exhaustion. Decision sets are visited in Gray-code order, so that consecutive sets differ in one item.

`Naive(; start, flip)` supplies an update for objectives that can be revised when one item changes, in place of calling the objective for every set: `start(z)` returns the value of the empty set and `flip(ℓ, included, z)` the value after item `ℓ` is added or removed.
"""
struct Naive{Z,I} <: CDCPSolver
	ids::Vector{Int}
	z::Z
	incremental::I
end
Naive(ids, z) = Naive(ids, z, nothing)
Naive(; start, flip, z=nothing) = Naive([0], z, (; start, flip))
