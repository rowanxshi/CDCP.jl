"""
    Naive{Z,I} <: CDCPSolver

A type for solving a [`CDCProblem`](@ref) with exhaustion. Decision sets are visited in Gray-code order, so that consecutive sets differ in one item.
"""
struct Naive{Z,I} <: CDCPSolver
	ids::Vector{Int}
	z::Z
	incremental::I
end
Naive(ids, z) = Naive(ids, z, nothing)
