"""
    BruteForce{Z} <: CDCPSolver

A type for solving a [`CDCProblem`](@ref) with a brute-force method.

!!! warning

This method exhaustively iterates through all possible choices,
which can be huge.
It is provided only for illustration.

# Usage
	solve(BruteForce, obj; z=nothing)
	solve!(cdcp::CDCProblem{<:BruteForce}; restart::Bool=false)

Pass the type `BruteForce` as the first argument to `solve` indicates the use of the brute-force method for the problem. Users are required to specify the objective function `obj` that returns the value evaluated at a choice vector `ℒ` with an optional parameter `z` that is typically a number. `obj` must have a method of either `obj(ℒ)` or `obj(ℒ, z)` with `ℒ` being a Boolean choice vector. `obj` must not restrict the specific type of `ℒ` but only assume `ℒ` is a vector with element type being `Bool`. Specifically, `obj` must *not* try to modify the elements in `ℒ` when it is called. It should only read from `ℒ` with `getindex`. 
## Keywords
- `z=nothing`: An optional parameter passed as the second argument to `obj` when evaluating `obj`.
"""
struct BruteForce{Z} <: CDCPSolver
	ids::Vector{Int}
	z::Z
end

# Borrowed from Combinatorics.jl just for generating the entire choice space
struct Combinations
	n::Int
	t::Int
end

# @inline is needed for avoid the allocation
@inline function Base.iterate(c::Combinations, s)
	if c.t == 0 # special case to generate 1 result for t==0
		isempty(s) && return (s, [1])
		return
	end
	for i in c.t:-1:1
		s[i] += 1
		if s[i] > (c.n - (c.t - i))
			continue
		end
		for j in i+1:c.t
			s[j] = s[j-1] + 1
		end
		break
	end
	s[1] > c.n - c.t + 1 && return
	(s, s)
end
