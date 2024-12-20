"""
    Naive{Z} <: CDCPSolver

A type for solving a [`CDCProblem`](@ref) with exhaustion.
"""
struct Naive{Z} <: CDCPSolver
	ids::Vector{Int}
	z::Z
end

# borrowed from Combinatorics.jl just for generating the entire choice space
struct Combinations
	n::Int
	t::Int
end

# @inline is needed to avoid allocation
@inline function Base.iterate(c::Combinations, s)
	if iszero(c.t) # special case to generate 1 result for t==0
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
