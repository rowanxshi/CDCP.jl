"""
	BruteForce{Z} <: CDCPSolver

A type for solving a [`CDCProblem`](@ref) with a brute-force method.

!!! warning

This method exhaustively iterates through all possible choices,
which can be huge.
It is provided only for illustration.

# Usage
	solve(BruteForce, obj; z=nothing)
	solve!(p::CDCProblem{<:BruteForce}; restart::Bool=false)

Pass the type `BruteForce` as the first argument to `solve`
indicates the use of the brute-force method for the problem.
Users are required to specify the objective function `obj`
that returns the value evaluated at a choice vector `x`
with an optional parameter `z` that is typically a number.
`obj` must have a method of either `obj(x)` or `obj(x, z)`
with `x` being a Boolean choice vector.
`obj` must not restrict the specific type of `x`
but only assume `x` is a vector with element type being `Bool`.
Specifically, `obj` must *not* try to modify the elements in `x` when it is called.
It should only read from `x` with `getindex`.

## Keywords
- `z=nothing`: An optional parameter passed as the second argument to `obj` when evaluating `obj`.
"""
struct BruteForce{Z} <: CDCPSolver
    ids::Vector{Int}
    z::Z
end

_init(::Type{BruteForce}, obj, args...; z=nothing, kwargs...) =
    BruteForce([0], z), copy(obj.x)

_setchoice(obj::Objective{<:Any,<:AbstractVector{Bool}}, ids::Vector{Int}) =
    (fill!(obj.x, false); fill!(view(obj.x, ids), true); obj)

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

function solve!(p::CDCProblem{<:BruteForce}; restart::Bool=true)
    obj, ids, z = p.obj, p.solver.ids, p.solver.z
    restart && (resize!(ids, 1); ids[1] = 0)
    S = length(p.x)
    # The case where no item is chosen
    fill!(obj.x, false)
    val, obj = value(obj, z)
    if val > p.fx
        p.fx = val
        p.x .= obj.x
    end
    # If restart = false, continue from the first combination of the same length as ids
    n1 = length(ids)
    for n in n1:S
        C = Combinations(S, n)
        # Set initial state; see how iterate is defined for Combinations
        resize!(ids, n)
        for i in 1:n
            ids[i] = min(n-1, i)
        end
        # iterate modifies ids in-place and returns (ids, ids)
        next = iterate(C, ids)
        while next !== nothing
            _setchoice(obj, ids)
            val, obj = value(obj, z)
            if val > p.fx
                p.fx = val
                p.x .= obj.x
            end
            next = iterate(C, ids)
        end
    end
    p.obj = obj
	return p
end
