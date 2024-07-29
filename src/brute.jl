struct BruteForce <: CDCPSolver
    ids::Vector{Int}
end

_init(::Type{BruteForce}, obj, args...; kwargs...) = BruteForce([0]), copy(obj.x)

_setchoice(obj::Objective{<:Any,Vector{Bool}}, ids::Vector{Int}) =
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

function _s(obj, C, ids, p)
    next = iterate(C, ids)
    while next !== nothing
        _setchoice(obj, ids)
        val, obj = value(obj, nothing)
        if val > p.fx
            p.fx = val
            p.x .= obj.x
        end
        next = iterate(C, ids)
    end
end

function solve!(p::CDCP{BruteForce}; restart::Bool=false)
    obj, ids = p.obj, p.solver.ids
    restart && (resize!(ids, 1); ids[1] = 0)
    S = length(p.x)
    # The case where no item is chosen
    fill!(obj.x, false)
    val, obj = value(obj, nothing)
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
            val, obj = value(obj, nothing)
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
