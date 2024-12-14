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

