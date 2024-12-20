function solve!(cdcp::CDCProblem{<:Naive}; restart::Bool=true)
	obj, ids, z = cdcp.obj, cdcp.solver.ids, cdcp.solver.z
	if restart
		resize!(ids, 1)
		ids[1] = 0
	end
	S = length(cdcp.x)
	# the case where no item is chosen
	fill!(obj.ℒ, false)
	value, obj = obj(z)
	if value > cdcp.value
		cdcp.value = value
		cdcp.x .= obj.ℒ
	end
	# if restart = false, continue from the first combination of the same length as ids
	n1 = length(ids)
	for n in n1:S
		C = Combinations(S, n)
		# Set initial state; see how iterate is defined for Combinations
		resize!(ids, n)
		for i in 1:n
			ids[i] = min(n - 1, i)
		end
		# iterate modifies ids in-place and returns (ids, ids)
		next = iterate(C, ids)
		while !isnothing(next)
			setℒ(obj, ids)
			value, obj = obj(z)
			if value > cdcp.value
				cdcp.value = value
				cdcp.x .= obj.ℒ
			end
			next = iterate(C, ids)
		end
	end
	cdcp.obj = obj
	return cdcp
end

