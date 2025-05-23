function branch!(cdcp::CDCProblem{<:Squeezing})
	(; branching) = cdcp.solver
	while !isempty(branching)
		itemstates = pop!(branching)
		itemstates, value, state = squeeze!(cdcp, itemstates)
		if (state == success) && (value > cdcp.value)
			cdcp.x = itemstates
			cdcp.value = value
		end
	end
	cdcp.state = success
	return cdcp
end

function branch(itemstates::AbstractVector{ItemState}, k::Int)
	# All aux are turned into undertermined
	k_in = squeeze(itemstates, included, k)
	# copy is here only in case itemstates is not SVector
	k_out = squeeze(copy(itemstates), excluded, k)
	return k_in, k_out
end
