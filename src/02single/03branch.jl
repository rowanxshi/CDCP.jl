function squeeze_branch!(cdcp::CDCProblem{<:Squeezing})
	branching = cdcp.solver.branching
	while !isempty(branching)
		itemstates = pop!(branching)
		itemstates, value, state = squeeze!(cdcp, itemstates)
		if state == success
			value > cdcp.value && begin
				cdcp.x = itemstates
				cdcp.value = value
			end
		elseif state == maxfcall_reached
			push!(branching, itemstates) # Put itemstates back
			return maxfcall_reached
		end
	end
	return success
end

function branch(itemstates::AbstractVector{ItemState}, k::Int)
	# All aux are turned into undertermined
	itemstates_in = _squeeze(itemstates, included, k)
	# copy is here only in case itemstates is not SVector
	itemstates_ex = _squeeze(copy(itemstates), excluded, k)
	return itemstates_in, itemstates_ex
end
