function squeeze_branch!(cdcp::CDCProblem{<:Squeezing})
	branching = cdcp.solver.branching
	while !isempty(branching)
		x = pop!(branching)
		x, fx, state = squeeze!(cdcp, x)
		if state == success
			fx > cdcp.fx && (cdcp.x = x; cdcp.fx = fx)
		elseif state == maxfcall_reached
			push!(branching, x) # Put x back
			return maxfcall_reached
		end
	end
	return success
end

function branch(x::AbstractVector{ItemState}, k::Int)
	# All aux are turned into undertermined
	xin = _squeeze(x, included, k)
	# copy is here only in case x is not SVector
	xex = _squeeze(copy(x), excluded, k)
	return xin, xex
end
