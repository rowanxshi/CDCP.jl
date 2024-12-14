function squeeze_branch!(p::CDCProblem{<:Squeezing})
	branching, tr = p.solver.branching, p.solver.trace
	while !isempty(branching)
		x = pop!(branching)
		tr === nothing || push!(tr, similar(tr[1], 0))
		x, fx, state = squeeze!(p, x)
		if state == success
			fx > p.fx && (p.x = x; p.fx = fx)
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
