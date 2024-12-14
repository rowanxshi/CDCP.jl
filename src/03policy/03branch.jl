function branching!(cdcp::CDCProblem{<:SqueezingPolicy})
	pool, branching, matcheds = cdcp.solver.pool, cdcp.solver.branching, cdcp.solver.matcheds
	ntasks = length(matcheds)
	for m in matcheds
		empty!(m) # Just to be safe
	end
	if ntasks == 1 # No multithreading
		for k in branching
			branching!(cdcp, k, 1)
		end
	else
		@sync for itask in 1:ntasks
			Threads.@spawn for ik in itask:ntasks:length(branching)
				branching!(cdcp, branching[ik], itask)
			end
		end
	end
	append!(pool, matcheds...)
	return inprogress
end

function branching!(cdcp::CDCProblem{<:SqueezingPolicy}, k::Int, itask::Int)
	pool, matched = cdcp.solver.pool, cdcp.solver.matcheds[itask]
	intervalchoice = pool[k]
	sp = cdcp.solver.singlesolvers[itask]
	xl = _solvesingle!(sp, intervalchoice.lb, intervalchoice.itemstates)
	xr = _solvesingle!(sp, intervalchoice.ub, intervalchoice.itemstates)
	if xl == xr
		pool[k] = IntervalChoice(intervalchoice.lb, intervalchoice.ub, xl)
	else
		search!(sp, matched, intervalchoice.lb, xl, intervalchoice.ub, xr, intervalchoice.itemstates, cdcp.solver.obj2, cdcp.solver.equal_obj)
		# Overwrite the old interval with aux
		pool[k] = pop!(matched)
	end
end

function _solvesingle!(cdcp::CDCProblem{<:Squeezing}, z, x0)
	_reset!(cdcp, z, x0)
	solve!(cdcp)
	cdcp.state == success || error("single-agenet solver fails with z = ", cdcp.solver.z)
	return cdcp.x
end

function search!(cdcp::CDCProblem{<:Squeezing}, matched, zl0, xl0, zr0, xr0, x0, obj2, equal_obj)
	zl, xl, zr, xr = zl0, xl0, zr0, xr0
	while true
		# xl and xr should always be different here
		obj1 = _setchoice(cdcp.obj, setsub(xl))
		obj2 = _setchoice(obj2, setsub(xr))
		z, obj = equal_obj(obj1, obj2, zl, zr)
		# ! TODO Rowan proof double-check
		if zl < z < zr # Additional cutoff points in between
			x0next = cdcp.solver.scdca ? x0 : setx0scdcb(x0, xl)
			xnew = _solvesingle!(cdcp, z, x0next)
			if xnew == xl || xnew == xr
				push!(matched, IntervalChoice(zl, z, xl), IntervalChoice(z, zr, xr))
				if zr != zr0 # Move to the interval on the right
					zl, xl = zr, xr
					zr, xr = zr0, xr0
				else
					return
				end
			else # More cutoff points on the left
				zr, xr = z, xnew
			end
		else # No split of interval
			if z == zl
				push!(matched, IntervalChoice(zl, zr, xr))
			elseif z == zr
				push!(matched, IntervalChoice(zl, zr, xl))
			else
				error("z is not in between")
			end
			if zr != zr0 # Move to the interval on the right
				zl, xl = zr, xr
				zr, xr = zr0, xr0
			else
				return
			end
		end
	end
end

# With SCD-C from below, choice for an item only switches once
function setx0scdcb(x0::SVector{S,ItemState}, xl::SVector{S,ItemState}) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(ifelse(xl[$i]==included, included, x0[$i])))
		end
		return :(SVector{S,ItemState}($ex))
	else
		return SVector{S,ItemState}(
			ntuple(i->ifelse(x[i]==included, included, x0[i]), S))
	end
end

function concat!(cdcp::CDCProblem{<:SqueezingPolicy})
	pool = cdcp.solver.pool
	sort!(pool, by=_lb)
	cutoffs = resize!(cdcp.x.cutoffs, 1)
	xs = resize!(cdcp.x.xs, 1)
	cutoffs[1] = pool[1].lb
	xs[1] = itemstates_last = pool[1].itemstates
	for intervalchoice in pool
		# Filter out potential singletons
		if intervalchoice.lb < intervalchoice.ub && intervalchoice.itemstates != itemstates_last
			push!(cutoffs, intervalchoice.lb)
			push!(xs, intervalchoice.itemstates)
			itemstates_last = intervalchoice.itemstates
		end
	end
	return success
end
