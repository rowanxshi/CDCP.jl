function branching!(cdcp::CDCProblem{<:SqueezingPolicy})
	intervalchoices, branching, matcheds = cdcp.solver.intervalchoices, cdcp.solver.branching, cdcp.solver.matcheds
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
	append!(intervalchoices, matcheds...)
	return inprogress
end

function branching!(cdcp::CDCProblem{<:SqueezingPolicy}, k::Int, itask::Int)
	intervalchoices, matched = cdcp.solver.intervalchoices, cdcp.solver.matcheds[itask]
	intervalchoice = intervalchoices[k]
	sp = cdcp.solver.singlesolvers[itask]
	itemstates_left = _solvesingle!(sp, intervalchoice.lb, intervalchoice.itemstates)
	itemstates_right = _solvesingle!(sp, intervalchoice.ub, intervalchoice.itemstates)
	if itemstates_left == itemstates_right
		intervalchoices[k] = IntervalChoice(intervalchoice.lb, intervalchoice.ub, itemstates_left)
	else
		search!(sp, matched, intervalchoice.lb, itemstates_left, intervalchoice.ub, itemstates_right, intervalchoice.itemstates, cdcp.solver.obj2, cdcp.solver.equal_obj)
		# Overwrite the old interval with aux
		intervalchoices[k] = pop!(matched)
	end
end

function _solvesingle!(cdcp::CDCProblem{<:Squeezing}, z, itemstates0)
	_reset!(cdcp, z, itemstates0)
	solve!(cdcp)
	cdcp.state == success || error("single-agenet solver fails with z = ", cdcp.solver.z)
	return cdcp.x
end

function search!(cdcp::CDCProblem{<:Squeezing}, matched, z_left0, itemstates_left0, z_right0, itemstates_right0, itemstates0, obj2, equal_obj)
	z_left, itemstates_left, z_right, itemstates_right = z_left0, itemstates_left0, z_right0, itemstates_right0
	while true
		# itemstates_left and itemstates_right should always be different here
		obj1 = _setchoice(cdcp.obj, setsub(itemstates_left))
		obj2 = _setchoice(obj2, setsub(itemstates_right))
		z_middle, _ = equal_obj(obj1, obj2, z_left, z_right)
		# ! TODO Rowan proof double-check
		if z_left < z_middle < z_right # Additional cutoff points in between
			itemstates0 = cdcp.solver.scdca ? itemstates0 : setitemstates_scdcb(itemstates0, itemstates_left)
			itemstates_new = _solvesingle!(cdcp, z_middle, itemstates0)
			if itemstates_new == itemstates_left || itemstates_new == itemstates_right
				push!(matched, IntervalChoice(z_left, z_middle, itemstates_left), IntervalChoice(z_middle, z_right, itemstates_right))
				if z_right != z_right0 # Move to the interval on the right
					z_left, itemstates_left = z_right, itemstates_right
					z_right, itemstates_right = z_right0, itemstates_right0
				else
					return
				end
			else # More cutoff points on the left
				z_right, itemstates_right = z_middle, itemstates_new
			end
		else # No split of interval
			if z_middle == z_left
				push!(matched, IntervalChoice(z_left, z_right, itemstates_right))
			elseif z_middle == z_right
				push!(matched, IntervalChoice(z_left, z_right, itemstates_left))
			else
				error("z is not in between")
			end
			if z_right != z_right0 # Move to the interval on the right
				z_left, itemstates_left = z_right, itemstates_right
				z_right, itemstates_right = z_right0, itemstates_right0
			else
				return
			end
		end
	end
end

# With SCD-C from below, choice for an item only switches once
function setitemstates_scdcb(itemstates::SVector{S,ItemState}, itemstates_left::SVector{S,ItemState}) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(ifelse(itemstates_left[$i]==included, included, itemstates[$i])))
		end
		return :(SVector{S,ItemState}($ex))
	else
		return SVector{S,ItemState}(
			ntuple(i->ifelse(itemstates_left[i]==included, included, itemstates[i]), S))
	end
end

function concat!(cdcp::CDCProblem{<:SqueezingPolicy})
	intervalchoices = cdcp.solver.intervalchoices
	sort!(intervalchoices, by=_lb)
	cutoffs = resize!(cdcp.x.cutoffs, 1)
	itemstates_s = resize!(cdcp.x.itemstates_s, 1)
	cutoffs[1] = intervalchoices[1].lb
	itemstates_s[1] = itemstates_last = intervalchoices[1].itemstates
	for intervalchoice in intervalchoices
		# Filter out potential singletons
		if intervalchoice.lb < intervalchoice.ub && intervalchoice.itemstates != itemstates_last
			push!(cutoffs, intervalchoice.lb)
			push!(itemstates_s, intervalchoice.itemstates)
			itemstates_last = intervalchoice.itemstates
		end
	end
	return success
end
