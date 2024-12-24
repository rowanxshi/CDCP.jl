# TODO maybe rename this file to search or something better than branch
# TODO remove matched and just put stuff onto the end of intervalchoices 
function search!(cdcp::CDCProblem{<:SqueezingPolicy})
	search!(cdcp.solver)
end
function search!(solver::SqueezingPolicy)
	(; branching_indices) = solver
	while !isempty(branching_indices)
		k = pop!(branching_indices)
		search!(solver, k)
	end
	return inprogress
end
function search!(solver::SqueezingPolicy, k::Int)
	(; intervalchoices) = solver
	intervalchoice = intervalchoices[k]
	(itemstates_left, itemstates_right) = solve_endpoints!(solver, intervalchoice)
	if itemstates_left == itemstates_right
		intervalchoices[k] = IntervalChoice(intervalchoice.zleft, intervalchoice.zright, itemstates_left)
	else
		search!(solver, intervalchoice.zleft, itemstates_left, intervalchoice.zright, itemstates_right, intervalchoice.itemstates)
		intervalchoices[k] = pop!(intervalchoices) # overwrite the old interval with aux
	end
end

function solve_endpoints!(solver::SqueezingPolicy, intervalchoice::IntervalChoice)
	(; singlecdcp) = solver
	itemstates_left = solve!(singlecdcp, intervalchoice.zleft, intervalchoice.itemstates).x
	itemstates_right = solve!(singlecdcp, intervalchoice.zright, intervalchoice.itemstates).x
	(itemstates_left, itemstates_right)
end
function solve!(cdcp::CDCProblem{<:Squeezing}, z, itemstates::AbstractVector{ItemState})
	reinit!(cdcp; z, fcall=false)
	cdcp.x = itemstates
	solve!(cdcp)
	(cdcp.state == success) || error("single-agent solver fails with z = ", cdcp.solver.z)
	return cdcp
end

function search!(solver::SqueezingPolicy, zleft0, itemstates_left0, zright0, itemstates_right0, itemstates0)
	(; intervalchoices, singlecdcp, obj2, equal_obj) = solver
	zleft, itemstates_left, zright, itemstates_right = zleft0, itemstates_left0, zright0, itemstates_right0
	while true
		obj1 = setℒ(singlecdcp.obj, to_sub(itemstates_left))
		obj2 = setℒ(obj2, to_sub(itemstates_right))
		zmiddle, _ = equal_obj(obj1, obj2, zleft, zright)
		(zleft <= zmiddle <= zright) || error("equal_obj returns z that is not between zleft and zright")
		if (zmiddle == zleft) || (zmiddle == zright)
			if zmiddle == zleft
				push!(intervalchoices, IntervalChoice(zleft, zright, itemstates_right))
			elseif zmiddle == zright
				push!(intervalchoices, IntervalChoice(zleft, zright, itemstates_left))
			end
			if zright != zright0 # move to the interval on the right
				zleft, itemstates_left = zright, itemstates_right
				zright, itemstates_right = zright0, itemstates_right0
			else
				return intervalchoices
			end
		else
			itemstates0 = singlecdcp.solver.scdca ? itemstates0 : setitemstates_scdcb(itemstates0, itemstates_left)
			itemstates_middle = solve!(singlecdcp, zmiddle, itemstates0).x
			if (itemstates_middle == itemstates_left) || (itemstates_middle == itemstates_right)
				push!(intervalchoices, IntervalChoice(zleft, zmiddle, itemstates_left), IntervalChoice(zmiddle, zright, itemstates_right))
				if zright != zright0 # move to the interval on the right
					zleft, itemstates_left = zright, itemstates_right
					zright, itemstates_right = zright0, itemstates_right0
				else
					return intervalchoices
				end
			else
				zright, itemstates_right = zmiddle, itemstates_middle
			end
		end
	end
end

# with scd-c from below, choice for an item only switches once
function setitemstates_scdcb(itemstates::SVector{S,ItemState}, itemstates_left::SVector{S,ItemState}) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(ifelse(itemstates_left[$i]==included, included, itemstates[$i])))
		end
		return :(SVector{S,ItemState}($ex))
	else
		return SVector{S,ItemState}(ntuple(i->ifelse(itemstates_left[i]==included, included, itemstates[i]), S))
	end
end

function concat!(cdcp::CDCProblem{<:SqueezingPolicy})
	(; intervalchoices) = cdcp.solver
	sort!(intervalchoices, by=(intervalchoice->getproperty(intervalchoice, :zleft)))
	cutoffs = resize!(cdcp.x.cutoffs, 1)
	itemstates_s = resize!(cdcp.x.itemstates_s, 1)
	cutoffs[1] = intervalchoices[1].zleft
	itemstates_s[1] = intervalchoices[1].itemstates
	itemstates_last = itemstates_s[1]
	for intervalchoice in intervalchoices
		if intervalchoice.itemstates != itemstates_last
			push!(cutoffs, intervalchoice.zleft)
			push!(itemstates_s, intervalchoice.itemstates)
			itemstates_last = intervalchoice.itemstates
		end
	end
	return success
end
