# TODO maybe rename this file to search or something better than branch
function search!(cdcp::CDCProblem{<:SqueezingPolicy})
	(; intervalchoices, branching_indices, matched) = cdcp.solver
	empty!(matched)
	for k in branching_indices
		search!(cdcp, k)
	end
	append!(intervalchoices, matched)
	return inprogress
end

function search!(cdcp::CDCProblem{<:SqueezingPolicy}, k::Int)
	matched = cdcp.solver.matched # vector of intervalchoices; TODO rename matcheds to be more descriptive, then remove this comment; but currently empty
	singlecdcp = cdcp.solver.singlecdcp

	(; intervalchoices)  = cdcp.solver
	intervalchoice = intervalchoices[k]
	itemstates_left = solve!(singlecdcp, intervalchoice.zleft, intervalchoice.itemstates).x
	itemstates_right = solve!(singlecdcp, intervalchoice.zright, intervalchoice.itemstates).x
	if itemstates_left == itemstates_right
		intervalchoices[k] = IntervalChoice(intervalchoice.zleft, intervalchoice.zright, itemstates_left)
	else
		search!(singlecdcp, matched, intervalchoice.zleft, itemstates_left, intervalchoice.zright, itemstates_right, intervalchoice.itemstates, cdcp.solver.obj2, cdcp.solver.equal_obj)
		intervalchoices[k] = pop!(matched) # overwrite the old interval with aux
	end
end

function solve!(cdcp::CDCProblem{<:Squeezing}, z, itemstates::AbstractVector{ItemState})
	reinit!(cdcp; z, fcall=false)
	cdcp.x = itemstates
	solve!(cdcp)
	(cdcp.state == success) || error("single-agent solver fails with z = ", cdcp.solver.z)
	return cdcp
end

function search!(cdcp::CDCProblem{<:Squeezing}, matched, zleft0, itemstates_left0, zright0, itemstates_right0, itemstates0, obj2, equal_obj)
	zleft, itemstates_left, zright, itemstates_right = zleft0, itemstates_left0, zright0, itemstates_right0
	while true
		# itemstates_left and itemstates_right should always be different here
		obj1 = setℒ(cdcp.obj, to_sub(itemstates_left))
		obj2 = setℒ(obj2, to_sub(itemstates_right))
		zmiddle, _ = equal_obj(obj1, obj2, zleft, zright)
		# ! TODO Rowan proof double-check
		if zleft < zmiddle < zright # Additional cutoff points in between
			itemstates0 = cdcp.solver.scdca ? itemstates0 : setitemstates_scdcb(itemstates0, itemstates_left)
			itemstates_new = solve!(cdcp, zmiddle, itemstates0).x
			if itemstates_new == itemstates_left || itemstates_new == itemstates_right
				push!(matched, IntervalChoice(zleft, zmiddle, itemstates_left), IntervalChoice(zmiddle, zright, itemstates_right))
				if zright != zright0 # Move to the interval on the right
					zleft, itemstates_left = zright, itemstates_right
					zright, itemstates_right = zright0, itemstates_right0
				else
					return
				end
			else # More cutoff points on the left
				zright, itemstates_right = zmiddle, itemstates_new
			end
		else # No split of interval
			if zmiddle == zleft
				push!(matched, IntervalChoice(zleft, zright, itemstates_right))
			elseif zmiddle == zright
				push!(matched, IntervalChoice(zleft, zright, itemstates_left))
			else
				error("z is not in between")
			end
			if zright != zright0 # Move to the interval on the right
				zleft, itemstates_left = zright, itemstates_right
				zright, itemstates_right = zright0, itemstates_right0
			else
				return
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
	intervalchoices = cdcp.solver.intervalchoices
	sort!(intervalchoices, by=(intervalchoice->getproperty(intervalchoice, :zleft)))
	cutoffs = resize!(cdcp.x.cutoffs, 1)
	itemstates_s = resize!(cdcp.x.itemstates_s, 1)
	cutoffs[1] = intervalchoices[1].zleft
	itemstates_s[1] = itemstates_last = intervalchoices[1].itemstates
	for intervalchoice in intervalchoices
		# filter out potential singletons
		if intervalchoice.zleft < intervalchoice.zright && intervalchoice.itemstates != itemstates_last
			push!(cutoffs, intervalchoice.zleft)
			push!(itemstates_s, intervalchoice.itemstates)
			itemstates_last = intervalchoice.itemstates
		end
	end
	return success
end
