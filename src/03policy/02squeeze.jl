function squeeze!(cdcp::CDCProblem{<:SqueezingPolicy})
	solver, squeezing, branching = cdcp.solver, cdcp.solver.squeezing, cdcp.solver.branching
	while !isempty(squeezing)
		k = pop!(squeezing)
		intervalchoice = solver.intervalchoices[k]
		i = next_undetermined(intervalchoice.itemstates)
		if i === nothing
			findfirst(==(aux), intervalchoice.itemstates) === nothing || push!(branching, k)
		else
			cdcp.obj.fcall < cdcp.solver.maxfcall || return maxfcall_reached
			intervalchoices = squeeze!(cdcp, intervalchoice, i)
			solver.intervalchoices[k] = intervalchoices[1]
			push!(squeezing, k)
			for j in 2:length(intervalchoices)
				push!(solver.intervalchoices, intervalchoices[j])
				push!(squeezing, length(solver.intervalchoices))
			end
		end
	end
	return inprogress
end

function squeeze!(cdcp::CDCProblem{<:SqueezingPolicy}, intervalchoice::IntervalChoice, i::Int)
	obj, scdca = cdcp.obj, cdcp.solver.scdca
	lookup = cdcp.solver.lookup_zero_margin
	obj = _setchoice(obj, scdca ? setsup(intervalchoice.itemstates) : setsub(intervalchoice.itemstates))
	key = (i, obj.ℒ)
	z = get(lookup, key, nothing)
	if isnothing(z)
		z, obj = cdcp.solver.zero_margin(obj, i, intervalchoice.zleft, intervalchoice.zright)
		lookup[key] = z
		cdcp.solver.zero_margin_call[] += 1
		cdcp.obj = obj
	end
	if z <= intervalchoice.zleft # Include the whole interval
		intervalchoice = _squeeze(intervalchoice, included, i)
		return (intervalchoice, )
	elseif z >= intervalchoice.zright # Cannot include any part of the interval
		return squeeze_exclude!(cdcp, intervalchoice, i)
	else
		intervalchoice_left = IntervalChoice(intervalchoice.zleft, z, intervalchoice.itemstates)
		intervalchoices = squeeze_exclude!(cdcp, intervalchoice_left, i)
		intervalchoice_right = IntervalChoice(z, intervalchoice.zright, _squeeze(intervalchoice.itemstates, included, i))
		return (intervalchoices..., intervalchoice_right)
	end
end

function squeeze_exclude!(cdcp::CDCProblem{<:SqueezingPolicy}, intervalchoice::IntervalChoice, i::Int)
	obj, scdca = cdcp.obj, cdcp.solver.scdca
	lookup = cdcp.solver.lookup_zero_margin
	obj = _setchoice(obj, scdca ? setsub(intervalchoice.itemstates) : setsup(intervalchoice.itemstates))
	key = (i, obj.ℒ)
	z = get(lookup, key, nothing)
	if z === nothing
		z, obj = cdcp.solver.zero_margin(obj, i, intervalchoice.zleft, intervalchoice.zright)
		lookup[key] = z
		cdcp.solver.zero_margin_call[] += 1
		cdcp.obj = obj
	end
	if z >= intervalchoice.zright # Exclude the whole interval
		intervalchoice = _squeeze(intervalchoice, excluded, i)
		return (intervalchoice, )
	elseif z <= intervalchoice.zleft # Cannot exclude any part of the interval
		intervalchoice = _setitemstate(intervalchoice, aux, i) # Don't use _squeeze
		return (intervalchoice, )
	else
		intervalchoice_left = IntervalChoice(intervalchoice.zleft, z, _squeeze(intervalchoice.itemstates, excluded, i))
		intervalchoice_right = IntervalChoice(z, intervalchoice.zright, setindex(intervalchoice.itemstates, aux, i))
		return intervalchoice_left, intervalchoice_right
	end
end

function _squeeze(intervalchoice::IntervalChoice, s::ItemState, i::Int)
	IntervalChoice(intervalchoice.zleft, intervalchoice.zright, _squeeze(intervalchoice.itemstates, s, i))
end

function _setitemstate(intervalchoice::IntervalChoice, s::ItemState, i::Int)
	IntervalChoice(intervalchoice.zleft, intervalchoice.zright, setindex(intervalchoice.itemstates, s, i))
end

