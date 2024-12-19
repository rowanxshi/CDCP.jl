function squeeze!(cdcp::CDCProblem{<:SqueezingPolicy})
	solver = cdcp.solver
	squeezing, branching = solver.squeezing, solver.branching
	while !isempty(squeezing)
		k = pop!(squeezing)
		intervalchoice = solver.intervalchoices[k]
		i = next_undetermined(intervalchoice.itemstates)
		if isnothing(i)
			isnothing(findfirst(==(aux), intervalchoice.itemstates)) || push!(branching, k)
		else
			(cdcp.obj.fcall < cdcp.solver.maxfcall) || return maxfcall_reached
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
	intervalchoices = squeeze_include!(cdcp, intervalchoice, i)
	if isone(length(intervalchoices))
		intervalchoice = first(intervalchoices)
		if (intervalchoice.itemstates[i] == included)
			return (intervalchoice, )
		else
			return squeeze_exclude!(cdcp, intervalchoice, i)
		end
	else
		intervalchoice_left, intervalchoice_right = intervalchoices
		intervalchoices = squeeze_exclude!(cdcp, intervalchoice_left, i)
		return (intervalchoices..., intervalchoice_right)
	end
end

function squeeze_include!(cdcp::CDCProblem{<:SqueezingPolicy}, intervalchoice::IntervalChoice, i::Int)
	obj, scdca = cdcp.obj, cdcp.solver.scdca
	lookup = cdcp.solver.lookup_zero_margin
	obj = _setchoice(obj, scdca ? setsup(intervalchoice.itemstates) : setsub(intervalchoice.itemstates))
	key = (i, obj.ℒ)
	z = get!(lookup, key) do
		z, obj = cdcp.solver.zero_margin(obj, i, intervalchoice.zleft, intervalchoice.zright)
		cdcp.solver.zero_margin_call[] += 1
		cdcp.obj = obj
		z
	end
	if z <= intervalchoice.zleft # include the whole interval
		intervalchoice = _squeeze(intervalchoice, included, i)
		return (intervalchoice, )
	elseif z >= intervalchoice.zright # cannot include any part of the interval
		return (intervalchoice, )
	else

		intervalchoice_left = IntervalChoice(intervalchoice.zleft, z, intervalchoice.itemstates)
		intervalchoice_right = IntervalChoice(z, intervalchoice.zright, _squeeze(intervalchoice.itemstates, included, i))
		return (intervalchoice_left, intervalchoice_right)
	end
end

function squeeze_exclude!(cdcp::CDCProblem{<:SqueezingPolicy}, intervalchoice::IntervalChoice, i::Int)
	obj, scdca = cdcp.obj, cdcp.solver.scdca
	lookup = cdcp.solver.lookup_zero_margin
	obj = _setchoice(obj, scdca ? setsub(intervalchoice.itemstates) : setsup(intervalchoice.itemstates))
	key = (i, obj.ℒ)
	z = get!(lookup, key) do
		z, obj = cdcp.solver.zero_margin(obj, i, intervalchoice.zleft, intervalchoice.zright)
		cdcp.solver.zero_margin_call[] += 1
		cdcp.obj = obj
		z
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

