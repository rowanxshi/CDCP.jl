function squeeze!(cdcp::CDCProblem{<:SqueezingPolicy}) pool, squeezing, branching = cdcp.solver.pool, cdcp.solver.squeezing, cdcp.solver.branching
	while !isempty(squeezing)
		k = pop!(squeezing)
		intervalchoice = pool[k]
		i = findfirst(==(undetermined), intervalchoice.itemstates)
		if i === nothing
			findfirst(==(aux), intervalchoice.itemstates) === nothing || push!(branching, k)
		else
			cdcp.obj.fcall < cdcp.maxfcall || return maxfcall_reached
			intervalchoices = squeeze!(cdcp, intervalchoice, i)
			pool[k] = intervalchoices[1]
			push!(squeezing, k)
			for j in 2:length(intervalchoices)
				push!(pool, intervalchoices[j])
				push!(squeezing, length(pool))
			end
		end
	end
	return inprogress
end

function squeeze!(cdcp::CDCProblem{<:SqueezingPolicy}, intervalchoice::IntervalChoice, i::Int)
	obj, scdca = cdcp.obj, cdcp.solver.scdca
	lookup = cdcp.solver.lookup_zero_margin
	obj = _setchoice(obj, scdca ? setsup(intervalchoice.itemstates) : setsub(intervalchoice.itemstates))
	key = (i, obj.x)
	z = get(lookup, key, nothing)
	if z === nothing
		z, obj = cdcp.solver.zero_margin(obj, i, intervalchoice.lb, intervalchoice.ub)
		lookup[key] = z
		cdcp.solver.zero_margin_call[] += 1
		cdcp.obj = obj
	end
	if z <= intervalchoice.lb # Include the whole interval
		intervalchoice = _squeeze(intervalchoice, included, i)
		return (intervalchoice, )
	elseif z >= intervalchoice.ub # Cannot include any part of the interval
		return squeeze_exclude!(cdcp, intervalchoice, i)
	else
		intervalchoice_left = IntervalChoice(intervalchoice.lb, z, intervalchoice.itemstates)
		intervalchoices = squeeze_exclude!(cdcp, intervalchoice_left, i)
		intervalchoice_right = IntervalChoice(z, intervalchoice.ub, intervalchoice.itemstates)
		return (intervalchoices..., intervalchoice_right)
	end
end

function squeeze_exclude!(cdcp::CDCProblem{<:SqueezingPolicy}, intervalchoice::IntervalChoice, i::Int)
	obj, scdca = cdcp.obj, cdcp.solver.scdca
	lookup = cdcp.solver.lookup_zero_margin
	obj = _setchoice(obj, scdca ? setsub(intervalchoice.itemstates) : setsup(intervalchoice.itemstates))
	key = (i, obj.x)
	z = get(lookup, key, nothing)
	if z === nothing
		z, obj = cdcp.solver.zero_margin(obj, i, intervalchoice.lb, intervalchoice.ub)
		lookup[key] = z
		cdcp.solver.zero_margin_call[] += 1
		cdcp.obj = obj
	end
	if z >= intervalchoice.ub # Exclude the whole interval
		intervalchoice = _squeeze(intervalchoice, excluded, i)
		return (intervalchoice, )
	elseif z <= intervalchoice.lb # Cannot exclude any part of the interval
		intervalchoice = _setitemstate(intervalchoice, aux, i) # Don't use _squeeze
		return (intervalchoice, )
	else
		intervalchoice_left = IntervalChoice(intervalchoice.lb, z, _squeeze(intervalchoice.itemstates, excluded, i))
		intervalchoice_right = IntervalChoice(z, intervalchoice.ub, _setitemstate(intervalchoice.itemstates, aux, i))
		return intervalchoice_left, intervalchoice_right
	end
end

function _squeeze(intervalchoice::IntervalChoice, s::ItemState, i::Int)
	IntervalChoice(intervalchoice.lb, intervalchoice.ub, _squeeze(intervalchoice.itemstates, s, i))
end

function _setitemstate(intervalchoice::IntervalChoice, s::ItemState, i::Int)
	IntervalChoice(intervalchoice.lb, intervalchoice.ub, _setitemstate(intervalchoice.itemstates, s, i))
end

