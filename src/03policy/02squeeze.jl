function squeeze!(cdcp::CDCProblem{<:SqueezingPolicy}) pool, squeezing, branching = cdcp.solver.pool, cdcp.solver.squeezing, cdcp.solver.branching
	while !isempty(squeezing)
		k = pop!(squeezing)
		intervalchoice = pool[k]
		i = findfirst(==(undetermined), intervalchoice.x)
		if i === nothing
			findfirst(==(aux), intervalchoice.x) === nothing || push!(branching, k)
		else
			cdcp.obj.fcall < cdcp.maxfcall || return maxfcall_reached
			xs = squeeze!(cdcp, intervalchoice, i)
			pool[k] = xs[1]
			push!(squeezing, k)
			for j in 2:length(xs)
				push!(pool, xs[j])
				push!(squeezing, length(pool))
			end
		end
	end
	return inprogress
end

function squeeze!(cdcp::CDCProblem{<:SqueezingPolicy}, intervalchoice::IntervalChoice, i::Int)
	obj, scdca, tr = cdcp.obj, cdcp.solver.scdca, cdcp.solver.trace
	lookup = cdcp.solver.lookup_zero_margin
	obj = _setchoice(obj, scdca ? setsup(intervalchoice.x) : setsub(intervalchoice.x))
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
		tr === nothing || push!(tr[end], SqueezingPolicyTrace(i, z, intervalchoice.lb, intervalchoice.ub, included))
		return (intervalchoice,)
	elseif z >= intervalchoice.ub # Cannot include any part of the interval
		tr === nothing ||
			push!(tr[end], SqueezingPolicyTrace(i, z, intervalchoice.lb, intervalchoice.ub, undetermined))
		return squeeze_exclude!(cdcp, intervalchoice, i)
	else
		tr === nothing ||
			push!(tr[end], SqueezingPolicyTrace(i, z, intervalchoice.lb, intervalchoice.ub, undetermined))
		x1 = IntervalChoice(intervalchoice.lb, z, intervalchoice.x)
		xs = squeeze_exclude!(cdcp, x1, i)
		x2 = IntervalChoice(z, intervalchoice.ub, intervalchoice.x)
		return xs..., x2
	end
end

function squeeze_exclude!(cdcp::CDCProblem{<:SqueezingPolicy}, intervalchoice::IntervalChoice, i::Int)
	obj, scdca, tr = cdcp.obj, cdcp.solver.scdca, cdcp.solver.trace
	lookup = cdcp.solver.lookup_zero_margin
	obj = _setchoice(obj, scdca ? setsub(intervalchoice.x) : setsup(intervalchoice.x))
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
		tr === nothing || push!(tr[end], SqueezingPolicyTrace(i, z, intervalchoice.lb, intervalchoice.ub, excluded))
		return (intervalchoice,)
	elseif z <= intervalchoice.lb # Cannot exclude any part of the interval
		intervalchoice = _setitemstate(intervalchoice, aux, i) # Don't use _squeeze
		tr === nothing || push!(tr[end], SqueezingPolicyTrace(i, z, intervalchoice.lb, intervalchoice.ub, aux))
		return (intervalchoice,)
	else
		tr === nothing || push!(tr[end], SqueezingPolicyTrace(i, z, intervalchoice.lb, intervalchoice.ub, excluded))
		x1 = IntervalChoice(intervalchoice.lb, z, _squeeze(intervalchoice.x, excluded, i))
		x2 = IntervalChoice(z, intervalchoice.ub, _setitemstate(intervalchoice.x, aux, i))
		return x1, x2
	end
end

function _squeeze(intervalchoice::IntervalChoice, s::ItemState, i::Int)
	IntervalChoice(intervalchoice.lb, intervalchoice.ub, _squeeze(intervalchoice.x, s, i))
end

function _setitemstate(intervalchoice::IntervalChoice, s::ItemState, i::Int)
	IntervalChoice(intervalchoice.lb, intervalchoice.ub, _setitemstate(intervalchoice.x, s, i))
end

