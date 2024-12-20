function squeeze!(cdcp::CDCProblem{<:SqueezingPolicy})
	solver = cdcp.solver
	squeezing, branching = solver.squeezing_indices, solver.branching_indices
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
	obj = setℒ(obj, scdca ? to_sup(intervalchoice.itemstates) : to_sub(intervalchoice.itemstates))
	key = (i, obj.ℒ)
	z, cdcp.obj = cdcp.solver.zero_margin(obj, i, intervalchoice.zleft, intervalchoice.zright)
	if z <= intervalchoice.zleft # include the whole interval
		intervalchoice = squeeze(intervalchoice, included, i)
		return (intervalchoice, )
	elseif z >= intervalchoice.zright # cannot include any part of the interval
		return (intervalchoice, )
	else
		intervalchoice_left = IntervalChoice(intervalchoice.zleft, z, intervalchoice.itemstates)
		intervalchoice_right = IntervalChoice(z, intervalchoice.zright, squeeze(intervalchoice.itemstates, included, i))
		return (intervalchoice_left, intervalchoice_right)
	end
end

function squeeze_exclude!(cdcp::CDCProblem{<:SqueezingPolicy}, intervalchoice::IntervalChoice, i::Int)
	obj, scdca = cdcp.obj, cdcp.solver.scdca
	obj = setℒ(obj, scdca ? to_sub(intervalchoice.itemstates) : to_sup(intervalchoice.itemstates))
	z, cdcp.obj = cdcp.solver.zero_margin(obj, i, intervalchoice.zleft, intervalchoice.zright)
	if z >= intervalchoice.zright # exclude the whole interval
		intervalchoice = squeeze(intervalchoice, excluded, i)
		return (intervalchoice, )
	elseif z <= intervalchoice.zleft # cannot exclude any part of the interval
		intervalchoice = setindex(intervalchoice, aux, i)
		return (intervalchoice, )
	else
		intervalchoice_left = IntervalChoice(intervalchoice.zleft, z, squeeze(intervalchoice.itemstates, excluded, i))
		intervalchoice_right = IntervalChoice(z, intervalchoice.zright, setindex(intervalchoice.itemstates, aux, i))
		return intervalchoice_left, intervalchoice_right
	end
end

function squeeze(intervalchoice::IntervalChoice, s::ItemState, i::Int)
	IntervalChoice(intervalchoice.zleft, intervalchoice.zright, squeeze(intervalchoice.itemstates, s, i))
end

function StaticArrays.setindex(intervalchoice::IntervalChoice, s::ItemState, i::Int)
	IntervalChoice(intervalchoice.zleft, intervalchoice.zright, setindex(intervalchoice.itemstates, s, i))
end

