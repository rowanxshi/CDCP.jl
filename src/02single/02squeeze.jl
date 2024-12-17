function squeeze!(cdcp::CDCProblem{<:Squeezing}, itemstates::AbstractVector{ItemState})
	S = length(itemstates)
	value = -Inf
	lastaux = nothing
	i = findfirst(==(undetermined), itemstates)

	if isnothing(i)
		lastaux = findfirst(==(aux), itemstates) # may find aux from initial value
		if isnothing(lastaux) # last value set by branching
			value, obj = cdcp.obj(cdcp.solver.z)
			cdcp.obj = obj
		end
	end

	while !isnothing(i) && (cdcp.obj.fcall < cdcp.maxfcall)
		itemstates, value, itemstate = squeeze!(cdcp, itemstates, i)
		if itemstate == aux
			lastaux = isnothing(lastaux) ? i : min(lastaux, i)
		else
			lastaux = nothing
		end
		if i < S
			i = findnext(==(undetermined), itemstates, i+1)
			isnothing(i) && (i = findfirst(==(undetermined), itemstates))
		else
			i = findfirst(==(undetermined), itemstates)
		end
	end

### TODO I FEEL LIKE THIS PART GOES INTO BRANCHING
	# whenever squeeze! makes progress, any aux is changed to undetermined
	# if aux is in itemstates, it can only be that lastaux !== nothing
	if isnothing(lastaux)
		state = cdcp.obj.fcall < cdcp.maxfcall ? success : maxfcall_reached
	else
		push!(cdcp.solver.branching, branch(itemstates, lastaux)...)
		state = inprogress
	end
	return itemstates, value, state
end

function squeeze!(cdcp::CDCProblem{<:Squeezing}, itemstates::AbstractVector{ItemState}, i::Int)
	# check exclude
	exclude, obj, value0 = isexcluded(cdcp, itemstates, i)	
	if exclude
		itemstates_new = _squeeze(itemstates, excluded, i)
		cdcp.obj = obj
		return itemstates_new, value0, excluded
	end
	# check include
	include, obj, value1 = isincluded(cdcp, itemstates, i)
	if include
		itemstates_new = _squeeze(itemstates, included, i)
		cdcp.obj = obj
		return itemstates_new, value1, included
	end
	# no progress
	value = convert(typeof(cdcp.value), -Inf)
	begin
		itemstates_new = _setitemstate(itemstates, aux, i)
		cdcp.obj = obj
		return itemstates_new, value, aux
	end
end

function isexcluded(cdcp::CDCProblem{<:Squeezing}, itemstates::AbstractVector{ItemState}, i::Int)
	obj, z = cdcp.obj, cdcp.solver.z
	if cdcp.solver.scdca
		sub = setsub(itemstates)
		obj = _setchoice(obj, sub)
		value1, value0, obj = margin(obj, i, z)
		exclude = (value1 - value0 <= 0)
	else
		sup = setsup(itemstates)
		obj = _setchoice(obj, sup)
		value1, value0, obj = margin(obj, i, z)
		exclude = (value1 - value0 < 0)
	end
	exclude, obj, value0
end

function isincluded(cdcp::CDCProblem{<:Squeezing}, itemstates::AbstractVector{ItemState}, i::Int)
	obj, z = cdcp.obj, cdcp.solver.z
	if cdcp.solver.scdca
		sup = setsup(itemstates)
		obj = _setchoice(obj, sup)
		value1, value0, obj = margin(obj, i, z)
		include = (value1 - value0 > 0)
	else
		sub = setsub(itemstates)
		obj = _setchoice(obj, sub)
		value1, value0, obj = margin(obj, i, z)
		include = (value1 - value0 >= 0)
	end
	return include, obj, value1
end

function _squeeze(itemstates::SVector{S,ItemState}, s::ItemState, k::Int) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(_squeeze(itemstates, s, k, $i)))
		end
		return :(SVector{S,ItemState}($ex))
	else
		return SVector{S,ItemState}(ntuple(i->_squeeze(itemstates, s, k, i), S))
	end
end

function _squeeze(itemstates::SVector{S,ItemState}, s::ItemState, k::Int, i::Int) where S
	if i == k
		return s
	else
		@inbounds itemstates_i = itemstates[i]
		if itemstates_i == aux
			return undetermined
		else
			return itemstates_i
		end
	end
end

