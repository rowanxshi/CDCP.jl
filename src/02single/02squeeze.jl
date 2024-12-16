function squeeze!(cdcp::CDCProblem{<:Squeezing}, itemstates::AbstractVector{ItemState})
	S = length(itemstates)
	value = -Inf
	lastaux = nothing
	i = findfirst(==(undetermined), itemstates)
	if i === nothing
		lastaux = findfirst(==(aux), itemstates) # May find aux from initial value
		if lastaux === nothing # Last value set by branching
			value, obj = cdcp.obj(cdcp.solver.z)
			cdcp.obj = obj
		end
	end
	while i !== nothing && cdcp.obj.fcall < cdcp.maxfcall
		itemstates, value, itemstate = squeeze!(cdcp, itemstates, i)
		if itemstate == aux
			lastaux = lastaux === nothing ? i : min(lastaux, i)
		else
			lastaux = nothing
		end
		#lastaux = ifelse(itemstate == aux, i, nothing)
		if i < S
			i = findnext(==(undetermined), itemstates, i+1)
			i === nothing && (i = findfirst(==(undetermined), itemstates))
		else
			i = findfirst(==(undetermined), itemstates)
		end
	end
	# Whenever squeeze! makes progress, any aux is changed to undetermined
	# If aux is in itemstates, it can only be that lastaux !== nothing
	if lastaux === nothing
		state = cdcp.obj.fcall < cdcp.maxfcall ? success : maxfcall_reached
	else
		push!(cdcp.solver.branching, branch(itemstates, lastaux)...)
		state = inprogress
	end
	return itemstates, value, state
end

function squeeze!(cdcp::CDCProblem{<:Squeezing}, itemstates::AbstractVector{ItemState}, i::Int)
	obj, scdca, z = cdcp.obj, cdcp.solver.scdca, cdcp.solver.z
	# For excluding: look at the best case
	if scdca
		obj = _setchoice(obj, setsub(itemstates))
		f1, f0, obj = margin(obj, i, z)
		exclude = f1 - f0 <= 0
	else
		obj = _setchoice(obj, setsup(itemstates))
		f1, f0, obj = margin(obj, i, z)
		exclude = f1 - f0 < 0
	end
	if exclude
		itemstates_new = _squeeze(itemstates, excluded, i)
		cdcp.obj = obj
		return itemstates_new, f0, excluded
	end

	# For including: look at the worst case
	if scdca
		obj = _setchoice(obj, setsup(itemstates))
		f1, f0, obj = margin(obj, i, z)
		include = f1 - f0 > 0
	else
		obj = _setchoice(obj, setsub(itemstates))
		f1, f0, obj = margin(obj, i, z)
		include = f1 - f0 >= 0
	end
	if include
		itemstates_new = _squeeze(itemstates, included, i)
		cdcp.obj = obj
		return itemstates_new, f1, included
	end

	itemstates_new = _setitemstate(itemstates, aux, i)
	value = convert(typeof(cdcp.value), -Inf)
	cdcp.obj = obj
	return itemstates_new, value, aux
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
		@inbounds itemstatesi = itemstates[i]
		if itemstatesi == aux
			return undetermined
		else
			return itemstatesi
		end
	end
end

