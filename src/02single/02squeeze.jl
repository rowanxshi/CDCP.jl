function squeeze!(cdcp::CDCProblem{<:Squeezing}, itemstates::AbstractVector{ItemState})
	value = -Inf
	lastaux = nothing
	i = next_undetermined(itemstates)
	if isnothing(i)
		lastaux = findfirst(==(aux), itemstates)
		if isnothing(lastaux)
			value, obj = cdcp.obj(cdcp.solver.z)
			cdcp.obj = obj
		end
	end
	while !isnothing(i)
		itemstates, value, itemstate = squeeze(cdcp, itemstates, i)
		if itemstate == aux
			lastaux = isnothing(lastaux) ? i : min(lastaux, i)
		else
			lastaux = nothing
		end
		i = next_undetermined(itemstates, i)
	end
	if isnothing(lastaux)
		cdcpstate = success
	else
		cdcpstate = inprogress
		push!(cdcp.solver.branching, branch(itemstates, lastaux)...)
	end
	return itemstates, value, cdcpstate
end

function squeeze(cdcp::CDCProblem{<:Squeezing}, itemstates::AbstractVector{ItemState}, i::Int)
	# check exclude
	exclude, obj, value0 = isexcluded(cdcp, itemstates, i)	
	if exclude
		itemstates_new = squeeze(itemstates, excluded, i)
		cdcp.obj = obj
		return itemstates_new, value0, excluded
	end
	# check include
	include, obj, value1 = isincluded(cdcp, itemstates, i)
	if include
		itemstates_new = squeeze(itemstates, included, i)
		cdcp.obj = obj
		return itemstates_new, value1, included
	end
	# no progress
	value = convert(typeof(cdcp.value), -Inf)
	begin
		itemstates_new = setindex(itemstates, aux, i)
		cdcp.obj = obj
		return itemstates_new, value, aux
	end
end

function isexcluded(cdcp::CDCProblem{<:Squeezing}, itemstates::AbstractVector{ItemState}, i::Int)
	obj, z = cdcp.obj, cdcp.solver.z
	if cdcp.solver.scdca
		sub = to_sub(itemstates)
		obj = setℒ(obj, sub)
		value1, value0, obj = margin(obj, i, z)
		exclude = (value1 - value0 <= 0)
	else
		sup = to_sup(itemstates)
		obj = setℒ(obj, sup)
		value1, value0, obj = margin(obj, i, z)
		exclude = (value1 - value0 < 0)
	end
	exclude, obj, value0
end

function isincluded(cdcp::CDCProblem{<:Squeezing}, itemstates::AbstractVector{ItemState}, i::Int)
	obj, z = cdcp.obj, cdcp.solver.z
	if cdcp.solver.scdca
		sup = to_sup(itemstates)
		obj = setℒ(obj, sup)
		value1, value0, obj = margin(obj, i, z)
		include = (value1 - value0 > 0)
	else
		sub = to_sub(itemstates)
		obj = setℒ(obj, sub)
		value1, value0, obj = margin(obj, i, z)
		include = (value1 - value0 >= 0)
	end
	return include, obj, value1
end

function squeeze(itemstates::SVector{S,ItemState}, s::ItemState, k::Int) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(squeeze(itemstates, s, k, $i)))
		end
		return :(SVector{S,ItemState}($ex))
	else
		return SVector{S,ItemState}(ntuple(i->squeeze(itemstates, s, k, i), S))
	end
end

function squeeze(itemstates::SVector{S,ItemState}, s::ItemState, k::Int, i::Int) where S
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

