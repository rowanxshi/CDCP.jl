function solve!(cdcp::CDCProblem{<:Squeezing}; restart::Bool=false, scdca=cdcp.solver.scdca, z=cdcp.solver.z)
	restart && (cdcp = _reinit!(cdcp; scdca=scdca, z=z))
	cdcp.x, cdcp.value, cdcp.state = squeeze!(cdcp, cdcp.x)
	cdcp.state == success || cdcp.state == maxfcall_reached && return cdcp
	cdcp.state = squeeze_branch!(cdcp)
	return cdcp
end

function _init(::Type{<:Squeezing}, obj, scdca::Bool; z=nothing, kwargs...)
	S = length(obj.ℒ)
	if obj.ℒ isa SVector
		itemstates = _fillstate(SVector{S,ItemState}, undetermined)
	else
		itemstates = fill(undetermined, S)
	end
	A = typeof(itemstates)
	return Squeezing(scdca, A[], z), itemstates
end

function _reinit!(cdcp::CDCProblem{<:Squeezing}; scdca=cdcp.solver.scdca, z=cdcp.solver.z)
	cdcp.obj = _clearfcall(cdcp.obj)
	empty!(cdcp.solver.branching)
	if cdcp.x isa SVector
		S = length(cdcp.x)
		cdcp.x = _fillstate(SVector{S,ItemState}, undetermined)
	else
		fill!(cdcp.x, undetermined)
	end
	cdcp.value = -Inf
	cdcp.solver = Squeezing(scdca, cdcp.solver.branching, z)
	return cdcp
end

function setsub(itemstates::SVector{S,ItemState}) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(itemstates[$i] == included))
		end
		return :(SVector{S,Bool}($ex))
	else
		return SVector{S,Bool}(ntuple(i->itemstates[i] == included, S))
	end
end

function setsup(itemstates::SVector{S,ItemState}) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(itemstates[$i] != excluded))
		end
		return :(SVector{S,Bool}($ex))
	else
		return SVector{S,Bool}(ntuple(i->itemstates[i] != excluded, S))
	end
end

function _setitemstate(itemstates::SVector{S,ItemState}, s::ItemState, i::Int) where S
	setindex(itemstates, s, i)
end

function _fillstate(::Type{<:SVector{S}}, s::ItemState) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(s))
		end
		return :(SVector{S,ItemState}($ex))
	else
		return SVector{S,ItemState}(ntuple(i->s, S))
	end
end

function _setchoice(obj::Objective{<:Any,A}, ℒ::A) where A
	Objective(obj.f, ℒ, obj.fcall)
end
