function solve!(cdcp::CDCProblem{<:Squeezing}; restart::Bool=false, scdca=cdcp.solver.scdca, z=cdcp.solver.z)
	restart && (cdcp = _reinit!(cdcp; scdca=scdca, z=z))
	cdcp.x, cdcp.fx, cdcp.state = squeeze!(cdcp, cdcp.x)
	cdcp.state == success || cdcp.state == maxfcall_reached && return cdcp
	cdcp.state = squeeze_branch!(cdcp)
	return cdcp
end

function _init(::Type{<:Squeezing}, obj, scdca::Bool; z=nothing, trace::Bool=false, valtype::Type=Float64, kwargs...)
	S = length(obj.x)
	if obj.x isa SVector
		x = _fillstate(SVector{S,ItemState}, undetermined)
	else
		x = fill(undetermined, S)
	end
	A = typeof(x)
	tr = trace ? [SqueezingTrace{typeof(x),valtype}[]] : nothing
	return Squeezing(scdca, A[], z, tr), x
end

function _reinit!(cdcp::CDCProblem{<:Squeezing}; scdca=cdcp.solver.scdca, z=cdcp.solver.z)
	cdcp.obj = _clearfcall(cdcp.obj)
	empty!(cdcp.solver.branching)
	if cdcp.solver.trace !== nothing
		empty!(cdcp.solver.trace)
		push!(eltype(cdcp.solver.trace)())
	end
	if cdcp.x isa SVector
		S = length(cdcp.x)
		cdcp.x = _fillstate(SVector{S,ItemState}, undetermined)
	else
		fill!(cdcp.x, undetermined)
	end
	cdcp.fx = -Inf
	cdcp.solver = Squeezing(scdca, cdcp.solver.branching, z, cdcp.solver.trace)
	return cdcp
end

function setsub(x::SVector{S,ItemState}) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(x[$i] == included))
		end
		return :(SVector{S,Bool}($ex))
	else
		return SVector{S,Bool}(ntuple(i->x[i] == included, S))
	end
end

function setsup(x::SVector{S,ItemState}) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(x[$i] != excluded))
		end
		return :(SVector{S,Bool}($ex))
	else
		return SVector{S,Bool}(ntuple(i->x[i] != excluded, S))
	end
end

function _setitemstate(x::SVector{S,ItemState}, s::ItemState, i::Int) where S
	setindex(x, s, i)
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

function _setchoice(obj::Objective{<:Any,A}, x::A) where A
	Objective(obj.f, x, obj.fcall)
end
