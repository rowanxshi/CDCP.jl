function solve!(cdcp::CDCProblem{<:Squeezing}; restart::Bool=false, scdca=cdcp.solver.scdca, z=cdcp.solver.z)
	restart && (cdcp = _reinit!(cdcp; scdca, z))
	cdcp.x, cdcp.value, cdcp.state = squeeze!(cdcp, cdcp.x)
	cdcp.state == success || cdcp.state == maxfcall_reached && return cdcp
	cdcp.state = squeeze_branch!(cdcp)
	return cdcp
end

function init_solverx(::Type{<:Squeezing}, obj, scdca::Bool; z=nothing, kwargs...)
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

