function solve!(cdcp::CDCProblem{<:Squeezing}; restart::Bool=false, scdca=cdcp.solver.scdca, z=cdcp.solver.z)
	restart && (cdcp = _reinit!(cdcp; scdca, z))
	cdcp.x, cdcp.value, cdcp.state = squeeze!(cdcp, cdcp.x)
	(cdcp.state != success) && squeeze_branch!(cdcp)
	cdcp
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

function next_undetermined(itemstates, i)
	if i < length(itemstates)
		i = findnext(==(undetermined), itemstates, i+1)
		isnothing(i) && (i = next_undetermined(itemstates))
	else
		i = next_undetermined(itemstates)
	end
	i
end

function next_undetermined(itemstates)
	findfirst(==(undetermined), itemstates)
end
