function solve!(cdcp::CDCProblem{<:Squeezing}; restart::Bool=false, scdca=cdcp.solver.scdca, z=cdcp.solver.z)
	restart && (cdcp = reinit!(cdcp; scdca, z))
	cdcp.x, cdcp.value, cdcp.state = squeeze!(cdcp, cdcp.x)
	(cdcp.state != success) && squeeze_branch!(cdcp)
	cdcp
end

function init_solverx(::Type{<:Squeezing}, obj, scdca::Bool; z=nothing, kwargs...)
	itemstates = allundetermined(obj)
	A = typeof(itemstates)
	return Squeezing(scdca, A[], z), itemstates
end

function reinit!(cdcp::CDCProblem{<:Squeezing}; scdca=cdcp.solver.scdca, z=cdcp.solver.z)
	cdcp.obj = clearfcall(cdcp.obj)
	empty!(cdcp.solver.branching)
	cdcp.x = allundetermined!(cdcp.x)
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
