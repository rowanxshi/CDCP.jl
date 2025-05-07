function solve!(cdcp::CDCProblem{<:Squeezing}; restart::Bool=false, scdca=cdcp.solver.scdca, z=cdcp.solver.z)
	restart && (cdcp = reinit!(cdcp; scdca, z))
	@debug "squeezing"
	cdcp.x, cdcp.value, cdcp.state = squeeze!(cdcp, cdcp.x)
	if (cdcp.state != success) && cdcp.solver.branch
		@debug("branching")
		branch!(cdcp)
	end
	return cdcp
end

function init_solverx(::Type{<:Squeezing}, obj, scdca::Bool; z=nothing, branch = true, kwargs...)
	itemstates = allundetermined(obj)
	V = typeof(itemstates)
	return Squeezing(scdca, V[], z, branch), itemstates
end

function reinit!(cdcp::CDCProblem{<:Squeezing}; obj=cdcp.obj, fcall=true, solverkw...)
	S = length(cdcp.x)
	cdcp = reinit!(cdcp, S; obj, fcall)
	cdcp = reinit_solverx!(cdcp; solverkw...)
	return cdcp
end
function reinit_solverx!(cdcp::CDCProblem{<:Squeezing}; solverkw...)
	cdcp.x = allundetermined!(cdcp.x)
	cdcp.solver = reinit!(cdcp.solver; solverkw...)
	cdcp
end
function reinit!(solver::Squeezing; scdca=solver.scdca, z=solver.z, branch=solver.branch)
	solver = Squeezing(scdca, empty!(solver.branching), z, branch)
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
