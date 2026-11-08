function solve!(cdcp::CDCProblem{<:Naive}; restart::Bool=true, ids=cdcp.solver.ids, z=cdcp.solver.z, incremental=cdcp.solver.incremental)
	(; obj) = cdcp
	S = length(cdcp.x)
	S < Sys.WORD_SIZE - 1 || throw(ArgumentError("exhaustion over $S items is not feasible"))
	if restart
		reinit!(cdcp; obj, ids, z, incremental)
	else
		startingset = ids[1] ⊻ (ids[1] >> 1)
		for ℓ in 1:S
			obj = setindex(obj, isodd(startingset >> (ℓ - 1)), ℓ)
		end
	end
	naive!(cdcp)
end

function init_solverx(::Type{Naive}, obj; z=nothing, incremental=nothing, kwargs...)
	solver = Naive([0], z, incremental)
	x = copy(obj.ℒ)
	solver, x
end
function reinit!(cdcp::CDCProblem{<:Naive}; obj=cdcp.obj, fcall=true, solverkw...)
	S = length(cdcp.x)
	cdcp = reinit!(cdcp, S; obj, fcall)
	cdcp = reinit_solverx!(cdcp; solverkw...)
	return cdcp
end
function reinit_solverx!(cdcp::CDCProblem{<:Naive}; solverkw...)
	cdcp.x = allundetermined!(cdcp.x)
	cdcp.solver = reinit!(cdcp.solver; solverkw...)
	cdcp
end
function reinit!(solver::Naive; ids=solver.ids, z=solver.z, incremental=solver.incremental)
	resize!(solver.ids, 1)
	solver.ids[1] = 0
	solver = Naive(ids, z, incremental)
end
