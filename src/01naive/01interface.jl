function init_solverx(::Type{Naive}, obj; z=nothing, kwargs...)
	init_solverx(Naive([0], z), obj)
end
function init_solverx(solver::Naive, obj; z=solver.z, kwargs...)
	Naive(solver.ids, z, solver.incremental), copy(obj.ℒ)
end
