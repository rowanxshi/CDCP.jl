function init_solverx(::Type{Naive}, obj, args...; z=nothing, kwargs...)
	solver = Naive([0], z)
	x = copy(obj.ℒ)
	solver, x
end

function setℒ(obj::Objective{<:Any,<:AbstractVector{Bool}}, ids::Vector{Int})
	fill!(obj.ℒ, false)
	fill!(view(obj.ℒ, ids), true)
	obj
end
