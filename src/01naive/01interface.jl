function _init(::Type{BruteForce}, obj, args...; z=nothing, kwargs...)
	BruteForce([0], z), copy(obj.x)
end

function _setchoice(obj::Objective{<:Any,<:AbstractVector{Bool}}, ids::Vector{Int})
	(fill!(obj.x, false); fill!(view(obj.x, ids), true); obj)
end
