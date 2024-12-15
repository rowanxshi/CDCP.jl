function _init(::Type{BruteForce}, obj, args...; z=nothing, kwargs...)
	BruteForce([0], z), copy(obj.ℒ)
end

function _setchoice(obj::Objective{<:Any,<:AbstractVector{Bool}}, ids::Vector{Int})
	(fill!(obj.ℒ, false); fill!(view(obj.ℒ, ids), true); obj)
end
