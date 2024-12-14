"""
    solve(SqueezingPolicy, obj, scdca::Bool, equal_obj, zbounds::Tuple{Z,Z}; kwargs...)
	solve(Squeezing, obj, scdca::Bool; kwargs...)
	solve(BruteForce, obj; kwargs...)

Solve a combinatorial discrete choice problem with a given solution algorithm
that can be [`SqueezingPolicy`](@ref), [`Squeezing`](@ref) or [`BruteForce`](@ref).
Results are returned as a [`CDCProblem`](@ref).
For details on usage,
see [`SqueezingPolicy`](@ref), [`Squeezing`](@ref) or [`BruteForce`](@ref) respectively.
An in-place version [`solve!`](@ref) can be used when a `CDCProblem` is preallocated.
"""
function solve(::Type{M}, args...; kwargs...) where M<:CDCPSolver
	solve!(init(M, args...; kwargs...))
end

function init(::Type{M}, obj, S::Integer, args...; maxfcall::Integer=1_000_000_000, valtype::Type=Float64, kwargs...) where M<:CDCPSolver
	S = Int(S)
	S > 0 || throw(ArgumentError("the number of items S must be a positive integer"))
	if obj isa Objective
		obj = _clearfcall(obj)
		length(obj.x) == S || throw(ArgumentError("length of obj.x is not $S"))
	else
		obj = Objective(obj, S < _static_threshold() ?
			SVector{S, Bool}(ntuple(i->false, S)) : Vector{Bool}(undef, S))
	end
	solver, x = _init(M, obj, args...; valtype=valtype, kwargs...)
	return CDCProblem{typeof(solver), typeof(obj), typeof(x), valtype}(
		solver, obj, x, convert(valtype,-Inf), maxfcall, inprogress)
end

# Default threshold for determining whether SVector is used for a choice
function _static_threshold()
	256
end
