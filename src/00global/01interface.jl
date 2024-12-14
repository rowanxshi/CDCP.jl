_setx(obj::Objective{<:Any,<:SVector}, v, i) =
    Objective(obj.f, setindex(obj.x, v, i), obj.fcall)

# Fallback method assumes x is mutable
_setx(obj::Objective, v, i) =
    Objective(obj.f, setindex!(obj.x, v, i), obj.fcall)

"""
    addfcall(obj::Objective, n=1)

Add `n` to the counter for function call for `obj`.
"""
addfcall(obj::Objective, n=1) = Objective(obj.f, obj.x, obj.fcall+n)
_clearfcall(obj::Objective) = Objective(obj.f, obj.x, 0)

"""
    value(obj::Objective, z)

Evaluate `obj` with parameter `z`.
The parameter is ignored if `z` is `nothing`.
"""
value(obj::Objective, z) = obj.f(obj.x, z), addfcall(obj)
value(obj::Objective, ::Nothing) = obj.f(obj.x), addfcall(obj)

"""
    margin(obj::Objective, i::Int, z)

Evalutate the change in `obj` with optional parameter `z`
when the `i`th item is included or not.
This corresponds to the `D_j` function in earlier implementation.
"""
function margin(obj::Objective, i::Int, z)
    obj = _setx(obj, true, i)
    f1, obj = value(obj, z)
    obj = _setx(obj, false, i)
    f0, obj = value(obj, z)
    return f1, f0, obj
end

# Default threshold for determining whether SVector is used for a choice
_static_threshold() = 256

function init(::Type{M}, obj, S::Integer, args...;
		maxfcall::Integer=1_000_000_000, valtype::Type=Float64,
        kwargs...) where M<:CDCPSolver
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
solve(::Type{M}, args...; kwargs...) where M<:CDCPSolver =
    solve!(init(M, args...; kwargs...))
