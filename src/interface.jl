struct Objective{F,A}
    f::F
    x::A
    fcall::Int
end

Objective(f, x) = Objective(f, x, 0)

_setx(obj::Objective{<:Any,<:SVector}, v, i) =
    Objective(obj.f, setindex(obj.x, v, i), obj.fcall)

# Fallback method assumes x is mutable
_setx(obj::Objective, v, i) =
    Objective(obj.f, setindex!(obj.x, v, i), obj.fcall)

addfcall(obj::Objective, n=1) = Objective(obj.f, obj.x, obj.fcall+n)
_clearfcall(obj::Objective) = Objective(obj.f, obj.x, 0)

value(obj::Objective, z) = obj.f(obj.x, z), addfcall(obj)
value(obj::Objective, ::Nothing) = obj.f(obj.x), addfcall(obj)

function margin(obj::Objective, i::Int, z)
    obj = _setx(obj, true, i)
    f1, obj = value(obj, z)
    obj = _setx(obj, false, i)
    f0, obj = value(obj, z)
    return f1, f0, obj
end

abstract type CDCPSolver end

@enum SolverState::Int8 begin
    inprogress
    success
    maxfcall_reached
end

@enum ItemState::Int8 begin
    undetermined
    included
    excluded
    aux
end

mutable struct CDCP{M<:CDCPSolver, O<:Objective, A, TF<:AbstractFloat}
    solver::M
	obj::O
    x::A
	fx::TF
	maxfcall::Int
    state::SolverState
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
    return CDCP{typeof(solver), typeof(obj), typeof(x), valtype}(
		solver, obj, x, convert(valtype,-Inf), maxfcall, inprogress)
end

solve(::Type{M}, args...; kwargs...) where M<:CDCPSolver =
    solve!(init(M, args...; kwargs...))
