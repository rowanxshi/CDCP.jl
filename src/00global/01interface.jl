"""
    solve(SqueezingPolicy, obj, scdca::Bool, equal_obj, zbounds::Tuple{Z,Z}; kwargs...)
    solve(Squeezing, obj, scdca::Bool; kwargs...)
    solve(Naive, obj; kwargs...)

Solve a combinatorial discrete choice problem with a given solution algorithm that can be [`SqueezingPolicy`](@ref), [`Squeezing`](@ref) or [`Naive`](@ref). Results are returned as a [`CDCProblem`](@ref).

For details on usage, see [`SqueezingPolicy`](@ref), [`Squeezing`](@ref) or [`Naive`](@ref) respectively.

An in-place version [`solve!`](@ref) can be used when a `CDCProblem` is preallocated.
"""
function solve(::Type{M}, args...; kwargs...) where M<:CDCPSolver
	solve!(init(M, args...; kwargs...))
end

function init(::Type{M}, obj, S::Integer, args...; valuetype::Type=Float64, kwargs...) where M<:CDCPSolver
	obj = init_Objective(obj, S)
	solver, x = init_solverx(M, obj, args...; kwargs...)
	return CDCProblem{typeof(solver), typeof(obj), typeof(x), valuetype}(solver, obj, x, convert(valuetype,-Inf), inprogress)
end

function allundetermined!(itemstates)
	if itemstates isa SVector
		S = length(itemstates)
		itemstates = fillstate(SVector{S,ItemState}, undetermined)
	else
		fill!(itemstates, undetermined)
	end
	itemstates
end
function allundetermined(obj::Objective)
	allundetermined(obj.ℒ)
end
function allundetermined(ℒ::AbstractVector)
	S = length(ℒ)
	if ℒ isa SVector
		itemstates = fillstate(SVector{S,ItemState}, undetermined)
	else
		itemstates = fill(undetermined, S)
	end
end

function fillstate(::Type{<:SVector{S}}, s::ItemState) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(s))
		end
		return :(SVector{S,ItemState}($ex))
	else
		return SVector{S,ItemState}(ntuple(i->s, S))
	end
end

function to_sub(itemstates::SVector{S,ItemState}) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(itemstates[$i] == included))
		end
		return :(SVector{S,Bool}($ex))
	else
		return SVector{S,Bool}(ntuple(i->itemstates[i] == included, S))
	end
end

function to_sup(itemstates::SVector{S,ItemState}) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(itemstates[$i] != excluded))
		end
		return :(SVector{S,Bool}($ex))
	else
		return SVector{S,Bool}(ntuple(i->itemstates[i] != excluded, S))
	end
end

