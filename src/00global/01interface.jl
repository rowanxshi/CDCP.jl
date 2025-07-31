"""
    solve(Algorithm::CDCPSolver, obj, S::Integer, args...; kwargs...) = cdcp::CDCProblem

Solve a combinatorial discrete choice problem over `S` possible items, using a given solution algorithm that can be [`SqueezingPolicy`](@ref), [`Squeezing`](@ref) or [`Naive`](@ref). Results are returned as a [`CDCProblem`](@ref). See the methods below.

The objective function `obj` should return the value evaluated at a choice vector `ℒ` with an optional parameter `z` that is typically a number (e.g. productivity).

`obj` must have a method of either `obj(ℒ)` (if no parameter is specified) or `obj(ℒ, z)` (if a parameter is specified). `obj` must not restrict the specific type of `ℒ` but only assume `ℒ` is a vector of length `S` with element type being `Bool`. Specifically, `obj` must *not* try to modify the elements in `ℒ` when it is called. It should only read from `ℒ` with `getindex`.

    solve(Squeezing, obj, scdca::Bool; z=nothing)

The problem should satisfy SCD-C from above if `scdca` is `true` and SCD-C from below if `scdca` is `false`.

Keywords
===
* `z=nothing`


    solve(SqueezingPolicy, obj, S::Integer, scdca::Bool, equal_obj, zbounds::Tuple{Z,Z}=(-Inf, Inf); kwargs...)

The problem should satisfy SCD-C from above if `scdca` is `true` and SCD-C from below if `scdca` is `false`.

The function `equal_obj` should returns the cutoff type `z` with `zleft <= z0 <= zright` that is indifferent between two decision sets `ℒ1` and `ℒ2`.

`equal_obj` can be defined with one of the two following methods:
* `equal_obj((ℒ1, ℒ2), zleft, zright)`: where the pair of input choices is accepted as a tuple
* `equal_obj(obj1::Objective, obj2::Objective, zleft, zright)`: where `obj1` and `obj2` are the same objective function attached with different input vectors `obj1.ℒ` and `obj2.ℒ` that correspond to `ℒ1` and `ℒ2` respectively

`zbounds` determines the domain of the parameter `z`, which defaults to `(-Inf, Inf)`.

Keywords
===
* `zero_margin=nothing`: a function that returns the type `z` that is indifferent to adding item `ℓ` to decision set `ℒ`; it should have the method `zero_margin(obj::Objective, ℓ::Int, zleft, zright)` with `ℒ` attached to `obj` and if not supplied, this function will be generated automatically from `equal_obj`
* `policy0::Policy=Policy(obj, zbounds)`: initial policy which on top of which optimisation occurs
* `skiprefinement::Bool=false`: skip the refinement step that searches resolves intervals for which squeezing doesn't pin down optimal policy
* `singlekw=NamedTuple()`: keyword arguments passed to the single-agent solver as a `NamedTuple` used in the branching stage.
* `maxfcall=1_000_000_000`: maximum calls to the objective function


    solve(Naive, obj, S::Integer; kwargs...)

Solve with exhaustion by checking all potential decision sets. This method should be used with extreme caution.

Keywords
===
* `z=nothing`: the type parameterising the problem
"""
function solve(::Type{Algorithm}, args...; kwargs...) where Algorithm<:CDCPSolver
	solve!(init(Algorithm, args...; kwargs...))
end

function init(::Type{Algorithm}, obj, S::Integer, args...; valuetype::Type=Float64, kwargs...) where Algorithm <: CDCPSolver
	obj = init_Objective(obj, S)
	solver, x = init_solverx(Algorithm, obj, args...; kwargs...)
	return CDCProblem{typeof(solver), typeof(obj), typeof(x), valuetype}(solver, obj, x, convert(valuetype,-Inf), inprogress)
end
function reinit!(cdcp::CDCProblem, S; obj=cdcp.obj, fcall=true)
	cdcp.value = convert(typeof(cdcp.value), -Inf)
	cdcp.state = inprogress
	cdcp.obj = reinit!(obj, S; fcall)
	cdcp
end
function reinit!(obj::Objective, S; fcall=true)
	if obj isa Objective
		fcall && (obj = clearfcall(obj))
	else
		ℒ = (S < static_threshold()) ? SVector{S, Bool}(ntuple(i->false, S)) : Vector{Bool}(undef, S)
		obj = Objective(obj, ℒ)
	end
	obj
end

function allundetermined!(itemstates::SVector)
	allundetermined(itemstates)
end
function allundetermined!(itemstates::AbstractVector)
	fill!(itemstates, undetermined)
end
function allundetermined(obj::Objective)
	allundetermined(obj.ℒ)
end
function allundetermined(::SVector{S}) where S
	fillstate(SVector{S,ItemState}, undetermined)
end
function allundetermined(ℒ::AbstractVector)
	fill(undetermined, length(ℒ))
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

