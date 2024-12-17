struct IntervalChoice{Z,A<:AbstractVector{ItemState}}# <: AbstractVector{ItemState}
	lb::Z
	ub::Z
	itemstates::A
end

function Base.size(intervalchoice::IntervalChoice)
	size(intervalchoice.itemstates)
end
Base.@propagate_inbounds function Base.getindex(intervalchoice::IntervalChoice, i::Int)
	intervalchoice.itemstates[i]
end
function _lb(intervalchoice::IntervalChoice)
	intervalchoice.lb
end

struct Policy{Z,A} <: AbstractVector{IntervalChoice{Z,A}}
	cutoffs::Vector{Z}
	itemstates_s::Vector{A}
	ub::Z
end

function Base.size(policy::Policy)
	size(policy.itemstates_s)
end
Base.@propagate_inbounds function Base.getindex(policy::Policy, i::Int)
	IntervalChoice(policy.cutoffs[i], i+1>length(policy.itemstates_s) ? policy.ub : policy.cutoffs[i+1], policy.itemstates_s[i])
end

struct DiffObj{Obj}
	obj1::Obj
	obj2::Obj
end

function (d::DiffObj)(z, fcall)
	# ℒ should have been set
	f1, obj1 = d.obj1(z)
	f2, obj2 = d.obj2(z)
	fcall[] += 2
	return f2 - f1
end

struct Equal_Obj{M,KW<:NamedTuple}
	m::M
	kwargs::KW
	fcall::RefValue{Int}
end

function Equal_Obj(m, kwargs::NamedTuple=NamedTuple())
	Equal_Obj(m, kwargs, Ref(0))
end

struct Default_Zero_Margin{F, O}
	equal_obj::F
	obj2::O
end

function (zm::Default_Zero_Margin)(obj::Objective, i, lb, ub)
	# the order of true and false matters in case there is no interior solution
	obj = _setℒ(obj, true, i)
	obj2 = _setℒ(_copyℒ(zm.obj2, obj.ℒ), false, i)
	return zm.equal_obj(obj, obj2, lb, ub)
end

"""
    SqueezingPolicy{Z,A,AO,F1,F2,O,S,TR} <: CDCPSolver

A type for solving a [`CDCProblem`](@ref) with a policy method as in Arkolakis, Eckert and Shi (2023).

# Usage
	solve(SqueezingPolicy, obj, scdca::Bool, equal_obj, zbounds::Tuple{Z,Z}=(-Inf, Inf); kwargs...)
	solve!(cdcp::CDCProblem{<:SqueezingPolicy}; restart::Bool=false)

Pass the type `SqueezingPolicy` as the first argument to `solve` indicates the use of the policy method for the problem. Users are required to specify the objective function `obj` that returns the value evaluated at a choice vector `ℒ` with a parameter `z` that is a number. `obj` must have a method of `obj(ℒ, z)` with `ℒ` being a Boolean choice vector. `obj` must not restrict the specific type of `ℒ` but only assume `ℒ` is a vector with element type being `Bool`. Specifically, `obj` must *not* try to modify the elements in `ℒ` when it is called. It should only read from `ℒ` with `getindex`. The problem should satisfy SCD-C from above if `scdca` is `true` and SCD-C from below if `scdca` is `false`. `zbounds` determines the range of the parameter `z`, which could be `(-Inf, Inf)` if `z` can be any real number.

`equal_obj` is a user-specified function that returns the cutoff point `z0` such that for a given pair of input choices `ℒ1` and `ℒ2`, `obj(ℒ1, z0)` equals to `obj(ℒ2, z0)` with `zl <= z0 <= zr`. It can be defined with one of the two alternative methods: - `equal_obj((ℒ1, ℒ2), zl, zr))` where the pair of input choices is accepted as a tuple.
- `equal_obj(obj1::Objective, obj2::Objective, zl, zr)` where `obj1` and `obj2`
are the same objective function attached with different input vectors
`obj1.ℒ` and `obj2.ℒ` that correspond to `ℒ1` and `ℒ2` respectively.

## Keywords
- `zero_margin=nothing`: An optionally specified function that returns `z0` such that the `i`th margin of the objective function is zero at `ℒ`. The function has a method `zero_margin(obj::Objective, i::Int, lb, ub)` with `ℒ` attached to `obj`.
- `ntasks=1`: Number of threads used in the branching process.
- `nobranching::Bool=false`: Skip the branching stage; only for inspecting the solver.
- `singlekw=NamedTuple()`: keyword arguments passed to the single-agent solver as a `NamedTuple`; a single-agent solver is used in the branching stage.

!!! info

	In case a cutoff point is not found,
	`equal_obj` or `zero_margin` should return `NaN` but not `nothing`.
	This requirement is a breaking change from earlier implementation.
"""
struct SqueezingPolicy{Z,A,AO,F1,F2,O,S} <: CDCPSolver
	scdca::Bool
	intervalchoices::Vector{IntervalChoice{Z,A}}
	squeezing::Vector{Int}
	branching::Vector{Int}
	lookup_zero_margin::Dict{Tuple{Int,AO},Z}
	zero_margin::F1
	matcheds::Vector{Vector{IntervalChoice{Z,A}}}
	singlesolvers::Vector{S}
	equal_obj::F2
	obj2::O
	zero_margin_call::RefValue{Int}
	equal_obj_call::RefValue{Int}
	nobranching::Bool
	maxfcall::Int
end
