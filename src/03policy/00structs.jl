struct IntervalChoice{Z,V<:AbstractVector{ItemState}}
	zleft::Z
	zright::Z
	itemstates::V
end
function Base.size(intervalchoice::IntervalChoice)
	size(intervalchoice.itemstates)
end
Base.@propagate_inbounds function Base.getindex(intervalchoice::IntervalChoice, i::Int)
	intervalchoice.itemstates[i]
end

"""
    Policy{Z,V<:AbstractVector{ItemState}} <: AbstractVector{IntervalChoice{Z,V}}

Policy function derived from [`SqueezingPolicy`](@ref).

Given a policy function `policy::Policy`, the optimal decision set of a type `z` can be retrieved as expected, i.e. `policy(z)`.

See also [`solve`](@ref), [`CDCProblem`](@ref).

Fields
===
* `cutoffs::Vector{Z}`: cutoff types at which optimal decision set switches
* `itemstates_s::Vector{V}`: a vector of decision sets, so that `itemstates_s[k]` corresponds to the types above `cutoffs[k]` and below `cutoffs[k+1]`
* `zright::Z`: the maximum type `z` (can be `Inf`).

"""
struct Policy{Z,V <: AbstractVector{ItemState}} <: AbstractVector{IntervalChoice{Z,V}}
	cutoffs::Vector{Z}
	itemstates_s::Vector{V}
	zright::Z
end
function (policy::Policy)(z)
	(z <= first(policy.cutoffs)) && error("provided type $z is less than the lowest cutoff $first(policy.cutoffs)")
	k = findlast(<=(z), policy.cutoffs)
	policy.itemstates_s[k]
end
function Policy(obj::Objective, zbounds; itemstates=allundetermined(obj))
	S = length(obj.ℒ)
	Policy([first(zbounds)], [itemstates], last(zbounds))
end
function Base.size(policy::Policy)
	size(policy.itemstates_s)
end
Base.@propagate_inbounds function Base.getindex(policy::Policy, i::Int)
	IntervalChoice(policy.cutoffs[i], i+1>length(policy.itemstates_s) ? policy.zright : policy.cutoffs[i+1], policy.itemstates_s[i])
end

struct DiffObj{Obj}
	obj1::Obj
	obj2::Obj
end
function (d::DiffObj)(z, fcall::RefValue)
	# ℒ should have been set
	value1, _ = d.obj1(z)
	value2, _ = d.obj2(z)
	fcall[] += 2
	(value2 - value1)
end

"""
	struct Default_Zero_Margin{F, O, D <: AbstractDict}

A type that acts as the _memoised_ function:

	`(::Default_Zero_Margin)(obj::Objective, i, x...)`

which uses `equal_obj` to calculate the type at which the marginal value of `i` in `obj.ℒ` is zero. Solver constructs and uses this type for `SqueezingPolicy.zero_margin` unless user supplies one; to avoid memoisation, supply `SqueezingPolicy.zero_margin`.

See also [`SqueezingPolicy`](@ref).

Fields
===
* `equal_obj::F`
* `obj2::O`
* `memo::D`
"""
struct Default_Zero_Margin{F, O, D <: AbstractDict}
	equal_obj::F
	obj2::O
	memo::D
end

function Default_Zero_Margin(equal_obj, obj2::Objective{<: Any, V}, ::Z) where {V, Z <: Real}
	Default_Zero_Margin(equal_obj, obj2, Dict{Tuple{Int,typeof(packℒ(obj2.ℒ))},Z}())
end
function (zm::Default_Zero_Margin)(obj::Objective, i, zleft, zright)
	# the order of true and false matters in case there is no interior solution
	key = (i, packℒ(obj.ℒ))
	z = get(zm.memo, key, nothing)
	if isnothing(z)
		obj = setindex(obj, true, i)
		obj2 = setindex(ℒ!(zm.obj2, obj.ℒ), false, i)
		z, obj = zm.equal_obj(obj, obj2, zleft, zright)
		zm.memo[key] = z
	end
	z, obj
end
function reinit!(zm::Default_Zero_Margin; memo = true)
	memo && empty!(zm.memo)
end
function packℒ(ℒ::AbstractVector{Bool})
	BitVector(ℒ).chunks
end
@generated function packℒ(ℒ::SVector{C,Bool}) where C
	len = cld(C, 64)
	words = map(1:len) do k
		word = :(zero(UInt64))
		for j in 1:min(64, C - 64(k-1))
			word = :($word | (UInt64(ℒ[$(64(k-1)+j)]) << $(j-1)))
		end
		word
	end
	Expr(:tuple, words...)
end

function ℒ!(obj::Objective{<:Any, V}, ℒ::V) where V<:SVector
	Objective(obj.f, ℒ, obj.fcall, obj.memo)
end
function ℒ!(obj::Objective{<:Any, V}, ℒ::V) where V
	Objective(obj.f, copyto!(obj.ℒ, ℒ), obj.fcall, obj.memo)
end

"""
    SqueezingPolicy{Z,V<:AbstractVector{ItemState},F1,F2,S,O} <: CDCPSolver

A type for solving a [`CDCProblem`](@ref) with policy function squeezing.

See also [`solve`](@ref).

Fields
===
* `scdca::Bool`
* `intervalchoices::Vector{IntervalChoice{Z,V}}`
* `squeezing_indices::Vector{Int}`
* `branching_indices::Vector{Int}`
* `zero_margin::F1`
* `equal_obj::F2`
* `singlecdcp::S`
* `obj2::O`
* `skiprefinement::Bool`
* `maxfcall::Int`
* `lastitemchecked::Vector{Int}`
"""
struct SqueezingPolicy{Z,V<:AbstractVector{ItemState},F1,F2,S,O} <: CDCPSolver
	scdca::Bool
	intervalchoices::Vector{IntervalChoice{Z,V}}
	squeezing_indices::Vector{Int}
	branching_indices::Vector{Int}
	zero_margin::F1
	equal_obj::F2
	singlecdcp::S
	obj2::O
	skiprefinement::Bool
	maxfcall::Int
	lastitemchecked::Vector{Int}
end
