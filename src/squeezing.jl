"""
	Squeezing{A,Z,TR} <: CDCPSolver

A type for solving a [`CDCProblem`](@ref) with a single-agent method
as in Arkolakis, Eckert and Shi (2023).

# Usage
	solve(Squeezing, obj, scdca::Bool; kwargs...)
	solve!(p::CDCProblem{<:Squeezing}; restart::Bool=false)

Pass the type `Squeezing` as the first argument to `solve`
indicates the use of the single-agent method for the problem.
Users are required to specify the objective function `obj`
that returns the value evaluated at a choice vector `x`
with an optional parameter `z` that is typically a number.
`obj` must have a method of either `obj(x)` or `obj(x, z)`
with `x` being a Boolean choice vector.
`obj` must not restrict the specific type of `x`
but only assume `x` is a vector with element type being `Bool`.
Specifically, `obj` must *not* try to modify the elements in `x` when it is called.
It should only read from `x` with `getindex`.
The problem should satisfy SCD-C from above if `scdca` is `true`
and SCD-C from below if `scdca` is `false`.

## Keywords
- `z=nothing`: An optional parameter passed as the second argument to `obj` when evaluating `obj`.
- `trace::Bool=false`: An experimental feature for tracing the solver progress; behavior may be changed at any time in future updates.
- `valtype::Type=Float64`: type of the value returned by `obj`.
"""
struct Squeezing{A,Z,TR} <: CDCPSolver
	scdca::Bool
    branching::Vector{A}
	z::Z
	trace::TR
end

struct SqueezingTrace{V,TF}
	k::Int
	x0::V
	newstate::ItemState
	fx::TF
end

function _fillstate(::Type{<:SVector{S}}, s::ItemState) where S
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

function _init(::Type{<:Squeezing}, obj, scdca::Bool;
		z=nothing, trace::Bool=false, valtype::Type=Float64, kwargs...)
	S = length(obj.x)
	if obj.x isa SVector
        x = _fillstate(SVector{S,ItemState}, undetermined)
    else
        x = fill(undetermined, S)
    end
	A = typeof(x)
	tr = trace ? [SqueezingTrace{typeof(x),valtype}[]] : nothing
	return Squeezing(scdca, A[], z, tr), x
end

function setsub(x::SVector{S,ItemState}) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(x[$i] == included))
		end
		return :(SVector{S,Bool}($ex))
	else
		return SVector{S,Bool}(ntuple(i->x[i] == included, S))
	end
end

function setsup(x::SVector{S,ItemState}) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(x[$i] != excluded))
		end
		return :(SVector{S,Bool}($ex))
	else
		return SVector{S,Bool}(ntuple(i->x[i] != excluded, S))
	end
end

_setchoice(obj::Objective{<:Any,A}, x::A) where A =
	Objective(obj.f, x, obj.fcall)

function _squeeze(x::SVector{S,ItemState}, s::ItemState, k::Int, i::Int) where S
	if i == k
		return s
	else
		@inbounds xi = x[i]
		if xi == aux
			return undetermined
		else
			return xi
		end
	end
end

function _squeeze(x::SVector{S,ItemState}, s::ItemState, k::Int) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(_squeeze(x, s, k, $i)))
		end
		return :(SVector{S,ItemState}($ex))
	else
		return SVector{S,ItemState}(ntuple(i->_squeeze(x, s, k, i), S))
	end
end

_setitemstate(x::SVector{S,ItemState}, s::ItemState, i::Int) where S = setindex(x, s, i)

function squeeze!(p::CDCProblem{<:Squeezing}, x::AbstractVector{ItemState}, i::Int)
	obj, scdca, z, tr = p.obj, p.solver.scdca, p.solver.z, p.solver.trace
    # For excluding: look at the best case
	if scdca
		obj = _setchoice(obj, setsub(x))
		f1, f0, obj = margin(obj, i, z)
		exclude = f1 - f0 <= 0
	else
		obj = _setchoice(obj, setsup(x))
		f1, f0, obj = margin(obj, i, z)
		exclude = f1 - f0 < 0
	end
	if exclude
		xnew = _squeeze(x, excluded, i)
		tr === nothing || push!(tr[end], SqueezingTrace(i, x, excluded, f0))
		p.obj = obj
		return xnew, f0, excluded
	end

	# For including: look at the worst case
	if scdca
		obj = _setchoice(obj, setsup(x))
		f1, f0, obj = margin(obj, i, z)
		include = f1 - f0 > 0
	else
		obj = _setchoice(obj, setsub(x))
		f1, f0, obj = margin(obj, i, z)
		include = f1 - f0 >= 0
	end
	if include
		xnew = _squeeze(x, included, i)
		tr === nothing || push!(tr[end], SqueezingTrace(i, x, included, f1))
		p.obj = obj
		return xnew, f1, included
	end

	xnew = _setitemstate(x, aux, i)
	fx = convert(typeof(p.fx), -Inf)
	tr === nothing || push!(tr[end], SqueezingTrace(i, x, aux, fx))
	p.obj = obj
	return xnew, fx, aux
end

function branch(x::AbstractVector{ItemState}, k::Int)
	# All aux are turned into undertermined
	xin = _squeeze(x, included, k)
	# copy is here only in case x is not SVector
	xex = _squeeze(copy(x), excluded, k)
	return xin, xex
end

function squeeze!(p::CDCProblem{<:Squeezing}, x::AbstractVector{ItemState})
	S = length(x)
	fx = -Inf
	lastaux = nothing
	i = findfirst(==(undetermined), x)
	if i === nothing
		lastaux = findfirst(==(aux), x) # May find aux from initial value
		if lastaux === nothing # Last value set by branching
			fx, obj = value(p.obj, p.solver.z)
			p.obj = obj
		end
	end
    while i !== nothing && p.obj.fcall < p.maxfcall
		x, fx, itemstate = squeeze!(p, x, i)
		if itemstate == aux
			lastaux = lastaux === nothing ? i : min(lastaux, i)
		else
			lastaux = nothing
		end
		#lastaux = ifelse(itemstate == aux, i, nothing)
		if i < S
			i = findnext(==(undetermined), x, i+1)
			i === nothing && (i = findfirst(==(undetermined), x))
		else
			i = findfirst(==(undetermined), x)
		end
	end
	# Whenever squeeze! makes progress, any aux is changed to undetermined
	# If aux is in x, it can only be that lastaux !== nothing
	if lastaux === nothing
		state = p.obj.fcall < p.maxfcall ? success : maxfcall_reached
	else
		push!(p.solver.branching, branch(x, lastaux)...)
		state = inprogress
	end
	return x, fx, state
end

function squeeze_branch!(p::CDCProblem{<:Squeezing})
	branching, tr = p.solver.branching, p.solver.trace
	while !isempty(branching)
		x = pop!(branching)
		tr === nothing || push!(tr, similar(tr[1], 0))
		x, fx, state = squeeze!(p, x)
		if state == success
			fx > p.fx && (p.x = x; p.fx = fx)
		elseif state == maxfcall_reached
			push!(branching, x) # Put x back
			return maxfcall_reached
		end
	end
	return success
end

function _reinit!(p::CDCProblem{<:Squeezing}; scdca=p.solver.scdca, z=p.solver.z)
	p.obj = _clearfcall(p.obj)
	empty!(p.solver.branching)
	if p.solver.trace !== nothing
		empty!(p.solver.trace)
		push!(eltype(p.solver.trace)())
	end
	if p.x isa SVector
		S = length(p.x)
		p.x = _fillstate(SVector{S,ItemState}, undetermined)
	else
		fill!(p.x, undetermined)
	end
	p.fx = -Inf
	p.solver = Squeezing(scdca, p.solver.branching, z, p.solver.trace)
	return p
end

function solve!(p::CDCProblem{<:Squeezing};
		restart::Bool=false, scdca=p.solver.scdca, z=p.solver.z)
	restart && (p = _reinit!(p; scdca=scdca, z=z))
	p.x, p.fx, p.state = squeeze!(p, p.x)
	p.state == success || p.state == maxfcall_reached && return p
	p.state = squeeze_branch!(p)
	return p
end
