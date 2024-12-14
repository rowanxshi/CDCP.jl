_setitemstate(x::SVector{S,ItemState}, s::ItemState, i::Int) where S = setindex(x, s, i)

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
