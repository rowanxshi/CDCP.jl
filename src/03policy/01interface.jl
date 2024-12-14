function solve!(p::CDCProblem{<:SqueezingPolicy}; restart::Bool=false,
		obj=p.obj, zero_margin=p.solver.zero_margin, equal_obj=p.solver.equal_obj,
		scdca=p.solver.scdca)
	restart && (p = _reinit!(p; obj=obj, zero_margin=zero_margin, equal_obj=equal_obj,
		scdca=scdca))
	p.state = squeeze!(p)
	if p.state == maxfcall_reached
		@warn "maxfcall is reached before convergence"
		return p
	end
	p.solver.nobranching || (p.state = branching!(p))
	concat!(p)
	p.solver.nobranching || (p.state = success)
	return p
end

function _init(::Type{<:SqueezingPolicy}, obj, scdca::Bool, equal_obj,
		zbounds::Tuple{Z,Z}=(-Inf, Inf);
		zero_margin=nothing, x0=nothing, ntasks=1, trace::Bool=false,
		nobranching::Bool=false, singlekw=NamedTuple(), kwargs...) where Z
	S = length(obj.x)
	if x0 === nothing
		if obj.x isa SVector
			allu = _fillstate(SVector{S,ItemState}, undetermined)
			x = Policy([zbounds[1]], [allu], zbounds[2])
		else
			allu = fill(undetermined, S)
			x = Policy([zbounds[1]], [allu], zbounds[2])
		end
	else
		x = x0
	end
	tr = trace ? [SqueezingPolicyTrace{Z}[]] : nothing
	obj2 = deepcopy(obj)
	# Harmonize user defined functions
	if !applicable(equal_obj, obj, obj2, zbounds...)
		# Assume equal_obj follows the old requirement for equalise_obj
		equal_obj = Wrapped_Equalise_Obj(equal_obj)
	end
	if zero_margin === nothing
		zero_margin = Default_Zero_Margin(equal_obj, obj2)
	elseif !applicable(zero_margin, obj, 1, zbounds...)
		@warn "Consider adapting `zero_margin` to the new method"
		zero_margin = Wrapped_Zero_D_j_Obj(zero_margin, fill(false, S))
	end
	pool = [x[i] for i in eachindex(x.xs)]
	A = eltype(x.xs)
	matcheds = [IntervalChoice{Z,A}[] for _ in 1:ntasks]
	ss = [init(Squeezing, obj, S, scdca; z=zero(Z), singlekw...) for _ in 1:ntasks]
	return SqueezingPolicy(scdca, pool, collect(1:length(x.xs)), Int[],
		Dict{Tuple{Int,typeof(obj.x)},Z}(), zero_margin, matcheds, ss,
		equal_obj, obj2, Ref(0), Ref(0), tr, nobranching), x
end

function _reset!(p::CDCProblem{<:Squeezing}, z, x)
	p.state = inprogress
	p.solver = Squeezing(p.solver.scdca, empty!(p.solver.branching), z, p.solver.trace)
	p.x = x
	p.fx = -Inf
	return p
end

function _reinit!(p::CDCProblem{<:SqueezingPolicy};
		obj=p.obj, zero_margin=p.solver.zero_margin, equal_obj=p.solver.equal_obj,
		scdca=p.solver.scdca)
	S = length(p.x.xs[1])
	if obj isa Objective
		p.obj = _clearfcall(obj)
	else
		p.obj = Objective(obj, S < _static_threshold() ?
			SVector{S, Bool}(ntuple(i->false, S)) : Vector{Bool}(undef, S))
	end
	zmin = p.x.cutoffs[1]
	zbounds = (zmin, p.x.ub)
	if p.obj.x isa SVector
		allu = _fillstate(SVector{S,ItemState}, undetermined)
		resize!(p.x.cutoffs, 1)
		p.x.cutoffs[1] = zmin
		resize!(p.x.xs, 1)
		p.x.xs[1] = allu
	else
		allu = fill(undetermined, S)
		resize!(p.x.cutoffs, 1)
		p.x.cutoffs[1] = zmin
		resize!(p.x.xs, 1)
		p.x.xs[1] = allu
	end
	p.fx = convert(typeof(p.fx), -Inf)
	p.state = inprogress
	sol = p.solver
	if sol.trace !== nothing
		resize!(sol.trace, 1)
		empty!(sol.trace[1])
	end
	obj2 = deepcopy(p.obj)
	# Harmonize user defined functions
	if !applicable(equal_obj, obj, obj2, zbounds...)
		# Assume equal_obj follows the old requirement for equalise_obj
		equal_obj = Wrapped_Equalise_Obj(equal_obj)
	end
	if zero_margin === nothing
		zero_margin = Default_Zero_Margin(equal_obj, obj2)
	elseif !applicable(zero_margin, obj, 1, zbounds...)
		@warn "Consider adapting `zero_margin` to the new method"
		zero_margin = Wrapped_Zero_D_j_Obj(zero_margin, fill(false, S))
	end
	resize!(sol.pool, 1)
	sol.pool[1] = p.x[1]
	resize!(sol.squeezing, 1)
	sol.squeezing[1] = 1
	empty!(sol.branching)
	empty!(sol.lookup_zero_margin)
	for m in sol.matcheds
		empty!(m)
	end
	for s in sol.singlesolvers
		s.obj = p.obj
		_reinit!(s; scdca=scdca)
	end
	sol.zero_margin_call[] = 0
	sol.equal_obj_call[] = 0
	p.solver = SqueezingPolicy(scdca, sol.pool, sol.squeezing, sol.branching,
		sol.lookup_zero_margin, zero_margin, sol.matcheds, sol.singlesolvers,
		equal_obj, obj2, sol.zero_margin_call, sol.equal_obj_call,
		sol.trace, sol.nobranching)
	return p
end

_copyx(obj::Objective{<:Any, A}, x::A) where A<:SVector =
	Objective(obj.f, x, obj.fcall)

# Fallback method assumes A is a mutable array
_copyx(obj::Objective{<:Any, A}, x::A) where A =
	Objective(obj.f, copyto!(obj.x, x), obj.fcall)
