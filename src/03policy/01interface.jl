# TODO: still nedd to review this whole file
function solve!(cdcp::CDCProblem{<:SqueezingPolicy}; restart::Bool=false, obj=cdcp.obj, zero_margin=cdcp.solver.zero_margin, equal_obj=cdcp.solver.equal_obj, scdca=cdcp.solver.scdca)
	restart && (cdcp = _reinit!(cdcp; obj, zero_margin, equal_obj, scdca))
	cdcp.state = squeeze!(cdcp)
	if cdcp.state == maxfcall_reached
		@warn "maxfcall is reached before convergence"
		return cdcp
	end
	cdcp.solver.nobranching || (cdcp.state = branching!(cdcp))
	concat!(cdcp)
	cdcp.solver.nobranching || (cdcp.state = success)
	return cdcp
end

function init_solverx(::Type{<:SqueezingPolicy}, obj, scdca::Bool, equal_obj, zbounds::Tuple{Z,Z}=(-Inf, Inf); zero_margin=nothing, policy0=nothing, ntasks=1, nobranching::Bool=false, singlekw=NamedTuple(), maxfcall=1_000_000_000, kwargs...) where Z
	S = length(obj.ℒ)
	if isnothing(policy0)
		policy = Policy(obj, zbounds)
	else
		policy = policy0
	end
	obj2 = deepcopy(obj)
	# Harmonize user defined functions
	if !applicable(equal_obj, obj, obj2, zbounds...)
		# Assume equal_obj follows the old requirement for equalise_obj
		equal_obj = Wrapped_Equalise_Obj(equal_obj)
	end
	if isnothing(zero_margin)
		zero_margin = Default_Zero_Margin(equal_obj, obj2)
	elseif !applicable(zero_margin, obj, 1, zbounds...)
		@warn "Consider adapting `zero_margin` to the new method"
		zero_margin = Wrapped_Zero_D_j_Obj(zero_margin, fill(false, S))
	end
	intervalchoices = [policy[i] for i in eachindex(policy.itemstates_s)]
	A = eltype(policy.itemstates_s)
	matcheds = [IntervalChoice{Z,A}[] for _ in 1:ntasks]
	singlesolvers = [init(Squeezing, obj, S, scdca; z=zero(Z), singlekw...) for _ in 1:ntasks]
	return SqueezingPolicy(scdca, intervalchoices, collect(1:length(policy.itemstates_s)), Int[], zero_margin, Ref(0), Dict{Tuple{Int,typeof(obj.ℒ)},Z}(), equal_obj, Ref(0), matcheds, singlesolvers, obj2, nobranching, maxfcall), policy
end

function _reset!(cdcp::CDCProblem{<:Squeezing}, z, itemstates)
	cdcp.state = inprogress
	cdcp.solver = Squeezing(cdcp.solver.scdca, empty!(cdcp.solver.branching), z)
	cdcp.x = itemstates
	cdcp.value = -Inf
	return cdcp
end

function _reinit!(cdcp::CDCProblem{<:SqueezingPolicy}; obj=cdcp.obj, zero_margin=cdcp.solver.zero_margin, equal_obj=cdcp.solver.equal_obj, scdca=cdcp.solver.scdca)
	S = length(cdcp.x.itemstates_s[1])
	if obj isa Objective
		cdcp.obj = _clearfcall(obj)
	else
		cdcp.obj = Objective(obj, S < _static_threshold() ? SVector{S, Bool}(ntuple(i->false, S)) : Vector{Bool}(undef, S))
	end
	zmin = cdcp.x.cutoffs[1]
	zbounds = (zmin, cdcp.x.zright)
	resize!(cdcp.x.cutoffs, 1)
	cdcp.x.cutoffs[1] = zmin
	resize!(cdcp.x.itemstates_s, 1)
	cdcp.x.itemstates_s[1] = allundetermined(obj)
	cdcp.value = convert(typeof(cdcp.value), -Inf)
	cdcp.state = inprogress
	solver = cdcp.solver
	obj2 = deepcopy(cdcp.obj)
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
	resize!(solver.intervalchoices, 1)
	solver.intervalchoices[1] = cdcp.x[1]
	resize!(solver.squeezing, 1)
	solver.squeezing[1] = 1
	empty!(solver.branching)
	empty!(solver.lookup_zero_margin)
	for m in solver.matcheds
		empty!(m)
	end
	for s in solver.singlesolvers
		s.obj = cdcp.obj
		_reinit!(s; scdca=scdca)
	end
	solver.zero_margin_call[] = 0
	solver.equal_obj_call[] = 0
	cdcp.solver = SqueezingPolicy(scdca, solver.intervalchoices, solver.squeezing, solver.branching, zero_margin, solver.zero_margin_call, solver.lookup_zero_margin, equal_obj, solver.equal_obj_call, solver.matcheds, solver.singlesolvers, obj2, solver.nobranching, solver.maxfcall)
	return cdcp
end

function _copyℒ(obj::Objective{<:Any, A}, ℒ::A) where A<:SVector
	Objective(obj.f, ℒ, obj.fcall)
end

# Fallback method assumes A is a mutable array
function _copyℒ(obj::Objective{<:Any, A}, ℒ::A) where A
	Objective(obj.f, copyto!(obj.ℒ, ℒ), obj.fcall)
end
