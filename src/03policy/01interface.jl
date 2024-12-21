function solve!(cdcp::CDCProblem{<:SqueezingPolicy}; restart::Bool=false, obj=cdcp.obj, zero_margin=cdcp.solver.zero_margin, equal_obj=cdcp.solver.equal_obj, scdca=cdcp.solver.scdca)
	restart && (cdcp = reinit!(cdcp; obj, zero_margin, equal_obj, scdca))
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

function init_solverx(::Type{<:SqueezingPolicy}, obj, scdca::Bool, equal_obj, zbounds::Tuple{Z,Z}=(-Inf, Inf); zero_margin=nothing, policy0::Policy=Policy(obj, zbounds), ntasks=1, nobranching::Bool=false, singlekw=NamedTuple(), maxfcall=1_000_000_000, kwargs...) where Z
	S = length(obj.ℒ)
	obj2 = deepcopy(obj)
	# harmonize user defined functions
	if !applicable(equal_obj, obj, obj2, zbounds...) # TODO: MOVE TO COMPAT SECTION
		# assume equal_obj follows the old requirement for equalise_obj
		equal_obj = Wrapped_Equalise_Obj(equal_obj)
	end
	if isnothing(zero_margin)
		zero_margin = Default_Zero_Margin(equal_obj, obj2, zero(Z))
	elseif !applicable(zero_margin, obj, 1, zbounds...)
		@warn "Consider adapting `zero_margin` to the new method"
		zero_margin = Wrapped_Zero_D_j_Obj(zero_margin, fill(false, S))
	end
	intervalchoices = [policy0[i] for i in eachindex(policy0.itemstates_s)]
	V = eltype(policy0.itemstates_s)
	matcheds = [IntervalChoice{Z,V}[] for _ in 1:ntasks]
	singlesolvers = [init(Squeezing, obj, S, scdca; z=zero(Z), singlekw...) for _ in 1:ntasks]
	return SqueezingPolicy(scdca, intervalchoices, collect(1:length(policy0.itemstates_s)), Int[], zero_margin, equal_obj, matcheds, singlesolvers, obj2, nobranching, maxfcall), policy0
end

function reinit!(cdcp::T; obj=cdcp.obj, zero_margin=cdcp.solver.zero_margin, equal_obj=cdcp.solver.equal_obj, scdca=cdcp.solver.scdca) where {T <: CDCProblem{<:SqueezingPolicy}}
	S = length(cdcp.x.itemstates_s[1])
	if obj isa Objective
		cdcp.obj = clearfcall(obj)
	else
		cdcp.obj = Objective(obj, S < static_threshold() ? SVector{S, Bool}(ntuple(i->false, S)) : Vector{Bool}(undef, S))
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
	# harmonize user defined functions
	if !applicable(equal_obj, obj, obj2, zbounds...)
		# assume equal_obj follows the old requirement for equalise_obj
		equal_obj = Wrapped_Equalise_Obj(equal_obj)
	end
	if isnothing(zero_margin)
		zero_margin = Default_Zero_Margin(equal_obj, obj2, zmin)
	elseif !applicable(zero_margin, obj, 1, zbounds...)
		@warn "Consider adapting `zero_margin` to the new method"
		zero_margin = Wrapped_Zero_D_j_Obj(zero_margin, fill(false, S))
	elseif zero_margin isa Default_Zero_Margin
		reinit!(zero_margin)
	end
	resize!(solver.intervalchoices, 1)
	solver.intervalchoices[1] = cdcp.x[1]
	resize!(solver.squeezing_indices, 1)
	solver.squeezing_indices[1] = 1
	empty!(solver.branching_indices)
	for m in solver.matcheds
		empty!(m)
	end
	for s in solver.singlesolvers
		s.obj = cdcp.obj
		reinit!(s; scdca=scdca)
	end
	cdcp.solver = SqueezingPolicy(scdca, solver.intervalchoices, solver.squeezing_indices, solver.branching_indices, zero_margin, equal_obj, solver.matcheds, solver.singlesolvers, obj2, solver.nobranching, solver.maxfcall)
	return cdcp
end