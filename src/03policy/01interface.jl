function solve!(cdcp::CDCProblem{<:SqueezingPolicy}; restart::Bool=false, obj=cdcp.obj, zero_margin=cdcp.solver.zero_margin, equal_obj=cdcp.solver.equal_obj, scdca=cdcp.solver.scdca)
	restart && (cdcp = reinit!(cdcp; obj, zero_margin, equal_obj, scdca))
	@debug "squeezing"
	squeeze!(cdcp)
	if cdcp.state == maxfcall_reached
		@warn "maxfcall is reached before convergence"
		return cdcp
	end
	cdcp.solver.skiprefinement || (@debug("refining"); cdcp.state = refine!(cdcp))
	concat!(cdcp)
	cdcp.solver.skiprefinement || (cdcp.state = success)
	return cdcp
end

function init_solverx(::Type{<:SqueezingPolicy}, obj, scdca::Bool, equal_obj, zbounds::Tuple{Z,Z}=(-Inf, Inf); zero_margin=nothing, policy0::Policy=Policy(obj, zbounds), skiprefinement::Bool=false, singlekw=NamedTuple(), maxfcall=1_000_000_000, kwargs...) where Z
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
	singlecdcp = init(Squeezing, obj, S, scdca; z=zero(Z), singlekw...)
	return SqueezingPolicy(scdca, intervalchoices, collect(eachindex(intervalchoices)), Int[], zero_margin, equal_obj, singlecdcp, obj2, skiprefinement, maxfcall), policy0
end

function reinit!(cdcp::CDCProblem{<:SqueezingPolicy}; obj=cdcp.obj, fcall=true, solverkw...)
	S = length(cdcp.x.itemstates_s[1])
	cdcp = reinit!(cdcp, S; obj, fcall)
	cdcp.obj = reinit!(obj, S)
	cdcp = reinit_solverx!(cdcp; obj, solverkw...)
	return cdcp
end

function reinit_solverx!(cdcp::CDCProblem{<:SqueezingPolicy}; obj=cdcp.obj, fcall=true, solverkw...)
	cdcp.x = reinit!(cdcp.x, obj)
	zbounds = (cdcp.x.cutoffs[1], cdcp.x.zright)
	cdcp.solver = reinit!(cdcp.solver, cdcp.x[1], zbounds; obj, fcall, solverkw...)
	cdcp
end
function reinit!(policy::Policy, obj)
	zmin = policy.cutoffs[1]
	resize!(policy.cutoffs, 1)
	policy.cutoffs[1] = zmin
	resize!(policy.itemstates_s, 1)
	policy.itemstates_s[1] = allundetermined(obj)
	policy
end
function reinit!(solver::T, intervalchoice::IntervalChoice, zbounds; obj, fcall=false, zero_margin=solver.zero_margin, equal_obj=solver.equal_obj, scdca=solver.scdca) where {T<:SqueezingPolicy}
	resize!(solver.intervalchoices, 1)
	solver.intervalchoices[1] = intervalchoice
	resize!(solver.squeezing_indices, 1)
	solver.squeezing_indices[1] = 1
	empty!(solver.branching_indices)
	obj2 = deepcopy(obj)
	zmin = first(zbounds)
	reinit!(solver.singlecdcp; obj, fcall, scdca)
	# harmonize user defined functions
	if isnothing(zero_margin)
		zero_margin = Default_Zero_Margin(equal_obj, obj2, zmin)
	elseif !applicable(zero_margin, obj, 1, zbounds...)
		@warn "Consider adapting `zero_margin` to the new method"
		zero_margin = Wrapped_Zero_D_j_Obj(zero_margin, fill(false, S))
	elseif zero_margin isa Default_Zero_Margin
		reinit!(zero_margin)
	end
	if !applicable(equal_obj, obj, obj2, zbounds...)
		# assume equal_obj follows the old requirement for equalise_obj
		equal_obj = Wrapped_Equalise_Obj(equal_obj)
	end
	solver = SqueezingPolicy(scdca, solver.intervalchoices, solver.squeezing_indices, solver.branching_indices, zero_margin, equal_obj, solver.singlecdcp, obj2, solver.skiprefinement, solver.maxfcall)
end
