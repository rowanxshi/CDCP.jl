function containers(Vs)
	working = [Vs; ]
	converged = similar(working)
	(working, converged)
end
function containers(C::Integer)
	sub = falses(C)
	sup = trues(C)
	aux = falses(C)
	return (sub, sup, aux)
end

function invert_triplet(itemstates::SVector{S,ItemState}, sub, sup, aux) where S
	checklength(itemstates, sub, sup, aux)
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(ItemState( @inbounds(sub[$i]), @inbounds(sup[$i]), @inbounds(aux[$i]))))
		end
		return :(SVector{S,ItemState}($ex))
	else
		return SVector{S,ItemState}(ntuple(i->ItemState(sub[i], sup[i], aux[i]), S))
	end
end
function ItemState(sub::Bool, sup::Bool, isaux::Bool)
	if isaux
		return aux
	elseif sub != sup
		return undetermined
	elseif sub
		return included
	else
		return excluded
	end
end
function checklength(::SVector{S,ItemState}, sub, sup, aux) where S
	length(sub) == length(sup) == length(aux) == S || throw(ArgumentError("length of sub, sup, aux are inconsistent with C"))
end

function invert_state!(sub::AbstractVector{Bool}, sup::AbstractVector{Bool}, isaux::AbstractVector{Bool}, itemstates::SVector{S,ItemState}) where S
	@inbounds for i in 1:S
		itemstate = itemstates[i]
		if itemstate == aux
			isaux[i] = true
		else
			isaux[i] = false
			if itemstate == undetermined
				sub[i] = false
				sup[i] = true
			elseif itemstate == included
				sub[i] = true
				sup[i] = true
			else
				sub[i] = false
				sup[i] = false
			end
		end
	end
end

function zero_D_j(equalise_obj, holder::AbstractVector{Bool})
	zero_D_j_obj(j::Integer, J::AbstractVector{Bool}, extras...) = let equalise_obj = equalise_obj, holder = holder
		J2 = copyto!(holder, J)
		J2 = setindex!(J2, !J[j], j)
		equalise_obj((J, J2), extras...)
	end
end

function restart!((working, converged, done), C)
	empty!.((working, converged, done))
	int = Interval(containers(C), -Inf, Inf)
	push!(working, int)
end

function invert_cutoffspolicies!(cdcp::CDCProblem{<: SqueezingPolicy}, containers)
	intervalchoices, squeezing = cdcp.solver.intervalchoices, cdcp.solver.squeezing_indices
	itemstates = first(intervalchoices).itemstates
	working = first(containers)
	if !isempty(working)
		resize!(intervalchoices, length(working))
		resize!(squeezing, length(working))
		for (i, interval) in enumerate(working)
			intervalchoices[i] = IntervalChoice(interval.l, interval.r, invert_triplet(itemstates, interval.sub, interval.sup, interval.aux))
			squeezing[i] = i
		end
	end
	# if cutoffspolicies !== nothing
	# 	if !initialized
	# 		cutoffs, policies = cutoffspolicies
	# 		resize!(intervalchoices, length(policies))
	# 		resize!(squeezing, length(policies))
	# 		for i in 1:length(policies)
	# 			intervalchoices[i] = IntervalChoice(cutoffs[i], cutoffs[i+1], policies[i])
	# 			squeezing[i] = i
	# 		end
	# 	end  
	# end
	cdcp
end
