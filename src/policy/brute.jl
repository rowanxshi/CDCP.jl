function brute!(policy_fn, (working, converged, done), π, equalise_π, scdca::Bool, memo...)
	endpoints = Vector{typeof(Inf)}(undef, 0)
	
	for int in done
		push!(endpoints, int.l)
		push!(endpoints, int.r)
	end
	unique!(sort!(endpoints))
	
	# go interval by interval
	@inbounds for n in 1:(length(endpoints)-1)
		# gather relevant decision sets into converged
		empty!(converged)
		for option in done
			if (option.l ≤ endpoints[n]) && (option.r ≥ endpoints[n+1])
				push!(converged, option)
			end
		end
		
		if isone(length(converged))
			patch!(policy_fn, pop!(converged))
		else
			# track the current subinterval's progress with working
			empty!(working)
			active_strategies = trues(length(converged))
			
			# subinterval is an interval struct, but the BitVectors track which strategies (in converged) are still in consideration
			push!(working, interval(active_strategies, endpoints[n], endpoints[n+1]))
			converge_brute!((policy_fn, working), converged, π, equalise_π, memo...)
		end
	end

	paste_adjacent!(policy_fn)
end

converge_brute!((policy_fn, working), converged, π, equalise_π, memo...) = while !isempty(working)
	subinterval = last(working)
	if sum(subinterval.sub) == 1
		i_option = findfirst(subinterval.sub)
		patch!(policy_fn, interval(converged[i_option].sub, subinterval.l, subinterval.r))
		pop!(working)
		continue
	end

	if !simple_filter!(subinterval, converged, π)
		i_J1 = findfirst(subinterval.sub)
		i_J2 = findlast(subinterval.sub)
		J1 = converged[i_J1].sub
		J2 = converged[i_J2].sub
		i_pair = (i_J1, i_J2)
		set_pair = (J1, J2)

		z_equal = if isempty(memo) # memoised if a dictionary is provided
			equalise_π(set_pair)
		else
			get!(first(memo), set_pair, equalise_π(set_pair))
		end
		left = subinterval.l
		right = subinterval.r
		if z_equal ≤ left
			simple_filter!(subinterval, i_pair, set_pair, π, left)
		elseif z_equal ≥ right
			simple_filter!(subinterval, i_pair, set_pair, π, right)
		else
			append!(working, brute_branch!(pop!(working), i_pair, set_pair, π, z_equal))
		end
	end
end

function brute_branch!(subinterval::interval, i_pair::NTuple{2, Int}, set_pair::NTuple{2, BitVector}, π, z_equal)
	subinterval_left = interval(subinterval.sub, subinterval.l, z_equal)
	subinterval_right = interval(copy(subinterval.sub), z_equal, subinterval.r)

	J1_better = simple_filter!(subinterval_left, i_pair, set_pair, π, (subinterval.l + z_equal)/2)
	simple_filter!(subinterval_right, i_pair, !J1_better)

	(subinterval_left, subinterval_right)
end

function simple_filter!(subinterval::interval, i_pair::NTuple{2, Int}, J1_better::Bool)
	subinterval.sub[first(i_pair)] = J1_better
	subinterval.sub[last(i_pair)] = !J1_better

	J1_better
end
function simple_filter!(subinterval::interval, i_pair::NTuple{2, Int}, set_pair::NTuple{2, BitVector}, π, z)
	J1_better = π(first(set_pair), z) > π(last(set_pair), z)
	simple_filter!(subinterval, i_pair, J1_better)

	J1_better
end

simple_filter!(subinterval::interval, converged::Vector{interval}, π) = @inbounds begin
	options = sum(subinterval.sub)
	
	max_left = maximum(enumerate(converged)) do (i, option)
		!subinterval.sub[i] && return -Inf
		
		π(option.sub, subinterval.l)
	end
	
	for (i, option) in enumerate(converged)
		!subinterval.sub[i] && continue
		
		(π(option.sub, subinterval.r) < max_left) && (subinterval.sub[i] = false)
	end
	
	return options > sum(subinterval.sub)
end

function paste_adjacent!((cutoffs, policies))
	i = 1
	for n in 1:(length(policies)-1)
		changed = false
		policies[i] == policies[i+1] && begin
			deleteat!(cutoffs, i+1)
			deleteat!(policies, i+1)
			changed = true
		end
		i += changed ? 0 : 1
	end
		
	cutoffs, policies
end
function patch!((cutoffs, policies), interval::interval)
	left_in = insorted(interval.l, cutoffs)
	right_in = insorted(interval.r, cutoffs)
	
	!any((left_in, right_in)) && begin
		place = searchsortedfirst(cutoffs, interval.l)
		cutoffs[place] < interval.r && error()
		!isnothing(policies[place-1]) && error()
		
		insert!(cutoffs, place, interval.r)
		insert!(policies, place, nothing)
		
		insert!(cutoffs, place, interval.l)
		insert!(policies, place, interval.sub)
		
		return place
	end
	
	all((left_in, right_in)) && begin
		place = searchsortedfirst(cutoffs, interval.l)
		cutoffs[place+1] != interval.r && error()
		!isnothing(policies[place]) && error()
		
		policies[place] = interval.sub
		
		return place
	end
	
	left_in && begin
		place = searchsortedfirst(cutoffs, interval.l)
		cutoffs[place+1] < interval.r && error()
		!isnothing(policies[place]) && error()
		
		insert!(cutoffs, place+1, interval.r)
		insert!(policies,  place, interval.sub)
		
		return place
	end
	
	right_in && begin
		place = searchsortedfirst(cutoffs, interval.l)
		cutoffs[place-1] > interval.l && error()
		!isnothing(policies[place-1]) && error()
		
		insert!(cutoffs, place, interval.l)
		insert!(policies, place, interval.sub)
		
		return place
	end
end
