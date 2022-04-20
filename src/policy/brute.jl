function brute!(policy, working, converged, done; cdcp...)
	cutoffs, policies = policy
	for int in done
		push!(cutoffs, int.l); push!(cutoffs, int.r)
	end
	unique!(sort!(cutoffs))
	
	empty!(policies)
	resize!(policies, length(cutoffs) - 1)
	fill!(policies, nothing)
	
	# go interval by interval
	while any(isnothing, policies)
		k = findfirst(isnothing, policies)
		
		# gather relevant decision sets into converged
		empty!(converged)
		for option in done
			if (option.l ≤ cutoffs[k]) && (option.r ≥ cutoffs[k+1])
				push!(converged, option)
			end
		end
##		
		if isone(length(converged))
			policies[k] = first(converged).sub
		else
#			track the current subinterval's progress with working
			empty!(working)
			active_strats = trues(length(converged))
			
			# subinterval is an interval struct, but the BitVectors track which strategies (in converged) are still in consideration
			push!(working, interval(active_strats, cutoffs[k], cutoffs[k+1]))
			converge_brute!(policy, working, converged; cdcp...)
		end
	end
#
	paste_adjacent!(policy)
end

converge_brute!(policy, working, converged; cdcp...) = while !isempty(working)
	subint = last(working)
	if sum(subint.sub) == 1
		i_option = findfirst(subint.sub)
		patch!(policy, interval(converged[i_option].sub, subint.l, subint.r))
		pop!(working)
		continue
	end

	if !simple_filter!(subint, converged; cdcp...)
		i_J1 = findfirst(subint.sub)
		i_J2 = findlast(subint.sub)
		J1 = converged[i_J1].sub
		J2 = converged[i_J2].sub
		i_pair = (i_J1, i_J2)
		pair = (J1, J2)

		z_equal = cdcp[:equalise_obj](pair, subint.l, subint.r)

		if isnothing(z_equal) || (z_equal ≤ subint.l) || (z_equal ≥ subint.r)
			simple_filter!(subint, converged, (subint.l + subint.r)/2, i_J1, i_J2; cdcp...)
		else
			append!(working, brute_branch!(pop!(working), converged, z_equal, i_J1, i_J2; cdcp...))
		end
	end
end

function brute_branch!(subint::interval, converged, z_equal, i_pair...; cdcp...)
	subint_left = interval(subint.sub, subint.l, z_equal)
	subint_right = interval(copy(subint.sub), z_equal, subint.r)

	simple_filter!(subint_left, converged, subint.l, i_pair...; cdcp...)
	simple_filter!(subint_right, converged, subint.r, i_pair...; cdcp...)

	(subint_left, subint_right)
end

function simple_filter!(subint::interval, converged, z, i_pair...; cdcp...)
	diff = cdcp[:obj](converged[first(i_pair)].sub, z) - cdcp[:obj](converged[last(i_pair)].sub, z)
	J1_worse = signbit(diff)
	subint.sub[J1_worse ? first(i_pair) : last(i_pair)] = false
	return true
end

@inbounds function simple_filter!(subint::interval, converged; cdcp...)
	options = sum(subint.sub)
	
	max_left = maximum(enumerate(converged)) do (i, option)
		!subint.sub[i] && return -Inf
		
		cdcp[:obj](option.sub, subint.l)
	end
	
	for (i, option) in enumerate(converged)
		!subint.sub[i] && continue
		
		(cdcp[:obj](option.sub, subint.r) < max_left) && (subint.sub[i] = false)
	end
	
	return options > sum(subint.sub)
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

function patch!((cutoffs, policies), int::interval)
	l_in = insorted(int.l, cutoffs)
	r_in = insorted(int.r, cutoffs)
	
	if (!l_in && !r_in)
		k = searchsortedfirst(cutoffs, int.l)
		cutoffs[k] < int.r && error()
		!isnothing(policies[k-1]) && error()
		
		insert!(cutoffs, k, int.r)
		insert!(policies, k, nothing)
		
		insert!(cutoffs, k, int.l)
		insert!(policies, k, int.sub)
		
		return k
	end
	
	if (l_in && r_in)
		k = searchsortedfirst(cutoffs, int.l)
		cutoffs[k+1] != int.r && error()
		!isnothing(policies[k]) && error("patching policy fn with $int but it already has strategy $(policies[k])")
		
		policies[k] = int.sub
		return k
	end
	
	l_in && begin
		k = searchsortedfirst(cutoffs, int.l)
		!isnothing(policies[k]) && error()
		cutoffs[k+1] < int.r && error()
		
		insert!(cutoffs, k+1, int.r)
		insert!(policies,  k, int.sub)
		
		return k
	end
	
	r_in && begin
		k = searchsortedfirst(cutoffs, int.l)
		cutoffs[k-1] > int.l && error()
		!isnothing(policies[k-1]) && error()
		
		insert!(cutoffs, k, int.l)
		insert!(policies, k, int.sub)
		
		return k
	end
end