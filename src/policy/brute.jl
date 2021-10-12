function brute!((working, converged, done), π, equalise_π, scdca::Bool)
	endpoints = Vector{typeof(Inf)}(undef, 0)
	within_endpoints = similar(endpoints)
	memo = Dict{NTuple{2, BitVector}, typeof(Inf)}()
	
	for int in done
		push!(endpoints, int.l)
		push!(endpoints, int.r)
	end
	unique!(sort!(endpoints))
	
	# go interval by interval
	@inbounds for n in 1:(length(endpoints)-1)
		# gather relevant decision sets into working
		empty!(working)
		for option in done
			if (option.l ≤ endpoints[n]) && (option.r ≥ endpoints[n+1])
				push!(working, option)
			end
		end
		if isone(length(working))
			# store final policy function intervals in converged
			push!(converged, pop!(working))
		else
			empty!(within_endpoints)
			# gather points where potentially policy changes
			push!(within_endpoints, endpoints[n])
			push!(within_endpoints, endpoints[n+1])
			for i in 1:length(working), j in (i+1):length(working)
				pair = (working[i].sub, working[j].sub)
				z_equal = get!(memo, pair, equalise_π(pair))
				(endpoints[n] < z_equal < endpoints[n+1]) && push!(within_endpoints, z_equal)
			end
			sort!(within_endpoints)
			
			# for each subinterval, evaluate at midpoint and find best
			for nn in 1:(length(within_endpoints)-1)
				z = if isinf(within_endpoints[nn])
					within_endpoints[nn+1] - one(within_endpoints[nn+1])
				elseif isinf(within_endpoints[nn+1])
					within_endpoints[nn] + one(within_endpoints[nn])
				else
					(within_endpoints[nn] + within_endpoints[nn+1])/2
				end
				max_π = -Inf
				i_max = 0
				for (i_option, option) in enumerate(working)
					option_π = π(option.sub, z)
					if option_π > max_π
						max_π = option_π
						i_max = i_option
					end
				end
				push!(converged, interval(working[i_max].sub, within_endpoints[nn], within_endpoints[nn+1]))
			end
		end
	end
	
	# paste together any adjacent intervals with the same policy
	concatenate!(sort!(converged, lt=int_isless), endpoints)
end
function concatenate!(converged, endpoints)
	empty!(endpoints)
	policies = Vector{BitVector}(undef, 0)
	push!(endpoints, first(converged).l)
	push!(policies, first(converged).sub)
	
	for int in converged
		isequal(int.sub, last(policies)) && continue
		push!(endpoints, int.l)
		push!(policies, int.sub)
	end
	
	push!(endpoints, last(converged).r)
	
	return endpoints, policies
end
