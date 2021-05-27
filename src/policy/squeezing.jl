function update!(int::interval, j::Int, zero_D_j_π, scdca, memo, emptyset)
	# above which j is in: worse case
	worst = scdca ? int.sup : int.sub
	z_in = get!(memo, (j, worst), zero_D_j_π(j, worst))
	
	# below which j is out: best case
	best = scdca ? int.sub : int.sup
	z_out = get!(memo, (j, best), zero_D_j_π(j, best))
	
	(z_in, z_out)
	
	# one new interval
	@inbounds if z_in ≤ int.l
		sub_new = int.sub[:]
		sub_new[j] = true
		return (interval(sub_new, int.sup, emptyset, int.l, int.r), )
	elseif z_out ≥ int.r
		sup_new = int.sup[:]
		sup_new[j] = false
		return (interval(int.sub, sup_new, emptyset, int.l, int.r), )
	end
	
	@inbounds begin
		aux_new = int.aux[:]
		aux_new[j] = true
	end
	# two new intervals
	@inbounds if (z_out ≤ int.l) && (int.l < z_in < int.r)
		sub_new = int.sub[:]
		sub_new[j] = true
		return interval(int.sub, int.sup, aux_new, int.l, z_in), interval(sub_new, int.sup, emptyset, z_in, int.r)
	elseif (z_in ≥ int.r) && (int.l < z_out < int.r)
		sup_new = int.sup[:]
		sup_new[j] = false
		return interval(int.sub, sup_new, emptyset, int.l, z_out), interval(int.sub, int.sup, aux_new, z_out, int.r)
	elseif int.l < z_in == z_out < int.r
		sub_new = int.sub[:]
		sub_new[j] = true
		sup_new = int.sup[:]
		sup_new[j] = false
		return interval(int.sub, sup_new, emptyset, int.l, z_out), interval(sub_new, int.sup, emptyset, z_in, int.r)
	end

	# three new intervals
	@inbounds if int.l < z_out < z_in < int.r
		sub_new = int.sub[:]
		sub_new[j] = true
		sup_new = int.sup[:]
		sup_new[j] = false
		return interval(int.sub, sup_new, emptyset, int.l, z_out), interval(int.sub, int.sup, aux_new, z_out, z_in), interval(sub_new, int.sup, emptyset, z_in, int.r)
	end	
	
	# otherwise, no update
	return (interval(int.sub, int.sup, aux_new, int.l, int.r), )
end
function converge!((working, converged)::NTuple{2, Vector{interval}}, zero_D_j_π, scdca, memo, emptyset)
	while !isempty(working)
		int = pop!(working)
		next = next_undetermined(int)
		if iszero(next)
			push!(converged, int)
		else
			append!(working, update!(int, next, zero_D_j_π, scdca, memo, emptyset))
		end
	end
	
	return converged
end
