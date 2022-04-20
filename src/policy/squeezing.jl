function update!(int::interval, j::Integer; cdcp...)
	# above which j is in: worse case
	worst = cdcp[:scdca] ? int.sup : int.sub
	
#	include on the whole interval: if D_j(worst, l) ≥ 0
	if isfinite(int.l) && (cdcp[:D_j_obj](j, worst, int.l) ≥ 0)
		j_in = setindex!(copy(int.sub), true, j)
		return (interval(j_in, int.sup, cdcp[:emptyset], int.l, int.r), )
	end

#	no conclusion on inclusion for the whole interval: z_in ≥ r
	if isfinite(int.r) && (cdcp[:D_j_obj](j, worst, int.r) ≤ 0) 
		return exclude_update!(int, j; cdcp...)
	end
	
#	otherwise, use zero_D_j_obj
	z_in = cdcp[:zero_D_j_obj](j, worst, int.l, int.r)
	if isnothing(z_in)
		add_j = cdcp[:D_j_obj](j, worse, middle(int)) ≥ 0
		!add_j && return exclude_update!(int, j; cdcp...)
		j_in = setindex!(copy(int.sub), true, j)
		return (interval(j_in, int.sup, cdcp[:emptyset], int.l, int.r), )
	end
	# include for some: l < z_in < r
	j_in = setindex!(copy(int.sub), true, j)
	return (exclude_update!(interval(int.sub, int.sup, int.aux, int.l, z_in), j; cdcp...)..., interval(j_in, int.sup, cdcp[:emptyset], z_in, int.r))
end

function exclude_update!(int::interval, j::Integer; cdcp...)
	best = cdcp[:scdca] ? int.sub : int.sup
	
#	exclude on the whole interval: z_out ≤ r
	if isfinite(int.r) && (cdcp[:D_j_obj](j, best, int.r) ≤ 0)
		j_out = setindex!(copy(int.sup), false, j)
		return (interval(int.sub, j_out, cdcp[:emptyset], int.l, int.r), )
	end
	
	j_aux = setindex!(copy(int.aux), true, j)
#	no conclusion on exclusion for the whole interval: z_out ≤ l
	if isfinite(int.l) && (cdcp[:D_j_obj](j, best, int.l) ≥ 0)
		return (interval(int.sub, int.sup, j_aux, int.l, int.r), )
	end
	
#	otherwise, use zero_D_j_obj
	z_out = cdcp[:zero_D_j_obj](j, best, int.l, int.r)
	if isnothing(z_out)
		drop_j = cdcp[:D_j_obj](j, best, middle(int)) ≤ 0
		!drop_j && return (interval(int.sub, int.sup, j_aux, int.l, int.r), )
		j_out = setindex!(copy(int.sup), false, j)
		return (interval(int.sub, j_out, cdcp[:emptyset], int.l, int.r), )
	end
	# exclude for some, then unsure for the rest
	j_out = setindex!(copy(int.sup), false, j)
	return (interval(int.sub, j_out, cdcp[:emptyset], int.l, z_out), interval(int.sub, int.sup, j_aux, z_out, int.r))
end

function converge!(working, converged; cdcp...)
	while !isempty(working)
		int = pop!(working)
		j = next_undetermined(int)
		if iszero(j)
			push!(converged, int)
		else
			append!(working, update!(int, j; cdcp...))
		end
	end
	
	return converged
end
