function update!(int::Interval, j::Integer; cdcp...)
	# above which j is in: worse case
	worst = cdcp[:scdca] ? int.sup : int.sub
	
#	include on the whole interval: if D_j(worst, l) ≥ 0
	if isfinite(int.l) && !signbit(cdcp[:D_j_obj](j, worst, int.l))
		j_in = setindex!(copy(int.sub), true, j)
		return (Interval(j_in, int.sup, cdcp[:emptyset], int.l, int.r), )
	end

#	no conclusion on inclusion for the whole interval: z_in ≥ r
	Djworst = cdcp[:D_j_obj](j, worst, int.r)
	if isfinite(int.r) && (signbit(Djworst) || iszero(Djworst)) 
		return exclude_update!(int, j; cdcp...)
	end
	
#	otherwise, use zero_D_j_obj
	z_in = cdcp[:zero_D_j_obj](j, worst, int.l, int.r)
	if isnothing(z_in)
		add_j = !signbit(cdcp[:D_j_obj](j, worst, middle(int)))
		!add_j && return exclude_update!(int, j; cdcp...)
		j_in = setindex!(copy(int.sub), true, j)
		return (Interval(j_in, int.sup, cdcp[:emptyset], int.l, int.r), )
	end
	# include for some: l < z_in < r
	j_in = setindex!(copy(int.sub), true, j)
	return (exclude_update!(Interval(int.sub, int.sup, int.aux, int.l, z_in), j; cdcp...)..., Interval(j_in, int.sup, cdcp[:emptyset], z_in, int.r))
end

function exclude_update!(int::Interval, j::Integer; cdcp...)
	best = cdcp[:scdca] ? int.sub : int.sup
	
#	exclude on the whole interval: z_out ≤ r
	Djbest = cdcp[:D_j_obj](j, best, int.r)
	if isfinite(int.r) && (signbit(Djbest) || iszero(Djbest))
		j_out = setindex!(copy(int.sup), false, j)
		return (Interval(int.sub, j_out, cdcp[:emptyset], int.l, int.r), )
	end
	
	j_aux = setindex!(copy(int.aux), true, j)
#	no conclusion on exclusion for the whole interval: z_out ≤ l
	if isfinite(int.l) && !signbit(cdcp[:D_j_obj](j, best, int.l))
		return (Interval(int.sub, int.sup, j_aux, int.l, int.r), )
	end
	
#	otherwise, use zero_D_j_obj
	z_out = cdcp[:zero_D_j_obj](j, best, int.l, int.r)
	if isnothing(z_out)
		Djbest = cdcp[:D_j_obj](j, best, middle(int))
		drop_j = signbit(Djbest) || iszero(Djbest)
		!drop_j && return (Interval(int.sub, int.sup, j_aux, int.l, int.r), )
		j_out = setindex!(copy(int.sup), false, j)
		return (Interval(int.sub, j_out, cdcp[:emptyset], int.l, int.r), )
	end
	# exclude for some, then unsure for the rest
	j_out = setindex!(copy(int.sup), false, j)
	return (Interval(int.sub, j_out, cdcp[:emptyset], int.l, z_out), Interval(int.sub, int.sup, j_aux, z_out, int.r))
end

function converge!(working, converged; cdcp...)
	while !isempty(working)
		int = pop!(working)
		isequal(int.l, int.r) && continue
		j = next_undetermined(int)
		if iszero(j)
			push!(converged, int)
		else
			append!(working, update!(int, j; cdcp...))
		end
	end
	
	return converged
end
