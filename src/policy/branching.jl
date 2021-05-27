function branch(int::interval)
	((a,b,c), (d,e,f)) = branch((int.sub, int.sup, int.aux))
	return interval(a, b, c, int.l, int.r), interval(d, e, f, int.l, int.r)
end
function converge_branches!((working, converged, done)::NTuple{3, Vector{interval}}, zero_D_j_π, scdca, memo, emptyset)
	while !isempty(converged)
		int = pop!(converged)
		if isequal(int.sub, int.sup)
			push!(done, int)
		else
			append!(working, branch(int))
			converge!((working, converged), zero_D_j_π, scdca, memo, emptyset)
		end
	end
	sort!(done, lt=int_isless)
	
	return working, converged, done
end
