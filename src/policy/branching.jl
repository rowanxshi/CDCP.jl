function converge_branches!(working, converged, done; cdcp...)
	while !isempty(converged)
		int = pop!(converged)
		if isequal(int.sub, int.sup)
			push!(done, int)
		else
			append!(working, branch(int))
			converge!(working, converged; cdcp...)
		end
	end
	sort!(done, lt=int_isless)
	
	return working, converged, done
end
