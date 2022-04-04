function _containers(Vs)
	working = [Vs; ]
	converged = similar(working)
	
	(working, converged)
end

function D_j(obj)
	D_j_obj(j::Integer, J::AbstractArray{Bool}, z::Real = 1) = let obj = obj
		bool_j = J[j]
		J = setindex!(J, true, j)
		marg = obj(J, z)
		J = setindex!(J, false, j)
		marg -= obj(J, z)
		J = setindex!(J, bool_j, j)
		return marg
	end
end
