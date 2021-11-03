function containers((sub, sup, aux)::NTuple{3, AbstractVector{Bool}})
	working = [(sub, sup, aux); ]
	converged = similar(working)
	
	(working, converged)
end

function D_j(π::Function)
	function D_j_π(j::Integer, J::A) where A <: AbstractArray{Bool}
		bool_j = J[j]
		J[j] = true
		marg = π(J)
		J[j] = false
		marg -= π(J)
		J[j] = bool_j
		return marg
	end
end
