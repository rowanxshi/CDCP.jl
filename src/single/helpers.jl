function containers((sub, sup, aux)::NTuple{3, AbstractVector{Bool}})
	working = [(sub, sup, aux); ]
	converged = similar(working)
	
	(working, converged)
end

function D_j(π::Function)
	D_j_π(j::Integer, J::AbstractArray{Bool}) = let π = π
		bool_j = J[j]
		J[j] = true
		marg = π(J)
		J[j] = false
		marg -= π(J)
		J[j] = bool_j
		return marg
	end
end
