import Random

function π(J::AbstractArray{Bool}, z::Float64, (C, var, scdca)::Tuple{Int, Vector{Float64}, Bool})
	δ = scdca ? 0.25 : 1.1
	f = range(0.1, length = C, step = 0.1)
	
	profits = z*sum(c -> J[c]*var[c], 1:C)^δ - sum(c -> J[c]*f[c], 1:C)
end
function zero_D_j_π(j::Integer, J::AbstractArray{Bool}, (C, var, scdca)::Tuple{Int, Vector{Float64}, Bool})
	δ = scdca ? 0.25 : 1.1
	f = range(0.1, length = C, step = 0.1)
	
	bool_j = J[j]
	J[j] = true
	z = sum(c -> J[c]*var[c], 1:C)^δ
	J[j] = false
	z -= sum(c -> J[c]*var[c], 1:C)^δ
	J[j] = bool_j
	z = f[j]/z
end
function equalise_π((J1, J2)::NTuple{2, AbstractVector{Bool}}, (C, var, scdca)::Tuple{Int, Vector{Float64}, Bool})
	δ = scdca ? 0.25 : 1.1
	f = range(0.1, length = C, step = 0.1)
	
	z = sum(c -> J1[c]*var[c], 1:C)^δ - sum(c -> J2[c]*var[c], 1:C)^δ
	z = sum(c -> (J1[c] - J2[c])*f[c], 1:C)/z
end
function initiate(C::Int = 10, scdca::Bool = false; seed::Int = 10)
	Random.seed!(seed)
	var = rand(C)
	π_params = (C, var, scdca)
	functions = (π = ((J, z) -> π(J, z, π_params)), equalise_π = (pair -> equalise_π(pair, π_params)), zero_D_j_π = ((j, J) -> zero_D_j_π(j, J, π_params)))
	return functions, π_params
end
