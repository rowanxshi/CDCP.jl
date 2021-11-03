import Test, Random, CDCP

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
	functions = ((J, z) -> π(J, z, π_params), pair -> equalise_π(pair, π_params), (j, J) -> zero_D_j_π(j, J, π_params))
	return functions, π_params
end

function coincide(C::Int, scdca::Bool; seed::Int = 10)
	(functions, ) = initiate(C, scdca, seed = seed)
	
	J = falses(C)
	Vs = CDCP.containers(C);
	
	(cutoffs, policies) = CDCP.policy(C, scdca, functions...)
	
	z_wrong = nothing
	coincide = all(0.01:0.1:50) do z
		CDCP.naive!(J, J -> first(functions)(J, z))	
		CDCP.solve!(Vs, scdca, J -> first(functions)(J, z))

		interval = searchsortedfirst(cutoffs, z)-1
		
		match = (J == first(Vs) == policies[interval])
		!match && (z_wrong = z)
		match
	end
	
	(coincide, z_wrong)
end


Test.@testset begin

for seed in 10:10:100
	Test.@test first(coincide(5, true, seed = seed))
	Test.@test first(coincide(5, false, seed = seed))
end

end