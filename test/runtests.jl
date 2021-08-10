import Test, Random, CDCP

function π(J::BitVector, z::Float64, (C, var, scdca)::Tuple{Int, Vector{Float64}, Bool})
    δ = scdca ? 0.25 : 1.1
    f = range(0.1, length = C, step = 0.1)
    
    profits = z*sum(c -> J[c]*var[c], 1:C)^δ - sum(c -> J[c]*f[c], 1:C)
end
function zero_D_j_π(j::Integer, J::AbstractVector{Bool}, (C, var, scdca)::Tuple{Int, Vector{Float64}, Bool})
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
function equalise_π((J1, J2)::Tuple{<: AbstractVector{Bool}, <: AbstractVector{Bool}}, (C, var, scdca)::Tuple{Int, Vector{Float64}, Bool})
    δ = scdca ? 0.25 : 1.1
    f = range(0.1, length = C, step = 0.1)
    
    z = sum(c -> J1[c]*var[c], 1:C)^δ - sum(c -> J2[c]*var[c], 1:C)^δ
    z = sum(c -> (J1[c] - J2[c])*f[c], 1:C)/z
end

Test.@testset begin

function test_single(C::Integer = 5, z::T = 1.; scdca::Bool = true) where T <: Real
	sub = falses(C)
	sup = trues(C)
	aux = falses(C)
	
	working = [(sub, sup, aux); ]
	converged = similar(working)
	Random.seed!(1)
	var = rand(C)
	
    π_params = (C, var, scdca)
    
	CDCP.solve!((sub, sup, aux), J -> π(J, z, π_params), scdca)
	CDCP.solve!((sub, sup, aux), J -> π(J, z, π_params), scdca, containers = (working, converged))

	CDCP.solve(C, J -> π(J, z, π_params), scdca)
	CDCP.solve(C, J -> π(J, z, π_params), scdca, containers = (working, converged))
end

Test.@test typeof(test_single(5)) <: AbstractVector{Bool}

function test_policy(C::Integer = 5, scdca::Bool = true)
	Random.seed!(1)
	var = rand(C)
	
    π_params = (C, var, scdca)
	
	CDCP.policy(C, (J, z) -> π(J, z, π_params), pair -> equalise_π(pair, π_params), scdca)
end

(cutoffs, policies) = test_policy(15)
Test.@test typeof(cutoffs) == Vector{Float64}
Test.@test typeof(policies) == Vector{BitVector}

end
