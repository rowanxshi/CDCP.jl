import Test, Random, CDCP

function π(J::V, z::Float64, (C, var, scdca)::Tuple{Int, Vector{Float64}, Bool}) where V <: AbstractArray{Bool}
	δ = scdca ? 0.25 : 1.1
	f = range(0.1, length = C, step = 0.1)
	
	profits = z*sum(c -> J[c]*var[c], 1:C)^δ - sum(c -> J[c]*f[c], 1:C)
end
function zero_D_j_π(j::Integer, J::V, (C, var, scdca)::Tuple{Int, Vector{Float64}, Bool}) where V <: AbstractVector{Bool}
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

function coincide(C::Int, scdca::Bool; seed::Int = 10)
	Random.seed!(seed)
	var = rand(C)
	π_params = (C, var, scdca)	
	π_test(J, z) = π(J, z, π_params)
	
	J = falses(C)
	sub = falses(C); sup = trues(C); aux = falses(C);
	
	(cutoffs, policies) = CDCP.policy(C, π_test, pair -> equalise_π(pair, π_params), scdca)
	
	z_wrong = nothing
	coincide = all(0.01:0.1:50) do z
		CDCP.naive!(J, J -> π_test(J, z))	
		CDCP.solve!((sub, sup, aux), J -> π_test(J, z), scdca)
		
		interval = searchsortedfirst(cutoffs, z)-1
		
		match = (J == sup == policies[interval])
		!match && (z_wrong = z)
		match
	end
	
	(coincide, z_wrong)
end

#println("Test solution methods coincide with brute force")
#Test.@test first(coincide(C, scdca, seed = seed))

function test_single(C::Integer = 5, z::T = 1.; scdca::Bool = true) where T <: Real
	Random.seed!(10)
	var = rand(C)
	π_params = (C, var, scdca)
	
	sub = falses(C); sup = trues(C); aux = falses(C)
	working = [(sub, sup, aux); ]; converged = similar(working)
	
	a = CDCP.solve!((sub, sup, aux), J -> π(J, z, π_params), scdca)
	b = CDCP.solve!((sub, sup, aux), J -> π(J, z, π_params), scdca, containers = (working, converged))
	c = CDCP.solve(C, J -> π(J, z, π_params), scdca)
	d = CDCP.solve(C, J -> π(J, z, π_params), scdca, containers = (working, converged))
	
	return a, b, c, d
end

println("Test all single-agent methods")
(a, b, c, d) = test_single(5)
Test.@test a == b
Test.@test a == c
Test.@test a == d
Test.@test typeof(a) <: AbstractVector{Bool}

function test_policy(C::Integer = 5, scdca::Bool = true)
	Random.seed!(10)
	var = rand(C)
	π_params = (C, var, scdca)
	memo = Dict{NTuple{2, BitVector}, Float64}()
	
	e, f = CDCP.policy(C, (J, z) -> π(J, z, π_params), (j, J) -> zero_D_j_π(j, J, π_params), pair -> equalise_π(pair, π_params), scdca)
	CDCP.policy(C, (J, z) -> π(J, z, π_params), (j, J) -> zero_D_j_π(j, J, π_params), pair -> equalise_π(pair, π_params), scdca, memo...)
	g, h = CDCP.policy(C, (J, z) -> π(J, z, π_params), pair -> equalise_π(pair, π_params), scdca)
	CDCP.policy(C, (J, z) -> π(J, z, π_params), pair -> equalise_π(pair, π_params), scdca, memo...)
	
	return e, f, g, h
end

println("Test all policy function methods")
(e, f, g, h) = test_policy(15)
Test.@test f == h
Test.@test e ≈ g

Test.@test typeof(e) == Vector{Float64}
Test.@test typeof(f) == Vector{Union{BitVector, Nothing}}
Test.@test length(first(f)) == 15

end
