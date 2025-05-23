function _obj(J::AbstractArray{Bool}, z::Real, (C, var, scdca)::Tuple{Int, Vector{<: Real}, Bool})
	δ = scdca ? 0.25 : 1.1
	f = range(0.1, length = C, step = 0.1)

	profits = z*sum(c -> J[c]*var[c], 1:C)^δ - sum(c -> J[c]*f[c], 1:C)
end
function _zero_D_j_obj(j::Integer, J::AbstractArray{Bool}, (C, var, scdca)::Tuple{Int, Vector{<: Real}, Bool})
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
function _equalise_obj((J1, J2), (C, var, scdca)::Tuple{Int, Vector{<: Real}, Bool})
	δ = scdca ? 0.25 : 1.1
	f = range(0.1, length = C, step = 0.1)

	z = sum(c -> J1[c]*var[c], 1:C)^δ - sum(c -> J2[c]*var[c], 1:C)^δ
	z = sum(c -> (J1[c] - J2[c])*f[c], 1:C)/z
end
function initiate(C::Int = 10, scdca::Bool = false; seed::Int = 10, z0::Real = 1)
	Random.seed!(seed)
	var = rand(C)
	params = (C, var, scdca)

	obj(J, z::Real = z0) = let params = params
		_obj(J, z, params)
	end
	equalise_obj(pair, extras...) = let params = params
		_equalise_obj(pair, params)
	end
	zero_D_j_obj(j, J, extras...) = let params = params
		_zero_D_j_obj(j, J, params)
	end

	cdcp = (; scdca, obj, equalise_obj, zero_D_j_obj)
	return cdcp, params
end

function coincide(C::Int, scdca::Bool = false; seed::Int = 10)
	(cdcp, ) = initiate(C, scdca, seed = seed)

	J = falses(C)
	Vs = CDCP.containers(C)
	
	(cutoffs, policies) = CDCP.policy(C; cdcp...)
	
	z_wrong = nothing
	coincide = all(0.01:0.1:50) do z
		CDCP.naive!(J; obj=(J->cdcp.obj(J, z)))
		CDCP.solve!(Vs; cdcp.scdca, obj=(J->cdcp.obj(J, z)))

		interval = searchsortedfirst(cutoffs, z)-1
		
		match = (J == first(Vs) == policies[interval])
		!match && (z_wrong = z)
		match
	end
	
	(coincide, z_wrong)
end

@testset "wrapper" begin
	C = 2
	cdcp = initiate(C)
	containers = (falses(C), trues(C), falses(C))
	@test CDCP.naive(C; cdcp[1].obj) isa AbstractVector{Bool}
	@test CDCP.naive(C; cdcp[1].obj, z = 2.0) isa AbstractVector{Bool}
	@test CDCP.naive!(falses(C); cdcp[1].obj) isa AbstractVector{Bool}
	@test CDCP.naive!(falses(C); cdcp[1].obj, z = 2.0) isa AbstractVector{Bool}
	@test solve(C; cdcp[1].scdca, cdcp[1].obj) isa AbstractVector{Bool}
	@test solve!(containers; cdcp[1].scdca, cdcp[1].obj) isa AbstractVector{Bool}
	cutoffs_strategies = policy(C; cdcp[1].scdca, cdcp[1].obj, cdcp[1].equalise_obj)
	@test cutoffs_strategies isa Tuple{<: Any, <: Any}
	@test cutoffs_strategies[1] isa AbstractVector{<: Real}
	@test cutoffs_strategies[2] isa AbstractVector{<: Union{Nothing, <: AbstractVector{Bool}}}
	for seed in 10:10:100
		@test first(coincide(10, true, seed = seed))
		@test first(coincide(10, false, seed = seed))
	end
end
