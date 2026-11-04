function solve!(cdcp::CDCProblem{<:Naive}; restart::Bool=true)
	(; obj) = cdcp
	(; ids, z, incremental) = cdcp.solver
	S = length(cdcp.x)
	S < Sys.WORD_SIZE - 1 || throw(ArgumentError("exhaustion over $S items is not feasible"))
	restart && resize!(ids, 1)
	restart && (ids[1] = 0)
	gray = ids[1] ⊻ (ids[1] >> 1)
	for ℓ in 1:S
		obj = setindex(obj, isodd(gray >> (ℓ - 1)), ℓ)
	end
	value, obj = firstvalue(incremental, obj, z)
	cdcp = improve!(cdcp, value, obj)
	for k in (ids[1] + 1):(2^S - 1)
		ℓ = trailing_zeros(k) + 1
		obj = setindex(obj, !obj.ℒ[ℓ], ℓ)
		value, obj = nextvalue(incremental, obj, ℓ, z)
		cdcp = improve!(cdcp, value, obj)
		ids[1] = k
	end
	cdcp.obj = obj
	return cdcp
end

function improve!(cdcp, value, obj)
	if value > cdcp.value
		cdcp.value = value
		cdcp.x .= obj.ℒ
	end
	cdcp
end

firstvalue(::Nothing, obj, z) = obj(z)
function firstvalue(incremental, obj, z)
	value = incremental.start(z)
	for ℓ in eachindex(obj.ℒ)
		obj.ℒ[ℓ] && (value = incremental.flip(ℓ, true, z))
	end
	value, addfcall(obj)
end

nextvalue(::Nothing, obj, ℓ, z) = obj(z)
nextvalue(incremental, obj, ℓ, z) = incremental.flip(ℓ, obj.ℒ[ℓ], z), addfcall(obj)
