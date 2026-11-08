function solve!(cdcp::CDCProblem{<:Naive}; restart::Bool=true)
	(; obj) = cdcp
	(; ids, z, incremental) = cdcp.solver
	S = length(cdcp.x)
	S < Sys.WORD_SIZE - 1 || throw(ArgumentError("exhaustion over $S items is not feasible"))
	if restart
		resize!(ids, 1)
		ids[1] = 0
	end
	gray = ids[1] ⊻ (ids[1] >> 1)
	for ℓ in 1:S
		obj = setindex(obj, isodd(gray >> (ℓ - 1)), ℓ)
	end
	best = cdcp.value
	value, obj = firstvalue(incremental, obj, z)
	if value > best
		best = value
		cdcp.value = value
		cdcp.x .= obj.ℒ
	end
	for k in (ids[1] + 1):(2^S - 1)
		ℓ = trailing_zeros(k) + 1
		obj = setindex(obj, !obj.ℒ[ℓ], ℓ)
		value, obj = nextvalue(incremental, obj, ℓ, z)
		if value > best
			best = value
			cdcp.value = value
			cdcp.x .= obj.ℒ
		end
		ids[1] = k
	end
	cdcp.obj = obj
	return cdcp
end

function firstvalue(::Nothing, obj, z)
	obj(z)
end
function firstvalue(incremental, obj, z)
	value = incremental.start(z)
	for ℓ in eachindex(obj.ℒ)
		obj.ℒ[ℓ] && (value = incremental.flip(ℓ, true, z))
	end
	value, addfcall(obj)
end

function nextvalue(::Nothing, obj, ℓ, z)
	obj(z)
end
function nextvalue(incremental, obj, ℓ, z)
	incremental.flip(ℓ, obj.ℒ[ℓ], z), addfcall(obj)
end
