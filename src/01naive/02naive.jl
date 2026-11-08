function naive!(cdcp)
	(; obj) = cdcp
	S = length(cdcp.x)
	(; ids, z, incremental) = cdcp.solver
	value, obj = firstvalue(incremental, obj, z)
	if value > cdcp.value
		cdcp.value = value
		cdcp.x .= obj.ℒ
	end
	for k in (ids[1] + 1):(2^S - 1)
		ℓ = trailing_zeros(k) + 1
		obj = setindex(obj, !obj.ℒ[ℓ], ℓ)
		value, obj = nextvalue(incremental, obj, ℓ, z)
		if value > cdcp.value
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
