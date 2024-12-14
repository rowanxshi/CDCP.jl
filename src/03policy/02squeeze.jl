function squeeze!(p::CDCProblem{<:SqueezingPolicy}) pool, squeezing, branching = p.solver.pool, p.solver.squeezing, p.solver.branching
	while !isempty(squeezing)
		k = pop!(squeezing)
		x = pool[k]
		i = findfirst(==(undetermined), x.x)
		if i === nothing
			findfirst(==(aux), x.x) === nothing || push!(branching, k)
		else
			p.obj.fcall < p.maxfcall || return maxfcall_reached
			xs = squeeze!(p, x, i)
			pool[k] = xs[1]
			push!(squeezing, k)
			for j in 2:length(xs)
				push!(pool, xs[j])
				push!(squeezing, length(pool))
			end
		end
	end
	return inprogress
end

function squeeze!(p::CDCProblem{<:SqueezingPolicy}, x::IntervalChoice, i::Int)
	obj, scdca, tr = p.obj, p.solver.scdca, p.solver.trace
	lookup = p.solver.lookup_zero_margin
	obj = _setchoice(obj, scdca ? setsup(x.x) : setsub(x.x))
	key = (i, obj.x)
	z = get(lookup, key, nothing)
	if z === nothing
		z, obj = p.solver.zero_margin(obj, i, x.lb, x.ub)
		lookup[key] = z
		p.solver.zero_margin_call[] += 1
		p.obj = obj
	end
	if z <= x.lb # Include the whole interval
		x = _squeeze(x, included, i)
		tr === nothing || push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, included))
		return (x,)
	elseif z >= x.ub # Cannot include any part of the interval
		tr === nothing ||
			push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, undetermined))
		return squeeze_exclude!(p, x, i)
	else
		tr === nothing ||
			push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, undetermined))
		x1 = IntervalChoice(x.lb, z, x.x)
		xs = squeeze_exclude!(p, x1, i)
		x2 = IntervalChoice(z, x.ub, x.x)
		return xs..., x2
	end
end

function squeeze_exclude!(p::CDCProblem{<:SqueezingPolicy}, x::IntervalChoice, i::Int)
	obj, scdca, tr = p.obj, p.solver.scdca, p.solver.trace
	lookup = p.solver.lookup_zero_margin
	obj = _setchoice(obj, scdca ? setsub(x.x) : setsup(x.x))
	key = (i, obj.x)
	z = get(lookup, key, nothing)
	if z === nothing
		z, obj = p.solver.zero_margin(obj, i, x.lb, x.ub)
		lookup[key] = z
		p.solver.zero_margin_call[] += 1
		p.obj = obj
	end
	if z >= x.ub # Exclude the whole interval
		x = _squeeze(x, excluded, i)
		tr === nothing || push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, excluded))
		return (x,)
	elseif z <= x.lb # Cannot exclude any part of the interval
		x = _setitemstate(x, aux, i) # Don't use _squeeze
		tr === nothing || push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, aux))
		return (x,)
	else
		tr === nothing || push!(tr[end], SqueezingPolicyTrace(i, z, x.lb, x.ub, excluded))
		x1 = IntervalChoice(x.lb, z, _squeeze(x.x, excluded, i))
		x2 = IntervalChoice(z, x.ub, _setitemstate(x.x, aux, i))
		return x1, x2
	end
end

function _squeeze(x::IntervalChoice, s::ItemState, i::Int)
	IntervalChoice(x.lb, x.ub, _squeeze(x.x, s, i))
end

function _setitemstate(x::IntervalChoice, s::ItemState, i::Int)
	IntervalChoice(x.lb, x.ub, _setitemstate(x.x, s, i))
end

