function squeeze!(p::CDCProblem{<:Squeezing}, x::AbstractVector{ItemState})
	S = length(x)
	fx = -Inf
	lastaux = nothing
	i = findfirst(==(undetermined), x)
	if i === nothing
		lastaux = findfirst(==(aux), x) # May find aux from initial value
		if lastaux === nothing # Last value set by branching
			fx, obj = value(p.obj, p.solver.z)
			p.obj = obj
		end
	end
	while i !== nothing && p.obj.fcall < p.maxfcall
		x, fx, itemstate = squeeze!(p, x, i)
		if itemstate == aux
			lastaux = lastaux === nothing ? i : min(lastaux, i)
		else
			lastaux = nothing
		end
		#lastaux = ifelse(itemstate == aux, i, nothing)
		if i < S
			i = findnext(==(undetermined), x, i+1)
			i === nothing && (i = findfirst(==(undetermined), x))
		else
			i = findfirst(==(undetermined), x)
		end
	end
	# Whenever squeeze! makes progress, any aux is changed to undetermined
	# If aux is in x, it can only be that lastaux !== nothing
	if lastaux === nothing
		state = p.obj.fcall < p.maxfcall ? success : maxfcall_reached
	else
		push!(p.solver.branching, branch(x, lastaux)...)
		state = inprogress
	end
	return x, fx, state
end

function squeeze!(p::CDCProblem{<:Squeezing}, x::AbstractVector{ItemState}, i::Int)
	obj, scdca, z, tr = p.obj, p.solver.scdca, p.solver.z, p.solver.trace
	# For excluding: look at the best case
	if scdca
		obj = _setchoice(obj, setsub(x))
		f1, f0, obj = margin(obj, i, z)
		exclude = f1 - f0 <= 0
	else
		obj = _setchoice(obj, setsup(x))
		f1, f0, obj = margin(obj, i, z)
		exclude = f1 - f0 < 0
	end
	if exclude
		xnew = _squeeze(x, excluded, i)
		tr === nothing || push!(tr[end], SqueezingTrace(i, x, excluded, f0))
		p.obj = obj
		return xnew, f0, excluded
	end

	# For including: look at the worst case
	if scdca
		obj = _setchoice(obj, setsup(x))
		f1, f0, obj = margin(obj, i, z)
		include = f1 - f0 > 0
	else
		obj = _setchoice(obj, setsub(x))
		f1, f0, obj = margin(obj, i, z)
		include = f1 - f0 >= 0
	end
	if include
		xnew = _squeeze(x, included, i)
		tr === nothing || push!(tr[end], SqueezingTrace(i, x, included, f1))
		p.obj = obj
		return xnew, f1, included
	end

	xnew = _setitemstate(x, aux, i)
	fx = convert(typeof(p.fx), -Inf)
	tr === nothing || push!(tr[end], SqueezingTrace(i, x, aux, fx))
	p.obj = obj
	return xnew, fx, aux
end

function _squeeze(x::SVector{S,ItemState}, s::ItemState, k::Int) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(_squeeze(x, s, k, $i)))
		end
		return :(SVector{S,ItemState}($ex))
	else
		return SVector{S,ItemState}(ntuple(i->_squeeze(x, s, k, i), S))
	end
end

function _squeeze(x::SVector{S,ItemState}, s::ItemState, k::Int, i::Int) where S
	if i == k
		return s
	else
		@inbounds xi = x[i]
		if xi == aux
			return undetermined
		else
			return xi
		end
	end
end

