function _containers(C::Integer)
	sub = falses(C)
	sup = trues(C)
	aux = falses(C)
	return (sub, sup, aux)
end

function _containers(Vs)
	working = [Vs; ]
	converged = similar(working)
	(working, converged)
end

function _checklength(::SVector{S,ItemState}, sub, sup, aux) where S
	length(sub) == length(sup) == length(aux) == S || throw(ArgumentError("length of sub, sup, aux are inconsistent with C"))
end

function _parse_state(sub::Bool, sup::Bool, isaux::Bool)
	if isaux
		return aux
	elseif sub != sup
		return undetermined
	elseif sub
		return included
	else
		return excluded
	end
end

function _parse_triplet(x::SVector{S,ItemState}, sub, sup, aux) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(_parse_state(
				@inbounds(sub[$i]), @inbounds(sup[$i]), @inbounds(aux[$i]))))
		end
		return :(_checklength(x, sub, sup, aux); SVector{S,ItemState}($ex))
	else
		_checklength(x, sub, sup, aux)
		return SVector{S,ItemState}(ntuple(i->_parse_state(sub[i], sup[i], aux[i]), S))
	end
end

function _invert_state!(sub::AbstractVector{Bool}, sup::AbstractVector{Bool},
		isaux::AbstractVector{Bool}, x::SVector{S,ItemState}) where S
	@inbounds for i in 1:S
		xi = x[i]
		if xi == aux
			isaux[i] = true
		else
			isaux[i] = false
			if xi == undetermined
				sub[i] = false
				sup[i] = true
			elseif xi == included
				sub[i] = true
				sup[i] = true
			else
				sub[i] = false
				sup[i] = false
			end
		end
	end
end

function zero_D_j(equalise_obj, holder::AbstractVector{Bool})
	zero_D_j_obj(j::Integer, J::AbstractVector{Bool}, extras...) = let equalise_obj = equalise_obj, holder = holder
		J2 = copyto!(holder, J)
		J2 = setindex!(J2, !J[j], j)
		equalise_obj((J, J2), extras...)
	end
end

function restart!((working, converged, done), C)
	empty!.((working, converged, done))
	int = Interval(_containers(C), -Inf, Inf)
	push!(working, int)
end

