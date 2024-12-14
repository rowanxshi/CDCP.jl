# With SCD-C from below, choice for an item only switches once
function setx0scdcb(x0::SVector{S,ItemState}, xl::SVector{S,ItemState}) where S
	if @generated
		ex = :(())
		for i in 1:S
			push!(ex.args, :(ifelse(xl[$i]==included, included, x0[$i])))
		end
		return :(SVector{S,ItemState}($ex))
	else
		return SVector{S,ItemState}(
            ntuple(i->ifelse(x[i]==included, included, x0[i]), S))
	end
end

function search!(p::CDCProblem{<:Squeezing}, matched, zl0, xl0, zr0, xr0, x0, obj2, equal_obj)
    zl, xl, zr, xr = zl0, xl0, zr0, xr0
    while true
        # xl and xr should always be different here
        obj1 = _setchoice(p.obj, setsub(xl))
        obj2 = _setchoice(obj2, setsub(xr))
        z, obj = equal_obj(obj1, obj2, zl, zr)
        # p.solver.equal_obj_call[] += 1
        # ! TODO Rowan proof double-check
        if zl < z < zr # Additional cutoff points in between
            x0next = p.solver.scdca ? x0 : setx0scdcb(x0, xl)
            xnew = _solvesingle!(p, z, x0next)
            if xnew == xl || xnew == xr
                push!(matched, IntervalChoice(zl, z, xl), IntervalChoice(z, zr, xr))
                if zr != zr0 # Move to the interval on the right
                    zl, xl = zr, xr
                    zr, xr = zr0, xr0
                else
                    return
                end
            else # More cutoff points on the left
                zr, xr = z, xnew
            end
        else # No split of interval
            if z == zl
                push!(matched, IntervalChoice(zl, zr, xr))
            elseif z == zr
                push!(matched, IntervalChoice(zl, zr, xl))
            else
                error("z is not in between")
            end
            if zr != zr0 # Move to the interval on the right
                zl, xl = zr, xr
                zr, xr = zr0, xr0
            else
                return
            end
        end
    end
end

function _solvesingle!(p::CDCProblem{<:Squeezing}, z, x0)
    _reset!(p, z, x0)
    solve!(p)
    p.state == success || error("single-agenet solver fails with z = ", p.solver.z)
    return p.x
end

function branching!(p::CDCProblem{<:SqueezingPolicy}, k::Int, itask::Int)
    pool, matched = p.solver.pool, p.solver.matcheds[itask]
    x = pool[k]
    sp = p.solver.singlesolvers[itask]
    xl = _solvesingle!(sp, x.lb, x.x)
    xr = _solvesingle!(sp, x.ub, x.x)
    if xl == xr
        pool[k] = IntervalChoice(x.lb, x.ub, xl)
    else
        search!(sp, matched, x.lb, xl, x.ub, xr, x.x, p.solver.obj2, p.solver.equal_obj)
        # Overwrite the old interval with aux
        pool[k] = pop!(matched)
    end
end

function branching!(p::CDCProblem{<:SqueezingPolicy})
    pool, branching, matcheds = p.solver.pool, p.solver.branching, p.solver.matcheds
    ntasks = length(matcheds)
    for m in matcheds
        empty!(m) # Just to be safe
    end
    if ntasks == 1 # No multithreading
        for k in branching
            branching!(p, k, 1)
        end
    else
        @sync for itask in 1:ntasks
            Threads.@spawn for ik in itask:ntasks:length(branching)
                branching!(p, branching[ik], itask)
            end
        end
    end
    append!(pool, matcheds...)
    return inprogress
end


function concat!(p::CDCProblem{<:SqueezingPolicy})
    pool = p.solver.pool
    sort!(pool, by=_lb)
    cutoffs = resize!(p.x.cutoffs, 1)
    xs = resize!(p.x.xs, 1)
    cutoffs[1] = pool[1].lb
    xs[1] = xlast = pool[1].x
    for x in pool
        # Filter out potential singletons
        if x.lb < x.ub && x.x != xlast
            push!(cutoffs, x.lb)
            push!(xs, x.x)
            xlast = x.x
        end
    end
    return success
end
