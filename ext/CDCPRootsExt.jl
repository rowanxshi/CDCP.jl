module CDCPRootsExt

using Roots
using CombinatorialDiscreteChoiceProblems

const CDCP = CombinatorialDiscreteChoiceProblems

function (eo::CDCP.Equal_Obj{<:Roots.AbstractUnivariateZeroMethod})(
        obj1::CDCP.Objective, obj2::CDCP.Objective, lb, ub)
    eo.fcall[] = 0
    f = CDCP.DiffObj(obj1, obj2)
    if isinf(lb) && isinf(ub)
        x0 = 0
    elseif isinf(lb)
        x0 = ub - 1
    elseif isinf(ub)
        x0 = lb + 1
    else
        x0 = (lb + ub) / 2
    end
    z = Roots.find_zero(f, x0, eo.m, eo.fcall; eo.kwargs...)
    return z, obj1
end

end # module
