module CDCPRootsExt

using Roots
using CombinatorialDiscreteChoiceProblems

const CDCP = CombinatorialDiscreteChoiceProblems

function (eo::CDCP.Equal_Obj{<:Roots.AbstractUnivariateZeroMethod})(
        obj1::CDCP.Objective, obj2::CDCP.Objective, zleft, zright)
    eo.fcall[] = 0
    f = CDCP.DiffObj(obj1, obj2)
    if isinf(zleft) && isinf(zright)
        x0 = 0
    elseif isinf(zleft)
        x0 = zright - 1
    elseif isinf(zright)
        x0 = zleft + 1
    else
        x0 = (zleft + zright) / 2
    end
    z = Roots.find_zero(f, x0, eo.m, eo.fcall; eo.kwargs...)
    return z, obj1
end

end # module
