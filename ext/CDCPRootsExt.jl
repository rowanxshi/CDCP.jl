module CDCPRootsExt

using Roots
using CombinatorialDiscreteChoice

const CDCP = CombinatorialDiscreteChoice

function (eo::CDCP.Equal_Obj)(obj1::CDCP.Objective, obj2::CDCP.Objective, lb, ub)
    eo.fcall[] = 0
    f = DiffObj(obj1, obj2)
    try
        z = Roots.find_zero(f, (lb, ub), eo.M, eo.fcall; eo.kwargs...)
    catch
        # Caller should recognize the implication for a value on the bound
        z = f(lb, eo.fcall) > 0 ? lb : ub
    end
    return z, CDCP.Objective(obj.f, obj.x, obj.fcall+zm.fcall[])
end

end # module
