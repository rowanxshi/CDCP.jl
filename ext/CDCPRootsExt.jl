module CDCPRootsExt

using Roots
using CombinatorialDiscreteChoice

const CDCP = CombinatorialDiscreteChoice

function (zm::CDCP.Zero_Margin)(obj::CDCP.Objective, i, lb, ub)
    zm.fcall[] = 0
    f = Margin_i(obj, i)
    try
        z = Roots.find_zero(f, (lb, ub), zm.M, zm.fcall; zm.kwargs...)
    catch
        z = f(lb, zm.fcall) > 0 ? lb : ub
    end
    return z, CDCP.Objective(obj.f, obj.x, obj.fcall+zm.fcall[])
end

function (eo::CDCP.Equal_Obj)(obj1::CDCP.Objective, obj2::CDCP.Objective, lb, ub)
    eo.fcall[] = 0
    f = DiffObj(obj1, obj2)
    # Should already know that a zero point exists between lb and ub
    z = Roots.find_zero(f, (lb, ub), eo.M, eo.fcall; eo.kwargs...)
    return z, CDCP.Objective(obj.f, obj.x, obj.fcall+zm.fcall[])
end

end # module
