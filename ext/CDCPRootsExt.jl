module CDCPRootsExt

using Roots
using CDCP

function (eo::CDCP.Equal_Obj{<:Roots.AbstractUnivariateZeroMethod})(
        obj1::CDCP.Objective, obj2::CDCP.Objective, zleft, zright)
    eo.fcall[] = 0
    f = CDCP.DiffObj(obj1, obj2)
    if isinf(zleft) && isinf(zright)
        z0 = 0
    elseif isinf(zleft)
        z0 = zright - 1
    elseif isinf(zright)
        z0 = zleft + 1
    else
        z0 = (zleft + zright) / 2
    end
    z = Roots.find_zero(f, z0, eo.m, eo.fcall; eo.kwargs...)
    return z, obj1
end

end # module
