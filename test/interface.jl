struct TestObj
    δ::Float64
    S::Int
    val::Vector{Float64}
end

function (obj::TestObj)(x::AbstractVector{Bool}, z::Real=1)
    length(x) == obj.S || throw(DimensionMismatch())
    f = range(0.1, length=obj.S, step=0.1)
    return z * sum(i->x[i]*obj.val[i], 1:obj.S)^obj.δ - sum(i->x[i]*f[i], 1:obj.S)
end

@testset "Objective" begin
    f = TestObj(0.25, 10, rand(10))
    obj = Objective(f, trues(10))
    @test obj.fcall == 0
    v, obj = obj(nothing)
    @test v == f(trues(10))
    @test obj.fcall == 1

    obj = _clearfcall(obj)
    @test obj.fcall == 0
    v, obj = obj(2)
    @test v == f(trues(10), 2)

    f1, f0, obj = margin(obj, 2, nothing)
    x1 = trues(10)
    x0 = copy(x1)
    x0[2] = false
    @test f1 - f0 == f(x1) - f(x0)
    @test obj.fcall == 3
end
