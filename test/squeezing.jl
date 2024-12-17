const v = [0.24952817563772145, 0.30744184685956255, 0.16527154238546204,
    0.21557003274986386, 0.3603164670276888, 0.8878175216561821, 0.005695423872065342,
    0.9097347266622412, 0.5856034888347053, 0.41973295503165087]

@testset "Squeezing Naive" begin
    f = TestObj(0.25, 10, v)
    obj0 = Objective(f, fill(true, 10))
    p0 = init(Naive, obj0, 10)
    @time solve!(p0)
    @test p0.obj.fcall == 2^10

    obj = Objective(f, SVector{10,Bool}(trues(10)))
    p = init(Squeezing, obj, 10, true)
    @time solve!(p)
    @test 100 < p.obj.fcall < 200
    @test (p.x .== included) == p0.x
    @test (p.x .!= excluded) == p0.x
    @test p.value == p0.value
end
