struct TestObj
	δ::Float64
	S::Int
	val::Vector{Float64}
end

function (obj::TestObj)(ℒ::AbstractVector{Bool}, z::Real=1)
	length(ℒ) == obj.S || throw(DimensionMismatch())
	f = range(0.1, length=obj.S, step=0.1)
	return z * sum(i->ℒ[i]*obj.val[i], 1:obj.S)^obj.δ - sum(i->ℒ[i]*f[i], 1:obj.S)
end

@testset "Objective" begin
	f = TestObj(0.25, 10, rand(10))
	obj = Objective(f, trues(10))
	@test obj.fcall == 0
	value, obj = obj(nothing)
	@test value == f(trues(10))
	@test isone(obj.fcall)

	obj = clearfcall(obj)
	@test iszero(obj.fcall)
	value, obj = obj(2)
	@test value == f(trues(10), 2)

	value1, value0, obj = margin(obj, 2, nothing)
	ℒ1 = trues(10)
	ℒ0 = copy(ℒ1)
	ℒ0[2] = false
	@test value1 - value0 == f(ℒ1) - f(ℒ0)
	@test obj.fcall == 3
end
