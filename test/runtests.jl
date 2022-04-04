import Test

include("example_obj.jl")

function coincide(C::Int, scdca::Bool = false; seed::Int = 10)
	(cdcp, ) = initiate(C, scdca, seed = seed)
	
	J = falses(C)
	Vs = CDCP._containers(C)
	
	(cutoffs, policies) = CDCP.policy(C; cdcp...)
	
	z_wrong = nothing
	coincide = all(0.01:0.1:50) do z
		CDCP.naive!(J; cdcp.obj)
		CDCP.solve!(Vs; cdcp.scdca, cdcp.obj)
#
		interval = searchsortedfirst(cutoffs, z)-1
#		
		match = (J == first(Vs) == policies[interval])
		!match && (z_wrong = z)
		match
	end
#	
	(coincide, z_wrong)
end

#Test.@testset begin
#
#for seed in 10:10:100
#	Test.@test first(coincide(5, true, seed = seed))
#	Test.@test first(coincide(5, false, seed = seed))
#end
#
#end