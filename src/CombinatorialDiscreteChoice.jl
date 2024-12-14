"""
Solve combinatorial discrete choice problems as in Arkolakis, Eckert and Shi (2023).
See [here](https://github.com/rowanxshi/CDCP.jl) for more.
"""
module CombinatorialDiscreteChoice

using Base: RefValue
using StaticArrays: SVector, setindex

import CommonSolve: init, solve!, solve

# CommonSolve
export init, solve!, solve

export Objective,
       addfcall,
       value,
       margin,
       CDCProblem,
       SolverState,
       ItemState,
       undetermined,
       included,
       excluded,
       aux,

       BruteForce,

       Squeezing,
       SqueezingTrace,

       Policy,
       SqueezingPolicy

include("00helpers/interface.jl")
include("00helpers/brute.jl")
include("01single/single.jl")
include("02policy/policy.jl")
include("03compat/wrapper.jl")

end
