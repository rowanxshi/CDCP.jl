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
       CDCP,
       SolverState,
       ItemState,
       undetermined,
       included,
       excluded,
       aux,

       BruteForce,

       Squeezing,
       SqueezingTrace,

       StateChoice,
       Policy,
       SqueezingPolicy

include("interface.jl")
include("brute.jl")
include("squeezing.jl")
include("squeezingpolicy.jl")

end
