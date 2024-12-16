"""
Solve combinatorial discrete choice problems as in Arkolakis, Eckert and Shi (2023).
See [here](https://github.com/rowanxshi/CDCP.jl) for more.
"""
module CombinatorialDiscreteChoiceProblems

using Base: RefValue
using StaticArrays: SVector, setindex

import CommonSolve: init, solve!, solve

# CommonSolve
export init, solve!, solve

export Objective,
       addfcall,
       margin,
       CDCProblem,
       SolverState,
       ItemState,
       undetermined,
       included,
       excluded,
       aux,

       Naive,

       Squeezing,

       Policy,
       SqueezingPolicy

# used for everything
include("00global/00structs.jl")
include("00global/01interface.jl")

# naive solution with exhaustion
include("01naive/00structs.jl")
include("01naive/01interface.jl")
include("01naive/02naive.jl")

# single-agent solution
include("02single/00structs.jl")
include("02single/01interface.jl")
include("02single/02squeeze.jl")
include("02single/03branch.jl")

# policy function solution
include("03policy/00structs.jl")
include("03policy/01interface.jl")
include("03policy/02squeeze.jl")
include("03policy/03branch.jl")

# for backwards compatibility
include("04compat/00structs.jl")
include("04compat/01interface.jl")
include("04compat/02helpers.jl")

end
