"""
Solve combinatorial discrete choice problems as in Arkolakis, Eckert and Shi (2024).

See also [`solve`](@ref), [`Squeezing`](@ref), [`SqueezingPolicy`](@ref).
"""
module CDCP

using Base: RefValue

import StaticArrays
using StaticArrays: SVector, setindex

import CommonSolve: init, solve!, solve

# CommonSolve
export init, solve!, solve

export Objective, margin, CDCProblem, CDCPSolver, SolverState, ItemState, undetermined, included, excluded, aux
export Naive, Squeezing, Policy, SqueezingPolicy
export policy, policy!

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
include("03policy/03refine.jl")

# for backwards compatibility
include("04compat/00structs.jl")
include("04compat/01interface.jl")
include("04compat/02helpers.jl")

end
