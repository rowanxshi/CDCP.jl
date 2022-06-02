"Solve CDCPs. Main functions provided are [`solve!`](@ref), [`solve`](@ref), and [`policy`](@ref). See [here](https://github.com/rowanxshi/CDCP.jl) for more."
module CDCP

include("helpers.jl")

include("./single/helpers.jl")
include("./single/squeezing.jl")
include("./single/branching.jl")
include("./single/solvers.jl")

include("./policy/helpers.jl")
include("./policy/squeezing.jl")
include("./policy/branching.jl")
include("./policy/brute.jl")
include("./policy/solvers.jl")

export solve!, solve, policy, policy!

end
