using Test
using CombinatorialDiscreteChoiceProblems

using Random
using Roots
using StaticArrays

using CombinatorialDiscreteChoiceProblems: Equal_Obj, clearfcall

const CDCP = CombinatorialDiscreteChoiceProblems

const tests = [
    "interface",
    "squeezing",
    "squeezingpolicy",
    "wrapper"
]

printstyled("Running tests:\n", color=:blue, bold=true)

@time for test in tests
    include("$test.jl")
    println("\033[1m\033[32mPASSED\033[0m: $(test)")
end
