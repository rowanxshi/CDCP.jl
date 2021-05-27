# CDCP

A Julia package to solve [combinatorial discrete choice problems](https://rowanxshi.github.io/papers/cdc.pdf), either for single agents or over a single dimension of heterogeneity.

## Installation

From command line:
``
julia -e 'using Pkg; Pkg.add(url="https://github.com/rowanxshi/CDCP.jl")
``

From within Julia:
``
import Pkg
Pkg.add(url="https://github.com/rowanxshi/CDCP.jl")
``

## Key functions

The package provides `solve!`, `solve`, and `policy`.

### Single agent problems

Either `solve!` or `solve` can handle single agent problems. The first solves the problem inplace, so it should be provided with three pre-allocated Boolean vectors `(sub, sup, aux)`, the objective function `π`, and whether the problem satisfies SCD-C from above.

``
solve!((sub, sup, aux), π, D_j_π, scdca::Bool; containers)
solve!((sub, sup, aux), π, scdca::Bool; containers)
``

The solver expects the objective function `π(J)` to accept a Boolean vector like `sub`.

There are two optional arguments: the marginal value function `D_j_π(j, J)` and `containers`. The marginal value function should accept an index `j` and a Boolean vector `J`, computing the marginal value of item `j` to the set `J`. If the marginal value function is omitted, the solver automatically generates one using the provided objective function `π`.

The `containers = (working, converged)` can be preallocated and passed to the solver. Both `working` and `converged` are both `Vectors` holding elements like `(sub, sup, aux)`. These are used for the branching step and will be automatically allocated if not supplied.

> **Examples**
> Suppose we have the objective function `π(J::BitVector)` and marginal value function `D_j_π(j::Integer, J::BitVector)` already defined and that the CDCP is over `C` items. We can preallocate everything, then call the solver as follows.
> ``
> sub = falses(C)
> sup = trues(C)
> aux = falses(C)
> working = [(sub, sup, aux); ]
> converged = similar(working)
>
> # if the objective obeys SCD-C from above
> solve!((sub, sup, aux), π, D_j_π, true, containers = (working, converged)) 
>
> # if the objective obeys SCD-C from below
> solve!((sub, sup, aux), π, D_j_π, false, containers = (working, converged)) 
> ``
> We could equally call the solver omitting the preallocated containers, the marginal value function, or both:
> ``
> solve!((sub, sup, aux), π, D_j_π, true)
> solve!((sub, sup, aux), π, true, containers = (working, converged)) 
> solve!((sub, sup, aux), π, true)
> ``

Pre-allocating is helpful if the problem will be solved many times because it avoids the solver allocating the vectors each time.

``
solve(C::Integer, π, D_j_π, scdca::Bool; containers)
solve(C::Integer, π, scdca::Bool; containers)
``

Non-inplace versions, where `C` is the number of items in the CDCP. ( `C` need not be specified for the inplace caller, since it can be inferred from the length of the Boolean vectors.) They are otherwise identical in usage to the inplace versions.

### Policy function

The package provides `policy` to identify the policy function for problems over a single dimension of heterogeneity (which we call _productivity_ from here).

``
policy(C::Integer, π, zero_D_j_π, equalise_π, scdca::Bool)
policy(C::Integer, π, equalise_π, scdca::Bool)
``

The solver now expects the objective function `π(J, z)` to accept a Boolean vector `J` and a real number `z` describing _productivity_. It also requires a function `equalise_π((J1, J2))`, which take a pair of Boolean vectors `(J1, J2)` and returns the _productivity_ of an agent indifferent between the two strategies.

A function `zero_D_j_π(j, J)` may optionally be provided, which accepts an index `j` and Boolean vector `J`. It should return the _productivity_ at which the marginal value of `j` to set `J` is zero. If the function is omitted, the solver automatically constructs one using the `equalise_π` function.

The policy function will be returned as a series of cutoff _productivities_ and the optimal decision sets for all types within each interval.

> **Examples**
> Suppose we have our three functions `π(J::BitVector, z::Float64)`, `zero_D_j_π(j::Integer, J::BitVector)` and `equalise_π((J1, J2))` defined. The CDCP is over `C = 3` items. We can call the solver as follows.
> ``
> # if the objective obeys SCD-C from above
> (cutoffs, policies) = policy(3, π, zero_D_j_π, equalise_π, true)
>
> # if the objective obeys SCD-C from below
> (cutoffs, policies) = policy(3, π, zero_D_j_π, equalise_π, false)
> ``
> We could equally omit the `zero_D_j_π` function, letting the solver generate one itself:
> ``
> (cutoffs, policies) = policy(3, π, equalise_π, false)
> ``
> The returned returned `cutoffs` will be a vector of cutoffs; say `cutoffs = [-Inf, 2, 4, 6, Inf]`. The returned `policies` will be a vector of Boolean arrays, say
> ``
> [
> 	[false; false; false];
> 	[false; false; true];
> 	[true; false; true];
> 	[true; true; true];
> ]
> ``
> The length of `cutoffs` will always exceed the length of `policies` by 1. This result means that for agents with _productivity_ in between `(-Inf, 2)`, the optimal strategy is `[false; false; false]` --- that is, to select none of the items. The optimal strategy for the next interval given by the `cutoffs`, `(2, 4)`, is `[false; false; true]`; and so on.  
