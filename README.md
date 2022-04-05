# CDCP

A Julia package to solve [combinatorial discrete choice problems](https://rowanxshi.github.io/papers/cdc.pdf), either for single agents or over a single dimension of heterogeneity.

## Installation

From command line:

```julia
julia -e 'import Pkg; Pkg.add(url="https://github.com/rowanxshi/CDCP.jl#main")'
```

From within Julia:

```julia
import Pkg
Pkg.add(url="https://github.com/rowanxshi/CDCP.jl#main")
```

## Key functions

The package provides `solve!`, `solve`, and `policy`.

### Single agent problems

Either `solve!` or `solve` can handle single agent problems. The first solves the problem inplace, so it should be provided with three pre-allocated Boolean vectors `(sub, sup, aux)`, the objective function `obj`, and whether the problem satisfies SCD-C from above.

```julia
solve!((sub, sup, aux); scdca::Bool, obj, [D_j_obj; containers])
```

The solver expects the objective function `obj(J)` to accept a Boolean vector like `sub`.

There are two optional arguments: the marginal value function `D_j_obj(j, J)` and `containers`. The marginal value function should accept an index `j` and a Boolean vector `J`, computing the marginal value of item `j` to the set `J`. If the marginal value function is omitted, the solver automatically generates one using the provided objective function `obj`.

The `containers = (working, converged)` can be preallocated and passed to the solver. Both `working` and `converged` are `Vectors` holding elements like `(sub, sup, aux)`. These are used for the branching step and will be automatically allocated if not supplied.

> **Examples**
>
> Suppose we have the objective function `obj(J)` and marginal value function `D_j_obj(j, J:)` already defined and that the CDCP is over `C` items. We can preallocate everything, then call the solver as follows.
>
> ```julia
> sub = falses(C)
> sup = trues(C)
> aux = falses(C)
> working = [(sub, sup, aux); ]
> converged = similar(working)
> ```
>
> If the objective obeys SCD-C from above:
>
> ```julia
> solve!((sub, sup, aux); scdca = true, obj, D_j_obj, containers = (working, converged)) 
> ```
>
> If the objective obeys SCD-C from below:
>
> ```julia
> solve!((sub, sup, aux); scdca = false, obj, D_j_obj, containers = (working, converged)) 
> ```
>
> We could equally call the solver omitting the preallocated containers, the marginal value function, or both:
>
> ```julia
> solve!((sub, sup, aux); scdca = true, obj, D_j_obj)
> solve!((sub, sup, aux); scdca = true, obj, containers = (working, converged)) 
> solve!((sub, sup, aux); scdca = true, obj)
> ```

Pre-allocating is helpful if the problem will be solved many times because it avoids the solver allocating the vectors each time.

```julia
solve(C::Integer; scdca::Bool, obj, [D_j_obj; containers])
```

Non-inplace versions, where `C` is the number of items in the CDCP. (`C` need not be specified for the inplace caller, since it can be inferred from the length of the Boolean vectors.) They are otherwise identical in usage to the inplace versions.

### Policy function

The package provides `policy` to identify the policy function for problems over a single dimension of heterogeneity (which we call _productivity_ from here).

```julia
policy(C::Integer; scdca::Bool, obj, equalise_obj, [D_j_obj, zero_D_j_obj])
```

The solver now expects the objective function `obj(J, z)` to accept a Boolean vector `J` and a real number `z` describing _productivity_. It also requires a function `equalise_obj((J1, J2), l, r)`, which takes a pair of Boolean vectors `(J1, J2)` and returns the _productivity_ of an agent indifferent between the two strategies. The solver provides the interval `[l, r]` within which the search should be performed; if the marginal type is not in the interval, it is sufficient to return `nothing`.

Similarly to the single-agent solver, the marginal value function `D_j_obj` is optional; it is automatically generated using the `obj` function if omitted.

A function `zero_D_j_obj(j, J, l, r)` may optionally be provided, which accepts an index `j`, Boolean vector `J`, and interval endpoints `l` and `r`. It should return the _productivity_ at which the marginal value of `j` to set `J` is zero. The solver provides the interval `[l, r]` within which this marginal type is located. If the function is omitted, the solver automatically constructs one using the `equalise_obj` function.

The policy function will be returned as a series of cutoff _productivities_ and the optimal decision sets for all types within each interval.

> **Examples**
>
> Suppose we have our three functions `obj(J, z)`, `zero_D_j_obj(j, J, l, r)` and `equalise_obj((J1, J2), l, r)` defined. The CDCP is over `C = 3` items. We can call the solver as follows.
>
> If the objective obeys SCD-C from above:
> ```julia
> (cutoffs, policies) = policy(3; scdca = true, obj, equalise_obj, zero_D_j_obj)
> ```
>
> If the objective obeys SCD-C from below
> ```julia
> (cutoffs, policies) = policy(3; scdca = false, obj, equalise_obj, zero_D_j_obj)
> ```
>
> We could equally omit the `zero_D_j_obj` function, letting the solver generate one itself:
>
> ```julia
> (cutoffs, policies) = policy(3; scdca = false, obj, equalise_obj)
> ```
>
> The returned returned `cutoffs` will be a vector of cutoffs; say `cutoffs = [-Inf, 2, 4, 6, Inf]`. The returned `policies` will be a vector of Boolean arrays, say
>
> ```julia
> [
> 	[false; false; false];
> 	[false; false; true];
> 	[true; false; true];
> 	[true; true; true];
> ]
> ```
>
> The length of `cutoffs` will always exceed the length of `policies` by 1. This result means that for agents with _productivity_ in between `(-Inf, 2)`, the optimal strategy is `[false; false; false]` --- that is, to select none of the items. The optimal strategy for the next interval given by the `cutoffs`, `(2, 4)`, is `[false; false; true]`; and so on.  
