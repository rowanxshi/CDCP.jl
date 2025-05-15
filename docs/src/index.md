# Combinatorial Discrete Choice Problems
```@docs
CDCP
```
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

## Examples

### Single agent problem

 Suppose we have the objective function `obj(J)` already defined and that the CDCP is over `C` items. We can call the solver as follows.


```julia
solve(Squeezing, obj, true) # If the objective obeys SCD-C from above
solve(Squeezing, obj, false) # If the objective obeys SCD-C from below
```

Suppose we have the objective function `obj(J, z)` already defined for parameter `z` (e.g. productivity). Then, we can call the solver similarly for type `z=5`:
```julia
solve(Squeezing, obj, true; z=5) # If the objective obeys SCD-C from above
solve(Squeezing, obj, false; z=5) # If the objective obeys SCD-C from below
```

### Policy function

The solver now expects the objective function `obj(J, z)` to accept a Boolean vector `J` and a real number `z` describing type (e.g. productivity). It also requires a function `equalise_obj((J1, J2), l, r)`, which takes a pair of Boolean vectors `(J1, J2)` and returns the type of an agent indifferent between the two strategies.

 Suppose we have our three functions `obj(J, z)` and `equalise_obj((J1, J2), l, r)` defined. We can call the solver as follows.

If the objective obeys SCD-C from above:
```julia
solve(SqueezingPolicy, obj, true, equal_obj) # If the objective obeys SCD-C from above
solve(SqueezingPolicy, obj, false, equal_obj) # If the objective obeys SCD-C from below
```
