# Backwards-compatible methods

```@docs
solve(C::Integer; scdca::Bool, obj)
solve!((sub, sup, aux); scdca::Bool, obj, restart::Bool=true, kwargs...)
policy(C::Integer; kwargs...)
policy!((cutoffs, policies), containers, C::Integer; scdca::Bool, obj, equalise_obj, zero_D_j_obj = zero_D_j(equalise_obj, falses(C)), show_time::Bool = false, restart::Bool = true, singlekw=NamedTuple(), kwargs...)
CDCP.naive(C::Integer; obj, z=nothing)
CDCP.naive!(J::AbstractVector{Bool}; obj)
```

