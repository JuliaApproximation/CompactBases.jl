# Inner products and norms

## Inner products

The inner product of two functions is simply computed as
```julia
julia> B = FiniteDifferences(10, 0.1)
Finite differences basis {Float64} on 0.0..1.1 with 10 points spaced by Δx = 0.1

julia> x = axes(B,1)
Inclusion(0.0..1.1)

julia> c = B \ sin.(2π*x)
10-element Array{Float64,1}:
  0.5877852522924731
  0.9510565162951535
  0.9510565162951536
  0.5877852522924732
  1.2246467991473532e-16
 -0.587785252292473
 -0.9510565162951535
 -0.9510565162951536
 -0.5877852522924734
 -2.4492935982947064e-16

julia> f = B*c
QuasiArrays.ApplyQuasiArray{Float64,1,typeof(*),Tuple{FiniteDifferences{Float64,Int64},Array{Float64,1}}}(*, (Finite differences basis {Float64} on 0.0..1.1 with 10 points spaced by Δx = 0.1, [0.5877852522924731, 0.9510565162951535, 0.9510565162951536, 0.5877852522924732, 1.2246467991473532e-16, -0.587785252292473, -0.9510565162951535, -0.9510565162951536, -0.5877852522924734, -2.4492935982947064e-16]))

julia> d = B \ cos.(2π*x)
10-element Array{Float64,1}:
  0.8090169943749475
  0.30901699437494745
 -0.30901699437494734
 -0.8090169943749473
 -1.0
 -0.8090169943749475
 -0.30901699437494756
  0.30901699437494723
  0.8090169943749473
  1.0

julia> g = B*d
QuasiArrays.ApplyQuasiArray{Float64,1,typeof(*),Tuple{FiniteDifferences{Float64,Int64},Array{Float64,1}}}(*, (Finite differences basis {Float64} on 0.0..1.1 with 10 points spaced by Δx = 0.1, [0.8090169943749475, 0.30901699437494745, -0.30901699437494734, -0.8090169943749473, -1.0, -0.8090169943749475, -0.30901699437494756, 0.30901699437494723, 0.8090169943749473, 1.0]))

julia> f'f
0.5000000000000001

julia> f'g
-3.0044051106072847e-17

julia> g'g
0.5
```

## Norms

Similarly, norms are easily computed as
```julia
julia> norm(f)
0.7071067811865476

julia> norm(g)
0.7071067811865476
```

For orthogonal bases, other ``p``-norms are also available:
```julia
julia> norm(g,1)
0.647213595499958

julia> norm(g,Inf)
1.0
```

## Caveats for non-orthogonal bases

Also for B-splines, the inner products and norms can be computed as
above:
```julia
c = ones(size(B,2))
f = B*c
```

For unrestricted bases, the metric is stored in the B-spline basis
object. Beware though, when working with [Restricted bases](@ref), the
metric needs to be recomputed every time an inner product is
computed.
```julia
@benchmark f'f
@benchmark norm(f)
```

When working with non-orthogonal basis functions, such as B-splines,
it is therefore advised to explicitly compute the metric once and use
that for all inner products and norms:
```julia
S = B'B
@benchmark dot(c,S,c)
```
