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
julia> R = BSpline(LinearKnotSet(7, 0.0, 1.0, 3))[:,2:end-1]
BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 7 on 0.0..1.0 (3 intervals), restricted to basis functions 2..8 ⊂ 1..9

julia> r = axes(R,1)
Inclusion(0.0..1.0)

julia> c = R \ sin.(2π*r)
7-element Array{Float64,1}:
  0.35236622516315186
  1.0329925571692269
  1.6856805930016943
 -5.528339670985388e-16
 -1.6856805930016923
 -1.0329925571692293
 -0.35236622516315186

julia> f = R*c
Spline on BSpline{Float64} basis with LinearKnotSet(Float64) of order
k = 7 on 0.0..1.0 (3 intervals)

julia> f'f
0.4999507480992042

julia> norm(f)
0.7070719539758342
```

For unrestricted bases, the metric is stored in the B-spline basis
object. Beware though, when working with [Restricted bases](@ref), the
metric needs to be recomputed every time an inner product is
computed.
```julia
julia> @benchmark f'f
BenchmarkTools.Trial:
  memory estimate:  928 bytes
  allocs estimate:  25
  --------------
  minimum time:     5.112 μs (0.00% GC)
  median time:      5.536 μs (0.00% GC)
  mean time:        6.227 μs (2.44% GC)
  maximum time:     1.539 ms (98.89% GC)
  --------------
  samples:          10000
  evals/sample:     6

julia> @benchmark norm(f)
BenchmarkTools.Trial:
  memory estimate:  944 bytes
  allocs estimate:  26
  --------------
  minimum time:     5.163 μs (0.00% GC)
  median time:      5.605 μs (0.00% GC)
  mean time:        6.209 μs (2.26% GC)
  maximum time:     1.415 ms (99.13% GC)
  --------------
  samples:          10000
  evals/sample:     6
```

When working with non-orthogonal basis functions, such as B-splines,
it is therefore advised to explicitly compute the metric once and use
that for all inner products and norms:
```julia
julia> S = R'R
7×7 BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:
 0.0368555    0.0285178   0.0102464    0.00267575  0.000472084  4.3428e-5   2.70996e-8
 0.0285178    0.0492836   0.0347641    0.0176553   0.00640938   0.00142346  4.3428e-5
 0.0102464    0.0347641   0.0402597    0.0319871   0.0180037    0.00640938  0.000472084
 0.00267575   0.0176553   0.0319871    0.0380685   0.0319871    0.0176553   0.00267575
 0.000472084  0.00640938  0.0180037    0.0319871   0.0402597    0.0347641   0.0102464
 4.3428e-5    0.00142346  0.00640938   0.0176553   0.0347641    0.0492836   0.0285178
 2.70996e-8   4.3428e-5   0.000472084  0.00267575  0.0102464    0.0285178   0.0368555

julia> @benchmark dot(c,S,c)
BenchmarkTools.Trial:
  memory estimate:  16 bytes
  allocs estimate:  1
  --------------
  minimum time:     80.383 ns (0.00% GC)
  median time:      88.495 ns (0.00% GC)
  mean time:        98.239 ns (1.42% GC)
  maximum time:     6.183 μs (98.38% GC)
  --------------
  samples:          10000
  evals/sample:     944
```
