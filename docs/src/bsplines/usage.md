# Usage

We first load the package
```julia
using CompactBases
```

Then we define a linearly spaced knot set between `0` and `1` with
five intervals and cubic splines. By default, full multiplicity of the
endpoints is assumed.
```julia
julia> t = LinearKnotSet(4, 0, 1, 5)
12-element LinearKnotSet{4,4,4,Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}:
 0.0
 0.0
 0.0
 0.0
 0.2
 0.4
 0.6
 0.8
 1.0
 1.0
 1.0
 1.0

julia> B = BSpline(t)
BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 4 (cubic) on 0.0..1.0 (5 intervals)

julia> size(B)
(ContinuumArrays.AlephInfinity{1}(), 8)
```
The last statement means that `B` is a quasimatrix with a continuous
first dimension which is spanned by 8 basis functions.

The overlap matrix is created simply:

```julia
julia> S = B'B
8×8 BandedMatrices.BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:
 0.0285714    0.0175      0.00369048  0.000238095  …   ⋅           ⋅           ⋅
 0.0175       0.0442857   0.03125     0.00690476       ⋅           ⋅           ⋅
 0.00369048   0.03125     0.0653571   0.0449206       3.96825e-5   ⋅           ⋅
 0.000238095  0.00690476  0.0449206   0.095873        0.00474206  5.95238e-5   ⋅
  ⋅           5.95238e-5  0.00474206  0.0472619       0.0449206   0.00690476  0.000238095
  ⋅            ⋅          3.96825e-5  0.00474206   …  0.0653571   0.03125     0.00369048
  ⋅            ⋅           ⋅          5.95238e-5      0.03125     0.0442857   0.0175
  ⋅            ⋅           ⋅           ⋅              0.00369048
  0.0175      0.0285714
```

## Number of quadrature points

As explained in the theory section on [Integrals](@ref), Gauß–Legendre
quadrature is used to approximate integrals between B-splines. It is
possible to decide how many quadrature points are used by passing an
argument to the [`BSpline`](@ref) constructor:

```julia
julia> B = BSpline(t,4)
BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 4 (cubic) on 0.0..1.0 (5 intervals)
```

If too few quadrature points are specified, a warning is issued:

```julia
julia> B = BSpline(t,3)
┌ Warning: N = 3 quadrature points not enough to calculate overlaps between polynomials of order k = 4
└ @ CompactBases ~/.julia/dev/CompactBases/src/quadrature.jl:48
BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 4 (cubic) on 0.0..1.0 (5 intervals)
```

## Dual bases

It is sometimes useful to work with dual bases, e.g. B-splines of
different orders. This is possible as long as both bases are resolved
on the same quadrature points:

```julia
julia> tl = LinearKnotSet(3, 0, 1, 2)
7-element LinearKnotSet{3,3,3,Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}:
 0.0
 0.0
 0.0
 0.5
 1.0
 1.0
 1.0

julia> tr = LinearKnotSet(4, 0, 1, 2)
9-element LinearKnotSet{4,4,4,Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}:
 0.0
 0.0
 0.0
 0.0
 0.5
 1.0
 1.0
 1.0
 1.0

julia> N = 4
4

julia> L = BSpline(tl, N)
BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 3 (parabolic) on 0.0..1.0 (2 intervals)

julia> R = BSpline(tr, N)
BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 4 (cubic) on 0.0..1.0 (2 intervals)

julia> S = L'R
4×5 BandedMatrices.BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:
 0.0833333   0.0645833   0.0166667  0.00208333   ⋅
 0.0375      0.129167    0.108333   0.0541667   0.00416667
 0.00416667  0.0541667   0.108333   0.129167    0.0375
 0.0         0.00208333  0.0166667  0.0645833   0.0833333
```


## Evaluation of B-splines

It is the possible to query the values of e.g. the first basis functions at some values of `x`:

```julia
julia> B[0:0.1:1,1]
11-element SparseArrays.SparseVector{Float64,Int64} with 2 stored entries:
  [1 ]  =  1.0
  [2 ]  =  0.125
```

We can also evaluate all the B-splines at the same time:

```julia
julia> B[0:0.25:1,:]
5×8 SparseArrays.SparseMatrixCSC{Float64,Int64} with 14 stored entries:
  [1, 1]  =  1.0
  [2, 2]  =  0.105469
  [2, 3]  =  0.576823
  [3, 3]  =  0.0208333
  [2, 4]  =  0.315104
  [3, 4]  =  0.479167
  [4, 4]  =  0.00260417
  [2, 5]  =  0.00260417
  [3, 5]  =  0.479167
  [4, 5]  =  0.315104
  [3, 6]  =  0.0208333
  [4, 6]  =  0.576823
  [4, 7]  =  0.105469
  [5, 8]  =  1.0
```

Since the B-splines have compact support, they are only locally
non-zero, hence the sparse storage.

Finally, we can compute the value of a single B-spline, a range, or
all of them at a single point `x`:

```julia
julia> B[0.5,4]
0.47916666666666674

julia> B[0.5,:]
8-element Array{Float64,1}:
 0.0
 0.0
 0.020833333333333325
 0.47916666666666674
 0.4791666666666666
 0.020833333333333322
 0.0
 0.0

julia> B[0.5,4:8]
5-element Array{Float64,1}:
 0.47916666666666674
 0.4791666666666666
 0.020833333333333322
 0.0
 0.0
```

## Reference

```@docs
BSpline
BSpline(t::CompactBases.AbstractKnotSet, N)
BSpline(t::CompactBases.AbstractKnotSet; k′=3)
```
