# CompactBases.jl

A package for representing various bases constructed from basis
functions with compact support as quasi-arrays.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaApproximation.github.io/CompactBases.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaApproximation.github.io/CompactBases.jl/dev)
[![Build Status](https://github.com/JuliaApproximation/CompactBases.jl/workflows/CI/badge.svg)](https://github.com/JuliaApproximation/CompactBases.jl/actions)
[![Codecov](https://codecov.io/gh/JuliaApproximation/CompactBases.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaApproximation/CompactBases.jl)

This package implements bases with compactly supported bases functions
as quasi-arrays where one one axes is continuous and the other axis is
discrete (countably infinite), as implemented in
[QuasiArrays.jl](https://github.com/JuliaApproximation/QuasiArrays.jl)
and
[ContinuumArrays.jl](https://github.com/JuliaApproximation/ContinuumArrays.jl). The
bases available initially are various finite-differences (FD),
finite-elements discrete-variable representation (FE-DVR), and
B-splines, with the possibility of adding more in the future.

Note that this package is written by a pragmatic physicist, so you may
find it lacking in mathematical rigour.

## Example

The big advantage of this framework is that the user code does not
need to be aware of the underlying details of the basis employed, at
least that is the goal. As an example, we look at how to construct the
mass matrices and second derivative matrices for a few different bases.

```julia
using CompactBases

function test_basis(B)
    println(repeat("-", 100))
    display(B)
    @info "Mass matrix"
    S = B'B
    display(S)

    # This is the continuous axis
    x = axes(B,1)

    # This corresponds to a operator L whose action on a function
    # f(x) is defined as Lf(x) = sin(2πx)*f(x). In physics this is a
    # potential.
    @info "Sine operator"
    f = QuasiDiagonal(sin.(2π*x))
    display(B'*f*B)

    @info "Laplacian"
    D = Derivative(x)
    display(B'*D'*D*B)
    println(repeat("-", 100))
    println()
end

a,b = 0,1 # Extents
N = 3 # Number of nodes
k = 5 # Order of FE-DVR/B-splines
```

### Finite-differences

The available finite-differences (as of present) are three-point
stencils, with the first grid point at either `Δx` (normal) or `Δx/2`
(staggered).

```julia
Δx = (b-a)/(N+1) # Grid spacing
# Standard, uniform finite-differences
test_basis(FiniteDifferences(N, Δx))
```

```
----------------------------------------------------------------------------------------------------
Finite differences basis {Float64} on 0.0..1.0 with 3 points spaced by Δx = 0.25
[ Info: Mass matrix
UniformScaling{Float64}
0.25*I
[ Info: Sine operator
3×3 Diagonal{Float64,Array{Float64,1}}:
 1.0   ⋅             ⋅
  ⋅   1.22465e-16    ⋅
  ⋅    ⋅           -1.0
[ Info: Laplacian
3×3 SymTridiagonal{Float64,Array{Float64,1}}:
 -32.0   16.0     ⋅
  16.0  -32.0   16.0
    ⋅    16.0  -32.0
----------------------------------------------------------------------------------------------------
```

```julia
# Staggered, uniform finite-differences
test_basis(StaggeredFiniteDifferences(N, Δx))
```

```
----------------------------------------------------------------------------------------------------
Staggered finite differences basis {Float64} on 0.0..0.875 with 3 points spaced by ρ = 0.25
[ Info: Mass matrix
UniformScaling{Float64}
0.25*I
[ Info: Sine operator
3×3 Diagonal{Float64,Array{Float64,1}}:
 0.707107   ⋅          ⋅
  ⋅        0.707107    ⋅
  ⋅         ⋅        -0.707107
[ Info: Laplacian
3×3 SymTridiagonal{Float64,Array{Float64,1}}:
 -65.25     21.3333     ⋅
  21.3333  -35.5556   17.0667
    ⋅       17.0667  -33.28
----------------------------------------------------------------------------------------------------
```

### FE-DVR

The FE-DVR implementation follows

- Rescigno, T. N., & McCurdy, C. W. (2000). Numerical Grid Methods for
  Quantum-Mechanical Scattering Problems. Physical Review A, 62(3),
  032706. http://dx.doi.org/10.1103/physreva.62.032706

The scalar operators are diagonal, whereas differential operators are
almost block-diagonal, with one-element overlaps.

```julia
# Finite-element boundaries
tf = range(a, stop=b, length=N+2)
# By indexing the second dimension, we can implement Dirichlet0
# boundary conditions.
test_basis(FEDVR(tf, max(2,k))[:,2:end-1])
```

```
----------------------------------------------------------------------------------------------------
FEDVR{Float64} basis with 4 elements on 0.0..1.0, restricted to elements 1:4, basis functions 2..16 ⊂ 1..17
[ Info: Mass matrix
UniformScaling{Bool}
true*I
[ Info: Sine operator
15×15 Diagonal{Float64,Array{Float64,1}}:
 0.267921   ⋅         ⋅         ⋅    ⋅         ⋅         ⋅         ⋅             ⋅          ⋅          ⋅          ⋅     ⋅          ⋅          ⋅
  ⋅        0.707107   ⋅         ⋅    ⋅         ⋅         ⋅         ⋅             ⋅          ⋅          ⋅          ⋅     ⋅          ⋅          ⋅
  ⋅         ⋅        0.963441   ⋅    ⋅         ⋅         ⋅         ⋅             ⋅          ⋅          ⋅          ⋅     ⋅          ⋅          ⋅
  ⋅         ⋅         ⋅        1.0   ⋅         ⋅         ⋅         ⋅             ⋅          ⋅          ⋅          ⋅     ⋅          ⋅          ⋅
  ⋅         ⋅         ⋅         ⋅   0.963441   ⋅         ⋅         ⋅             ⋅          ⋅          ⋅          ⋅     ⋅          ⋅          ⋅
  ⋅         ⋅         ⋅         ⋅    ⋅        0.707107   ⋅         ⋅             ⋅          ⋅          ⋅          ⋅     ⋅          ⋅          ⋅
  ⋅         ⋅         ⋅         ⋅    ⋅         ⋅        0.267921   ⋅             ⋅          ⋅          ⋅          ⋅     ⋅          ⋅          ⋅
  ⋅         ⋅         ⋅         ⋅    ⋅         ⋅         ⋅        1.22465e-16    ⋅          ⋅          ⋅          ⋅     ⋅          ⋅          ⋅
  ⋅         ⋅         ⋅         ⋅    ⋅         ⋅         ⋅         ⋅           -0.267921    ⋅          ⋅          ⋅     ⋅          ⋅          ⋅
  ⋅         ⋅         ⋅         ⋅    ⋅         ⋅         ⋅         ⋅             ⋅        -0.707107    ⋅          ⋅     ⋅          ⋅          ⋅
  ⋅         ⋅         ⋅         ⋅    ⋅         ⋅         ⋅         ⋅             ⋅          ⋅        -0.963441    ⋅     ⋅          ⋅          ⋅
  ⋅         ⋅         ⋅         ⋅    ⋅         ⋅         ⋅         ⋅             ⋅          ⋅          ⋅        -1.0    ⋅          ⋅          ⋅
  ⋅         ⋅         ⋅         ⋅    ⋅         ⋅         ⋅         ⋅             ⋅          ⋅          ⋅          ⋅   -0.963441    ⋅          ⋅
  ⋅         ⋅         ⋅         ⋅    ⋅         ⋅         ⋅         ⋅             ⋅          ⋅          ⋅          ⋅     ⋅        -0.707107    ⋅
  ⋅         ⋅         ⋅         ⋅    ⋅         ⋅         ⋅         ⋅             ⋅          ⋅          ⋅          ⋅     ⋅          ⋅        -0.267921
[ Info: Laplacian
7×7-blocked 15×15 BlockBandedMatrices.BlockSkylineMatrix{Float64,Array{Float64,1},BlockBandedMatrices.BlockSkylineSizes{Tuple{BlockArrays.BlockedUnitRange{Array{Int64,1}},BlockArrays.BlockedUnitRange{Array{Int64,1}}},Array{Int64,1},Array{Int64,1},BandedMatrix{Int64,Array{Int64,2},Base.OneTo{Int64}},Array{Int64,1}}}:
 -746.667    298.667    -74.6667  │     33.0583  │      ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅
  298.667   -426.667    298.667   │    -90.5097  │      ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅
  -74.6667   298.667   -746.667   │    758.901   │      ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅
 ─────────────────────────────────┼──────────────┼───────────────────────────────────┼──────────────┼───────────────────────────────────┼──────────────┼─────────────────────────────────
   33.0583   -90.5097   758.901   │  -2240.0     │   758.901    -90.5097    33.0583  │    -16.0     │      ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅
 ─────────────────────────────────┼──────────────┼───────────────────────────────────┼──────────────┼───────────────────────────────────┼──────────────┼─────────────────────────────────
     ⋅          ⋅          ⋅      │    758.901   │  -746.667    298.667    -74.6667  │     33.0583  │      ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅
     ⋅          ⋅          ⋅      │    -90.5097  │   298.667   -426.667    298.667   │    -90.5097  │      ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅
     ⋅          ⋅          ⋅      │     33.0583  │   -74.6667   298.667   -746.667   │    758.901   │      ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅
 ─────────────────────────────────┼──────────────┼───────────────────────────────────┼──────────────┼───────────────────────────────────┼──────────────┼─────────────────────────────────
     ⋅          ⋅          ⋅      │    -16.0     │    33.0583   -90.5097   758.901   │  -2240.0     │   758.901    -90.5097    33.0583  │    -16.0     │      ⋅          ⋅          ⋅
 ─────────────────────────────────┼──────────────┼───────────────────────────────────┼──────────────┼───────────────────────────────────┼──────────────┼─────────────────────────────────
     ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅      │    758.901   │  -746.667    298.667    -74.6667  │     33.0583  │      ⋅          ⋅          ⋅
     ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅      │    -90.5097  │   298.667   -426.667    298.667   │    -90.5097  │      ⋅          ⋅          ⋅
     ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅      │     33.0583  │   -74.6667   298.667   -746.667   │    758.901   │      ⋅          ⋅          ⋅
 ─────────────────────────────────┼──────────────┼───────────────────────────────────┼──────────────┼───────────────────────────────────┼──────────────┼─────────────────────────────────
     ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅      │    -16.0     │    33.0583   -90.5097   758.901   │  -2240.0     │   758.901    -90.5097    33.0583
 ─────────────────────────────────┼──────────────┼───────────────────────────────────┼──────────────┼───────────────────────────────────┼──────────────┼─────────────────────────────────
     ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅      │    758.901   │  -746.667    298.667    -74.6667
     ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅      │    -90.5097  │   298.667   -426.667    298.667
     ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅      │       ⋅      │      ⋅          ⋅          ⋅      │     33.0583  │   -74.6667   298.667   -746.667
----------------------------------------------------------------------------------------------------
```

### B-splines

All operators become banded when using B-splines, including the mass
matrix, which leads to generalized eigenvalue problems, among other
things.

```julia
tb = LinearKnotSet(k, a, b, N+1)
test_basis(BSpline(tb)[:,2:end-1])
```

```
----------------------------------------------------------------------------------------------------
BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 5 (quartic) on 0.0..1.0 (4 intervals), restricted to basis functions 2..7 ⊂ 1..8
[ Info: Mass matrix
6×6 BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:
 0.0414683   0.0307567    0.00934744  0.00101411  2.75573e-6    ⋅
 0.0307567   0.0581184    0.0437077   0.0124663   0.000610548  2.75573e-6
 0.00934744  0.0437077    0.0786449   0.0543455   0.0124663    0.00101411
 0.00101411  0.0124663    0.0543455   0.0786449   0.0437077    0.00934744
 2.75573e-6  0.000610548  0.0124663   0.0437077   0.0581184    0.0307567
  ⋅          2.75573e-6   0.00101411  0.00934744  0.0307567    0.0414683
[ Info: Sine operator
6×6 BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:
 0.0243786     0.0237789    0.00817377    0.000914686   1.89458e-6     ⋅
 0.0237789     0.0508203    0.0345991     0.00707979    7.42166e-20  -1.89458e-6
 0.00817377    0.0345991    0.034669      5.37956e-18  -0.00707979   -0.000914686
 0.000914686   0.00707979   6.40955e-18  -0.034669     -0.0345991    -0.00817377
 1.89458e-6    5.7276e-20  -0.00707979   -0.0345991    -0.0508203    -0.0237789
  ⋅           -1.89458e-6  -0.000914686  -0.00817377   -0.0237789    -0.0243786
[ Info: Laplacian
6×6 BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:
 -7.08571    -0.530159    1.01587    0.253968   0.0031746    ⋅
 -0.530159   -2.53333    -0.26455    0.756966   0.161552    0.0031746
  1.01587    -0.26455    -1.8455    -0.310406   0.756966    0.253968
  0.253968    0.756966   -0.310406  -1.8455    -0.26455     1.01587
  0.0031746   0.161552    0.756966  -0.26455   -2.53333    -0.530159
   ⋅          0.0031746   0.253968   1.01587   -0.530159   -7.08571
----------------------------------------------------------------------------------------------------
```
