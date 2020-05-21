# Densities

The _mutual density_ between two functions ``f(x)`` and ``g(x)`` is
defined as

```math
(f\cdot g)(x) \defd \conj{f}(x)g(x)
```

where the term density is used in analogy with [density matrices in
quantum mechanics](https://en.wikipedia.org/wiki/Density_matrix), and
for which reason the first function is conjugated. Without the
conjugate, it is equivalent to the [pointwise
product](https://en.wikipedia.org/wiki/Pointwise_product).

The mutual density ``h(x)`` of two functions expanded over two sets of basis
functions

```math
\begin{aligned}
f(x) &= \sum_i f_i \ket{A_i}, &
g(x) &= \sum_i g_i \ket{B_i},
\end{aligned}
```

can be expressed as one expansion over one set of basis functions,

```math
h(x) = \sum_i h_i \ket{C_i},
```

and this can be computed using [`Density`](@ref), as long as
``\ket{A_i}`` and ``\ket{B_i}`` are compatible. The space spanned by
``\ket{C_i}`` must be able to express products of the original spaces,
i.e. if ``f(x)`` and ``g(x)`` are both quadratic functions, formally
``h(x)`` will be a quartic function. There exist rather complicated
algorithms for producing the expansion coefficients for the product
function expanded over B-splines of higher order than the constituent
    factors (see e.g. [^moerken1991]), however this is rarely needed in
practice since we either 
1. use enough intervals to successfully approximate the function as
   piecewise polynomials of lower order,
2. use high enough polynomial order of our basis functions that we can
   represent the factors and their product satisfactorily,
3. or both of the above.

Therefore, we instead turn to the pragmatic approach which is based on
the [Vandermonde
matrix](https://en.wikipedia.org/wiki/Vandermonde_matrix): if we want
to find the expansion coefficients of ``h(x)`` on the basis
``\ket{C_i}``, we compute the Vandermonde matrix ``\mat{V}_C`` (using
[`vandermonde`](@ref)) for that basis and solve

```math
\mat{V}_C \vec{h} = h(\vec{x})
```

where the right-hand side ``h(\vec{x})`` denotes the function values
of ``h(x)`` on all the interpolation points (returned by
[`locs`](@ref)), but these are in turn given by ``(f\cdot
g)(\vec{x})``, which we can reconstruct from the expansion
coefficients of ``f(x)`` and ``g(x)`` by multiplying those by _their_
respective Vandermonde matrices ``\mat{V}_A`` and ``\mat{V}_B``. In
total, the expansion coefficients for the product function are thus
given by

```math
\vec{h} =
\mat{V}_C^+
[\mat{V}_A
\conj{\vec{f}}
\odot
\mat{V}_B
\vec{g}]
```

where ``\mat{V}_C^+`` is the Moore–Penrose inverse of ``\mat{V}_C``
which can easily be [computed using its singular-value
decomposition](https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Singular_value_decomposition_(SVD)),
and ``\odot`` denotes the [elementwise or Hadamard
product](https://en.wikipedia.org/wiki/Hadamard_product_(matrices)),
i.e. Julia's `.*` product. For this to work, naturally the bases
``\ket{A_i}``, ``\ket{B_i}``, and ``\ket{C_i}`` must share
interpolation points ``\vec{x}`` and thus their Vandermonde matrices
must have the same number of rows. For simplicity, CompactBases.jl
assumes that ``\ket{B_i}`` and ``\ket{C_i}`` are the same and that
``\ket{A_i}`` is "compatible" (checked by
[`assert_compatible_bases`](@ref)).

For orthogonal bases such as finite-difference and FE-DVR, where the
Vandermonde matrix is diagonal and same for both bases, the above
formula simplifies to

```math
\vec{h} =
\mat{V}
[\conj{\vec{f}}
\odot
\vec{g}],
```

and in the case of finite-differences on a uniform grid, where the
expansion coefficients coincide with the function values at the
interpolation points, the Vandermonde matrix reduces to the identity
matrix and the product formula is simpler still:

```math
\vec{h} =
\conj{\vec{f}}
\odot
\vec{g}.
```

## Example

```jldoctest
julia> import CompactBases: applied

julia> f = x -> sin(2π*x);

julia> g = x -> x*exp(-x);

julia> h = x -> f(x)*g(x);

julia> rmax,k,N = 10.0,7,71
(10.0, 7, 71)

julia> R = BSpline(LinearKnotSet(k, 0, rmax, N))
BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 7 on 0.0..10.0 (71 intervals)

julia> r = axes(R,1)
Inclusion(0.0..10.0)

julia> cf = R \ f.(r);

julia> cg = R \ g.(r);

julia> ch = R \ h.(r); # Direct expansion of density for comparison

julia> ρ = Density(applied(*,R,cf), applied(*,R,cg));

julia> ρ.ρ # Expansion coefficients computed via Vandermonde interpolation
77-element Array{Float64,1}:
  2.5990302560741227e-6
 -1.2089008321819097e-5
  0.016650515757241868
  0.08609474140758361
  0.22264566607626868
  0.28428278922020694
  0.04487154524021238
 -0.3097583765428373
 -0.45136208173491543
 -0.23816476137564302
  0.16360184909444878
  0.43464800730290326
  0.3730492731330008
  0.046047307534972166
 -0.2846933885764034
 -0.3797631540817138
  ⋮
  0.002725515878867043
  0.0015277308449457526
 -0.00041784201453291337
 -0.0016586000070113436
 -0.0015285004379182573
 -0.00041624819494116155
  0.0007250143682327532
  0.0011338505005407673
  0.0007023771068706905
 -9.73144173106993e-5
 -0.0006544921666713287
 -0.0006548469749281168
 -0.0004243163194033977
 -0.00020768615644842965
 -6.69583797265351e-5
 -5.669216664398177e-10

julia> norm(ρ.ρ - ch)
1.5656487551918337e-6
```

![Mutual densities](figures/mutual_densities.svg)

## Bibliography

[^moerken1991]:   K. Mørken (1991). Some Identities for Products and Degree Raising of Splines. Constructive Approximation, 7(1), 195–208. <http://dx.doi.org/10.1007/bf01888153>
