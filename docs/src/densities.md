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

## Bibliography

[^moerken1991]:   K. Mørken (1991). Some Identities for Products and Degree Raising of Splines. Constructive Approximation, 7(1), 195–208. http://dx.doi.org/10.1007/bf01888153
