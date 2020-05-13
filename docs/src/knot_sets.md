# Knot sets

Three built-in knot sets are provided out-of-the-box:

### Arbitrary knot set

Arbitrary knot sets allow custom placement of the individual knots, as
well as custom multiplicity of interior knots. This can be used for
problems where discontinuity at a certain location is desired.

```@docs
ArbitraryKnotSet
ArbitraryKnotSet(k::Integer, t::AbstractVector, ml::Integer=k, mr::Integer=k)
```

### Linear knot set

```@docs
LinearKnotSet
LinearKnotSet(k::Integer, a, b, N::Integer, ml::Integer=k, mr::Integer=k)
```

### Exponential knot set

Exponential knot sets are useful for approximating exponentially
varying functions (e.g. bound states of atoms). The quadrature points
of each interval are distributed as were the knot set piecewise
linear.

```@docs
ExpKnotSet
ExpKnotSet(k::Integer, a::T, b::T, N::Integer, ml::Integer=k, mr::Integer=k; base::T=T(10), include0::Bool=true) where T
```

## Reference

```@docs
order
numintervals
numfunctions
first
last
length
getindex
nonempty_intervals
```

## Internals

```@docs
CompactBases.AbstractKnotSet
CompactBases.assert_multiplicities
CompactBases.leftmultiplicity
CompactBases.rightmultiplicity
CompactBases.find_interval
CompactBases.within_interval
CompactBases.within_support
```

### Quadrature functions

[Gauß–Legendre
quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss–Legendre_quadrature)
can be used to exactly calculate integrals of polynomials, and very
accurately the integrals of smoothly varying functions. The approximation is given by

$$\begin{equation}
\int\limits_a^b \diff{x}f(x)\approx
\frac{b-a}{2}\sum_{i=1}^n w_i
f\left(\frac{b-a}{2}x_i+\frac{a+b}{2}\right),
\end{equation}$$

where $x_i$ are the roots of the quadrature and $w_i$ the
corresponding weights, given on the elementary interval $[-1,1]$.

```@docs
CompactBases.num_quadrature_points
CompactBases.lgwt!
CompactBases.lgwt
```

