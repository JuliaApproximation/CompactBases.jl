# Non-uniform grids

By employing non-uniform grids, more effort can be concentrated to
those parts where the function is expected to vary more, hopefully
leading to a better
approximation. [`StaggeredFiniteDifferences`](@ref) supports both
uniform and non-uniform node distributions; in the former case, the
node locations are stored as an `AbstractRange`, in the latter, as any
other kind of `AbstractVector`.

## Integrals

Since the quadrature weights in the non-uniform case are
location-dependent, they are baked into to the expansion coefficients
(as in the [FE-DVR case](@ref FE-DVR quadrature)), such that
```math
\int\diff{x} f(x) \approx \sum_i f_i,
```
where ``f_i`` is the ``i``th expansion coefficient of ``f(x)``. In the uniform
case, the expansion coefficients instead equal the function value at
the nodes, which means
```math
\int\diff{x} f(x) \approx \Delta x \sum_i f_i.
```
To be sure that you always get the right integration weight, use the
metric to perform integrals, e.g.
```math
\int\diff{x} f(x)
```
is computed using
```julia
integral = sum(S*f)
```
where the metric matrix can be found using
```julia
S = B'B
```
Similarly, integrals of the kind
```math
\int\diff{x} \conj{f}(x)g(x)
```
are easily computed as
```julia
integral = dot(f, S, g)
```

## Example

We start by constructing two [`StaggeredFiniteDifferences`](@ref)
grids, one uniform and one log–linear, of approximately the same
extent (at least including ``r=30``):

```jldoctest
julia> N = 50
50

julia> rmax = 30.0
30.0

julia> B1 = StaggeredFiniteDifferences(rmax, N)
Staggered finite differences basis {Float64} on 0.0..30.606060606060606 with 50 points spaced by ρ = 0.6060606060606061

julia> ρmin = 0.001 # Min step-size
0.001

julia> ρmax = 1.0 # Max step-size
1.0

julia> α = 0.1 # Step-size change rate
0.1

julia> B2 = StaggeredFiniteDifferences(ρmin, ρmax, α, rmax)
Staggered finite differences basis {Float64} on 0.0..31.475628188861677 with 103 points
```

As mentioned above, the metrics are different in the uniform and
non-uniform cases:

```jldoctest
julia> S1 = B1'B1
UniformScaling{Float64}
0.6060606060606061*I

julia> S2 = B2'B2
UniformScaling{Bool}
true*I
```

Since the step-size is changing for non-uniform grids,
[`step`](@ref) is ill-defined, but returns unity to simplify the
implementation of the differential operators (the coefficients are
divided by `step(B)^k`, where `k` is the order of the differential;
since that is already accounted for in the case of non-uniform grids
by scaling the coefficients, `1^k` does not change the result):

```jldoctest
julia> step(B1)
0.6060606060606061

julia> step(B2)
1.0
```

The node locations and weights are shown in the figure below:
![Non-uniform grid](figures/fd/compare_staggered_non_uniform_grid.svg)

The basis functions become asymmetric in the non-uniform case:
```jldoctest
julia> # Compute basis functions on a dense grid
       ξ = 10.0 .^ range(-3, stop=log10(25), length=1000);

julia> χ1 = B1[ξ, :];

julia> χ2 = B2[ξ, :];
```

![Non-uniform basis functions](figures/fd/compare_staggered_non_uniform_basis.svg)

We then expand the function
```math
f(x) = \sin(2\pi x/30)\exp(-4x/30)
```
on both grids and compare the error:

```jldoctest
julia> f = x -> sin(2π*x/rmax)*exp(-4x/rmax)
#47 (generic function with 1 method)

julia> xx1 = axes(B1,1)
Inclusion(0.0..30.606060606060606)

julia> c1 = B1 \ f.(xx1);


julia> xx2 = axes(B2,1)
Inclusion(0.0..31.475628188861677)

julia> c2 = B2 \ f.(xx2);

julia> f1 = χ1*c1;

julia> f2 = χ2*c2;

julia> fe = f.(ξ);
```
![Non-uniform basis functions](figures/fd/compare_staggered_non_uniform_reconstruction.svg)
