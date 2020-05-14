# CompactBases.jl

CompactBases.jl is a package for various bases for function
approximation with compact support in the framework of
[ContinuumArrays](https://github.com/JuliaApproximation/ContinuumArrays.jl)

![CompactBases logo](assets/logo.svg)

## Example usage

```julia-repl
julia> a,b = 0.0,1.0 # Extents
(0.0, 1.0)

julia> N = 13 # Number of nodes
13

julia> k = 5 # Order of FE-DVR/B-splines
5

julia> x = range(a, stop=b, length=1000)
0.0:0.001001001001001001:1.0

julia> Δx = (b-a)/(N+1) # Grid spacing
0.07142857142857142

julia> # Standard, uniform finite-differences
           fd = FiniteDifferences(N, Δx)
Finite differences basis {Float64} on 0.0..1.0 with 13 points spaced by Δx = 0.07142857142857142

julia> # Staggered, uniform finite-differences
           sfd_uni = StaggeredFiniteDifferences(N, Δx)
Staggered finite differences basis {Float64} on 0.0..0.9642857142857144 with 13 points spaced by ρ = 0.07142857142857142

julia> # Staggered, non-uniform finite-differences
           sfd_nonuni = StaggeredFiniteDifferences(0.01, 0.5, 0.1, b)
Staggered finite differences basis {Float64} on 0.0..1.1112312795594823 with 39 points

julia> # Finite-element boundaries
           tf = range(a, stop=b, length=N+2)
0.0:0.07142857142857142:1.0

julia> # We can vary the polynomial order in each element
           forder = vcat(7, fill(4,length(tf)-2))
14-element Array{Int64,1}:
 7
 4
 4
 4
 4
 4
 4
 4
 4
 4
 4
 4
 4
 4

julia> # By indexing the second dimension, we can implement Dirichlet0
           # boundary conditions.
           fem = FEDVR(tf, forder)[:,2:end-1]
FEDVR{Float64} basis with 14 elements on 0.0..1.0, restricted to elements 1:14, basis functions 2..45 ⊂ 1..46

julia> tb = ExpKnotSet(k, -2.0, log10(b), N+1)
23-element ExpKnotSet{5,5,5,Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},Array{Float64,1}}:
 0.0
 0.0
 0.0
 0.0
 0.0
 0.01
 0.014251026703029978
 0.020309176209047358
 0.028942661247167503
 0.04124626382901352
 0.05878016072274912
 0.0837677640068292
 0.11937766417144363
 0.17012542798525887
 0.24244620170823283
 0.3455107294592219
 0.4923882631706739
 0.7017038286703828
 1.0
 1.0
 1.0
 1.0
 1.0

julia> splines = BSpline(tb)[:,2:end-1]
BSpline{Float64} basis with ExpKnotSet(Float64) of  on order k = 5 (quartic) on 0,0.01..1.0 (14 intervals), restricted to basis functions 2..17 ⊂ 1..18
```

![Simple example](figures/simple_example.svg)

## Reference

```@index
```

```@autodocs
Modules = [CompactBases]
```
