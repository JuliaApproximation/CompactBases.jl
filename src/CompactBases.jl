module CompactBases

import Base: eltype, axes, size, ==, getindex, checkbounds, copyto!, similar, show, step
import Base.Broadcast: materialize

using ContinuumArrays
import ContinuumArrays: Basis, ℵ₁, Derivative, Inclusion, @simplify,
    AbstractQuasiArrayApplyStyle

using QuasiArrays
import QuasiArrays: AbstractQuasiMatrix, QuasiAdjoint, MulQuasiArray,
    PInvQuasiMatrix, InvQuasiMatrix, QuasiDiagonal,
    BroadcastQuasiArray, SubQuasiArray

using IntervalSets

using LazyArrays

using LinearAlgebra
import LinearAlgebra: Matrix, dot
using SparseArrays
using FillArrays

using OffsetArrays
using BandedMatrices
using BlockBandedMatrices

using FastGaussQuadrature

using Formatting
using RecipesBase

vandermonde(B::Basis) = B[locs(B),:]

include("distributions.jl")
include("restricted_bases.jl")
include("types.jl")
include("unit_vectors.jl")
include("intervals.jl")
include("knot_sets.jl")
include("quadrature.jl")
include("skewtridiag.jl")

include("materialize_dsl.jl")

include("fornberg.jl")
include("finite_differences.jl")
include("fd_operators.jl")
include("fd_derivatives.jl")

include("fedvr.jl")
include("fedvr_operators.jl")
include("fedvr_derivatives.jl")

include("bsplines.jl")
include("bspline_derivatives.jl")

include("inner_products.jl")
include("densities.jl")
include("linear_operators.jl")


"""
    centers(B)

Return the locations of the mass centers of all basis functions of
`B`; for orthogonal bases such as finite-differences and FE-DVR, this
is simply [`locs`](@ref), i.e. the location of the quadrature nodes.
"""
centers(B::BasisOrRestricted) = locs(B)

export AbstractFiniteDifferences,
    FiniteDifferences, StaggeredFiniteDifferences, ImplicitFiniteDifferences,
    Derivative, dot, QuasiDiagonal, Inclusion, .., distribution,
    FEDVR, Derivative, @elem, dot,
    BSpline,
    ArbitraryKnotSet, LinearKnotSet, ExpKnotSet,
    order, numintervals, numfunctions, nonempty_intervals,
    centers, vandermonde, Density,
    LinearOperator, DiagonalOperator, ShiftAndInvert

end # module