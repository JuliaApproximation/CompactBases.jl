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

using BandedMatrices
using BlockBandedMatrices

using FastGaussQuadrature

using Printf
using RecipesBase

include("restricted_bases.jl")
include("unit_vectors.jl")
include("intervals.jl")
include("knot_sets.jl")
include("quadrature.jl")
include("skewtridiag.jl")

include("materialize_dsl.jl")

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

export AbstractFiniteDifferences, FiniteDifferences, RadialDifferences, ImplicitFiniteDifferences,
    Derivative, dot, QuasiDiagonal, Inclusion, ..

export FEDVR, Derivative, @elem, dot

export BSpline

export ArbitraryKnotSet, LinearKnotSet, ExpKnotSet, order, numintervals, numfunctions, nonempty_intervals

end # module
