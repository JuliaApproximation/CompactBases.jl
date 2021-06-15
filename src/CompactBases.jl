module CompactBases
# The package does support precompilation, but doing `using CompactBases` with
# precompilation on can lead to a huge spike in memory. That does not happen if disable
# precompilation for CompactBases with __precompile__(false).
#
# So, if this is a problem for the user, it is possible to disable precompilation by setting
# the COMPACTBASES_NO_PRECOMPILE environment variable. E.g. the user can set it when
# launching Julia by doing
#
#     COMPACTBASES_NO_PRECOMPILE= julia
#
# or set it in the ~/.julia/config/startup.jl file by setting
#
#     ENV["COMPACTBASES_NO_PRECOMPILE"] = ""
#
# Note that any value will do, even an empty string, as only the presence of the variable is
# checked.
haskey(ENV, "COMPACTBASES_NO_PRECOMPILE") && __precompile__(false)

import Base: eltype, axes, size, ==, getindex, checkbounds, copyto!, similar, show, step
import Base.Broadcast: materialize

using ContinuumArrays
import ContinuumArrays: Basis, ℵ₁, Derivative, Inclusion, @simplify, simplifiable, _simplify,
    AbstractQuasiArrayApplyStyle

using QuasiArrays
import QuasiArrays: AbstractQuasiMatrix, QuasiAdjoint, MulQuasiArray,
    PInvQuasiMatrix, InvQuasiMatrix, QuasiDiagonal,
    BroadcastQuasiArray, SubQuasiArray, LazyQuasiMatrix

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

include("coulomb_derivative.jl")

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

"""
    metric(B)

Returns the metric or mass matrix of the basis `B`, equivalent to
`S=B'B`.
"""
metric(B::BasisOrRestricted) = B'B

"""
    metric_shape(B)

Returns the shape of the metric or mass matrix of the basis `B`,
i.e. `Diagonal`, `BandedMatrix`, etc.
"""
metric_shape(B::BasisOrRestricted) = typeof(metric(B))

function matrix_element_metric(B::BasisOrRestricted)
    # For non-orthogonal bases, we need to apply the metric inverse,
    # after computing the action of a matrix representation of a
    # linear operator on a ket. However, when computing the inner
    # product with the bra afterwards, we again use the metric. Here,
    # we short-circuit this by combining them into them identity
    # operator I, if the operator metric equals the metric. For
    # orthogonal bases, with the integration weights not built into the
    # coefficients, and an identity matrix operator metric, we /do/
    # need to use the metric for the inner product.
    S = metric(B)
    oS = operator_metric(B)
    oS == S ? I : (oS \ S)
end

export AbstractFiniteDifferences,
    FiniteDifferences, StaggeredFiniteDifferences, ImplicitFiniteDifferences,
    Derivative, dot, QuasiDiagonal, Inclusion, .., distribution,
    FEDVR, Derivative, @elem, dot,
    BSpline,
    ArbitraryKnotSet, LinearKnotSet, ExpKnotSet,
    order, numintervals, numfunctions, nonempty_intervals,
    centers, vandermonde,
    metric, metric_shape, matrix_element_metric,
    FunctionProduct, Density,
    LinearOperator, DiagonalOperator, ShiftAndInvert

end # module
