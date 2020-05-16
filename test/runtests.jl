using CompactBases
import CompactBases: locs,
    unrestricted_basis, restriction_extents

using IntervalSets
using QuasiArrays
using ContinuumArrays
import ContinuumArrays: ℵ₁, Inclusion

using LinearAlgebra
using BandedMatrices
using BlockBandedMatrices
using LazyArrays
import LazyArrays: materialize
using FillArrays

using Test

include("derivative_accuracy_utils.jl")
include("fedvr/runtests.jl")
include("inner_products.jl")
