using CompactBases
import CompactBases: locs,
    unrestricted_basis, restriction_extents

using IntervalSets
using QuasiArrays
import QuasiArrays: AbstractQuasiArray,  AbstractQuasiMatrix, MulQuasiArray
using ContinuumArrays
import ContinuumArrays: ℵ₁, Inclusion

using LinearAlgebra
using BandedMatrices
using BlockBandedMatrices
using SparseArrays

using LazyArrays
import LazyArrays: materialize
using FillArrays

using ArnoldiMethod
using Random

using Test

function vecdist(a::AbstractVector, b::AbstractVector,
                 ϵ = eps(eltype(a)))
    δ = √(sum(abs2, a-b))
    δ, δ/√(sum(abs2, a .+ ϵ))
end

include("derivative_accuracy_utils.jl")

@testset "CompactBases" begin
    include("fd/runtests.jl")
    include("fedvr/runtests.jl")
    include("bsplines/runtests.jl")
    include("inner_products.jl")
end
