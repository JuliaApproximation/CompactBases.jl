abstract type AbstractFiniteDifferences{T} <: Basis{T} end
const RestrictedFiniteDifferences{T,B<:AbstractFiniteDifferences{T}} =
    RestrictedQuasiArray{T,2,B}
const FiniteDifferencesOrRestricted = BasisOrRestricted{<:AbstractFiniteDifferences}
const AdjointFiniteDifferencesOrRestricted = AdjointBasisOrRestricted{<:AbstractFiniteDifferences}

eltype(::AbstractFiniteDifferences{T}) where T = T

axes(B::AbstractFiniteDifferences{T}) where T = (Inclusion(leftendpoint(B)..rightendpoint(B)), Base.OneTo(length(locs(B))))
size(B::AbstractFiniteDifferences) = (ℵ₁, length(locs(B)))

distribution(B::AbstractFiniteDifferences) = distribution(locs(B))

# step is only well-defined for uniform grids
step(::Uniform, B::AbstractFiniteDifferences) = step(locs(B))
# All-the-same, we define step for non-uniform grids to make the mass
# matrix and derivative stencils work properly.
step(::NonUniform, B::AbstractFiniteDifferences{T}) where T = one(T)
step(B::AbstractFiniteDifferences) = step(distribution(B), B)

"""
    local_step(B, i)

The step size around grid point `i`. For uniform grids, this is
equivalent to [`step`](@ref).
"""
local_step(B::AbstractFiniteDifferences, _) = step(B)

weight(::Uniform, B::AbstractFiniteDifferences{T}, _) where T = one(T)
weight(::NonUniform, B::AbstractFiniteDifferences, j) = 1/√(local_step(B, j))
weight(B::AbstractFiniteDifferences, j) = weight(distribution(B), B, j)

weights( ::Uniform, B::AbstractFiniteDifferences{T}) where T = ones(T,size(B,2))
weights(d::NonUniform, B::AbstractFiniteDifferences) = weight.(Ref(d), Ref(B), 1:size(B,2))
weights(B::AbstractFiniteDifferences) = weights(distribution(B), B)

inverse_weight(::Uniform, B::AbstractFiniteDifferences{T}, _) where T = one(T)
inverse_weight(::NonUniform, B::AbstractFiniteDifferences, j) = √(local_step(B, j))
inverse_weight(B::AbstractFiniteDifferences, j) = inverse_weight(distribution(B), B, j)

inverse_weights( ::Uniform, B::AbstractFiniteDifferences{T}) where T = ones(T,size(B,2))
inverse_weights(d::NonUniform, B::AbstractFiniteDifferences) = inverse_weight.(Ref(d), Ref(B), 1:size(B,2))
inverse_weights(B::AbstractFiniteDifferences) = inverse_weights(distribution(B), B)

==(A::AbstractFiniteDifferences,B::AbstractFiniteDifferences) = locs(A) == locs(B)
Base.hash(R::FD, h::UInt) where {FD<:AbstractFiniteDifferences} = hash(locs(R), hash(FD, h))

assert_compatible_bases(A::FiniteDifferencesOrRestricted, B::FiniteDifferencesOrRestricted) =
    locs(A) == locs(B) ||
        throw(ArgumentError("Can only multiply finite-differences sharing the same nodes"))

function loc(B::AbstractFiniteDifferences, i)
    N = size(B,2)
    if i < 1
        loc(B,1) - local_step(B,1)
    elseif i > N
        loc(B,N) + local_step(B,N)
    else
        locs(B)[i]
    end
end

tent(x, a, b, c) =
    1 - clamp((b-x)/(b-a), 0, 1) - clamp((b-x)/(b-c), 0, 1)

getindex(B::AbstractFiniteDifferences, x::Real, i::Integer) =
    weight(B,i)*tent(x, loc(B, i-1), loc(B, i), loc(B, i+1))

step(B::RestrictedFiniteDifferences) = step(parent(B))
local_step(B::RestrictedFiniteDifferences, i) = local_step(parent(B), j+first(indices(B,2))-1)

function getindex(B::AbstractFiniteDifferences{T}, x::AbstractVector, sel::AbstractVector) where T
    l = locs(B)
    χ = spzeros(T, length(x), length(sel))
    o = sel[1] - 1

    n = length(l)
    j = sel[1]
    b = l[j]
    a = j > 1 ? l[j-1] : 2b - l[j+1]
    for j in sel
        c = j < n ? l[j+1] : 2b - l[j-1]
        x₀ = b
        J = j - o
        δx = local_step(B, j)
        w = weight(B, j)
        for i ∈ findall(in(a..c), x)
            χ[i,J] = w*tent(x[i], a, b, c)
        end
        a,b = b,c
    end
    χ
end

getindex(B::FiniteDifferencesOrRestricted, x::AbstractVector, ::Colon) =
    getindex(parent(B), x, indices(B,2))

getindex(B::RestrictedFiniteDifferences, x::AbstractVector, sel::AbstractVector) =
    getindex(parent(B), x, indices(B,2)[sel])

vandermonde(::Uniform, B::AbstractFiniteDifferences{T}) where T =
    BandedMatrices._BandedMatrix(Ones{T}(1,size(B,2)), axes(B,2), 0, 0)

vandermonde(::NonUniform, B::AbstractFiniteDifferences) =
    BandedMatrix(0 => weights(B))

vandermonde(B::AbstractFiniteDifferences) = vandermonde(distribution(B), B)

metric_shape(B::FiniteDifferencesOrRestricted) = Diagonal

# * Types

const FDArray{T,N,B<:FiniteDifferencesOrRestricted} = FuncArray{T,N,B}
const FDVector{T,B<:FiniteDifferencesOrRestricted} = FuncVector{T,B}
const FDMatrix{T,B<:FiniteDifferencesOrRestricted} = FuncMatrix{T,B}
const FDVecOrMat{T,B<:FiniteDifferencesOrRestricted} = FuncVecOrMat{T,B}

const AdjointFDArray{T,N,B<:FiniteDifferencesOrRestricted} = AdjointFuncArray{T,N,B}
const AdjointFDVector{T,B<:FiniteDifferencesOrRestricted} = AdjointFuncVector{T,B}
const AdjointFDMatrix{T,B<:FiniteDifferencesOrRestricted} = AdjointFuncMatrix{T,B}
const AdjointFDVecOrMat{T,B<:FiniteDifferencesOrRestricted} = AdjointFuncVecOrMat{T,B}

# This is an operator on the form O = B*M*B⁻¹, where M is a matrix
# acting on the expansion coefficients. When applied to e.g. a
# FDVector, v = B*c, we get O*v = B*M*B⁻¹*B*c = B*M*c. Acting on an
# AdjointFDVector is also trivial: v'O' = (B*c)'*(B*M*B⁻¹)' =
# B'c'M'. It is less clear if it is possible to form Hermitian
# operator that does _not_ need be transposed before application to
# adjoint vectors.
const FDOperator{T,B<:FiniteDifferencesOrRestricted,M<:AbstractMatrix} = MulQuasiArray{T,2,<:Tuple{<:B,<:M,
                                                                                                   <:PInvQuasiMatrix{<:Any,<:Tuple{<:B}}}}

const FDMatrixElement{T,B<:FiniteDifferencesOrRestricted,M<:AbstractMatrix,V<:AbstractVector} =
    MulQuasiArray{T,0,<:Mul{<:Any,
                            <:Tuple{<:Adjoint{<:Any,<:V},<:QuasiAdjoint{<:Any,<:B},
                                    <:B,<:M,<:QuasiAdjoint{<:Any,<:B},
                                    <:B,<:V}}}

# * Mass matrix

@simplify function *(Ac::QuasiAdjoint{<:Any,<:FiniteDifferencesOrRestricted}, B::FiniteDifferencesOrRestricted)
    step(B)*combined_restriction(parent(Ac),B)
end

# * Basis inverses

@simplify function *(A⁻¹::PInvQuasiMatrix{<:Any,<:Tuple{<:AbstractFiniteDifferences}}, B::AbstractFiniteDifferences)
    A = parent(A⁻¹)
    A == B || throw(ArgumentError("Cannot multiply basis with inverse of other basis"))
    one(eltype(A))*I
end

@simplify function *(A::AbstractFiniteDifferences, B⁻¹::PInvQuasiMatrix{<:Any,<:Tuple{<:AbstractFiniteDifferences}})
    B = parent(B⁻¹)
    A == B || throw(ArgumentError("Cannot multiply basis with inverse of other basis"))
    one(eltype(A))*I
end

@simplify function *(A⁻¹::PInvQuasiMatrix{<:Any,<:Tuple{<:AbstractFiniteDifferences}}, v::FDArray)
    A = parent(A⁻¹)
    B,c = v.args
    A == B || throw(ArgumentError("Cannot multiply basis with inverse of other basis"))
    c
end

@simplify function *(v::AdjointFDArray, B⁻¹::PInvQuasiMatrix{<:Any,<:Tuple{<:AbstractFiniteDifferences}})
    c,Ac = v.args
    B = parent(B⁻¹)
    parent(Ac) == B || throw(ArgumentError("Cannot multiply basis with inverse of other basis"))
    c
end

# * Various finite differences
# ** Cartesian finite differences
struct FiniteDifferences{T,R<:AbstractRange{T}} <: AbstractFiniteDifferences{T}
    r::R
end
FiniteDifferences(n::Integer, Δx::Real) = FiniteDifferences((1:n)*Δx)

locs(B::FiniteDifferences) = B.r
step(B::FiniteDifferences) = step(B.r)

IntervalSets.leftendpoint(B::FiniteDifferences) = B.r[1] - step(B)
IntervalSets.rightendpoint(B::FiniteDifferences) = B.r[end] + step(B)

show(io::IO, B::FiniteDifferences{T}) where {T} =
    write(io, "Finite differences basis {$T} on $(axes(B,1).domain) with $(size(B,2)) points spaced by Δx = $(step(B))")

# ** Staggered finite differences

function log_lin_grid!(r, ρₘᵢₙ, ρₘₐₓ, α, js)
    0 < ρₘᵢₙ < Inf && 0 < ρₘₐₓ < Inf ||
        throw(DomainError((ρₘᵢₙ, ρₘₐₓ), "Log–linear grid steps need to be finite and larger than zero"))
    δρ = ρₘₐₓ-ρₘᵢₙ
    for j = js
        r[j] = r[j-1] + ρₘᵢₙ + (1-exp(-α*r[j-1]))*δρ
    end
end

function log_lin_grid(ρₘᵢₙ, ρₘₐₓ, α, n::Integer)
    r = zeros(n)
    r[1] = ρₘᵢₙ/2
    log_lin_grid!(r, ρₘᵢₙ, ρₘₐₓ, α, 2:n)
    r
end

function log_lin_grid(ρₘᵢₙ, ρₘₐₓ, α, rₘₐₓ)
    # Guess n
    n = ceil(Int, rₘₐₓ/ρₘₐₓ)
    r = zeros(n)
    r[1] = ρₘᵢₙ/2
    log_lin_grid!(r, ρₘᵢₙ, ρₘₐₓ, α, 2:n)
    nprev = n
    while r[end] < rₘₐₓ
        n *= 2
        resize!(r, n)
        log_lin_grid!(r, ρₘᵢₙ, ρₘₐₓ, α, nprev+1:n)
        nprev = n
    end
    j = findlast(<(rₘₐₓ), r)
    r[1:j+1]
end

"""
    StaggeredFiniteDifferences

Staggered finite differences with variationally derived three-points
stencils for the case where there is Dirichlet0 boundary condition at
r = 0. Supports non-uniform grids, c.f.

- Krause, J. L., & Schafer, K. J. (1999). Control of THz Emission from
  Stark Wave Packets. The Journal of Physical Chemistry A, 103(49),
  10118–10125. http://dx.doi.org/10.1021/jp992144

"""
struct StaggeredFiniteDifferences{T,V<:AbstractVector{T}} <: AbstractFiniteDifferences{T}
    r::V

    # Uniform case
    function StaggeredFiniteDifferences(r::V) where {T,V<:AbstractRange{T}}
        issorted(r) ||
            throw(ArgumentError("Node locations must be non-decreasing"))
        2r[1] ≈ step(r) || throw(ArgumentError("First node location must be half grid step"))

        new{T,V}(r)
    end

    # Arbitrary case
    function StaggeredFiniteDifferences(r::V) where {T,V<:AbstractVector{T}}
        issorted(r) ||
            throw(ArgumentError("Node locations must be non-decreasing"))

        new{T,V}(r)
    end

    # Constructor of uniform grid
    StaggeredFiniteDifferences(n::Integer, ρ) =
        StaggeredFiniteDifferences(ρ*(1:n) .- ρ/2)

    """
    StaggeredFiniteDifferences(rₘₐₓ::T, n::I, args...)

Convenience constructor for [`StaggeredFiniteDifferences`](@ref) covering the
open interval `(0,rₘₐₓ)` with `n` grid points.
"""
    StaggeredFiniteDifferences(rₘₐₓ::T, n::I) where {I<:Integer, T} =
        StaggeredFiniteDifferences(n, rₘₐₓ/(n-one(T)/2))

    # Constructor of log-linear grid
    StaggeredFiniteDifferences(ρₘᵢₙ, ρₘₐₓ, α, n::Integer) =
        StaggeredFiniteDifferences(log_lin_grid(ρₘᵢₙ, ρₘₐₓ, α, n::Integer))

    StaggeredFiniteDifferences(ρₘᵢₙ::T, ρₘₐₓ::T, α::T, rₘₐₓ::T) where {T<:Real} =
        StaggeredFiniteDifferences(log_lin_grid(ρₘᵢₙ, ρₘₐₓ, α, rₘₐₓ))
end

locs(B::StaggeredFiniteDifferences) = B.r
local_step(::Uniform, B::StaggeredFiniteDifferences, _) = step(B.r)
function local_step(::NonUniform, B::StaggeredFiniteDifferences, j)
    r = B.r
    if j == 1
        (r[2] - r[1])
    elseif j == length(r)
        r[j] - r[j-1]
    else
        (r[j+1]-r[j-1])/2
    end
end
local_step(B::StaggeredFiniteDifferences, j) = local_step(distribution(B), B, j)

IntervalSets.leftendpoint(B::StaggeredFiniteDifferences{T}) where T = zero(T)
IntervalSets.rightendpoint(B::StaggeredFiniteDifferences{T}) where T = 2B.r[end]-B.r[end-1]

function show(io::IO, B::StaggeredFiniteDifferences{T}) where T
    ρ = Base.diff(B.r)
    mi,ma = extrema(ρ)
    write(io, "Staggered finite differences basis {$T} on $(axes(B,1).domain) with $(size(B,2)) points with spacing varying from $(mi) to $(ma)")
end

show(io::IO, B::StaggeredFiniteDifferences{T,<:AbstractRange}) where T =
    write(io, "Staggered finite differences basis {$T} on $(axes(B,1).domain) with $(size(B,2)) points spaced by ρ = $(step(B))")

# ** Implicit finite differences
struct ImplicitFiniteDifferences{T,V<:AbstractVector{T}} <: AbstractFiniteDifferences{T}
    r::V

    # Uniform case
    function ImplicitFiniteDifferences(r::V) where {T,V<:AbstractRange{T}}
        issorted(r) ||
            throw(ArgumentError("Node locations must be non-decreasing"))

        new{T,V}(r)
    end

    # Arbitrary case
    function ImplicitFiniteDifferences(r::V) where {T,V<:AbstractVector{T}}
        issorted(r) ||
            throw(ArgumentError("Node locations must be non-decreasing"))

        new{T,V}(r)
    end

    # Constructor of uniform grid
    ImplicitFiniteDifferences(n::Integer, ρ) =
        ImplicitFiniteDifferences(ρ*(1:n))

    """
    ImplicitFiniteDifferences(rₘₐₓ::T, n::I, args...)

Convenience constructor for [`ImplicitFiniteDifferences`](@ref) covering the
open interval `(0,rₘₐₓ)` with `n` grid points.
"""
    ImplicitFiniteDifferences(rₘₐₓ::T, n::I) where {I<:Integer, T} =
        ImplicitFiniteDifferences(n, rₘₐₓ/(n+1))

    # Constructor of log-linear grid
    ImplicitFiniteDifferences(ρₘᵢₙ, ρₘₐₓ, α, n::Integer) =
        ImplicitFiniteDifferences(log_lin_grid(ρₘᵢₙ, ρₘₐₓ, α, n::Integer))

    ImplicitFiniteDifferences(ρₘᵢₙ::T, ρₘₐₓ::T, α::T, rₘₐₓ::T) where {T<:Real} =
        ImplicitFiniteDifferences(log_lin_grid(ρₘᵢₙ, ρₘₐₓ, α, rₘₐₓ))
end

locs(B::ImplicitFiniteDifferences) = B.r
IntervalSets.leftendpoint(B::ImplicitFiniteDifferences) = B.r[1] - local_step(B,1)
IntervalSets.rightendpoint(B::ImplicitFiniteDifferences) = B.r[end] + local_step(B,length(B.r))

ContinuumArrays.MemoryLayout(::Type{<:BasisOrRestricted{<:ImplicitFiniteDifferences}}) = ContinuumArrays.BasisLayout()
ContinuumArrays.MemoryLayout(::Type{<:AdjointBasisOrRestricted{<:ImplicitFiniteDifferences}}) = ContinuumArrays.AdjointBasisLayout()

show(io::IO, B::ImplicitFiniteDifferences{T}) where {T} =
    write(io, "Implicit finite differences basis {$T} on $(axes(B,1).domain) with $(size(B,2)) points spaced by Δx = $(step(B))")

# *** Implicit finite-differences derivatives

mutable struct ImplicitDerivative{T,Tri,Mat,MatFact} <: AbstractMatrix{T}
    Δ::Tri
    M::Mat
    M⁻¹::MatFact
    c::T
end
ImplicitDerivative(Δ::Tri, M::Mat, c::T=true) where {T,Tri,Mat} =
    ImplicitDerivative(Δ, M, factorize(M), c)

Base.show(io::IO, ∂::ND) where {ND<:ImplicitDerivative} =
    write(io, "$(size(∂,1))×$(size(∂,2)) $(ND)")

Base.show(io::IO, ::MIME"text/plain", ∂::ND) where {ND<:ImplicitDerivative} =
    show(io, ∂)

Base.size(∂::ND, args...) where {ND<:ImplicitDerivative} = size(∂.Δ, args...)
Base.eltype(∂::ND) where {ND<:ImplicitDerivative} = eltype(∂.Δ)
Base.axes(∂::ND, args...) where {ND<:ImplicitDerivative} = axes(∂.Δ, args...)

Base.:(*)(∂::ND,a::T) where {T<:Number,ND<:ImplicitDerivative} =
    ImplicitDerivative(∂.Δ, ∂.M, ∂.M⁻¹, ∂.c*a)

Base.:(*)(a::T,∂::ND) where {T<:Number,ND<:ImplicitDerivative} =
    ∂ * a

Base.:(/)(∂::ND,a::T) where {T<:Number,ND<:ImplicitDerivative} =
    ∂ * inv(a)

Base.complex(∂::ImplicitDerivative) =
    ImplicitDerivative(complex(∂.Δ), complex(∂.M), complex(∂.c))

function LinearAlgebra.mul!(y::Y, ∂::ND, x::X,
                            α::Number=true, β::Number=false) where {Y<:AbstractVector,
                                                                    ND<:ImplicitDerivative,
                                                                    X<:AbstractVector}
    mul!(y, ∂.Δ, x, α, β)
    ldiv!(∂.M⁻¹, y)
    lmul!(∂.c, y)
    y
end

function Base.copyto!(y::Y, ∂::Mul{<:Any,Tuple{<:ImplicitDerivative, X}}) where {X<:AbstractVector,Y<:AbstractVector}
    C = ∂.C
    mul!(C, ∂.A, ∂.B)
    lmul!(∂.α, C)
    C
end

for op in [:(+), :(-)]
    for Mat in [:Diagonal, :UniformScaling]
        @eval begin
            function Base.$op(∂::ND, B::$Mat) where {ND<:ImplicitDerivative}
                B̃ = inv(∂.c)*∂.M*B
                ImplicitDerivative($op(∂.Δ, B̃), ∂.M, ∂.M⁻¹, ∂.c)
            end
        end
    end
end

function Base.:(*)(∂::ImplicitDerivative{T}, x::AbstractVector{U}) where {T,U}
    size(∂, 2) == size(x,1) ||
        throw(DimensionMismatch("Number of rows of x $(size(x,1)) does not agree with number of columns of ∂ $(size(∂,2))"))
    mul!(similar(x, promote_type(T, U)), ∂, x)
end

function Base.:(*)(x::Adjoint{T,<:AbstractVector{T}}, ∂::ImplicitDerivative{U}) where {T,U}
    size(∂, 1) == size(x,2) ||
        throw(DimensionMismatch("Number of columns of x $(size(x,2)) does not agree with number of rows of ∂ $(size(∂,1))"))
    lmul!(∂.c, rdiv!(copy(x), ∂.M⁻¹)*∂.Δ)
end

# *** Implicit finite-differences derivative factorizations

struct ImplicitFactorization{TriFact,Mat}
    Δ⁻¹::TriFact
    M::Mat
end

Base.size(∂⁻¹::NF, args...) where {NF<:ImplicitFactorization} = size(∂⁻¹.M, args...)
Base.eltype(∂⁻¹::NF) where {NF<:ImplicitFactorization} = eltype(∂⁻¹.M)

LinearAlgebra.factorize(∂::ND) where {ND<:ImplicitDerivative} =
    ImplicitFactorization(factorize(∂.c*∂.Δ), ∂.M)

function LinearAlgebra.ldiv!(y, ∂⁻¹::NF, x) where {NF<:ImplicitFactorization}
    mul!(y, ∂⁻¹.M, x)
    ldiv!(∂⁻¹.Δ⁻¹, y)
    y
end

# * Projections

# Vandermonde interpolation for finite differences is equivalent to
# evaluating the function on the grid points, since the basis
# functions are orthogonal (in the sense of a quadrature) and there is
# no overlap between adjacent basis functions.
function Base.:(\ )(B::FiniteDifferencesOrRestricted, f::BroadcastQuasiArray)
    axes(f,1) == axes(B,1) ||
        throw(DimensionMismatch("Function on $(axes(f,1).domain) cannot be interpolated over basis on $(axes(B,1).domain)"))

    x = locs(B)
    getindex.(Ref(f), x) ./ weights(B)
end
