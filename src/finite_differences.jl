abstract type AbstractFiniteDifferences{T,I<:Integer} <: Basis{T} end
const RestrictedFiniteDifferences{T,B<:AbstractFiniteDifferences{T}} =
    RestrictedQuasiArray{T,2,B}
const FiniteDifferencesOrRestricted = BasisOrRestricted{<:AbstractFiniteDifferences}

eltype(::AbstractFiniteDifferences{T}) where T = T

axes(B::AbstractFiniteDifferences{T}) where T = (Inclusion(leftendpoint(B)..rightendpoint(B)), Base.OneTo(length(locs(B))))
size(B::AbstractFiniteDifferences) = (ℵ₁, length(locs(B)))

locs(B::RestrictedFiniteDifferences) = locs(parent(B))[indices(B,2)]

==(A::AbstractFiniteDifferences,B::AbstractFiniteDifferences) = locs(A) == locs(B)

tent(x, x₀::T, δx::T) where T =
    one(T) - abs(clamp((x-x₀)/δx, -one(T), one(T)))

getindex(B::AbstractFiniteDifferences, x::Real, i::Integer) =
    tent(x, locs(B)[i], step(B))

step(B::RestrictedFiniteDifferences) = step(parent(B))

"""
    within_interval(x, interval)

Return the indices of the elements of `x` that lie within the given
closed `interval`.
"""
function within_interval(x::AbstractRange, interval::ClosedInterval)
    a = leftendpoint(interval)
    b = rightendpoint(interval)
    δx = step(x)
    max(ceil(Int, (a-x[1])/δx),1):min(floor(Int, (b-x[1])/δx)+1,length(x))
end

function getindex(B::AbstractFiniteDifferences{T}, x::AbstractRange, sel::AbstractVector) where T
    l = locs(B)
    δx = step(B)
    χ = spzeros(T, length(x), length(sel))
    o = sel[1] - 1
    for j in sel
        x₀ = l[j]
        J = j - o
        for i ∈ within_interval(x, x₀ - δx..x₀ + δx)
            χ[i,J] = tent(x[i], x₀, δx)
        end
    end
    χ
end

getindex(B::FiniteDifferencesOrRestricted, x::AbstractRange, ::Colon) =
    getindex(parent(B), x, indices(B,2))

getindex(B::RestrictedFiniteDifferences, x::AbstractRange, sel::AbstractVector) =
    getindex(parent(B), x, indices(B,2)[sel])

# * Types

const FDArray{T,N,B<:FiniteDifferencesOrRestricted} = MulQuasiArray{T,N,<:Tuple{<:B,<:AbstractArray{T,N}}}
const FDVector{T,B} = FDArray{T,1,B}
const FDMatrix{T,B} = FDArray{T,2,B}
const FDVecOrMat{T,B} = Union{FDVector{T,B},FDMatrix{T,B}}

const AdjointFDArray{T,N,B<:FiniteDifferencesOrRestricted} = MulQuasiArray{T,<:Any,<:Tuple{<:Adjoint{T,<:AbstractArray{T,N}},
                                                                                           <:QuasiAdjoint{T,<:B}}}
const AdjointFDVector{T,B} = AdjointFDArray{T,1,B}
const AdjointFDMatrix{T,B} = AdjointFDArray{T,2,B}
const AdjointFDVecOrMat{T,B} = Union{AdjointFDVector{T,B},AdjointFDMatrix{T,B}}

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

const FDInnerProduct{T,U,B<:FiniteDifferencesOrRestricted{U},V1<:AbstractVector{T},V2<:AbstractVector{T}} =
    Mul{<:Any, <:Tuple{<:Adjoint{<:Any,<:V1},<:QuasiAdjoint{<:Any,<:B},<:B,<:V2}}

const LazyFDInnerProduct{FD<:FiniteDifferencesOrRestricted} = Mul{<:Any,<:Tuple{
    <:Mul{<:Any, <:Tuple{
        <:Adjoint{<:Any,<:AbstractVector},
        <:QuasiAdjoint{<:Any,<:FD}}},
    <:Mul{<:Any, <:Tuple{
        <:FD,
        <:AbstractVector}}}}

# * Mass matrix

@simplify function *(Ac::QuasiAdjoint{<:Any,<:AbstractFiniteDifferences}, B::AbstractFiniteDifferences)
    A = parent(Ac)
    A == B || throw(ArgumentError("Cannot multiply functions on different grids"))
    I*step(B)
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
struct FiniteDifferences{T,I} <: AbstractFiniteDifferences{T,I}
    j::UnitRange{I}
    Δx::T
end
FiniteDifferences(j::UnitRange{I}, Δx::T) where {T<:Real,I<:Integer} =
    FiniteDifferences{T,I}(j, Δx)

FiniteDifferences(n::I, Δx::T) where {T<:Real,I<:Integer} =
    FiniteDifferences(1:n, Δx)

locs(B::FiniteDifferences) = B.j*B.Δx
step(B::FiniteDifferences{T}) where {T} = B.Δx

IntervalSets.leftendpoint(B::FiniteDifferences) = (B.j[1]-1)*step(B)
IntervalSets.rightendpoint(B::FiniteDifferences) = (B.j[end]+1)*step(B)

show(io::IO, B::FiniteDifferences{T}) where {T} =
    write(io, "Finite differences basis {$T} on $(axes(B,1).domain) with $(size(B,2)) points spaced by Δx = $(B.Δx)")

# ** Radial finite differences
# Specialized finite differences for the case where there is
# Dirichlet0 boundary condition at r = 0.
struct RadialDifferences{T,I} <: AbstractFiniteDifferences{T,I}
    j::Base.OneTo{I}
    ρ::T
    Z::T
    δβ₁::T # Correction used for bare Coulomb potentials, Eq. (22) Schafer2009

    RadialDifferences(n::I, ρ::T, Z=one(T),
                      δβ₁=Z*ρ/8 * (one(T) + Z*ρ)) where {I<:Integer, T} =
                          new{T,I}(Base.OneTo(n), ρ, convert(T,Z), convert(T,δβ₁))
end

"""
    RadialDifferences(rₘₐₓ::T, n::I, args...)

Convenience constructor for [`RadialDifferences`](@ref) covering the
open interval `(0,rₘₐₓ)` with `n` grid points.
"""
RadialDifferences(rₘₐₓ::T, n::I, args...) where {I<:Integer, T} =
    RadialDifferences(n, rₘₐₓ/(n+one(T)/2), args...)

locs(B::RadialDifferences{T}) where T = (B.j .- one(T)/2)*B.ρ
step(B::RadialDifferences{T}) where {T} = B.ρ

IntervalSets.leftendpoint(B::RadialDifferences{T}) where T = zero(T)
IntervalSets.rightendpoint(B::RadialDifferences{T}) where T = (B.j[end]+one(T)/2)*step(B)

==(A::RadialDifferences,B::RadialDifferences) = locs(A) == locs(B) && A.Z == B.Z && A.δβ₁ == B.δβ₁

show(io::IO, B::RadialDifferences{T}) where T =
    write(io, "Radial finite differences basis {$T} on $(axes(B,1).domain) with $(size(B,2)) points spaced by ρ = $(B.ρ)")


# ** Numerov finite differences
#=
This is an implementation of finite difference scheme described in

- Muller, H. G. (1999). An Efficient Propagation Scheme for the
  Time-Dependent Schrödinger equation in the Velocity Gauge. Laser
  Physics, 9(1), 138–148.

where the first derivative is approximated as

\[\partial f =
\left(1+\frac{h^2}{6}\Delta_2\right)^{-1}
\Delta_1 f
\equiv
M_1^{-1}\tilde{\Delta}_1 f,\]
where
\[M_1 \equiv
\frac{1}{6}
\begin{bmatrix}
4+\lambda' & 1 &\\
1 & 4 & 1\\
& 1 & 4 & \ddots\\
&&\ddots&\ddots\\
\end{bmatrix},\]
and
\[\tilde{\Delta}_1 \equiv
\frac{1}{2h}
\begin{bmatrix}
\lambda & 1 &\\
-1 &  & 1\\
& -1 &  & \ddots\\
&&\ddots&\ddots\\
\end{bmatrix},\]

where \(\lambda=\lambda'=\sqrt{3}-2\) for problems with a singularity
at the boundary \(r=0\) and zero otherwise; and the second derivative
as

\[\partial^2 f =
\left(1+\frac{h^2}{12}\Delta_2\right)^{-1}
\Delta_2 f
\equiv
-2M_2^{-1}\Delta_2 f,\]
where
\[M_2 \equiv
-\frac{1}{6}
\begin{bmatrix}
10-2\delta\beta_1 & 1 &\\
1 & 10 & 1\\
& 1 & 10 & \ddots\\
&&\ddots&\ddots\\
\end{bmatrix},\]
and
\[\Delta_2 \equiv
\frac{1}{h^2}
\begin{bmatrix}
-2(1+\delta\beta_1) & 1 &\\
1 & -2 & 1\\
& 1 & -2 & \ddots\\
&&\ddots&\ddots\\
\end{bmatrix},\]

where, again, \(\delta\beta_1 = -Zh[12-10Zh]^{-1}\) is a correction
introduced for problems singular at the origin.

=#

struct NumerovFiniteDifferences{T,I} <: AbstractFiniteDifferences{T,I}
    j::UnitRange{I}
    Δx::T
    # Used for radial problems with a singularity at r = 0.
    λ::T
    δβ₁::T
end

function NumerovFiniteDifferences(j::UnitRange{I}, Δx::T, singular_origin::Bool=false, Z=zero(T)) where {T<:Real,I<:Integer}
    λ,δβ₁ = if singular_origin
        first(j) == 1 ||
            throw(ArgumentError("Singular origin correction only valid when grid starts at Δx (i.e. `j[1] == 1`)"))
        # Eqs. (20,17), Muller (1999)
        (√3 - 2),(-Z*Δx/(12 - 10Z*Δx))
    else
        zero(T),zero(T)
    end
    NumerovFiniteDifferences{T,I}(j, Δx, λ, δβ₁)
end

NumerovFiniteDifferences(n::I, Δx::T, args...) where {T<:Real,I<:Integer} =
    NumerovFiniteDifferences(1:n, Δx, args...)

locs(B::NumerovFiniteDifferences) = B.j*B.Δx
step(B::NumerovFiniteDifferences{T}) where {T} = B.Δx

IntervalSets.leftendpoint(B::NumerovFiniteDifferences) = (B.j[1]-1)*step(B)
IntervalSets.rightendpoint(B::NumerovFiniteDifferences) = (B.j[end]+1)*step(B)

show(io::IO, B::NumerovFiniteDifferences{T}) where {T} =
    write(io, "Numerov finite differences basis {$T} on $(axes(B,1).domain) with $(size(B,2)) points spaced by Δx = $(B.Δx)")

mutable struct NumerovDerivative{T,Tri,Mat,MatFact} <: AbstractMatrix{T}
    Δ::Tri
    M::Mat
    M⁻¹::MatFact
    c::T
end
NumerovDerivative(Δ::Tri, M::Mat, c::T) where {T,Tri,Mat} =
    NumerovDerivative(Δ, M, factorize(M), c)

Base.show(io::IO, ∂::ND) where {ND<:NumerovDerivative} =
    write(io, "$(size(∂,1))×$(size(∂,2)) $(ND)")

Base.show(io::IO, ::MIME"text/plain", ∂::ND) where {ND<:NumerovDerivative} =
    show(io, ∂)

Base.size(∂::ND, args...) where {ND<:NumerovDerivative} = size(∂.Δ, args...)
Base.eltype(∂::ND) where {ND<:NumerovDerivative} = eltype(∂.Δ)
Base.axes(∂::ND, args...) where {ND<:NumerovDerivative} = axes(∂.Δ, args...)

Base.:(*)(∂::ND,a::T) where {T<:Number,ND<:NumerovDerivative} =
    NumerovDerivative(∂.Δ, ∂.M, ∂.M⁻¹, ∂.c*a)

Base.:(*)(a::T,∂::ND) where {T<:Number,ND<:NumerovDerivative} =
    ∂ * a

Base.:(/)(∂::ND,a::T) where {T<:Number,ND<:NumerovDerivative} =
    ∂ * inv(a)

function LinearAlgebra.mul!(y::Y, ∂::ND, x::X) where {Y<:AbstractVector,
                                                      ND<:NumerovDerivative,
                                                      X<:AbstractVector}
    mul!(y, ∂.Δ, x)
    ldiv!(∂.M⁻¹, y)
    lmul!(∂.c, y)
    y
end

function Base.copyto!(y::Y, ∂::Mul{<:Any,Tuple{<:NumerovDerivative, X}}) where {X<:AbstractVector,Y<:AbstractVector}
    C = ∂.C
    mul!(C, ∂.A, ∂.B)
    lmul!(∂.α, C)
    C
end

for op in [:(+), :(-)]
    for Mat in [:Diagonal, :Tridiagonal, :SymTridiagonal, :UniformScaling]
        @eval begin
            function Base.$op(∂::ND, B::$Mat) where {ND<:NumerovDerivative}
                B̃ = inv(∂.c)*∂.M*B
                NumerovDerivative($op(∂.Δ, B̃), ∂.M, ∂.M⁻¹, ∂.c)
            end
        end
    end
end

struct NumerovFactorization{TriFact,Mat}
    Δ⁻¹::TriFact
    M::Mat
end

Base.size(∂⁻¹::NF, args...) where {NF<:NumerovFactorization} = size(∂⁻¹.M, args...)
Base.eltype(∂⁻¹::NF) where {NF<:NumerovFactorization} = eltype(∂⁻¹.M)

LinearAlgebra.factorize(∂::ND) where {ND<:NumerovDerivative} =
    NumerovFactorization(factorize(∂.c*∂.Δ), ∂.M)

function LinearAlgebra.ldiv!(y, ∂⁻¹::NF, x) where {NF<:NumerovFactorization}
    mul!(y, ∂⁻¹.M, x)
    ldiv!(∂⁻¹.Δ⁻¹, y)
    y
end

# * Projections

# Vandermonde interpolation for finite differences is equivalent to
# evaluating the function on the grid points, since the basis
# functions are orthogonal (in the sense of a quadrature) and there is
# no overlap between adjacent basis functions.
function Base.:(\ )(B::FD, f::BroadcastQuasiArray) where {T,FD<:AbstractFiniteDifferences{T}}
    axes(f,1) == axes(B,1) ||
        throw(DimensionMismatch("Function on $(axes(f,1).domain) cannot be interpolated over basis on $(axes(B,1).domain)"))
    getindex.(Ref(f), locs(B))
end
