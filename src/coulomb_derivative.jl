"""
    CoulombDerivative(D, Z, ℓ)

Helper operator wrapping the `Derivative` operator in the case of a
Coulomb potential of charge `Z` and centrifugal potential
``ℓ(ℓ+1)/2r²``. Its main use is to correctly apply boundary conditions
at ``r=0``, when materializing derivatives.
"""
struct CoulombDerivative{T,Der<:Derivative{T},U<:Number} <: LazyQuasiMatrix{T}
    D::Der
    Z::U
    ℓ::Int
end

axes(D::CoulombDerivative) = axes(D.D)
==(a::CoulombDerivative, b::CoulombDerivative) = a.D == b.D && a.Z == b.Z && a.ℓ == b.ℓ
copy(D::CoulombDerivative) = CoulombDerivative(copy(D.D), D.Z, D.ℓ)

Base.show(io::IO, D::CoulombDerivative) =
    write(io, "CoulombDerivative($(D.D), Z = $(D.Z), ℓ = $(D.ℓ))")

export CoulombDerivative
