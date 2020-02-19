# * Derivatives

"""
    lagrangeder!(xⁱ, m, L′)

Calculate the derivative of the Lagrange interpolating polynomial
Lⁱₘ(x), given the roots `xⁱ`, *at* the roots, and storing the result
in `L′`.

∂ₓ Lⁱₘ(xⁱₘ,) = (xⁱₘ-xⁱₘ,)⁻¹ ∏(k≠m,m′) (xⁱₘ,-xⁱₖ)/(xⁱₘ-xⁱₖ), m≠m′,
                [δ(m,n) - δ(m,1)]/2wⁱₘ,                    m=m′

Eq. (20) Rescigno2000
"""
function lagrangeder!(xⁱ::AbstractVector, wⁱ::AbstractVector,
                      L′::AbstractMatrix)
    δ(a,b) = a == b ? 1 : 0
    n = length(xⁱ)
    for m in 1:n
        L′[m,m] = (δ(m,n)-δ(m,1))/2wⁱ[m]
        for m′ in 1:n
            m′ == m && continue
            f = 1
            for k = 1:n
                k in [m,m′] && continue
                f *= (xⁱ[m′]-xⁱ[k])/(xⁱ[m]-xⁱ[k])
            end
            L′[m,m′] = f/(xⁱ[m]-xⁱ[m′])
        end
    end
end

function diff!(D, L̃, L′, wⁱ, nⁱ)
    n = size(L′,1)
    for m = 1:n
        for m′ = 1:n
            D[m,m′] = nⁱ[m]*nⁱ[m′]*dot(L̃[m,:],wⁱ.*L′[m′,:])
        end
    end
    D
end

function diff(B::FEDVR{T}, n::Integer, i::Integer) where T
    o = B.order[i]

    # L′ contains derivatives of un-normalized basis functions at
    # the quadrature roots.
    L′ = Matrix{T}(undef, o, o)
    lagrangeder!(@elem(B,x,i), B.wⁱ[i],L′)

    # D contains ⟨ξᵢ|χⱼ′⟩ where ξᵢ = χᵢ⁽ⁿ⁻¹⁾
    D = similar(L′)
    L̃ = n == 1 ? Matrix{T}(I, size(L′)...) : -L′ # ∂ᴴ = -∂
    diff!(D, L̃, L′, B.wⁱ[i], @elem(B,n,i))
    i ≥ B.i₀ && (D ./= √(B.eiϕ)) # TODO Check if correct
    D
end

difffun(B::FEDVR, n::Integer) = i -> diff(B,n,i)
difffun(B::RestrictedFEDVR, n::Integer) = i -> diff(parent(B),n,i)

derop!(A, B::FEDVROrRestricted, n::Integer) =
    set_elements!(difffun(B,n), A, B)

const FlatFirstDerivative = Mul{<:Any, <:Tuple{
    <:QuasiAdjoint{<:Any, <:FEDVR},
    <:Derivative,
    <:FEDVR}}
const FlatRestrictedFirstDerivative = Mul{<:Any, <:Tuple{
    <:AdjointRestrictedFEDVR,
    <:Derivative,
    <:RestrictedFEDVR}}

const LazyFirstDerivative = Mul{<:Any, <:Tuple{
    <:Mul{<:Any, <:Tuple{
        <:QuasiAdjoint{<:Any, <:FEDVR},
        <:Derivative}},
    <:FEDVR}}

const LazyRestrictedFirstDerivative = Mul{<:Any, <:Tuple{
    <:Mul{<:Any,<:Tuple{
        <:AdjointRestrictedFEDVR,
        <:Derivative}},
    <:RestrictedFEDVR}}

const FirstDerivative = Union{FlatFirstDerivative, FlatRestrictedFirstDerivative,
                              LazyFirstDerivative, LazyRestrictedFirstDerivative}

const FlatSecondDerivative = Mul{<:Any, <:Tuple{
    <:QuasiAdjoint{<:Any, <:FEDVR},
    <:QuasiAdjoint{<:Any, <:Derivative},
    <:Derivative,
    <:FEDVR}}
const FlatRestrictedSecondDerivative = Mul{<:Any, <:Tuple{
    <:AdjointRestrictedFEDVR,
    <:QuasiAdjoint{<:Any, <:Derivative},
    <:Derivative,
    <:RestrictedFEDVR}}

const LazySecondDerivative = Mul{<:Any, <:Tuple{
    <:Mul{<:Any, <:Tuple{
        <:Mul{<:Any, <:Tuple{
            <:QuasiAdjoint{<:Any, <:FEDVR}, <:QuasiAdjoint{<:Any, <:Derivative}}},
        <:Derivative}},
    <:FEDVR}}

const LazyRestrictedSecondDerivative = Mul{<:Any, <:Tuple{
    <:Mul{<:Any,<:Tuple{
        <:Mul{<:Any,<:Tuple{
            <:AdjointRestrictedFEDVR,
            <:QuasiAdjoint{<:Any,<:Derivative}}},
        <:Derivative}},
    <:RestrictedFEDVR}}

const SecondDerivative = Union{FlatSecondDerivative,FlatRestrictedSecondDerivative,
                               LazySecondDerivative,LazyRestrictedSecondDerivative}

const FirstOrSecondDerivative = Union{FirstDerivative,SecondDerivative}

difforder(::FirstDerivative) = 1
difforder(::SecondDerivative) = 2

function copyto!(dest::Union{Tridiagonal,BlockSkylineMatrix}, M::FirstOrSecondDerivative)
    axes(dest) == axes(M) || throw(DimensionMismatch("axes must be same"))
    derop!(dest, basis(M), difforder(M))
    dest
end

basis(M::FirstOrSecondDerivative) = last(M.args)

similar(M::FirstOrSecondDerivative, ::Type{T}) where T = Matrix(undef, basis(M))
materialize(M::FirstOrSecondDerivative) = copyto!(similar(M, eltype(M)), M)
