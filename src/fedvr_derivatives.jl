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

ContinuumArrays.MemoryLayout(::Type{<:BasisOrRestricted{<:FEDVR}}) = ContinuumArrays.BasisLayout()
ContinuumArrays.MemoryLayout(::Type{<:AdjointBasisOrRestricted{<:FEDVR}}) = ContinuumArrays.AdjointBasisLayout()

@materialize function *(Ac::AdjointBasisOrRestricted{<:FEDVR},
                        D::Derivative,
                        B::BasisOrRestricted{<:FEDVR})
    T -> begin
        Matrix(undef, B, T)
    end
    dest::AbstractMatrix{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))

        derop!(dest, B, 1)
    end
end

@materialize function *(Ac::AdjointBasisOrRestricted{<:FEDVR},
                        Dc::QuasiAdjoint{<:Any,<:Derivative},
                        D::Derivative,
                        B::BasisOrRestricted{<:FEDVR})
    T -> begin
        Matrix(undef, B, T)
    end
    dest::AbstractMatrix{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))

        derop!(dest, B, 2)
    end
end
