# * Derivatives
# ** Three-point stencils
α(::FiniteDifferences{T}, ::Integer) where T = one(T)
β(::FiniteDifferences{T}, ::Integer) where T = one(T)
γ(::FiniteDifferences{T}, ::Integer) where T = zero(T)

# Three-point stencil with variationally derived coefficients for
# singularity at r=0.
α( ::RadialDifferences,    j::Integer)         = j^2/(j^2 - 1/4)
β(B::RadialDifferences{T}, j::Integer) where T = (j^2 - j + 1/2)/(j^2 - j + 1/4) + (j == 1 ? B.δβ₁ : zero(T))
γ( ::RadialDifferences{T},  ::Integer) where T = zero(T)

α( ::ImplicitFiniteDifferences{T},  ::Integer) where T = one(T)
β(B::ImplicitFiniteDifferences{T}, j::Integer) where T = one(T) + (j == 1 ? B.δβ₁ : zero(T))
γ(B::ImplicitFiniteDifferences{T}, j::Integer) where T = (j == 1 ? B.λ : zero(T))

α(B::AbstractFiniteDifferences) = α.(Ref(B), B.j[1:end-1])
β(B::AbstractFiniteDifferences) = β.(Ref(B), B.j)
γ(B::AbstractFiniteDifferences) = γ.(Ref(B), B.j)

function α(B̃::RestrictedFiniteDifferences)
    B = parent(B̃)
    j = indices(B̃,2)
    α.(Ref(B), B.j[j[1]:(j[end]-1)])
end

for f in [:β, :γ]
    @eval begin
        function $f(B̃::RestrictedFiniteDifferences)
            B = parent(B̃)
            $f.(Ref(B), B.j[indices(B̃,2)])
        end
    end
end

for f in [:α, :β, :γ]
    @eval begin
        function $f(Ã::BasisOrRestricted{<:AbstractFiniteDifferences},
                    B̃::BasisOrRestricted{<:AbstractFiniteDifferences},
                    band=0)
            i,i′ = extrema(indices(Ã,2))
            j,j′ = extrema(indices(B̃,2))
            k = max(i+min(0,band),j-max(0,band)):min(i′+min(0,band),j′-max(0,band))
            B = parent(B̃)
            $f.(Ref(B), B.j[k])
        end
    end
end

# ** Materialization

function derivative_matrix(::Type{T}, Ac, B, diff_order) where T
    A = parent(Ac)
    parent(A) == parent(B) ||
        throw(ArgumentError("Cannot multiply functions on different grids"))
    Ai,Bi = last(indices(A)), last(indices(B))
    if Ai == Bi
        n = length(Ai)
        if diff_order == 1
            Tridiagonal(Vector{T}(undef, n-1),
                        Vector{T}(undef, n),
                        Vector{T}(undef, n-1))
        elseif diff_order == 2
            SymTridiagonal(Vector{T}(undef, n),
                           Vector{T}(undef, n-1))
        end
    else
        m,n = length(Ai),length(Bi)
        offset = Ai[1]-Bi[1]
        BandedMatrix{T}(undef, (m,n), (1-offset,1+offset))
    end
end

@materialize function *(Ac::AdjointBasisOrRestricted{<:AbstractFiniteDifferences},
                        D::Derivative,
                        B::BasisOrRestricted{<:AbstractFiniteDifferences})
    FiniteDifferencesStyle
    T -> begin
        derivative_matrix(T, Ac, B, 1)
    end
    dest::Tridiagonal{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))

        # Central difference approximation
        μ = 2step(B)
        a = α(B)/μ
        dest.dl .= -a
        dest.d .= γ(B)/μ
        dest.du .= a
    end
    dest::BandedMatrix{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))

        dl,d,du = bandrange(dest)
        # Central difference approximation
        μ = 2step(B)
        dest[Band(dl)] .= -α(A,B,-1)/μ
        dest[Band(d)] .= γ(A,B)/μ
        dest[Band(du)] .= α(A,B,1)/μ
    end
end

@materialize function *(Ac::AdjointBasisOrRestricted{<:AbstractFiniteDifferences},
                        Dc::QuasiAdjoint{<:Any,<:Derivative},
                        D::Derivative,
                        B::BasisOrRestricted{<:AbstractFiniteDifferences})
    FiniteDifferencesStyle
    T -> begin
        derivative_matrix(T, Ac, B, 2)
    end
    dest::SymTridiagonal{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))

        δ² = step(B)^2
        dest.dv .= -2β(B)/δ²
        dest.ev .= α(B)/δ²
    end
    dest::BandedMatrix{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))

        dl,d,du = bandrange(dest)
        δ² = step(B)^2
        dest[Band(dl)] .= α(A,B,-1)/δ²
        dest[Band(d)] .= -2β(A,B)/δ²
        dest[Band(du)] .= α(A,B,1)/δ²
    end
end

# # *** Implicit derivatives

function implicit_lhs(B::BasisOrRestricted{<:ImplicitFiniteDifferences{T}}, difforder) where T
    m = size(B,2)
    M = SymTridiagonal(Vector{T}(undef, m), Vector{T}(undef, m-1))

    j₁ = first(indices(B, 2))

    if difforder == 1
        M.dv .= 4
        M.ev .= 1

        # Eq. (20ff), Muller (1999), λ = λ′ = √3 - 2, but only if
        # `singular_origin==true`.
        j₁ == 1 && (M.dv[1] += parent(B).λ)

        M /= 6
    elseif difforder == 2
        M.dv .= 10
        M.ev .= 1

        # Eq. (17), Muller (1999)
        j₁ == 1 && (M.dv[1] -= 2parent(B).δβ₁)

        M /= -6
    end
end

# f'  =   M₁⁻¹*Δ₁*f   Eq. (19)
# f'' = -2M₂⁻¹*Δ₂*f   Eq. (13)
implicit_coeff(::Type{T}, difforder) where T =
    [one(T), -2one(T)][difforder]

implicit_derivative(::Type{T}, M, difforder) where T =
    ImplicitDerivative(copyto!(similar(M, T), M),
                       implicit_lhs(last(M.args), difforder),
                       implicit_coeff(T, difforder))

materialize(M::Mul{FiniteDifferencesStyle, <:Tuple{
    <:AdjointBasisOrRestricted{B},
    <:Derivative,
    <:BasisOrRestricted{B}}}) where {T,B<:ImplicitFiniteDifferences{T}} =
        implicit_derivative(T, M, 1)

materialize(M::Mul{FiniteDifferencesStyle, <:Tuple{
    <:AdjointBasisOrRestricted{B},
    <:QuasiAdjoint{T,<:Derivative},
    <:Derivative,
    <:BasisOrRestricted{B}}}) where {T,B<:ImplicitFiniteDifferences{T}} =
        implicit_derivative(T, M, 2)
