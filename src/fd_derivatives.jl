# * Derivatives
# ** Three-point stencils

# α = super-/subdiagonal of second derivative
# β = -1/2 diagonal of second derivative
# γ = diagonal of first derivative
# δ = super-/subdiagonal of first derivative

α(::FiniteDifferences{T}, ::Integer) where T = one(T)
β(::FiniteDifferences{T}, ::Integer) where T = one(T)
γ(::FiniteDifferences{T}, ::Integer) where T = zero(T)
δ(::FiniteDifferences{T}, ::Integer) where T = one(T)

α(B::StaggeredFiniteDifferences,    j::Integer)         = B.α[j]
β(B::StaggeredFiniteDifferences,    j::Integer)         = B.β[j]
γ( ::StaggeredFiniteDifferences{T},  ::Integer) where T = zero(T)
δ(B::StaggeredFiniteDifferences,    j::Integer)         = B.δ[j]

α(B::StaggeredFiniteDifferences) = B.α
β(B::StaggeredFiniteDifferences) = B.β
γ(B::StaggeredFiniteDifferences{T}) where T = zeros(T, length(B.r))
δ(B::StaggeredFiniteDifferences) = B.δ

α( ::ImplicitFiniteDifferences{T},  ::Integer) where T = one(T)
β(B::ImplicitFiniteDifferences{T}, j::Integer) where T = one(T) + (j == 1 ? B.δβ₁ : zero(T))
γ(B::ImplicitFiniteDifferences{T}, j::Integer) where T = (j == 1 ? B.λ : zero(T))
δ( ::ImplicitFiniteDifferences{T},  ::Integer) where T = one(T)

α(B::AbstractFiniteDifferences) = α.(Ref(B), B.j[1:end-1])
β(B::AbstractFiniteDifferences) = β.(Ref(B), B.j)
γ(B::AbstractFiniteDifferences) = γ.(Ref(B), B.j)
δ(B::AbstractFiniteDifferences) = δ.(Ref(B), B.j[1:end-1])

for f in [:α, :δ]
    @eval begin
        function $f(B̃::RestrictedFiniteDifferences)
            j = indices(B̃,2)
            $f(parent(B̃))[j[1]:(j[end]-1)]
        end
    end
end

for f in [:β, :γ]
    @eval begin
        function $f(B̃::RestrictedFiniteDifferences)
            $f(parent(B̃))[indices(B̃,2)]
        end
    end
end

for f in [:α, :β, :γ, :δ]
    @eval begin
        function $f(Ã::BasisOrRestricted{<:AbstractFiniteDifferences},
                    B̃::BasisOrRestricted{<:AbstractFiniteDifferences},
                    band=0)
            i,i′ = extrema(indices(Ã,2))
            j,j′ = extrema(indices(B̃,2))
            k = max(i+min(0,band),j-max(0,band)):min(i′+min(0,band),j′-max(0,band))
            $f(parent(B̃))[k]
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
    T -> begin
        derivative_matrix(T, Ac, B, 1)
    end
    dest::Tridiagonal{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))

        # Central difference approximation
        μ = 2step(B)
        d = δ(B)/μ
        dest.dl .= -d
        dest.d .= γ(B)/μ
        dest.du .= d
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

LazyArrays.simplifiable(::typeof(*), ::AdjointBasisOrRestricted{B}, ::Derivative, ::BasisOrRestricted{B}) where {T,B<:ImplicitFiniteDifferences{T}} = Val(true)
LazyArrays._simplify(::typeof(*), A::AdjointBasisOrRestricted{B}, D::Derivative, C::BasisOrRestricted{B})  where {T,B<:ImplicitFiniteDifferences{T}} =
        implicit_derivative(T, ApplyQuasiArray(*,A,D,C), 1)

LazyArrays.simplifiable(::typeof(*), ::AdjointBasisOrRestricted{B}, ::QuasiAdjoint{T,<:Derivative}, ::Derivative, ::BasisOrRestricted{B}) where {T,B<:ImplicitFiniteDifferences{T}} = Val(true)
LazyArrays._simplify(::typeof(*), Ac::AdjointBasisOrRestricted{B}, Dc::QuasiAdjoint{T,<:Derivative}, D::Derivative, A::BasisOrRestricted{B}) where {T,B<:ImplicitFiniteDifferences{T}} = 
        implicit_derivative(T, ApplyQuasiArray(*,Ac,Dc,D,A), 2)
