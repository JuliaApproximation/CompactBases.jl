# * Derivatives

function diffop!(dest::AbstractMatrix,
                 L::BSplineOrRestricted, R::BSplineOrRestricted,
                 o)
    assert_compatible_BSplines(L,R)

    k = max(order(L),order(R))
    (bandwidth(dest,1) ≤ k-1) && (bandwidth(dest,2) ≤ k-1) ||
        throw(DimensionMismatch("Insufficient bandwidths of destination matrix"))

    d,r = divrem(o,2)
    a,b = d,d+r

    # ∂' = -∂ ⇒ for weak Laplacians (in case of Dirichlet0 boundary
    # conditions), we negate the left basis.
    χ = (iseven(a) ? 1 : -1)*basis_functions(L, a)
    ξ = basis_functions(R, b)

    overlap_matrix!(dest, χ, ξ, weights(parent(R)))

    dest
end

diffop!(dest::BandedMatrix, B::BSplineOrRestricted, o) =
    diffop!(dest, B, B, o)

@materialize function *(Ac::AdjointBasisOrRestricted{<:BSpline},
                        D::Derivative,
                        B::BasisOrRestricted{<:BSpline})
    BSplineStyle
    T -> begin
        Matrix(undef, parent(Ac), B, T)
    end
    dest::AbstractMatrix{T} -> begin
        diffop!(dest, parent(Ac), B, 1)
    end
end

@materialize function *(Ac::AdjointBasisOrRestricted{<:BSpline},
                        Dc::QuasiAdjoint{<:Any,<:Derivative},
                        D::Derivative,
                        B::BasisOrRestricted{<:BSpline})
    BSplineStyle
    T -> begin
        Matrix(undef, parent(Ac), B, T)
    end
    dest::AbstractMatrix{T} -> begin
        diffop!(dest, parent(Ac), B, 2)
    end
end
