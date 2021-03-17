# * Derivatives

function diffop!(dest::AbstractMatrix,
                 L::BSplineOrRestricted, R::BSplineOrRestricted,
                 a, b, op=I)
    assert_compatible_bases(L,R)

    k = max(order(L),order(R))
    (bandwidth(dest,1) ≤ k-1) && (bandwidth(dest,2) ≤ k-1) ||
        throw(DimensionMismatch("Insufficient bandwidths of destination matrix"))

    # ∂' = -∂ ⇒ for weak Laplacians (in case of Dirichlet0 boundary
    # conditions), we negate the left basis.
    χ = (iseven(a) ? 1 : -1)*basis_functions(L, a)
    ξ = basis_functions(R, b)

    overlap_matrix!(dest, χ, op*ξ, weights(parent(R)))

    dest
end

diffop!(dest::BandedMatrix, B::BSplineOrRestricted, o) =
    diffop!(dest, B, B, o)

@materialize function *(Ac::AdjointBasisOrRestricted{<:BSpline},
                        D::Derivative,
                        B::BasisOrRestricted{<:BSpline})
    T -> begin
        Matrix(undef, parent(Ac), B, T)
    end
    dest::AbstractMatrix{T} -> begin
        diffop!(dest, parent(Ac), B, 0, 1)
    end
end

@materialize function *(Ac::AdjointBasisOrRestricted{<:BSpline},
                        Dc::QuasiAdjoint{<:Any,<:Derivative},
                        D::Derivative,
                        B::BasisOrRestricted{<:BSpline})
    T -> begin
        Matrix(undef, parent(Ac), B, T)
    end
    dest::AbstractMatrix{T} -> begin
        diffop!(dest, parent(Ac), B, 1, 1)
    end
end

@materialize function *(Ac::AdjointBasisOrRestricted{<:BSpline},
                        Dc::QuasiAdjoint{<:Any,<:Derivative},
                        O::QuasiDiagonal,
                        D::Derivative,
                        B::BasisOrRestricted{<:BSpline})
    T -> begin
        Matrix(undef, parent(Ac), B, T)
    end
    dest::AbstractMatrix{T} -> begin
        # Evaluate the quasi-diagonal operator O on the quadrature roots
        # of B.
        op = Diagonal(getindex.(Ref(O.diag), locs(parent(B))))
        diffop!(dest, parent(Ac), B, 1, 1, op)
    end
end

@materialize function *(Ac::AdjointBasisOrRestricted{<:BSpline},
                        Dc::QuasiAdjoint{<:Any,<:Derivative},
                        O::QuasiDiagonal,
                        P::QuasiDiagonal,
                        D::Derivative,
                        B::BasisOrRestricted{<:BSpline})
    T -> begin
        Matrix(undef, parent(Ac), B, T)
    end
    dest::AbstractMatrix{T} -> begin
        # Evaluate the quasi-diagonal operators O & P on the
        # quadrature roots of B.
        op = Diagonal(getindex.(Ref(O.diag), locs(parent(B))) .*
                      getindex.(Ref(P.diag), locs(parent(B))))
        diffop!(dest, parent(Ac), B, 1, 1, op)
    end
end
