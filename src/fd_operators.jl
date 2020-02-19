# * Scalar operators

import ContinuumArrays: AbstractQuasiArrayApplyStyle
struct FiniteDifferencesStyle <: AbstractQuasiArrayApplyStyle end

@materialize function *(Ac::AdjointBasisOrRestricted{<:AbstractFiniteDifferences},
                        D::QuasiDiagonal,
                        B::BasisOrRestricted{<:AbstractFiniteDifferences})
    FiniteDifferencesStyle
    T -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))
        Ai,Bi = indices(A,2), indices(B,2)
        if Ai == Bi
            Diagonal(Vector{T}(undef, length(Ai)))
        else
            m,n = length(Ai),length(Bi)
            offset = Ai[1]-Bi[1]
            BandedMatrix{T}(undef, (m,n), (-offset,offset))
        end
    end
    dest::Diagonal{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))
        for (i,d) in enumerate(locs(B))
            dest.diag[i] = D.diag[d]
        end
    end
    dest::BandedMatrix{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))
        for (i,d) in enumerate(locs(B))
            dest.data[i] = D.diag[d]
        end
    end
end
