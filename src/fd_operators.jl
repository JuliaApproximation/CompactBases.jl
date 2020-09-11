# * Scalar operators

function similar(ADB::MulQuasiArray{<:Any,2,
                                    <:Tuple{<:AdjointBasisOrRestricted{<:AbstractFiniteDifferences},
                                            <:QuasiDiagonal,
                                            <:BasisOrRestricted{<:AbstractFiniteDifferences}}}, ::Type{T}=eltype(ADB)) where T
    Ac,D,B = ADB.args

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

function copyto!(dest::Diagonal{T}, ADB::MulQuasiArray{<:Any,2,
                                                       <:Tuple{<:AdjointBasisOrRestricted{<:AbstractFiniteDifferences},
                                                               <:QuasiDiagonal,
                                                               <:BasisOrRestricted{<:AbstractFiniteDifferences}}}) where T
    axes(dest) == axes(ADB) || throw(DimensionMismatch("axes must be same"))
    Ac,D,B = ADB.args

    A = parent(Ac)
    parent(A) == parent(B) ||
        throw(ArgumentError("Cannot multiply functions on different grids"))
    for (i,d) in enumerate(locs(B))
        dest.diag[i] = D.diag[d]
    end

    dest
end

function copyto!(dest::BandedMatrix{T}, ADB::MulQuasiArray{<:Any,2,
                                                           <:Tuple{<:AdjointBasisOrRestricted{<:AbstractFiniteDifferences},
                                                                   <:QuasiDiagonal,
                                                                   <:BasisOrRestricted{<:AbstractFiniteDifferences}}}) where T
    axes(dest) == axes(ADB) || throw(DimensionMismatch("axes must be same"))
    Ac,D,B = ADB.args

    A = parent(Ac)
    parent(A) == parent(B) ||
        throw(ArgumentError("Cannot multiply functions on different grids"))
    for (i,d) in enumerate(locs(B))
        dest.data[i] = D.diag[d]
    end

    dest
end

@simplify function *(Ac::AdjointBasisOrRestricted{<:AbstractFiniteDifferences},
                     D::QuasiDiagonal,
                     B::BasisOrRestricted{<:AbstractFiniteDifferences})
    M = ApplyQuasiArray(*, Ac, D, B)
    copyto!(similar(M), M)
end
