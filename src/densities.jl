# * FD Densities

function Base.Broadcast.broadcasted(::typeof(*), a::FDArray{T,N}, b::FDArray{T,N}) where {T,N}
    axes(a) == axes(b) || throw(DimensionMismatch("Incompatible axes"))
    A,ca = a.args
    B,cb = b.args
    A == B || throw(DimensionMismatch("Incompatible bases"))
    # We want the first MulQuasiArray to be conjugated, if complex
    c = conj.(ca) .* cb
    A*c
end

struct FDDensity{T,B<:FiniteDifferencesOrRestricted,
                 U<:AbstractVecOrMat{T},V<:AbstractVecOrMat{T}}
    R::B
    u::U
    v::V
end

function _FDDensity(Ra::FiniteDifferencesOrRestricted, ca::AbstractVecOrMat,
                    Rb::FiniteDifferencesOrRestricted, cb::AbstractVecOrMat)
    # Ra == Rb || throw(DimensionMismatch("Incompatible bases"))
    FDDensity(Ra, ca, cb)
end

function Base.copyto!(cρ::AbstractVecOrMat{T}, ld::FDDensity{T,R}, Rρ::R) where {T,R}
    # Rρ == ld.R || throw(DimensionMismatch("Incompatible bases"))
    # size(cρ) == size(ld.u) || throw(DimensionMismatch("Incompatible sizes"))
    # We want the first MulQuasiArray to be conjugated, if complex
    cρ .= conj.(ld.u) .* ld.v
    cρ
end

function Base.Broadcast.broadcasted(::typeof(⋆), a::V₁, b::V₂) where {T,B<:FiniteDifferencesOrRestricted,V₁<:FDVecOrMat{T,B},V₂<:FDVecOrMat{T,B}}
    # axes(a) == axes(b) || throw(DimensionMismatch("Incompatible axes"))
    _FDDensity(a.args..., b.args...)
end

function Base.copyto!(ρ::FDVecOrMat{T,R}, ld::FDDensity{T,R}) where {T,R}
    copyto!(ρ.args[2], ld, ρ.args[1])
    ρ
end

function Base.Broadcast.broadcasted(::typeof(⋆), a::V₁, b::V₂) where {T,B<:FiniteDifferencesOrRestricted,
                                                                      V₁<:Mul{<:Any, <:Tuple{B,<:AbstractVector{T}}},
                                                                      V₂<:Mul{<:Any, <:Tuple{B,<:AbstractVector{T}}}}
    # axes(a) == axes(b) || throw(DimensionMismatch("Incompatible axes"))
    _FDDensity(a.args..., b.args...)
end

function Base.copyto!(ρ::Mul{<:Any, <:Tuple{R,<:AbstractVector{T}}}, ld::FDDensity{T,R}) where {T,R}
    copyto!(ρ.args[2], ld, ρ.args[1])
    ρ
end

# * FE-DVR Densities

function Base.Broadcast.broadcasted(::typeof(*), a::M, b::M) where {T,N,M<:FEDVRArray{T,N}}
    axes(a) == axes(b) || throw(DimensionMismatch("Incompatible axes"))
    A,ca = a.args
    B,cb = b.args
    A == B || throw(DimensionMismatch("Incompatible bases"))
    c = similar(ca)
    a,b = restriction_extents(A)
    n = unrestricted_basis(A).n
    # We want the first MulQuasiArray to be conjugated, if complex
    @. c = conj(ca) * cb * @view(n[1+a:end-b])
    A*c
end

struct FEDVRDensity{T,B<:FEDVROrRestricted,
                    U<:AbstractVecOrMat{T},V<:AbstractVecOrMat{T}}
    R::B
    u::U
    v::V
end

function _FEDVRDensity(Ra::FEDVROrRestricted, ca::AbstractVecOrMat,
                       Rb::FEDVROrRestricted, cb::AbstractVecOrMat)
    Ra == Rb || throw(DimensionMismatch("Incompatible bases"))
    FEDVRDensity(Ra, ca, cb)
end

function Base.copyto!(cρ::AbstractVecOrMat{T}, ld::FEDVRDensity{T,R}, Rρ::R) where {T,R}
    Rρ == ld.R || throw(DimensionMismatch("Incompatible bases"))
    size(cρ) == size(ld.u) || throw(DimensionMismatch("Incompatible sizes"))
    a,b = restriction_extents(Rρ)
    n = unrestricted_basis(Rρ).n
    # We want the first MulQuasiArray to be conjugated, if complex
    @. cρ = conj(ld.u) * ld.v * @view(n[1+a:end-b])
end

function Base.Broadcast.broadcasted(::typeof(⋆), a::V, b::V) where {T,B<:FEDVR,V<:FEDVRVecOrMat{T,B}}
    axes(a) == axes(b) || throw(DimensionMismatch("Incompatible axes"))
    _FEDVRDensity(a.args..., b.args...)
end

function Base.copyto!(ρ::FEDVRVecOrMat{T,R}, ld::FEDVRDensity{T,R}) where {T,R}
    copyto!(ρ.args[2], ld, ρ.args[1])
    ρ
end

function Base.Broadcast.broadcasted(::typeof(⋆), a::V₁, b::V₂) where {T,B<:FEDVROrRestricted,
                                                                      V₁<:Mul{<:Any, <:Tuple{B,<:AbstractVector{T}}},
                                                                      V₂<:Mul{<:Any, <:Tuple{B,<:AbstractVector{T}}}}
    axes(a) == axes(b) || throw(DimensionMismatch("Incompatible axes"))
    _FEDVRDensity(a.args..., b.args...)
end

function Base.copyto!(ρ::Mul{<:Any, <:Tuple{R,<:AbstractVector{T}}}, ld::FEDVRDensity{T,R}) where {T,R}
    copyto!(ρ.args[2], ld, ρ.args[1])
    ρ
end
