@doc raw"""
    FunctionProduct{Conjugated}

Helper object to compute the expansion coefficients of ``ρ(x) \defd
f^\circ(x)g(x)``, where ``f^\circ`` denotes that ``f`` may be
conjugated, if so desired.
"""
struct FunctionProduct{Conjugated,T,A,B,V,M,FT,GT}
    ρ::Vector{T}
    L::A
    LV::V
    R::B
    RV::V
    C::M
    ftmp::FT
    gtmp::GT
    FunctionProduct{Conjugated}(ρ::Vector{T}, L::A, LV::V, R::B, RV::V, C::M,
                                ftmp::FT, gtmp::GT) where {Conjugated,T,A,B,V,M,FT,GT} =
                                    new{Conjugated,T,A,B,V,M,FT,GT}(ρ, L, LV, R, RV, C, ftmp, gtmp)
end

"""
    Density

Type-alias for [`FunctionProduct`](@ref) where the first function is
conjugated, as is necessary in complex linear algebra, when computing
mutual densities.
"""
const Density = FunctionProduct{true}

"""
    FunctionProduct{Conjugated}(R, L[, T]) where {Conjugated,T}

Construct a [`FunctionProduct`](@ref) for computing the product of two
functions expanded over `R` and `L`, respectively.
"""
function FunctionProduct{Conjugated}(R::BasisOrRestricted, L::BasisOrRestricted,
                                     ::Type{T}=promote_type(eltype(R), eltype(L))) where {Conjugated,T}
    assert_compatible_bases(L,R)

    LV = vandermonde(L)
    RV = vandermonde(R)

    ρ = if LV == RV && sum(bandwidths(LV)) == 0
        # Orthogonal case
        RV = RV[indices(R,2),:]
        uniform = all(isone, RV.data)
        FunctionProduct{Conjugated}(Vector{T}(undef, size(R,2)),
                                    L, I, R, I,
                                    uniform ? I : RV,
                                    uniform ? nothing : Vector{T}(undef, size(RV, 1)),
                                    nothing)
    else
        # General, non-orthogonal case
        FunctionProduct{Conjugated}(Vector{T}(undef, size(R,2)),
                                    L, LV, R, RV,
                                    pinv(Matrix(RV)),
                                    Vector{T}(undef, size(LV, 1)),
                                    Vector{T}(undef, size(RV, 1)))
    end
end

"""
    FunctionProduct{Conjugated}(f, g) where Conjugated

Construct a [`FunctionProduct`](@ref) for computing the product of the
two functions `f` and `g`.
"""
function FunctionProduct{Conjugated}(f, g) where Conjugated
    T = promote_type(eltype(f),eltype(g))
    L,cf = f.args
    R,cg = f.args

    ρ = FunctionProduct{Conjugated}(L, R, T)

    copyto!(ρ, f, g)

    ρ
end

Base.eltype(ρ::FunctionProduct{<:Any,T}) where T = T

"""
    copyto!(ρ::FunctionProduct{Conjugated}, f::AbstractVector, g::AbstractVector) where Conjugated

Update the [`FunctionProduct`](@ref) `ρ` from the vectors of expansion
coefficients, `f` and `g`.
"""
function Base.copyto!(ρ::FunctionProduct{Conjugated}, f::AbstractVector, g::AbstractVector) where Conjugated
    # Non-orthogonal
    mul!(ρ.ftmp, ρ.RV, f)
    Conjugated && conj!(ρ.ftmp)
    mul!(ρ.gtmp, ρ.LV, g)
    ρ.ftmp .*= ρ.gtmp

    mul!(ρ.ρ, ρ.C, ρ.ftmp)
end

function Base.copyto!(ρ::FunctionProduct{Conjugated,<:Any,<:Any,<:Any,UniformScaling{Bool},<:AbstractMatrix},
                      f::AbstractVector, g::AbstractVector) where Conjugated
    # Orthogonal, non-uniform
    if Conjugated
        ρ.ftmp .= conj.(f) .* g
    else
        ρ.ftmp .= f .* g
    end
    mul!(ρ.ρ, ρ.C, ρ.ftmp)
end

function Base.copyto!(ρ::FunctionProduct{Conjugated,<:Any,<:Any,<:Any,UniformScaling{Bool},UniformScaling{Bool}},
                      f::AbstractVector, g::AbstractVector) where Conjugated
    # Orthogonal, uniform
    if Conjugated
        ρ.ρ .= conj.(f) .* g
    else
        ρ.ρ .= f .* g
    end
end

"""
    Base.copyto!(ρ::FunctionProduct, f, g)

Update the [`FunctionProduct`](@ref) `ρ` from the functions `f` and `g`.
"""
function Base.copyto!(ρ::FunctionProduct, f, g)
    L,cf = f.args
    @assert L == ρ.L
    R,cg = g.args
    @assert R == ρ.R

    copyto!(ρ, cf, cg)
end

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
    cρ .= conj.(ld.u) .* ld.v .* weights(ld.R)
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
