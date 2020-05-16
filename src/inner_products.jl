# * Norms

const FDOrFEDVRVecOrMat = Union{FDVecOrMat,
                                Mul{<:Any, <:Tuple{<:FiniteDifferencesOrRestricted, <:AbstractArray}},
                                FEDVRVecOrMat,
                                Mul{<:Any, <:Tuple{<:FEDVROrRestricted, <:AbstractArray}}}

_norm(B::FiniteDifferencesOrRestricted, c::AbstractArray, p::Real=2) =
    norm(c, p)*(step(B)^(inv(p)))

_norm(R::FEDVROrRestricted, ϕ::AbstractArray, p::Real=2) = norm(ϕ, p)

LinearAlgebra.norm(v::FDOrFEDVRVecOrMat, p::Real=2) = _norm(v.args..., p)

function LinearAlgebra.normalize!(v::FDOrFEDVRVecOrMat, p::Real=2)
    v.args[2][:] /= norm(v,p)
    v
end

# * FD Inner products

function _inner_product(a::Adjoint{<:Any,<:AbstractVector}, A::QuasiAdjoint{<:Any,<:FD},
                        B::FD, b::AbstractVector) where {FD<:FiniteDifferencesOrRestricted}
    axes(A.parent) == axes(B) || throw(ArgumentError("Incompatible axes"))
    a*b*step(B)
end

LazyArrays.materialize(inner_product::FDInnerProduct{T,U,FD}) where {T,U,FD<:FiniteDifferencesOrRestricted{U}} =
    _inner_product(inner_product.args...)

function LazyArrays.materialize(inner_product::LazyFDInnerProduct{FD}) where {FD<:FiniteDifferencesOrRestricted}
    aA,Bb = inner_product.args
    _inner_product(aA.args..., Bb.args...)
end

function LazyArrays.materialize(s::Mul{<:Any, <:Tuple{
    <:Mul{<:Any, <:Tuple{
        <:Adjoint{<:Any,<:AbstractVector},
        <:ContinuumArrays.QuasiArrays.QuasiAdjoint{<:Any, <:FD}}},
    <:Mul{<:Any, <:Tuple{
        <:FD,
        <:Diagonal,
        <:ContinuumArrays.QuasiArrays.QuasiAdjoint{<:Any, <:FD}}},
    <:Mul{<:Any, <:Tuple{
        <:FD,
        <:AbstractVector}}}}) where {FD<:FiniteDifferencesOrRestricted}
    a,o,b = s.args
    axes(a.args[2].parent) == axes(o.args[1]) &&
        axes(o.args[3].parent) == axes(b.args[1]) ||
        throw(ArgumentError("Incompatible axes"))
    av = first(a.args)
    ov = o.args[2].diag
    bv = last(b.args)
    v = zero(promote_type(eltype(av),eltype(ov),eltype(bv)))
    @inbounds for i in eachindex(bv)
        v += av[i]*ov[i]*bv[i]
    end
    v*step(first(b.args))
end


# * FE-DVR Inner products

const FEDVRInnerProduct =
    Mul{<:Any, <:Tuple{<:Adjoint{<:Any,<:AbstractVecOrMat},
                       <:AdjointFEDVROrRestricted,
                       <:FEDVROrRestricted,
                       <:AbstractVecOrMat}}

const LazyFEDVRInnerProduct{B₁<:AdjointFEDVROrRestricted,B₂<:FEDVROrRestricted} = Mul{<:Any,<:Tuple{
    <:Mul{<:Any, <:Tuple{
        <:Adjoint{<:Any,<:AbstractVecOrMat},
        <:B₁}},
    <:Mul{<:Any, <:Tuple{
        <:B₂,
        <:AbstractVecOrMat}}}}

# No restrictions
function _inner_product(u::Adjoint{<:Any,<:AbstractArray}, A::QuasiAdjoint{<:Any,<:FEDVR},
                        B::FEDVR, v::AbstractArray)
    A' == B || throw(DimensionMismatch("Incompatible bases"))
    u*v
end

selview(v::AbstractVector, sel) = view(v, sel)
# Using ' would conjugate the elements, which we do not want since
# they are already conjugated.
selview(v::Adjoint{<:Any,<:AbstractVector}, sel) = transpose(view(v, sel))

selview(A::AbstractMatrix, sel) = view(A, sel, :)
selview(A::Adjoint{<:Any,<:AbstractMatrix}, sel) = view(A, :, sel)

# Implementation for various restrictions
function _inner_product(u::Adjoint{<:Any,<:AbstractVecOrMat}, a₁, b₁,
                        A::FEDVR, B::FEDVR,
                        v::AbstractVecOrMat, a₂, b₂)
    A == B || throw(DimensionMismatch("Incompatible bases"))
    n = A.n
    N = length(n)
    sel = (1+max(a₁,a₂)):(N-max(b₁,b₂))

    selview(u, sel .- a₁) * selview(v, sel .- a₂)
end

# A & B restricted
function _inner_product(u::Adjoint{<:Any,<:AbstractVecOrMat}, A::AdjointRestrictedFEDVR,
                        B::RestrictedFEDVR, v::AbstractVecOrMat)
    A′ = unrestricted_basis(A')
    B′ = unrestricted_basis(B)
    a₁,b₁ = restriction_extents(A')
    a₂,b₂ = restriction_extents(B)
    _inner_product(u, a₁, b₁, A′, B′, v, a₂, b₂)
end

# A unrestricted, B restricted
function _inner_product(u::Adjoint{<:Any,<:AbstractVecOrMat}, A::QuasiAdjoint{<:Any,<:FEDVR},
                        B::RestrictedFEDVR, v::AbstractVecOrMat)
    B′ = unrestricted_basis(B)
    a,b = restriction_extents(B)
    _inner_product(u, 0, 0, A', B′, v, a, b)
end

# A restricted, B unrestricted
function _inner_product(u::Adjoint{<:Any,<:AbstractVecOrMat}, A::AdjointRestrictedFEDVR,
                        B::FEDVR, v::AbstractVecOrMat)
    A′ = unrestricted_basis(A')
    a,b = restriction_extents(A')
    _inner_product(u, a, b, A′, B, v, 0, 0)
end

LazyArrays.materialize(inner_product::FEDVRInnerProduct) =
    _inner_product(inner_product.args...)

function LazyArrays.materialize(inner_product::LazyFEDVRInnerProduct{<:AdjointFEDVROrRestricted,<:FEDVROrRestricted})
    aA,Bb = inner_product.args
    _inner_product(aA.args..., Bb.args...)
end
