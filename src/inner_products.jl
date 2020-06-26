# * Inner products

_inner_product(a, B::AbstractFiniteDifferences, b) = step(B)*(a*b)
_inner_product(a,  ::FEDVR, b) = (a*b)

# u or v is a matrix and the result should be an array
function _inner_product(u, Ac, B, v)
    A = parent(Ac)
    lsel,rsel = combined_restriction_selection(A, B)
    _inner_product(view(u, :, lsel), B, view(v, rsel, :))
end

# Left and right functions are vectors and their inner product should
# yield a scalar
function _inner_product(u::Adjoint{<:Any,<:AbstractVector},
                        Ac, B,
                        v::AbstractVector)
    A = parent(Ac)
    lsel,rsel = combined_restriction_selection(A, B)
    _inner_product(transpose(view(u, lsel)), B, view(v, rsel))
end

# vector'vector -> scalar
function _inner_product(u::Adjoint{<:Any, <:AbstractVector},
                        Ac::QuasiAdjoint{T,BB}, B::BB,
                        v::AbstractVector) where {T,BB<:BSpline{T}}
    A = parent(Ac)
    S = A === B ? B.S : (Ac*B)

    dot(u', S, v)
end

# vector'vector -> scalar
_inner_product(u::Adjoint{<:Any, <:AbstractVector},
               Ac::AdjointBSplineOrRestricted, B::BSplineOrRestricted,
               v::AbstractVector) =
    dot(u', Ac*B, v)

# u or v is a matrix and the result should be an array
function _inner_product(u::Adjoint{<:Any, <:AbstractVecOrMat},
                        Ac::QuasiAdjoint{T,BB}, B::BB,
                        v::AbstractVecOrMat) where {T,BB<:BSpline{T}}
    A = parent(Ac)
    S = A === B ? B.S : (Ac*B)

    u*(S*v)
end

const FlatFuncInnerProduct{T,B} = Tuple{<:Adjoint{T,<:AbstractVecOrMat{T}},
                                        <:AdjointBasisOrRestricted{B},
                                        <:BasisOrRestricted{B},
                                        <:AbstractVecOrMat{T}}

materialize(M::Mul{<:Any, <:FlatFuncInnerProduct}) = _inner_product(M.args...)

materialize(M::Mul{<:Any, Tuple{AdjointFuncVector{<:Any,BB},
                                FuncVector{<:Any,BB}}}) where BB =
                                    _inner_product(M.args[1].args..., M.args[2].args...)

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

LinearAlgebra.norm(v::SplineVector) = √(v'v)
function LinearAlgebra.normalize!(v::SplineVector)
    v.args[2][:] /= norm(v)
    v
end
