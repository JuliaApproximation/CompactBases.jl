import LinearAlgebra: Matrix, AbstractMatrix, Array,
    transpose, adjoint, diag, mul!, _modify!, MulAddMul
import Base: size, similar,
    +, -, *, /, ==,
    conj, copy, real, imag,
    require_one_based_indexing, getindex

struct SkewTridiagonal{T,V<:AbstractVector{T}} <: AbstractMatrix{T}
    ev::V                        # subdiagonal
    function SkewTridiagonal{T,V}(ev) where {T,V<:AbstractVector{T}}
        require_one_based_indexing(ev)
        new{T,V}(ev)
    end
end

SkewTridiagonal(ev::V) where {T,V<:AbstractVector{T}} = SkewTridiagonal{T}(ev)
SkewTridiagonal{T}(ev::V) where {T,V<:AbstractVector{T}} = SkewTridiagonal{T,V}(ev)
SkewTridiagonal{T}(ev::AbstractVector) where {T} =
    SkewTridiagonal(convert(AbstractVector{T}, ev)::AbstractVector{T})

function SkewTridiagonal(A::AbstractMatrix)
    if all(iszero.(diag(A,0))) && diag(A,1) == -diag(A,-1)
        SkewTridiagonal(diag(A,-1))
    else
        throw(ArgumentError("matrix is not skew-symmetric; cannot convert to SkewTridiagonal"))
    end
end

SkewTridiagonal{T,V}(S::SkewTridiagonal{T,V}) where {T,V<:AbstractVector{T}} = S
SkewTridiagonal{T,V}(S::SkewTridiagonal) where {T,V<:AbstractVector{T}} =
    SkewTridiagonal(convert(V, S.ev)::V)
SkewTridiagonal{T}(S::SkewTridiagonal{T}) where {T} = S
SkewTridiagonal{T}(S::SkewTridiagonal) where {T} =
    SkewTridiagonal(convert(AbstractVector{T}, S.ev)::AbstractVector{T})
SkewTridiagonal(S::SkewTridiagonal) = S

AbstractMatrix{T}(S::SkewTridiagonal) where {T} =
    SkewTridiagonal(convert(AbstractVector{T}, S.ev)::AbstractVector{T})

function Matrix{T}(M::SkewTridiagonal) where T
    n = size(M, 1)
    Mf = zeros(T, n, n)
    @inbounds begin
        @simd for i = 1:n-1
            Mf[i+1,i] = M.ev[i]
            Mf[i,i+1] = -M.ev[i]
        end
    end
    return Mf
end
Matrix(M::SkewTridiagonal{T}) where {T} = Matrix{T}(M)
Array(M::SkewTridiagonal) = Matrix(M)

size(A::SkewTridiagonal) = (length(A.ev)+1, length(A.ev)+1)
function size(A::SkewTridiagonal, d::Integer)
    if d < 1
        throw(ArgumentError("dimension must be ≥ 1, got $d"))
    elseif d<=2
        return length(A.ev)+1
    else
        return 1
    end
end

similar(S::SkewTridiagonal, ::Type{T}) where {T} = SkewTridiagonal(similar(S.ev, T))

#Elementary operations
for func in (:conj, :copy, :real, :imag)
    @eval ($func)(M::SkewTridiagonal) = SkewTridiagonal(($func)(M.ev))
end

transpose(S::SkewTridiagonal) = -S
adjoint(S::SkewTridiagonal{<:Real}) = -S
adjoint(S::SkewTridiagonal) = Adjoint(-S)
Base.copy(S::Adjoint{<:Any,<:SkewTridiagonal}) = SkewTridiagonal(-copy.(adjoint.(S.parent.ev)))
Base.copy(S::Transpose{<:Any,<:SkewTridiagonal}) = SkewTridiagonal(-copy.(transpose.(S.parent.ev)))

function diag(M::SkewTridiagonal, n::Integer=0)
    # every branch call similar(..., ::Int) to make sure the
    # same vector type is returned independent of n
    absn = abs(n)
    if absn==1
        v = copyto!(similar(M.ev, length(M.ev)), M.ev)
        n == 1 && lmul!(-1, v)
        return v
    elseif absn <= size(M,1)
        return fill!(similar(M.ev, size(M,1)-absn), 0)
    else
        throw(ArgumentError(string("requested diagonal, $n, must be at least $(-size(M, 1)) ",
            "and at most $(size(M, 2)) for an $(size(M, 1))-by-$(size(M, 2)) matrix")))
    end
end

+(A::SkewTridiagonal, B::SkewTridiagonal) = SkewTridiagonal(A.ev+B.ev)
-(A::SkewTridiagonal, B::SkewTridiagonal) = SkewTridiagonal(A.ev-B.ev)
-(A::SkewTridiagonal) = SkewTridiagonal(-A.ev)
*(A::SkewTridiagonal, B::Number) = SkewTridiagonal(A.ev*B)
*(B::Number, A::SkewTridiagonal) = A*B
/(A::SkewTridiagonal, B::Number) = SkewTridiagonal(A.ev/B)
==(A::SkewTridiagonal, B::SkewTridiagonal) = (A.ev==B.ev)

@inline mul!(A::StridedVecOrMat, B::SkewTridiagonal, C::StridedVecOrMat,
             alpha::Number, beta::Number) =
                 _mul!(A, B, C, MulAddMul(alpha, beta))

@inline function _mul!(C::StridedVecOrMat, S::SkewTridiagonal, B::StridedVecOrMat,
                       _add::MulAddMul)
    m, n = size(B, 1), size(B, 2)
    if !(m == size(S, 1) == size(C, 1))
        throw(DimensionMismatch("A has first dimension $(size(S,1)), B has $(size(B,1)), C has $(size(C,1)) but all must match"))
    end
    if n != size(C, 2)
        throw(DimensionMismatch("second dimension of B, $n, doesn't match second dimension of C, $(size(C,2))"))
    end

    if m == 0
        return C
    elseif iszero(_add.alpha)
        return _rmul_or_fill!(C, _add.beta)
    end

    β = S.ev
    @inbounds begin
        for j = 1:n
            x₊ = B[1, j]
            x₀ = zero(x₊)
            # If m == 1 then β[1] is out of bounds
            β₀ = m > 1 ? zero(β[1]) : zero(eltype(β))
            for i = 1:m - 1
                x₋, x₀, x₊ = x₀, x₊, B[i + 1, j]
                β₋, β₀ = β₀, β[i]
                _modify!(_add, β₋*x₋ - β₀*x₊, C, (i, j))
            end
            _modify!(_add, β₀*x₀, C, (m, j))
        end
    end

    return C
end

###################
# Generic methods #
###################

## structured matrix methods ##
function Base.replace_in_print_matrix(A::SkewTridiagonal, i::Integer, j::Integer, s::AbstractString)
    i==j-1||i==j+1 ? s : Base.replace_with_centered_mark(s)
end

function getindex(A::SkewTridiagonal{T}, i::Integer, j::Integer) where T
    if !(1 <= i <= size(A,2) && 1 <= j <= size(A,2))
        throw(BoundsError(A, (i,j)))
    end
    if i == j + 1
        return A.ev[j]
    elseif i + 1 == j
        return -A.ev[i]
    else
        return zero(T)
    end
end
