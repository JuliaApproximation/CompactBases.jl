# For generalized eigenvalue problems, and application of linear
# operators, we need to also apply the metric (or its inverse) to the
# coefficients, to stay in the correct space.
operator_metric(B) = B'B

# However, for uniform finite-differences, the coefficients coincide
# with the function value at the nodes, hence we should /not/ apply
# the metric.
operator_metric(::Uniform,     ::AbstractFiniteDifferences) = I
# For non-uniform finite-differences, we should, however.
operator_metric(::NonUniform, B::AbstractFiniteDifferences) = B'B

operator_metric(B::AbstractFiniteDifferences) =
    operator_metric(distribution(B), B)

# * Linear operator
struct LinearOperator{TA,TB,TT}
    A::TA
    B⁻¹::TB
    temp::TT
end

Base.size(S::LinearOperator, args...) = size(S.A, args...)
Base.eltype(S::LinearOperator) = eltype(S.A)

function LinearAlgebra.mul!(y, M::LinearOperator, x,
                            α::Number=true, β::Number=false)
    mul!(M.temp, M.A, x, α, false)
    if iszero(β)
        ldiv!(y, M.B⁻¹, M.temp)
    else
        y *= β
        ldiv!(M.B⁻¹, M.temp)
        y += M.temp
    end
end

LinearAlgebra.mul!(y, M::LinearOperator{<:Any,Nothing}, x,
                   α::Number=true, β::Number=false) =
    mul!(y, M.A, x, α, β)

Base.:(*)(M::LinearOperator, x) =
    LinearAlgebra.mul!(similar(x), M, x)

LinearOperator(A, B) =
    LinearOperator(A, factorize(B), Vector{eltype(A)}(undef, size(A,1)))

LinearOperator(A, ::UniformScaling) =
    LinearOperator(A, nothing, nothing)

LinearOperator(A, R::BasisOrRestricted) =
    LinearOperator(A, operator_metric(R))

# * Diagonal operator

struct DiagonalOperator{T,Diag,Vandermonde,B,Mat<:AbstractMatrix{T},Tmp}
    diag::Diag
    V::Vandermonde
    R::B
    A::Mat
    tmp::Tmp
end

function DiagonalOperator(f)
    T = eltype(f)
    R,c = f.args

    V = vandermonde(R)

    L = if sum(bandwidths(V)) == 0
        # Orthogonal case
        V = V[indices(R,2),:]
        uniform = all(isone, V.data)
        if uniform
            DiagonalOperator(c, I, R, Diagonal(c), nothing)
        else
            DiagonalOperator(c, Diagonal(V), R, Diagonal(similar(c, T)), nothing)
        end
    elseif R isa BSplineOrRestricted
        # Non-orthogonal case
        A = Matrix(undef, R, R, T)
        DiagonalOperator(c, V, R, A, V*c)
    else
        throw(ArgumentError("Basis $(R) not yet supported"))
    end
    copyto!(L, c)
    L
end

Base.size(S::DiagonalOperator, args...) = size(S.A, args...)
Base.eltype(S::DiagonalOperator) = eltype(S.A)

function Base.copyto!(o::DiagonalOperator{<:Any,<:Any,<:UniformScaling}, diag::AbstractVector)
    # Orthogonal, uniform
    copyto!(o.diag, diag)
    o
end

function Base.copyto!(o::DiagonalOperator{<:Any,<:Any,<:Diagonal}, diag::AbstractVector)
    # Orthogonal, non-uniform
    copyto!(o.diag, diag)
    mul!(o.A.diag, o.V, o.diag)
    o
end

function Base.copyto!(o::DiagonalOperator{<:Any,<:Any,<:Any,<:BSplineOrRestricted},
                 diag::AbstractVector)
    # Non-orthogonal (but not general case, only B-splines so far)
    copyto!(o.diag, diag)
    # Reconstruct function using Vandermonde matrix
    mul!(o.tmp, o.V, diag)
    # Compute matrix elements via quadrature (the basis function are
    # resolved on the quadrature nodes, as is the function
    # reconstruction).
    overlap_matrix!(o.A, o.R', o.R, Diagonal(o.tmp))
    o
end

LinearAlgebra.mul!(y, L::DiagonalOperator, x, α::Number=true, β::Number=false) =
    mul!(y, L.A, x, α, β)

Base.:(*)(L::DiagonalOperator, x) =
    LinearAlgebra.mul!(similar(x), L, x)

# * Shift-and-invert

struct ShiftAndInvert{TA,TB,TT}
    A⁻¹::TA
    B::TB
    temp::TT

    function ShiftAndInvert(A, B, σ::Number=0)
        A⁻¹ = factorize(A-σ*B)
        temp = Vector{eltype(A)}(undef, size(A,1))

        new{typeof(A⁻¹),typeof(B),typeof(temp)}(A⁻¹,B,temp)
    end

    function ShiftAndInvert(A, ::UniformScaling, σ::Number=0)
        A⁻¹ = factorize(A-σ*I)
        new{typeof(A⁻¹), Nothing, Nothing}(A⁻¹, nothing, nothing)
    end
end

ShiftAndInvert(A, R::BasisOrRestricted, σ::Number=0) =
    ShiftAndInvert(A, operator_metric(R), σ)

Base.size(S::ShiftAndInvert, args...) = size(S.A⁻¹, args...)
Base.eltype(S::ShiftAndInvert) = eltype(S.A⁻¹)

function LinearAlgebra.mul!(y,M::ShiftAndInvert,x)
    mul!(M.temp, M.B, x)
    ldiv!(y, M.A⁻¹, M.temp)
end

LinearAlgebra.mul!(y,M::ShiftAndInvert{<:Any,Nothing},x) =
    ldiv!(y, M.A⁻¹, x)
