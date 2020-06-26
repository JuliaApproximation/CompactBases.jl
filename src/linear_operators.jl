# For generalized eigenvalue problems, and application of linear
# operators, we need to also apply the metric (or its inverse) to the
# coefficients, to stay in the correct space.
operator_metric(B) = B'B

# However, for uniform finite-differences, the coefficients coincide
# with the function value at the nodes, hence we should /not/ apply
# the metric.
operator_metric(::Uniform,     ::AbstractFiniteDifferences) = I
# For non-uniform finite-differences, we should, however.
operator_metric(::NonUniform, B::AbstractFiniteDifferences) = I

operator_metric(B::AbstractFiniteDifferences) =
    operator_metric(distribution(B), B)

operator_metric(B::FEDVROrRestricted) = I

# * Linear operator
"""
    LinearOperator(A, B⁻¹, temp)

A helper object used to apply the action of a linear operator, whose
matrix representation is given by `A`, in the basis whose metric
inverse is `B⁻¹`. For orthogonal bases (such as
e.g. [`FiniteDifferences`](@ref) and [`FEDVR`](@ref)), only
multiplication by `A` is necessary, but non-orthonal bases (such as
[`BSpline`](@ref)) necessitates the application of the metric inverse as well.
"""
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
        y .*= β
        ldiv!(M.B⁻¹, M.temp)
        y .+= M.temp
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

"""
    LinearOperator(A, R)

Construct the linear operator whose matrix representation in the basis
`R` is `A`, automatically deducing if application of the metric
inverse is also necessary.
"""
LinearOperator(A, R::BasisOrRestricted) =
    LinearOperator(A, operator_metric(R))

# * Diagonal operator

"""
    DiagonalOperator

Helper structure that when acting on a function ``f(x)`` produces the
product ``f(x)g(x)``, where ``g(x)`` is e.g. a potential. This is
similar to `R'*QuasiDiagonal(g.(x))*R` in intent, but simplifies the
case when ``g(x)`` needs to be updated, without computing a lot
overlap integrals. This is accomplished by computing the
[`FunctionProduct`](@ref) of the two functions directly.
"""
struct DiagonalOperator{Diag,FunctionProduct,B}
    diag::Diag
    ρ::FunctionProduct
    R::B
end

"""
    DiagonalOperator(f)

Construct [`DiagonalOperator`](@ref) from a function expansion `f`.
"""
function DiagonalOperator(f)
    T = eltype(f)
    R,c = f.args

    DiagonalOperator(c, FunctionProduct{false}(R, R), R)
end

Base.size(S::DiagonalOperator) = (length(S.ρ.ρ),length(S.ρ.ρ))
Base.size(S::DiagonalOperator, i) = size(S)[i]
Base.eltype(S::DiagonalOperator) = eltype(S.ρ)

function copyto_if_different!(dst, src)
    dst === src && return dst
    copyto!(dst, src)
end

"""
    copyto!(o::DiagonalOperator, diag::AbstractVector)

Update the [`DiagonalOperator`](@ref) to represent multiplication by
the function whose expansion coefficients are `diag`.
"""
function Base.copyto!(o::DiagonalOperator, diag::AbstractVector)
    copyto_if_different!(o.diag, diag)
    o
end

"""
    mul!(y, L::DiagonalOperator, x[, α=1, β=0])

Compute the action of the [`DiagonalOperator`](@ref) `L` on `x` and
store the result in `y`.
"""
function LinearAlgebra.mul!(y, L::DiagonalOperator, x,
                            α::Number=true, β::Number=false)
    copyto!(L.ρ, L.diag, x)
    if iszero(β)
        y .= false
    elseif !isone(β)
        lmul!(β, y)
    end
    y .+= α .* L.ρ.ρ
end

Base.:(*)(L::DiagonalOperator, x) =
    LinearAlgebra.mul!(similar(x), L, x)

# * Shift-and-invert

"""
    ShiftAndInvert(A⁻¹, B, temp)

Represents a shifted-and-inverted operator, suitable for Krylov
iterations when targetting an interior point of the eigenspectrum.
`A⁻¹` is a factorization of the shifted operator and `B` is the metric
matrix.
"""
struct ShiftAndInvert{TA,TB,TT}
    A⁻¹::TA
    B::TB
    temp::TT

@doc raw"""
    ShiftAndInvert(A, B[, σ=0])

Construct the shifted-and-inverted operator corresponding to the
eigenproblem ``(\mat{A}-σ\mat{B})\vec{x} = \lambda\mat{B}\vec{x}``.
"""
    function ShiftAndInvert(A, B, σ::Number=0)
        A⁻¹ = factorize(A-σ*B)
        temp = Vector{eltype(A)}(undef, size(A,1))

        new{typeof(A⁻¹),typeof(B),typeof(temp)}(A⁻¹,B,temp)
    end

@doc raw"""
    ShiftAndInvert(A, ::UniformScaling[, σ=0])

Construct the shifted-and-inverted operator corresponding to the
eigenproblem ``(\mat{A}-σ\mat{I})\vec{x} = \lambda\vec{x}``.
"""
    function ShiftAndInvert(A, ::UniformScaling, σ::Number=0)
        A⁻¹ = factorize(A-σ*I)
        new{typeof(A⁻¹), Nothing, Nothing}(A⁻¹, nothing, nothing)
    end
end


@doc raw"""
    ShiftAndInvert(A, R[, σ=0])

Construct the shifted-and-inverted operator corresponding to the
eigenproblem ``(\mat{A}-σ\mat{B})\vec{x} = \lambda\mat{B}\vec{x}``,
where ``\mat{B}`` is the suitable operator metric, depending on the
`AbstractQuasiMatrix` `R` employed.
"""
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
