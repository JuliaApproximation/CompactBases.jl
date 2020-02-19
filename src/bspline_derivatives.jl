# * Derivative types

# ** First derivatives

const FlatFirstDerivative = Mul{<:Any, <:Tuple{
    <:QuasiAdjoint{<:Any, <:BSpline},
    <:Derivative,
    <:BSpline}}
const FlatRestrictedFirstDerivative = Mul{<:Any, <:Tuple{
    <:Adjoint{<:Any,<:RestrictionMatrix},
    <:QuasiAdjoint{<:Any, <:BSpline},
    <:Derivative,
    <:BSpline,
    <:RestrictionMatrix}}

const LazyFirstDerivative = Mul{<:Any, <:Tuple{
    <:Mul{<:Any, <:Tuple{
        <:QuasiAdjoint{<:Any, <:BSpline},
        <:Derivative}},
    <:BSpline}}

const LazyRestrictedFirstDerivative = Mul{<:Any, <:Tuple{
    <:Mul{<:Any,<:Tuple{
        <:MulQuasiArray{<:Any, 2, <:Mul{<:Any, <:Tuple{
            <:Adjoint{<:Any,<:RestrictionMatrix},
            <:QuasiAdjoint{<:Any,<:BSpline}}}},
        <:Derivative}},
    <:RestrictedQuasiArray{<:Any,2,<:BSpline}}}

const FirstDerivative = Union{FlatFirstDerivative, FlatRestrictedFirstDerivative,
                              LazyFirstDerivative, LazyRestrictedFirstDerivative}

# ** Second derivatives

const FlatSecondDerivative = Mul{<:Any, <:Tuple{
    <:QuasiAdjoint{<:Any, <:BSpline},
    <:QuasiAdjoint{<:Any, <:Derivative},
    <:Derivative,
    <:BSpline}}
const FlatRestrictedSecondDerivative = Mul{<:Any, <:Tuple{
    <:Adjoint{<:Any,<:RestrictionMatrix},
    <:QuasiAdjoint{<:Any, <:BSpline},
    <:QuasiAdjoint{<:Any, <:Derivative},
    <:Derivative,
    <:BSpline,
    <:RestrictionMatrix}}

const LazySecondDerivative = Mul{<:Any, <:Tuple{
    <:Mul{<:Any, <:Tuple{
        <:Mul{<:Any, <:Tuple{
            <:QuasiAdjoint{<:Any, <:BSpline}, <:QuasiAdjoint{<:Any, <:Derivative}}},
        <:Derivative}},
    <:BSpline}}

const LazyRestrictedSecondDerivative = Mul{<:Any, <:Tuple{
    <:Mul{<:Any,<:Tuple{
        <:Mul{<:Any,<:Tuple{
            <:MulQuasiArray{<:Any, 2, <:Mul{<:Any, <:Tuple{
                <:Adjoint{<:Any,<:RestrictionMatrix},
                <:QuasiAdjoint{<:Any,<:BSpline}}}},
            <:QuasiAdjoint{<:Any,<:Derivative}}},
        <:Derivative}},
    <:RestrictedQuasiArray{<:Any,2,<:BSpline}}}

const SecondDerivative = Union{FlatSecondDerivative,FlatRestrictedSecondDerivative,
                               LazySecondDerivative,LazyRestrictedSecondDerivative}

const FirstOrSecondDerivative = Union{FirstDerivative,SecondDerivative}

difforder(::FirstDerivative) = 1
difforder(::SecondDerivative) = 2

# * Materialization

function diffop!(dest::BandedMatrix,
                 L::BSplineOrRestricted, R::BSplineOrRestricted,
                 o)
    k = max(order(L),order(R))
    (bandwidth(dest,1) ≤ k-1) && (bandwidth(dest,2) ≤ k-1) ||
        throw(DimensionMismatch("Insufficient bandwidths of destination matrix"))

    d,r = divrem(o,2)
    a,b = d,d+r

    # ∂' = -∂ ⇒ for weak Laplacians, we negate the left basis.
    χ = (iseven(a) ? 1 : -1)*basis_functions(L, a)
    ξ = basis_functions(R, b)

    overlap_matrix!(dest, χ, ξ, weights(R))

    dest
end

diffop!(dest::BandedMatrix, B::BSplineOrRestricted, o) =
    diffop!(dest, B, B, o)

function copyto!(dest::BandedMatrix, M::FirstOrSecondDerivative)
    axes(dest) == axes(M) || throw(DimensionMismatch("axes must be same"))

    diffop!(dest, leftbasis(M), rightbasis(M), difforder(M))
end

leftbasis(M::Union{FlatFirstDerivative,LazyFirstDerivative,
                   FlatSecondDerivative,LazySecondDerivative}) =
                       adjoint(first(M.args))

leftbasis(M::Union{FlatRestrictedFirstDerivative, FlatRestrictedSecondDerivative}) =
    adjoint(M.args[1]*M.args[2])

rightbasis(M::Union{FlatFirstDerivative,LazyFirstDerivative,
                    FlatSecondDerivative,LazySecondDerivative}) = last(M.args)

rightbasis(M::Union{FlatRestrictedFirstDerivative, FlatRestrictedSecondDerivative}) =
    M.args[end-1]*M.args[end]

function similar(M::FirstOrSecondDerivative, ::Type{T}) where T
    k = max(order(leftbasis(M)),order(rightbasis(M)))
    BandedMatrix(Zeros{T}(size(M)), (k-1,k-1))
end
materialize(M::FirstOrSecondDerivative) = copyto!(similar(M, eltype(M)), M)
