@doc raw"""
    FunctionProduct{Conjugated}

Helper object to compute the expansion coefficients of ``ρ(x) \defd
f^\circ(x)g(x)``, where ``f^\circ`` denotes that ``f`` may be
conjugated, if so desired.
"""
struct FunctionProduct{Conjugated,T,P,A,B,V,M,FT,GT,Tmp}
    ρ::P
    L::A
    LV::V
    R::B
    RV::V
    C::M
    ftmp::FT
    gtmp::GT
    tmp::Tmp
    FunctionProduct{Conjugated}(ρ::P, L::A, LV::V, R::B, RV::V, C::M,
                                ftmp::FT, gtmp::GT, tmp::Tmp) where {Conjugated,T,P<:AbstractVecOrMat{T},A,B,V,M,FT,GT,Tmp} =
                                    new{Conjugated,T,P,A,B,V,M,FT,GT,Tmp}(ρ, L, LV, R, RV, C, ftmp, gtmp, tmp)
end

"""
    Density

Type-alias for [`FunctionProduct`](@ref) where the first function is
conjugated, as is necessary in complex linear algebra, when computing
mutual densities.
"""
const Density = FunctionProduct{true}

tmpvec(::Type{T}, m, ::Nothing) where T = Vector{T}(undef, m)
tmpvec(::Type{T}, m, n) where T = Matrix{T}(undef, m, n)

@doc raw"""
    FunctionProduct{Conjugated}(L, R[, T; w=one]) where {Conjugated,T}

Construct a [`FunctionProduct`](@ref) for computing the product of two
functions expanded over `L` and `R`, respectively:

```math
h(x) = f^\circ(x)g(x)w(x)
```

where ``w(x)`` is an optional weight function, over the basis of
``g(x)``, via Vandermonde interpolation.
"""
function FunctionProduct{Conjugated}(L::BasisOrRestricted, R::BasisOrRestricted,
                                     ::Type{T}=promote_type(eltype(L), eltype(R)),
                                     nlhs=nothing, nrhs=nothing;
                                     w=one) where {Conjugated,T}
    assert_compatible_bases(L,R)

    LV = vandermonde(L)
    RV = vandermonde(R)

    rnlhs = something(nlhs,1)
    rnrhs = something(nrhs,1)
    rnlhs == rnrhs || isone(rnlhs) || isone(rnrhs) ||
        throw(DimensionMismatch("Cannot broadcast $(rnlhs) LHS and $(rnrhs) RHS to common size"))

    w = if w == one
        I
    else
        Diagonal(w.(locs(parent(R))))
    end

    nρ = (isone(rnlhs) && isone(rnrhs)
          ? nothing
          : max(rnlhs,rnrhs))

    ρρ = tmpvec(T, size(R,2), nρ)

    if LV == RV && sum(bandwidths(LV)) == 0
        # Orthogonal case
        RV = RV[indices(R,2),:]
        uniform = all(isone, RV.data)
        C = if uniform
            w
        else
            RV*(w isa AbstractMatrix ? w[indices(R,2),indices(R,2)] : w)
        end
        tmp = C === I ? nothing : tmpvec(T, size(RV, 1), nρ)
        FunctionProduct{Conjugated}(ρρ,
                                    L, I, R, I,
                                    C,
                                    nothing, nothing, tmp)
    else
        # General, non-orthogonal case
        FunctionProduct{Conjugated}(ρρ,
                                    L, LV, R, RV,
                                    pinv(Matrix(RV))*w,
                                    tmpvec(T, size(LV, 1), nlhs),
                                    tmpvec(T, size(RV, 1), nrhs),
                                    tmpvec(T, size(RV, 1), nρ))
    end
end

padding( ::AbstractVector) = nothing
padding(A::AbstractMatrix) = size(A,2)

"""
    FunctionProduct{Conjugated}(f, g) where Conjugated

Construct a [`FunctionProduct`](@ref) for computing the product of the
two functions `f` and `g`.
"""
function FunctionProduct{Conjugated}(f, g; kwargs...) where Conjugated
    T = promote_type(eltype(f),eltype(g))
    L,cf = f.args
    R,cg = g.args

    ρ = FunctionProduct{Conjugated}(L, R, T, padding(cf), padding(cg); kwargs...)

    copyto!(ρ, f, g)

    ρ
end

Base.eltype(ρ::FunctionProduct{<:Any,T}) where T = T

function Base.show(io::IO, ρ::FunctionProduct{Conjugated, T}) where {Conjugated,T}
    mL = size(ρ.L,2)
    isnothing(ρ.ftmp) || (ρ.ftmp isa AbstractMatrix) && (mL = (mL,size(ρ.ftmp,2)))
    mR = size(ρ.R,2)
    isnothing(ρ.gtmp) || (ρ.gtmp isa AbstractMatrix) && (mR = (mR,size(ρ.gtmp,2)))
    mρ = if ρ.ρ isa AbstractVector
        length(ρ.ρ)
    else
        isnothing(ρ.ftmp) && (mL = "($(mL), :)")
        isnothing(ρ.gtmp) && (mR = "($(mR), :)")
        size(ρ.ρ)
    end
    write(io, "$(mL) .* $(mR) -> $(mρ) FunctionProduct $T")
    Conjugated && write(io, ", conjugated (<=> Density)")
end

function Base.show(io::IO, ::MIME"text/plain", ρ::FunctionProduct)
    show(io, ρ)
    write(io, "; L .* R -> R, with\n  L: ")
    show(io, ρ.L)
    write(io, "\n  R: ")
    show(io, ρ.R)
end

"""
    copyto!(ρ::FunctionProduct{Conjugated}, f::AbstractVector, g::AbstractVector) where Conjugated

Update the [`FunctionProduct`](@ref) `ρ` from the vectors of expansion
coefficients, `f` and `g`.
"""
function Base.copyto!(ρ::FunctionProduct{Conjugated}, f::AbstractVecOrMat, g::AbstractVecOrMat) where Conjugated
    # Non-orthogonal
    mul!(ρ.ftmp, ρ.RV, f)
    Conjugated && conj!(ρ.ftmp)
    mul!(ρ.gtmp, ρ.LV, g)
    ρ.tmp .= ρ.ftmp .* ρ.gtmp

    mul!(ρ.ρ, ρ.C, ρ.tmp)
end

function Base.copyto!(ρ::FunctionProduct{Conjugated,<:Any,<:Any,<:Any,<:Any,UniformScaling{Bool},<:AbstractMatrix},
                      f::AbstractVecOrMat, g::AbstractVecOrMat) where Conjugated
    # Orthogonal, non-uniform
    if Conjugated
        ρ.tmp .= conj.(f) .* g
    else
        ρ.tmp .= f .* g
    end
    mul!(ρ.ρ, ρ.C, ρ.tmp)
end

function Base.copyto!(ρ::FunctionProduct{Conjugated,<:Any,<:Any,<:Any,<:Any,UniformScaling{Bool},UniformScaling{Bool}},
                      f::AbstractVecOrMat, g::AbstractVecOrMat) where Conjugated
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

# Orthogonal, uniform or non-uniform
Base.copyto!(M::Diagonal, ρ::FunctionProduct{<:Any,<:Any,<:Any,<:Any,<:Any,UniformScaling{Bool},<:Any}) =
    mul!(M.diag, ρ.C, ρ.ρ)

# Non-orthogonal (but not general case, only B-splines so far)
Base.copyto!(M::AbstractMatrix,
             ρ::FunctionProduct{<:Any,<:Any,<:Any,<:BSplineOrRestricted,<:BSplineOrRestricted}) =
                 overlap_matrix!(M, ρ.L', ρ.R, Diagonal(ρ.tmp))
