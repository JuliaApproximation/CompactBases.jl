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

function Base.show(io::IO, ρ::FunctionProduct{Conjugated, T}) where {Conjugated,T}
    mL = size(ρ.L,2)
    mR = size(ρ.R,2)
    mρ = length(ρ.ρ)
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
