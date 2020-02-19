"""
    UnitVector{T}(N, k)

Helper vector type of length `N` where the `k`th element is `one(T)`
and all the others `zero(T)`.
"""
struct UnitVector{T} <: AbstractVector{T}
    N::Int
    k::Int
end

Base.size(e::UnitVector) = (e.N,)
Base.getindex(e::UnitVector{T}, i::Int) where T = i == e.k ? one(T) : zero(T)
Base.show(io::IO, e::UnitVector{T}) where T = write(io, "ê{$T}($(e.k)|$(e.N))")
