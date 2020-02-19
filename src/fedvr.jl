# * Basis construction

struct FEDVR{T,R<:Real,O<:AbstractVector} <: Basis{T}
    t::AbstractVector{R}
    order::O
    i₀::Int
    t₀::R
    eiϕ::T
    x::Vector{T}
    wⁱ::Vector{Vector{R}}
    n::Vector{T}
    elems::Vector{UnitRange{Int}}
    function FEDVR(t::AbstractVector{R}, order::O; t₀::R=zero(R), ϕ::R=zero(R)) where {R<:Real,O<:AbstractVector}
        @assert length(order) == length(t)-1
        @assert all(order .> 1)

        i₀,eiϕ,T = if ϕ ≠ zero(R)
            findfirst(tt -> tt ≥ t₀, t),exp(im*ϕ),complex(R)
        else
            1,one(R),R
        end
        i₀ === nothing && throw(ArgumentError("Complex scaling starting point outside grid $(t)"))
        t₀ = t[i₀]

        x = Vector{Vector{T}}()
        w = Vector{Vector{T}}()
        for i in eachindex(order)
            xw = element_grid(order[i], t[i], t[i+1], t₀, i ≥ i₀ ? eiϕ : one(T))
            push!(x, xw[1])
            push!(w, xw[2])
        end

        rot = (i,v) -> i ≥ i₀ ? eiϕ*v : v
        n = [one(T) ./ .√(rot(i,wⁱ)) for (i,wⁱ) in enumerate(w)]
        for i in 1:length(order)-1
            n[i][end] = n[i+1][1] = one(T) ./ √(rot(i,w[i][end]) + rot(i+1,w[i+1][1]))
        end

        X = vcat(x[1], [x[i][2:end] for i in 2:length(order)]...)
        N = vcat(n[1], [n[i][2:end] for i in 2:length(order)]...)

        elems = [1:order[1]]

        l = order[1]
        for i = 2:length(order)
            l′ = l+order[i]-1
            push!(elems, l:l′)
            l = l′
        end

        new{T,R,O}(t,order,i₀,t₀,eiϕ,X,[w...],N,elems)
    end
end

FEDVR(t::AbstractVector{T},order::Integer; kwargs...) where T =
    FEDVR(t, Fill(order,length(t)-1); kwargs...)

const RestrictedFEDVR{T} = RestrictedQuasiArray{T,2,<:FEDVR{T}}
const AdjointRestrictedFEDVR{T} = AdjointRestrictedQuasiArray{T,2,<:FEDVR{T}}

const FEDVROrRestricted{T} = BasisOrRestricted{<:FEDVR{T}}
const AdjointFEDVROrRestricted{T} = AdjointBasisOrRestricted{<:FEDVR{T}}

# * Properties

axes(B::FEDVR) = (Inclusion(first(B.t)..last(B.t)), Base.OneTo(length(B.x)))
size(B::FEDVR) = (ℵ₁, length(B.x))
==(A::FEDVR,B::FEDVR) = A.t == B.t && A.order == B.order

order(B::FEDVR) = B.order
order(B::RestrictedFEDVR) = order(parent(B))[elements(B)]
order(B::FEDVROrRestricted, i) = order(parent(B))[i]

nel(B::FEDVR) = length(order(B))
element_boundaries(B::FEDVR) = vcat(1,last.(B.elems))

elements(B::FEDVR) = 1:nel(B)
function elements(B::RestrictedFEDVR)
    is=indices(B,2)
    a,b = first(is),last(is)
    elems = parent(B).elems
    i = findfirst(e -> a ∈ e, elems)
    j = findlast(e -> b ∈ e, elems)
    (isnothing(i) || isnothing(j)) && throw(ArgumentError("Invalid restriction $(a)..$(b)"))
    i:j
end

complex_rotate(x,B::FEDVR{T}) where {T<:Real} = x
complex_rotate(x,B::FEDVR{T}) where {T<:Complex} = x < B.t₀ ? x : B.t₀ + (x-B.t₀)*B.eiϕ

macro elem(B,v,i)
    :(@view($(esc(B)).$v[$(esc(B)).elems[$(esc(i))]]))
end

function show(io::IO, B::FEDVR{T}) where T
    write(io, "FEDVR{$(T)} basis with $(nel(B)) elements on $(axes(B,1).domain)")
    if T <: Complex
        rot = @printf(io, " with %s @ %.2f°", B.t₀ <= first(B.t) ? "ICS" : "ECS", rad2deg(angle(B.eiϕ)))
        B.t₀ > first(B.t) && @printf(io, " starting at %.2g", B.t₀)
    end
end

function show(io::IO, B::RestrictedFEDVR{T}) where T
    B′ = parent(B)
    is = indices(B,2)
    a,b = first(is),last(is)
    els = elements(B)
    N = length(B′.x)
    show(io, B′)
    write(io, ", restricted to elements $(els), basis functions $(a)..$(b) $(a>1 || b>N ? "⊂" : "⊆") 1..$(N)")
end

function block_orders(B::FEDVROrRestricted)
    is = indices(B,2)
    elb = element_boundaries(parent(B))
    els = elements(B)
    a,b = first(is),last(is)

    # If the restriction only includes the bridge function of the
    # first/last element covered, compute the block structure for one
    # element further in, since it simplifies.
    first(els) ≠ nel(parent(B)) && a == elb[first(els)+1] && (els = els[2:end])
    last(els) > first(els) && b == elb[last(els)] && (els = els[1:end-1])
    # Find out how much we need to "shave off" from the first and last
    # elements (inside the restriction), respectively.
    s = a-elb[first(els)]
    e = elb[last(els)+1]-b

    o = order(parent(B))[els]
    o,s,e
end

function block_structure(B::FEDVROrRestricted)
    o,s,e = block_orders(B)
    if length(o) > 1
        # First compute interior block structure; for each element,
        # the interior block and the block for the bridge function
        # connecting to /next/ element is added.
        b = o -> o > 2 ? [o-2,1] : [1]
        bs = reduce(vcat, map(b, o[2:end-1]), init=Int[])
        # Add first element.
        prepend!(bs, [o[1]-1-s,1])
        # Add last element.
        push!(bs, o[end]-1-e)
        bs
    else
        [first(o)-(s+e)]
    end
end

function block_bandwidths(B::FEDVROrRestricted, rows::Vector{<:Integer}=block_structure(B))
    nrows = length(rows)
    if nrows > 1
        bw = o -> o > 2 ? [2,1] : [1]
        # This actually generates one element "too much", but it does
        # not matter, since for the second last block, a lower
        # bandwidth of two extends outside the matrix, and that is
        # discarded.
        bws = bw.(first(block_orders(B)))
        l = vcat(1,bws[2:end]...)
        u = vcat(reverse.(bws[1:end-1])...,1)
        length(l) < nrows && (l = vcat(l,0))
        length(u) < nrows && (u = vcat(0,u))
        l,u
    else
        [0],[0]
    end
end

locs(B::FEDVR) = B.x

rlocs(B::FEDVR{<:Real}) = locs(B)

function rlocs(B::FEDVR{T}) where T
    R = real(T)
    nx = length(B.x)
    x = Vector{R}(undef, nx)
    ii = 1
    for i in eachindex(B.order)
        xw = element_grid(B.order[i], B.t[i], B.t[i+1], B.t[1], one(R))
        copyto!(view(x, ii:nx), xw[1])
        ii += length(xw[1])-1
    end
    x
end

function locs(B::RestrictedFEDVR)
    B′ = parent(B)
    a,b = restriction_extents(B)
    B′.x[1+a:end-b]
end

function rlocs(B::RestrictedFEDVR)
    B′ = parent(B)
    a,b = restriction_extents(B)
    rlocs(B′)[1+a:end-b]
end

IntervalSets.leftendpoint(B::FEDVR) = B.x[1]
IntervalSets.rightendpoint(B::FEDVR) = B.x[end]

IntervalSets.leftendpoint(B::RestrictedQuasiArray{<:Any,2,<:FEDVR}) =
    leftendpoint(parent(B))
IntervalSets.rightendpoint(B::RestrictedQuasiArray{<:Any,2,<:FEDVR}) =
    rightendpoint(parent(B))

# * Basis functions

function getindex(B::FEDVR{T}, x::Real, i::Integer, m::Integer) where T
    (x < B.t[i] || x > B.t[i+1]) && return zero(T)
    xⁱ = @elem(B,x,i)
    χ = @elem(B,n,i)[m]
    x′ = complex_rotate(x, B)
    for j = 1:B.order[i]
        j == m && continue
        χ *= (x′ - xⁱ[j])/(xⁱ[m]-xⁱ[j])
    end
    χ
end

checkbounds(B::FEDVR{T}, x::Real, k::Integer) where T =
    x ∈ axes(B,1) && 1 ≤ k ≤ size(B,2) || throw(BoundsError())

@inline function getindex(B::FEDVR{T}, x::Real, k::Integer) where T
    # @boundscheck checkbounds(B, x, k) # Slow
    i = 1
    m = k
    while i < length(B.order) && m > B.order[i]
        m -= B.order[i] - 1
        i += 1
    end
    x < B.t[i] && return zero(T)
    if x > B.t[i+1]
        if i < length(B.t)-1 && m == B.order[i]
            i += 1
            m = 1
        else
            return zero(T)
        end
    end
    B[x,i,m]
end

# * Types

const FEDVRArray{T,N,B<:FEDVROrRestricted} = MulQuasiArray{T,N,<:Tuple{B,<:AbstractArray{T,N}}}
const FEDVRVector{T,B<:FEDVROrRestricted} = FEDVRArray{T,1,B}
const FEDVRMatrix{T,B<:FEDVROrRestricted} = FEDVRArray{T,2,B}
const FEDVRVecOrMat{T,B<:FEDVROrRestricted} = Union{FEDVRVector{T,B},FEDVRMatrix{T,B}}

const AdjointFEDVRArray{T,N,B<:FEDVROrRestricted} = MulQuasiArray{T,<:Any,<:Tuple{<:Adjoint{T,<:AbstractArray{T,N}},
                                                                                  <:QuasiAdjoint{T,<:B}}}
const AdjointFEDVRVector{T,B<:FEDVROrRestricted} = AdjointFEDVRArray{T,1,B}
const AdjointFEDVRMatrix{T,B<:FEDVROrRestricted} = AdjointFEDVRArray{T,2,B}
const AdjointFEDVRVecOrMat{T,B<:FEDVROrRestricted} = Union{AdjointFEDVRVector{T,B},AdjointFEDVRMatrix{T,B}}

# * Mass matrix

@simplify function *(Ac::QuasiAdjoint{<:Any,<:FEDVR}, B::FEDVR)
    A = parent(Ac)
    A == B || throw(ArgumentError("Cannot multiply incompatible FEDVR expansions"))
    I
end

# A & B restricted
@simplify function *(Ac::AdjointRestrictedFEDVR, B::RestrictedFEDVR)
    # This is mainly for type-stability; it would be trivial to
    # generate the proper restriction matrix from the combination of
    # two differently restricted bases, but we would like to have
    # UniformScaling as the result if they are equal, and this has
    # higher priority. On the other hand, you typically only compute
    # the mass matrix in the beginning of the calculation, and thus
    # type-instability is not a big problem, so this behaviour may
    # change in the future.
    reverse(axes(Ac)) == axes(B) || throw(DimensionMismatch("axes must be same"))
    A = parent(Ac)
    parent(A) == parent(B) || throw(ArgumentError("Cannot multiply incompatible FEDVR expansions"))

    I
end

function materialize(M::Mul{<:Any,<:Tuple{<:QuasiAdjoint{T,<:FEDVR{T}},
                                          <:RestrictedFEDVR{T}}}) where T
    Ac,B = M.args
    axes(Ac,2) == axes(B,1) || throw(DimensionMismatch("axes must be same"))
    A = parent(Ac)
    A == parent(B) || throw(ArgumentError("Cannot multiply incompatible FEDVR expansions"))
    restriction(B)
end

function materialize(M::Mul{<:Any,<:Tuple{<:AdjointRestrictedFEDVR{T},
                                          <:FEDVR{T}}}) where T
    Ac,B = M.args
    axes(Ac,2) == axes(B,1) || throw(DimensionMismatch("axes must be same"))
    A = parent(Ac)
    parent(A) == B || throw(ArgumentError("Cannot multiply incompatible FEDVR expansions"))
    restriction(A)'
end

# * Basis inverses

@simplify function *(A⁻¹::PInvQuasiMatrix{<:Any,<:Tuple{<:FEDVR}},
                     B::FEDVR)
    A = parent(A⁻¹)
    A == B || throw(ArgumentError("Cannot multiply basis with inverse of other basis"))
    I
end

@simplify function *(A::FEDVR,
                     B⁻¹::PInvQuasiMatrix{<:Any,<:Tuple{<:FEDVR}})
    B = parent(B⁻¹)
    A == B || throw(ArgumentError("Cannot multiply basis with inverse of other basis"))
    I
end

@simplify function *(A⁻¹::PInvQuasiMatrix{<:Any,<:Tuple{<:FEDVR}},
                     v::FEDVRArray)
    A = parent(A⁻¹)
    B,c = v.args
    A == B || throw(ArgumentError("Cannot multiply basis with inverse of other basis"))
    c
end

@simplify function *(v::AdjointFEDVRArray,
                     B⁻¹::PInvQuasiMatrix{<:Any,<:Tuple{<:FEDVR}})
    c,Ac = v.args
    B = parent(B⁻¹)
    parent(Ac) == B || throw(ArgumentError("Cannot multiply basis with inverse of other basis"))
    c
end

# * Function interpolation

function Base.:(\ )(B::FEDVR{T}, f::BroadcastQuasiArray) where T
    axes(f,1) == axes(B,1) ||
        throw(DimensionMismatch("Function on $(axes(f,1).domain) cannot be interpolated over basis on $(axes(B,1).domain)"))

    v = zeros(T, size(B,2))
    for i ∈ 1:nel(B)
        @. v[B.elems[i]] += B.wⁱ[i]*f[@elem(B,x,i)]
    end
    v .*= B.n
    v
end

function Base.:(\ )(B::RestrictedFEDVR{T}, f::BroadcastQuasiArray) where T
    axes(f,1) == axes(B,1) ||
        throw(DimensionMismatch("Function on $(axes(f,1).domain) cannot be interpolated over basis on $(axes(B,1).domain)"))

    B′ = parent(B)
    is = indices(B,2)
    a,b = first(is),last(is)

    n = size(B′,2)
    v = zeros(T, size(B,2))
    for i ∈ elements(B)
        sel = B′.elems[i]
        # Find which basis functions of finite-element `i` should be
        # evaluated.
        subsel = if a<sel[1] && b > sel[end]
            # Element is completely within the restricted basis.
            Colon()
        else
            # Element straddles restriction; find subset of functions
            # that are within the restriction.
            s = min(max(a,sel[1]),sel[end])
            e = max(min(b,sel[end]),sel[1])
            findfirst(isequal(s),sel):findfirst(isequal(e),sel)
        end

        @. v[sel[subsel] .- (a-1)] += @view((B′.wⁱ[i]*f[@elem(B′,x,i)])[subsel])
    end
    v .*= @view(B′.n[a:b])
    v
end
