# * Diagonal matrices
DiagonalBlockDiagonal(A::AbstractMatrix, rows, cols) =
    BandedBlockBandedMatrix(A, rows,cols, (0,0), (0,0))

DiagonalBlockDiagonal(A::AbstractMatrix, rows) =
    DiagonalBlockDiagonal(A, rows, rows)

function (B::FEDVROrRestricted)(D::Diagonal)
    n = size(B,2)
    @assert size(D) == (n,n)
    all(order(B) .== 2) ? D : DiagonalBlockDiagonal(D, block_structure(B))
end

# * Dense operators
# ** Matrix construction

function Matrix(::UndefInitializer, B::Union{FEDVR{T},RestrictedFEDVR{T}}, ::Type{U}=T) where {T,U}
    if all(order(B) .== 2)
        n = size(B,2)
        dl = Vector{U}(undef, n-1)
        d = Vector{U}(undef, n)
        du = Vector{U}(undef, n-1)
        Tridiagonal(dl, d, du)
    else
        # @show B
        # @show
        rows = block_structure(B)
        # @show
        l,u = block_bandwidths(B,rows)

        BlockSkylineMatrix{U}(undef, rows, rows, (l,u))
    end
end

function Base.zeros(B::Union{FEDVR{T},RestrictedFEDVR{T}}, ::Type{U}=T) where {T,U}
    M = Matrix(undef, B, U)
    M .= zero(U)
    M
end

# ** Individual elements access

function set_element_common!(fun::Function, A, B::FEDVROrRestricted, i)
        els = elements(B)
    i ∉ els && throw(ArgumentError("Invalid element $(i) ∉ $(els)"))
    el1,eln = first(els),last(els)

    is = indices(B,2)

    elb = element_boundaries(parent(B))
    # If the first element only contributes through the bridge
    # function, also consider the second element as the "first"
    # element, since it will be exactly in the north-west corner of
    # the matrix.
    el1′ = el1 + (first(is) == elb[el1+1] ? 1 : 0)
    isfirst = i == el1 || i == el1′
    # If the last element only contributes through the bridge
    # function, also consider the next-to-last element as the "last"
    # element, since it will be exactly in the south-east corner of
    # the matrix.
    eln′ = eln - (last(is) == elb[eln] ? 1 : 0)
    islast = i == eln || i == eln′

    b = fun(i)
    o = order(B,i)
    @assert size(b,1) == size(b,2) == o

    # Find start and end of the interior block; the first and last
    # elements (within the restriction) have their interior blocks as
    # their outer blocks, which are cropped as the restriction
    # prescribes.
    s = 1 + (i == el1 ? first(is)-elb[el1] : (i == el1′ ? first(is)-elb[el1′] : 1))
    e = o - (i == eln ? elb[eln+1]-last(is) : (i == eln′ ? elb[eln]-last(is) : 1))

    # First element only contributes through the bridge function
    if isfirst && s == o
        A[1,1] += b[end,end]
        return
    end

    # Last element only contributes through the bridge function
    if islast && e == 1
        A[end,end] += b[1,1]
        return
    end

    b,o,s:e,isfirst,islast,el1′
end

function set_element!(fun::Function, A::BlockSkylineMatrix, B::FEDVROrRestricted, i)
    res = set_element_common!(fun, A, B, i)
    isnothing(res) && return
    b,o,sel,isfirst,islast,el1′ = res

    # Find block index of the north-west block of element i; this is a
    # linear search which is not optimal for large number of elements.
    j = isfirst ? 1 : 0
    if !isfirst
        for el in el1′:i-1
            j += 1 + (order(B,el) > 2 || el == el1′ ? 1 : 0)
        end
    end

    if o > 2
        ji = j+!isfirst
        # Interior block
        A[Block(ji,ji)] += b[sel,sel]

        if !isfirst
            # North-west corner
            A[Block(j,j)] .+= b[1,1]
            # West edge
            A[Block(j+1,j)] = b[sel,1:1]
            # North edge
            A[Block(j,ji)] = b[1:1,sel]
            if !islast
                # Nort-east corner
                A[Block(j,ji+1)] = b[1:1,end:end]
                # South-west corner
                A[Block(j+2,j)] = b[end:end,1:1]
            end
        end
        if !islast
            # South-east corner
            A[Block(ji+1,ji+1)] .+= b[end,end]
            # East edge
            A[Block(ji,ji+1)] = b[sel,end:end]
            # South edge
            A[Block(ji+1,ji)] = b[end:end, sel]
        end
    else
        A[Block(j,j)] .+= b[1,1]
        A[Block(j+1,j+1)] .+= b[2,2]
        A[Block(j,j+1)] .= b[1,2]
        A[Block(j+1,j)] .= b[1,2]
    end
end

function set_element!(fun::Function, A::Tridiagonal, B::FEDVROrRestricted, i)
    res = set_element_common!(fun, A, B, i)
    isnothing(res) && return
    b,o,sel,isfirst,islast,el1′ = res

    ii = i - el1′ + 1

    A.du[ii] = b[1,2]
    A.d[ii:ii+1] .+= diag(b)
    A.dl[ii] = b[2,1]
end

function set_elements!(fun::Function, A, B::FEDVROrRestricted)
    A .= false
    for i in elements(B)
        set_element!(fun, A, B, i)
    end
    A
end

# * Scalar operators

@simplify function *(Ac::QuasiAdjoint{<:Any,<:FEDVR},
                     D::QuasiDiagonal,
                     B::FEDVR)
    A = parent(Ac)
    A == B || throw(ArgumentError("Cannot multiply functions on different grids"))

    Diagonal(getindex.(Ref(D.diag), B.x))
end

# A & B restricted
@simplify function *(Ac::AdjointRestrictedFEDVR,
                     D::QuasiDiagonal,
                     B::RestrictedFEDVR)
    # This is mainly for type-stability; it would be trivial to
    # generate the proper banded matrix with one off-diagonal, from
    # the combination of two differently restricted bases, but we
    # would like to have Diagonal as the result if they are equal, and
    # this has higher priority. On the other hand, you typically only
    # compute scalar operators in the beginning of the calculation,
    # and thus type-instability is not a big problem, so this
    # behaviour may change in the future.
    reverse(axes(Ac)) == axes(B) || throw(DimensionMismatch("axes must be same"))
    A = parent(Ac)
    parent(A) == parent(B) || throw(ArgumentError("Cannot multiply incompatible FEDVR expansions"))

    a,b = restriction_extents(B)
    Diagonal(getindex.(Ref(D.diag), parent(B).x[1+a:end-b]))
end
