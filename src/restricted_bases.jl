# * Auxilliary type definitions for restricted bases

const RestrictionMatrix = BandedMatrix{<:Int, <:FillArrays.Ones}
const RestrictedQuasiArray{T,N,B<:Basis} = SubQuasiArray{T,N,B}
const AdjointRestrictedQuasiArray{T,N,B<:Basis} = QuasiAdjoint{T,<:RestrictedQuasiArray{T,N,B}}

const BasisOrRestricted{B<:Basis} = Union{B,<:RestrictedQuasiArray{<:Any,<:Any,<:B}}
const AdjointBasisOrRestricted{B<:Basis} = Union{<:QuasiAdjoint{<:Any,B},<:AdjointRestrictedQuasiArray{<:Any,<:Any,<:B}}

indices(B::Basis) = axes(B)
indices(B::RestrictedQuasiArray) = B.indices
indices(B, i) = indices(B)[i]

unrestricted_basis(R::AbstractQuasiMatrix) = R
unrestricted_basis(R::RestrictedQuasiArray) = parent(R)

==(A::BasisOrRestricted, B::BasisOrRestricted) =
    unrestricted_basis(A) == unrestricted_basis(B)

restriction_extents(::Basis) = 0,0
function restriction_extents(B̃::RestrictedQuasiArray)
    B = parent(B̃)
    a,b = B̃.indices[2][[1,end]]
    a-1,size(B,2)-b
end

restriction(B) = Diagonal(Ones{Int}(size(B,2)))
restriction(B̃::RestrictedQuasiArray) = last(LazyArrays.arguments(B̃))

function combined_restriction(A,B)
    parent(A) == parent(B) ||
        throw(ArgumentError("Cannot multiply functions on different grids"))
    restriction(A)'*restriction(B)
end

function combined_restriction_selection(A,B)
    parent(A) == parent(B) ||
        throw(ArgumentError("Cannot multiply functions on different grids"))

    la,lb = restriction_extents(A)
    ra,rb = restriction_extents(B)

    lsel = 1+ra:(size(A,2)-max(0,rb-lb))
    rsel = 1+la:(size(B,2)-max(0,lb-rb))
    lsel, rsel
end

function show(io::IO, B̃::RestrictedQuasiArray{<:Any,2})
    B = parent(B̃)
    a,b = B̃.indices[2][[1,end]]
    N = size(B,2)
    show(io, B)
    write(io, ", restricted to basis functions $(a)..$(b) $(a>1 || b<N ? "⊂" : "⊆") 1..$(N)")
end

IntervalSets.leftendpoint(B::RestrictedQuasiArray) =
    leftendpoint(parent(B))
IntervalSets.rightendpoint(B::RestrictedQuasiArray) =
    rightendpoint(parent(B))

locs(B::RestrictedQuasiArray) = locs(parent(B))[indices(B,2)]
real_locs(B::RestrictedQuasiArray) = real_locs(parent(B))[indices(B,2)]

distribution(B::RestrictedQuasiArray) = distribution(parent(B))
