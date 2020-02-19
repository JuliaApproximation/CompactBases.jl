# * Diagonal operators

# """
#     B'f*B

# Generate the matrix representing the quasi-diagonal operator `f` in
# the B-spline space `B`.

# # Examples

# ```jldoctest
# julia> B = BSpline(LinearKnotSet(3, 0, 1, 3))
# BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 3 (parabolic) on 0.0..1.0 (3 intervals)

# julia> x = QuasiDiagonal(axes(B,1))
# QuasiDiagonal{Float64,QuasiArrays.Inclusion{Float64,IntervalSets.Interval{:closed,:closed,Float64}}}(QuasiArrays.Inclusion{Float64,IntervalSets.Interval{:closed,:closed,Float64}}(0.0..1.0))

# julia> B'x*B
# 5×5 BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:
#  0.0037037    0.00462963  0.000925926   ⋅           ⋅
#  0.00462963   0.0259259   0.0236111    0.00138889   ⋅
#  0.000925926  0.0236111   0.0916667    0.0458333   0.00462963
#   ⋅           0.00138889  0.0458333    0.0851852   0.0342593
#   ⋅            ⋅          0.00462963   0.0342593   0.062963
# ```
# """
@simplify function *(Lc::QuasiAdjoint{<:Any,<:BSpline},
                     D::QuasiDiagonal,
                     R::BSpline)
    L = parent(Lc)
    χ = L.B
    ξ = Diagonal(getindex.(Ref(D.diag), R.x))*R.B
    k = max(order(L),order(R))
    m,n = size(L,2),size(R,2)

    T = promote_type(eltype(Lc),eltype(ξ),eltype(R))
    S = BandedMatrix(Zeros{T}(m,n), (k-1,k-1))
    overlap_matrix!(S, χ, ξ, weights(R))
end

# A & B restricted
function materialize(M::Mul{<:Any,<:Tuple{<:Adjoint{<:Any,<:RestrictionMatrix},
                                          <:QuasiAdjoint{<:Any,<:BSpline{T}},
                                          <:QuasiDiagonal,
                                          <:BSpline{T},
                                          <:RestrictionMatrix}}) where T
    restLc,Lc,D,R,restR = M.args
    axes(Lc,2) == axes(R,1) || throw(DimensionMismatch("axes must be same"))
    La,Lb = restriction_extents(parent(restLc))
    Ra,Rb = restriction_extents(restR)
    M = Lc*D*R
    M[1+La:end-Lb,1+Ra:end-Rb]
end

export QuasiDiagonal
