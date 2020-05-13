# Inner products and norms

The inner product between two splines is given by

$$\begin{equation}
\label{eqn:inner-product}
\braket{s}{\tilde{s}} \defd
\int\diff{x} \conj{s}(x) \tilde{s}(x) \approx
\sum_l w_l \conj{s}(x_l) \tilde{s}(x_l).
\end{equation}$$

From the [definition of a spline](splines.md), we have

$$\begin{equation}
\tag{\ref{eqn:inner-product}*}
\braket{s}{\tilde{s}} = (B\vec{c})^HB\vec{\tilde{c}} = \vec{c}^H \mat{S} \vec{\tilde{c}},
\end{equation}$$

where the overlap matrix is given by

$$\begin{equation}
\label{eqn:overlap-matrix}
\mat{S}_{ij} \defd
\braket{\B{i}{k}}{\B{j}{k}} \approx
\sum_l w_l \conj{\B{i}{k}(x_l)}{\B{j}{k}(x_l)}.
\end{equation}$$

When the B-spline basis is created, the B-splines are evaluated on the
quadrature roots and stored in a sparse matrix, and the banded overlap
matrix is formed:

```julia
julia> B = BSpline(LinearKnotSet(7, 0, 1, 10))
BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 7 on 0.0..1.0 (10 intervals)

julia> B.B
80×16 SparseArrays.SparseMatrixCSC{Float64,Int64} with 560 stored entries:
  [1 ,  1]  =  0.886629
  [2 ,  1]  =  0.525563
  [3 ,  1]  =  0.196947
  [4 ,  1]  =  0.0429226
  [5 ,  1]  =  0.00463197
  [6 ,  1]  =  0.000178262
  [7 ,  1]  =  1.10427e-6
  [8 ,  1]  =  6.12673e-11
  [1 ,  2]  =  0.11053
  [2 ,  2]  =  0.411336
  [3 ,  2]  =  0.543708
  [4 ,  2]  =  0.422368
  [5 ,  2]  =  0.234511
  [6 ,  2]  =  0.111732
  [7 ,  2]  =  0.0558643
  [8 ,  2]  =  0.0351626
  ⋮
  [73, 15]  =  0.0351626
  [74, 15]  =  0.0558643
  [75, 15]  =  0.111732
  [76, 15]  =  0.234511
  [77, 15]  =  0.422368
  [78, 15]  =  0.543708
  [79, 15]  =  0.411336
  [80, 15]  =  0.11053
  [73, 16]  =  6.12673e-11
  [74, 16]  =  1.10427e-6
  [75, 16]  =  0.000178262
  [76, 16]  =  0.00463197
  [77, 16]  =  0.0429226
  [78, 16]  =  0.196947
  [79, 16]  =  0.525563
  [80, 16]  =  0.886629

julia> B.S
16×16 BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:
 0.00769231   0.00492814   0.00142795   0.000218646  1.79179e-5   7.31907e-7   …   ⋅            ⋅            ⋅            ⋅            ⋅
 0.00492814   0.0110567    0.00855534   0.00328365   0.000674239  7.04773e-5       ⋅            ⋅            ⋅            ⋅            ⋅
 0.00142795   0.00855534   0.0147851    0.0118803    0.00501486   0.00109218       ⋅            ⋅            ⋅            ⋅            ⋅
 0.000218646  0.00328365   0.0118803    0.0187451    0.015232     0.00647614       ⋅            ⋅            ⋅            ⋅            ⋅
 1.79179e-5   0.000674239  0.00501486   0.015232     0.0232999    0.0188319        ⋅            ⋅            ⋅            ⋅            ⋅
 7.31907e-7   7.04773e-5   0.00109218   0.00647614   0.0188319    0.0289487    …  1.92709e-11   ⋅            ⋅            ⋅            ⋅
 1.15625e-8   2.93435e-6   0.000100846  0.00126285   0.0074747    0.0228901       1.57581e-7   3.01107e-11   ⋅            ⋅            ⋅
  ⋅           3.61328e-10  5.83108e-7   4.39448e-5   0.000854471  0.00665362      2.83876e-5   2.46171e-7   7.13735e-11   ⋅            ⋅
  ⋅            ⋅           7.13735e-11  2.46171e-7   2.83876e-5   0.000726484     0.000854471  4.39448e-5   5.83108e-7   3.61328e-10   ⋅
  ⋅            ⋅            ⋅           3.01107e-11  1.57581e-7   2.37367e-5      0.0074747    0.00126285   0.000100846  2.93435e-6   1.15625e-8
  ⋅            ⋅            ⋅            ⋅           1.92709e-11  1.31324e-7   …  0.0188319    0.00647614   0.00109218   7.04773e-5   7.31907e-7
  ⋅            ⋅            ⋅            ⋅            ⋅           1.92709e-11     0.0232999    0.015232     0.00501486   0.000674239  1.79179e-5
  ⋅            ⋅            ⋅            ⋅            ⋅            ⋅              0.015232     0.0187451    0.0118803    0.00328365   0.000218646
  ⋅            ⋅            ⋅            ⋅            ⋅            ⋅              0.00501486   0.0118803    0.0147851    0.00855534   0.00142795
  ⋅            ⋅            ⋅            ⋅            ⋅            ⋅              0.000674239  0.00328365   0.00855534   0.0110567    0.00492814
  ⋅            ⋅            ⋅            ⋅            ⋅            ⋅           …  1.79179e-5   0.000218646  0.00142795   0.00492814   0.00769231
```

If we thus have to vectors of expansion coefficients, we can calculate
the inner product of the correspond splines very simply:

```julia
julia> c = B \ x -> sin(2π*x)
16-element Array{Float64,1}:
  1.5675074227379815e-7
  0.10471902453028906
  0.3141611768331591
  0.6159137789819404
  0.9438219608666177
  1.1224640274531983
  0.9080930336327867
  0.34686036375515594
 -0.3468603637551552
 -0.9080930336327886
 -1.1224640274531954
 -0.9438219608666232
 -0.6159137789819348
 -0.3141611768331642
 -0.10471902453028659
 -1.567507426314564e-7

julia> c̃ = B \ x -> cos(2π*x)
16-element Array{Float64,1}:
  1.0000000051659867
  1.0000000496779038
  0.9736805198549477
  0.8552480183377358
  0.5498061181490912
  3.468447866195647e-7
 -0.6597681895565798
 -1.0675270371816847
 -1.0675270371816876
 -0.6597681895565795
  3.4684478615396205e-7
  0.5498061181490913
  0.855248018337735
  0.9736805198549482
  1.0000000496779042
  1.0000000051659863

julia> c'*B.S*c̃
1-element Array{Float64,1}:
 9.71445146547012e-17
```

[orthogonal, as expected, since $\int_0^1\diff{x}\sin(2\pi x)\cos(2\pi
x)=0$]. However, instead of manually supplying the overlap matrix
$\mat{S}$, it is more convenient to first construct the splines and
then compute the inner product:

```julia
julia> s = B*c
Spline on BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 7 on 0.0..1.0 (10 intervals)

julia> s̃ = B*c̃
Spline on BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 7 on 0.0..1.0 (10 intervals)

julia> s's̃
ContinuumArrays.QuasiArrays.ApplyQuasiArray{Float64,1,LazyArrays.Applied{ContinuumArrays.QuasiArrays.LazyQuasiArrayApplyStyle,typeof(*),Tuple{Adjoint{Float64,Array{Float64,1}},ContinuumArrays.QuasiArrays.QuasiAdjoint{Float64,BSpline{Float64,Float64,LinearKnotSet{7,7,7,Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},Array{Float64,1},Array{Float64,1},SparseArrays.SparseMatrixCSC{Float64,Int64},BandedMatrices.BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}}},BSpline{Float64,Float64,LinearKnotSet{7,7,7,Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},Array{Float64,1},Array{Float64,1},SparseArrays.SparseMatrixCSC{Float64,Int64},BandedMatrices.BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}},Array{Float64,1}}}}([1.56751e-7 0.104719 … -0.104719 -1.56751e-7]⋆ContinuumArrays.QuasiArrays.QuasiAdjoint{Float64,BSpline{Float64,Float64,LinearKnotSet{7,7,7,Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},Array{Float64,1},Array{Float64,1},SparseArrays.SparseMatrixCSC{Float64,Int64},BandedMatrices.BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}}}(BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 7 on 0.0..1.0 (10 intervals))⋆BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 7 on 0.0..1.0 (10 intervals)⋆[1.0, 1.0, 0.973681, 0.855248, 0.549806, 3.46845e-7, -0.659768, -1.06753, -1.06753, -0.659768, 3.46845e-7, 0.549806, 0.855248, 0.973681, 1.0, 1.0])
```