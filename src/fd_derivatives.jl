# * Derivatives
# ** Three-point stencils

# α = super-/subdiagonal of second derivative
# α̃ = subdiagonal of second derivative, in case not same as superdiagonal
# β = -1/2 diagonal of second derivative
# γ = diagonal of first derivative
# δ = super-/subdiagonal of first derivative

# *** Cartesian finite-differences

α(::FiniteDifferences{T}, ::Integer) where T = one(T)
β(::FiniteDifferences{T}, ::Integer) where T = one(T)
γ(::FiniteDifferences{T}, ::Integer) where T = zero(T)
δ(::FiniteDifferences{T}, ::Integer) where T = one(T)

# *** Staggered finite-differences

# Eq. (20), Schafer (2000)
α(::Uniform,  ::StaggeredFiniteDifferences{T}, j::Integer) where T = j^2/(j^2 - one(T)/4)
β(::Uniform, B::StaggeredFiniteDifferences{T}, j::Integer) where T = (j^2 - j + one(T)/2)/(j^2 - j + one(T)/4)

α(d::Uniform, B::StaggeredFiniteDifferences) = α.(Ref(d), Ref(B), 1:length(B.r)-1)
β(d::Uniform, B::StaggeredFiniteDifferences) = β.(Ref(d), Ref(B), 1:length(B.r))
δ(d::Uniform, B::StaggeredFiniteDifferences) = α(d, B)

function staggered_stencil_common(fun, B::StaggeredFiniteDifferences{T}, δN=0) where T
    r = B.r
    N = length(r)
    v = zeros(T, N+δN)

    r̃ = vcat(r, 2r[end]-r[end-1])
    a = 2r[1]-r[2]
    b,c,d = r[1],r[2],r[3]

    for j = 1:N+δN
        if j < N
            d = r̃[j+2]
        end
        v[j] = fun(a,b,c,d)

        a,b,c = b,c,d
    end
    v
end

function α(::NonUniform, B::StaggeredFiniteDifferences)
    staggered_stencil_common(B, -1) do a,b,c,d
        # Eq. (A13) Krause 1999
        2/((c-b)*√((d-b)*(c-a)))*((b+c)/2)^2/(b*c)
    end
end

function α(::NonUniform, B::StaggeredFiniteDifferences, f::Function)
    staggered_stencil_common(B, -1) do a,b,c,d
        rⱼ₊₁₂ = (b+c)/2
        # Eq. (A13) Krause 1999
        2f(rⱼ₊₁₂)/((c-b)*√((d-b)*(c-a)))*rⱼ₊₁₂^2/(b*c)
    end
end

function β(::NonUniform, B::StaggeredFiniteDifferences)
    staggered_stencil_common(B) do a,b,c,d
        # Eq. (A14) Krause 1999
        fp = ((c+b)/2b)^2
        fm = ((b+a)/2b)^2
        1/(c-a)*(1/(c-b)*fp + 1/(b-a)*fm)
    end
end

function β(::NonUniform, B::StaggeredFiniteDifferences, f::Function)
    staggered_stencil_common(B) do a,b,c,d
        # Eq. (A14) Krause 1999
        b⁻² = inv(b^2)
        rⱼ₊₁₂ = (b+c)/2
        rⱼ₋₁₂ = (a+b)/2
        fp = rⱼ₊₁₂^2 * b⁻²
        fm = rⱼ₋₁₂^2 * b⁻²
        1/(c-a)*(f(rⱼ₊₁₂)/(c-b)*fp + f(rⱼ₋₁₂)/(b-a)*fm)
    end
end

function δ(::NonUniform, B::StaggeredFiniteDifferences)
    staggered_stencil_common(B, -1) do a,b,c,d
        2/(√((d-b)*(c-a)))*((b+c)/2)^2/(b*c)
    end
end

α(B::StaggeredFiniteDifferences, args...)   = α(distribution(B), B, args...)
β(B::StaggeredFiniteDifferences, args...)   = β(distribution(B), B, args...)
γ(B::StaggeredFiniteDifferences{T}) where T = zeros(T, length(B.r))
δ(B::StaggeredFiniteDifferences)            = δ(distribution(B), B)

# *** Implicit finite-differences

α( ::Uniform, B::ImplicitFiniteDifferences{T}) where T = ones(T, length(B.r)-1)
α̃(d::Uniform, B::ImplicitFiniteDifferences) = α(d, B)
β( ::Uniform, B::ImplicitFiniteDifferences{T}) where T = ones(T, length(B.r))
γ( ::Uniform, B::ImplicitFiniteDifferences{T}) where T = zeros(T, length(B.r))
δ( ::Uniform, B::ImplicitFiniteDifferences{T}) where T = ones(T, length(B.r)-1)
δ̃(d::Uniform, B::ImplicitFiniteDifferences) = δ(d, B)

# Non-uniform, implicit finite differences, as described by
# Patchkovskii (2016)

function implicit_stencil_common(fun, B::ImplicitFiniteDifferences{T}, s=1, δN=0) where T
    r = B.r
    N = length(r)
    v = zeros(T, N+δN)

    r̃ = vcat(r[s:end], 2r[end]-r[end-1])
    a = r̃[1]-(s > 1 ? r[s-1] : 0)
    b = r̃[2]-r̃[1]
    c = zero(T)
    for j = 1:N+δN
        if j < N - (s-1)
            c = r̃[j+2]-r̃[j+1]
        end

        v[j] = fun(a,b,c)

        a,b = b,c
    end
    v
end

function α̃(::NonUniform, B::ImplicitFiniteDifferences)
    implicit_stencil_common(B, 2, -1) do a,b,c
        2/(a*(a+b)) # Eq. (33)
    end
end

function α(::NonUniform, B::ImplicitFiniteDifferences)
    implicit_stencil_common(B, 1, -1) do a,b,c
        2/(b*(a+b)) # Eq. (33)
    end
end

function β(::NonUniform, B::ImplicitFiniteDifferences)
    implicit_stencil_common(B) do a,b,c
        1/(a*b) # Eq. (33), remember β is multiplied by -2
    end
end

α(B::ImplicitFiniteDifferences) = α(distribution(B), B)
α̃(B::ImplicitFiniteDifferences) = α̃(distribution(B), B)
β(B::ImplicitFiniteDifferences) = β(distribution(B), B)
γ(B::ImplicitFiniteDifferences) = γ(distribution(B), B)
δ(B::ImplicitFiniteDifferences) = δ(distribution(B), B)
δ̃(B::ImplicitFiniteDifferences) = δ̃(distribution(B), B)

# *** Common

α(B::AbstractFiniteDifferences) = α.(Ref(B), B.j[1:end-1])
α̃(B::AbstractFiniteDifferences) = α(B)
β(B::AbstractFiniteDifferences) = β.(Ref(B), B.j)
γ(B::AbstractFiniteDifferences) = γ.(Ref(B), B.j)
δ(B::AbstractFiniteDifferences) = δ.(Ref(B), B.j[1:end-1])
δ̃(B::AbstractFiniteDifferences) = δ(B)

# *** Restrictions

for f in [:α, :α̃, :δ, :δ̃]
    @eval begin
        function $f(B̃::RestrictedFiniteDifferences)
            j = indices(B̃,2)
            $f(parent(B̃))[j[1]:(j[end]-1)]
        end
    end
end

for f in [:β, :γ]
    @eval begin
        function $f(B̃::RestrictedFiniteDifferences)
            $f(parent(B̃))[indices(B̃,2)]
        end
    end
end

for f in [:α, :α̃, :β, :γ, :δ, :δ̃]
    @eval begin
        function $f(Ã::BasisOrRestricted{<:AbstractFiniteDifferences},
                    B̃::BasisOrRestricted{<:AbstractFiniteDifferences},
                    band=0)
            i,i′ = extrema(indices(Ã,2))
            j,j′ = extrema(indices(B̃,2))
            k = max(i+min(0,band),j-max(0,band)):min(i′+min(0,band),j′-max(0,band))
            $f(parent(B̃))[k]
        end
    end
end

# ** Materialization

non_symmetric_laplacian(::Uniform, B::ImplicitFiniteDifferences) = false
non_symmetric_laplacian(::NonUniform, B::ImplicitFiniteDifferences) = true
non_symmetric_laplacian(B::ImplicitFiniteDifferences) = non_symmetric_laplacian(distribution(B), B)
non_symmetric_laplacian(::AbstractFiniteDifferences) = false

function derivative_matrix(::Type{T}, Ac, B, diff_order) where T
    A = parent(Ac)
    parent(A) == parent(B) ||
        throw(ArgumentError("Cannot multiply functions on different grids"))
    Ai,Bi = last(indices(A)), last(indices(B))
    if Ai == Bi
        n = length(Ai)
        if diff_order == 1 || non_symmetric_laplacian(parent(B))
            Tridiagonal(Vector{T}(undef, n-1),
                        Vector{T}(undef, n),
                        Vector{T}(undef, n-1))
        elseif diff_order == 2
            SymTridiagonal(Vector{T}(undef, n),
                           Vector{T}(undef, n-1))
        end
    else
        m,n = length(Ai),length(Bi)
        offset = Ai[1]-Bi[1]
        BandedMatrix{T}(undef, (m,n), (1-offset,1+offset))
    end
end

@materialize function *(Ac::AdjointBasisOrRestricted{<:AbstractFiniteDifferences},
                        D::Derivative,
                        B::BasisOrRestricted{<:AbstractFiniteDifferences})
    T -> begin
        derivative_matrix(T, Ac, B, 1)
    end
    dest::Tridiagonal{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))

        # Central difference approximation
        μ = 2step(B)
        dest.dl .= -δ̃(B)/μ
        dest.d .= γ(B)/μ
        dest.du .= δ(B)/μ
    end
    dest::BandedMatrix{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))

        dl,d,du = bandrange(dest)
        # Central difference approximation
        μ = 2step(B)
        dest[Band(dl)] .= -δ̃(A,B,-1)/μ
        dest[Band(d)] .= γ(A,B)/μ
        dest[Band(du)] .= δ(A,B,1)/μ
    end
end

@materialize function *(Ac::AdjointBasisOrRestricted{<:AbstractFiniteDifferences},
                        Dc::QuasiAdjoint{<:Any,<:Derivative},
                        D::Derivative,
                        B::BasisOrRestricted{<:AbstractFiniteDifferences})
    T -> begin
        derivative_matrix(T, Ac, B, 2)
    end
    dest::SymTridiagonal{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))

        δ² = step(B)^2
        dest.dv .= -2β(B)/δ²
        dest.ev .= α(B)/δ²
    end
    dest::Tridiagonal{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))

        δ² = step(B)^2
        dest.dl .= α̃(B)/δ²
        dest.d .= -2β(B)/δ²
        dest.du .= α(B)/δ²
    end
    dest::BandedMatrix{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))

        dl,d,du = bandrange(dest)
        δ² = step(B)^2
        dest[Band(dl)] .= α̃(A,B,-1)/δ²
        dest[Band(d)] .= -2β(A,B)/δ²
        dest[Band(du)] .= α(A,B,1)/δ²
    end
end

@materialize function *(Ac::AdjointBasisOrRestricted{<:StaggeredFiniteDifferences},
                        Dc::QuasiAdjoint{<:Any,<:Derivative},
                        O::QuasiDiagonal,
                        D::Derivative,
                        B::BasisOrRestricted{<:StaggeredFiniteDifferences})
    T -> begin
        derivative_matrix(T, Ac, B, 2)
    end
    dest::SymTridiagonal{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))

        δ² = step(B)^2
        f = r -> O.diag[max(zero(T),r)]
        dest.dv .= -2β(B,f)/δ²
        dest.ev .= α(B,f)/δ²
    end
    dest::Tridiagonal{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))

        δ² = step(B)^2
        dest.dl .= α̃(B)/δ²
        dest.d .= -2β(B)/δ²
        dest.du .= α(B)/δ²
    end
    dest::BandedMatrix{T} -> begin
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Cannot multiply functions on different grids"))

        dl,d,du = bandrange(dest)
        δ² = step(B)^2
        dest[Band(dl)] .= α̃(A,B,-1)/δ²
        dest[Band(d)] .= -2β(A,B)/δ²
        dest[Band(du)] .= α(A,B,1)/δ²
    end
end

# *** Implicit derivatives

#=
This is an implementation of finite difference scheme described in

- Muller, H. G. (1999). An Efficient Propagation Scheme for the
  Time-Dependent Schrödinger equation in the Velocity Gauge. Laser
  Physics, 9(1), 138–148.

where the first derivative is approximated as

\[\partial f =
\left(1+\frac{h^2}{6}\Delta_2\right)^{-1}
\Delta_1 f
\equiv
M_1^{-1}\tilde{\Delta}_1 f,\]
where
\[M_1 \equiv
\frac{1}{6}
\begin{bmatrix}
4+\lambda' & 1 &\\
1 & 4 & 1\\
& 1 & 4 & \ddots\\
&&\ddots&\ddots\\
\end{bmatrix},\]
and
\[\tilde{\Delta}_1 \equiv
\frac{1}{2h}
\begin{bmatrix}
\lambda & 1 &\\
-1 &  & 1\\
& -1 &  & \ddots\\
&&\ddots&\ddots\\
\end{bmatrix},\]

where \(\lambda=\lambda'=\sqrt{3}-2\) for problems with a singularity
at the boundary \(r=0\) and zero otherwise; and the second derivative
as

\[\partial^2 f =
\left(1+\frac{h^2}{12}\Delta_2\right)^{-1}
\Delta_2 f
\equiv
-2M_2^{-1}\Delta_2 f,\]
where
\[M_2 \equiv
-\frac{1}{6}
\begin{bmatrix}
10-2\delta\beta_1 & 1 &\\
1 & 10 & 1\\
& 1 & 10 & \ddots\\
&&\ddots&\ddots\\
\end{bmatrix},\]
and
\[\Delta_2 \equiv
\frac{1}{h^2}
\begin{bmatrix}
-2(1+\delta\beta_1) & 1 &\\
1 & -2 & 1\\
& 1 & -2 & \ddots\\
&&\ddots&\ddots\\
\end{bmatrix},\]

where, again, \(\delta\beta_1 = -Zh[12-10Zh]^{-1}\) is a correction
introduced for problems singular at the origin.

=#

function implicit_lhs(::Uniform, B::BasisOrRestricted{<:ImplicitFiniteDifferences{T}}, difforder) where T
    m = size(B,2)
    M = SymTridiagonal(Vector{T}(undef, m), Vector{T}(undef, m-1))

    j₁ = first(indices(B, 2))

    if difforder == 1
        M.dv .= 4
        M.ev .= 1

        # # Eq. (20ff), Muller (1999).
        # j₁ == 1 && !iszero(Z) && (M.dv[1] += √3 - 2)

        # f'  =   M₁⁻¹*Δ₁*f   Eq. (19)
        M /= 6
    elseif difforder == 2
        M.dv .= 10
        M.ev .= 1

        # f'' = -2M₂⁻¹*Δ₂*f   Eq. (13)
        M /= (-6)*(-2)
    end
end

function implicit_lhs(::NonUniform, B::BasisOrRestricted{<:ImplicitFiniteDifferences{T}}, difforder) where T
    dl, d, du = if difforder == 1
    elseif difforder == 2
        # Eq. (33), Patchkovskii (2016)
        dl = implicit_stencil_common(B, 2, -1) do a,b,c
            (a^2 + a*b - b^2)/(6*a*(a+b))
        end
        d = implicit_stencil_common(B) do a,b,c
            (a^2 + 3a*b + b^2)/(6*a*b)
        end
        du = implicit_stencil_common(B, 1, -1) do a,b,c
            (-a^2 + a*b + b^2)/(6*b*(a+b))
        end
        dl, d, du
    end

    Tridiagonal(dl, d, du)
end

implicit_lhs(B::BasisOrRestricted{<:ImplicitFiniteDifferences}, args...) =
    implicit_lhs(distribution(B), B, args...)

implicit_derivative(::Type{T}, M, difforder) where T =
    ImplicitDerivative(copyto!(similar(M, T), M),
                       implicit_lhs(last(M.args), difforder))

LazyArrays.simplifiable(::typeof(*), ::AdjointBasisOrRestricted{B}, ::Derivative, ::BasisOrRestricted{B}) where {T,B<:ImplicitFiniteDifferences{T}} = Val(true)
LazyArrays._simplify(::typeof(*), A::AdjointBasisOrRestricted{B}, D::Derivative, C::BasisOrRestricted{B})  where {T,B<:ImplicitFiniteDifferences{T}} =
        implicit_derivative(T, ApplyQuasiArray(*,A,D,C), 1)

LazyArrays.simplifiable(::typeof(*), ::AdjointBasisOrRestricted{B}, ::QuasiAdjoint{T,<:Derivative}, ::Derivative, ::BasisOrRestricted{B}) where {T,B<:ImplicitFiniteDifferences{T}} = Val(true)
LazyArrays._simplify(::typeof(*), Ac::AdjointBasisOrRestricted{B}, Dc::QuasiAdjoint{T,<:Derivative}, D::Derivative, A::BasisOrRestricted{B}) where {T,B<:ImplicitFiniteDifferences{T}} = 
        implicit_derivative(T, ApplyQuasiArray(*,Ac,Dc,D,A), 2)

# *** Coulomb derivatives

# Below, we implement some corrections to the derivative matrices in
# case of an attractive Coulomb problem [-∂²/2 - Z/r + ℓ(ℓ+1)/2r^2].
#
# For StaggeredFiniteDifferences, the correction is just a heuristic
# quoted by
#
# - Schafer, K. J. (2009). Numerical Methods in Strong Field Physics. In
#   T. Brabec (Eds.), (pp. 111–145). Springer.
#
# whereas for ImplicitFiniteDifferences, it is derived via Taylor
# expansions, for the uniform case, see
#
# - Muller, H. G. (1999). An Efficient Propagation Scheme for the
#   Time-Dependent Schrödinger equation in the Velocity Gauge. Laser
#     Physics, 9(1), 138–148.
#
# and for the non-uniform case, see
#
# - Patchkovskii, S., & Muller, H. (2016). Simple, Accurate, and
#   Efficient Implementation of 1-Electron Atomic Time-Dependent
#   SchröDinger Equation in Spherical Coordinates. Computer Physics
#   Communications, 199,
#   153–169. http://dx.doi.org/10.1016/j.cpc.2015.10.014

gradient_correction!(∂, B::AbstractFiniteDifferences, Z, ℓ) = ∂

function materialize(M::Mul{<:Any, <:Tuple{
    <:AdjointBasisOrRestricted{<:AbstractFiniteDifferences},
    <:CoulombDerivative,
    <:BasisOrRestricted{<:AbstractFiniteDifferences}}}) where T
    Ac,CD,B = M.args
    D = CD.D
    ∂ = Ac*D*B
    if first(indices(parent(Ac),2)) == first(indices(B,2)) == 1
        gradient_correction!(∂, parent(B), CD.Z, CD.ℓ)
    end
    ∂
end

function laplacian_correction!(∂²::AbstractMatrix, B::StaggeredFiniteDifferences, Z, ℓ)
    if ℓ == 0 && Z ≠ 0
        ρ = local_step(B,1)
        # Eq. (22) Schafer (2009)
        δβ₁ = 1/ρ^2*Z*ρ/8*(1 + Z*ρ)
        ∂²[1,1] -= 2δβ₁
    end

    ∂²
end

function laplacian_correction!(∂²::ImplicitDerivative, ::Uniform, B::ImplicitFiniteDifferences, Z, ℓ)
    if ℓ == 0 && Z ≠ 0
        ρ = step(B)
        # Eq. (17), Muller (1999)
        δβ₁ = -Z*ρ/(12 - 10Z*ρ)
        ∂².Δ[1,1] -= 2δβ₁/ρ^2
        ∂².M[1,1] -= 2δβ₁/((-6)*(-2))

        ∂².M⁻¹ = factorize(∂².M)
    end
    ∂²
end

function laplacian_correction!(∂²::ImplicitDerivative, ::NonUniform, B::ImplicitFiniteDifferences, Z, ℓ)
    Z == 0 && return ∂²
    a = B.r[1]
    b = B.r[2] - a
    if ℓ == 0
        # Eq. (34), Patchkovskii (2016)
        num = b*(a*(3a+4b)*Z - 6*(a+b))
        ∂².Δ[1,1] = 2*(a+b)*(6a - (3a - b)*(a+b)*Z)/(a^2*num) # /local_step(B, 1)
        ∂².M[1,1] = (a+b)*(a*(a+b)*(a+3b)*Z - 3*(a^2 + 3b*a + b^2))/(3a*num)
        ∂².M[1,2] = (a^3*(-Z) + a^2*(b*Z+3) + b*a*(2*b*Z-3) - 3*b^2)/(3*num)
    else
        a = B.r[1]
        b = B.r[2] - a

        # Eq. (35)
        num = b*(3a+4b)
        ∂².Δ[1,1] = (3a - b)*(a+b)^2/(a^3*numb)
        ∂².M[1,1] = (a+b)^2*(a+3b)/(3a*num)
        ∂².M[1,2] = -(a-2b)*(a+b)/(3num)
    end
    ∂².M⁻¹ = factorize(∂².M)
    ∂²
end

laplacian_correction!(∂²::ImplicitDerivative, B::ImplicitFiniteDifferences, Z, ℓ) =
    laplacian_correction!(∂², distribution(B), B, Z, ℓ)

laplacian_correction!(∂², B::AbstractFiniteDifferences, Z, ℓ) = ∂²

function materialize(M::Mul{<:Any, <:Tuple{
    <:AdjointBasisOrRestricted{<:AbstractFiniteDifferences},
    <:QuasiAdjoint{T,<:CoulombDerivative},
    <:CoulombDerivative,
    <:BasisOrRestricted{<:AbstractFiniteDifferences}}}) where T
    Ac,CDc,CD,B = M.args
    @assert parent(CDc).Z == CD.Z
    Dc = adjoint(parent(CDc).D)
    D = CD.D
    ∂² = Ac*Dc*D*B
    if first(indices(parent(Ac),2)) == first(indices(B,2)) == 1
        laplacian_correction!(∂², parent(B), CD.Z, CD.ℓ)
    end
    ∂²
end

function materialize(M::Mul{<:Any, <:Tuple{
    <:AdjointBasisOrRestricted{<:AbstractFiniteDifferences},
    <:QuasiAdjoint{T,<:CoulombDerivative},
    <:QuasiDiagonal,
    <:CoulombDerivative,
    <:BasisOrRestricted{<:AbstractFiniteDifferences}}}) where T
    Ac,CDc,O,CD,B = M.args
    @assert parent(CDc).Z == CD.Z
    Dc = adjoint(parent(CDc).D)
    D = CD.D
    ∂² = Ac*Dc*O*D*B
    if first(indices(parent(Ac),2)) == first(indices(B,2)) == 1
        laplacian_correction!(∂², parent(B), CD.Z, CD.ℓ)
    end
    ∂²
end
