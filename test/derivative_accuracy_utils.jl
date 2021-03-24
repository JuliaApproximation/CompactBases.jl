using LinearAlgebra

using Test
using PrettyTables
using ProgressMeter
using UnicodePlots

using RollingFunctions

function test_derivatives(R::AbstractQuasiMatrix,
                          f::Function, g::Function, h::Function)
    r = axes(R,1)

    D = Derivative(r)

    ∇ = R' * D * R
    ∇² = R' * D' * D * R

    S = R'R
    nonortho = bandwidths(S) != (0,0)
    S⁻¹ = nonortho ? qr(S) : I

    fv = R \ f.(r)

    gv = similar(fv)
    mul!(gv, ∇, fv)
    nonortho && ldiv!(S⁻¹, gv)
    δgv = R \ g.(r) - gv

    hv = similar(gv)
    mul!(hv, ∇, gv)
    nonortho && ldiv!(S⁻¹, hv)
    δhv = R \ h.(r) - hv

    hv′ = similar(fv)
    mul!(hv′, ∇², fv)
    nonortho && ldiv!(S⁻¹, hv′)
    δhv′ = R \ h.(r) - hv′

    R,fv,gv,hv,hv′,δgv,δhv,δhv′
end

function estimate_convergence_rate(ρ, ϵ)
    i = findall(!iszero, ϵ)
    x = log10.(ρ[i])
    y = log10.(abs.(ϵ[i]))

    f = (x,y) -> ([x ones(length(x))] \ y)[1]

    # Find all slopes for consecutive groups of three errors
    slopes = rolling(f, x, y, 3)
    e = findfirst(s -> s < 0, slopes)
    if e !== nothing && e > 5
        slopes = slopes[1:max(1,e-1)]
    end

    maximum(slopes)
end

function error_slopes(hs, ϵ, names, verbosity=0)
    loghs = log10.(hs)
    logϵ = log10.(abs.(ϵ))
    n = size(logϵ,2)

    slopes = zeros(n)
    for j = 1:n
        slopes[j] = estimate_convergence_rate(hs, ϵ[:,j])
    end

    if verbosity > 0
        pretty_table([hs ϵ], vcat("h", ["$(name) [$(slope)]"
                                        for (name,slope) in zip(names,slopes)]))
        plt = lineplot(loghs, logϵ[:,1], name=names[1], xlabel="log10(h)", ylabel="log10(error)")
        for j = 2:n
            lineplot!(plt, loghs, logϵ[:,j], name=names[j])
        end
        # for (q,o) in enumerate(["linear", "quadratic", "cubic", "quartic", "quintic"])
        #     lineplot!(plt, loghs, q*loghs, name=o)
        # end
        display(plt)
        println()
        println()
    end

    slopes
end

function derivative_test_functions(d)
    a,b = √d*[-1,1]

    # The functions vanish at the boundaries; Dirichlet0 boundary
    # conditions. They /cannot/ be evaluated at the boundaries, since
    # that results in NaNs.
    f = x -> exp(-1/(d-x^2))
    g = x -> -2*exp(-1/(d-x^2))*x/((d-x^2)^2)
    h = x -> -2*exp(-1/(d-x^2))*(d^2 + 2*(d-1)*x^2-3x^4)/((d-x^2)^4)

    f,g,h,a,b
end

function compute_derivative_errors(fun::Function, Ns,
                                   f::Function, g::Function, h::Function,
                                   verbosity=0)
    p = Progress(length(Ns))
    errors = map(Ns) do N
        R = fun(N)
        R,fv,gv,hv,hv′,δgv,δhv,δhv′ = test_derivatives(R, f, g, h)
        S = R'R
        ProgressMeter.next!(p)
        [√(dot(δgv, S, δgv)) √(dot(δhv, S, δhv)) √(dot(δhv′, S, δhv′))]
    end |> e -> vcat(e...)

    ϵg = errors[:,1]
    ϵh = errors[:,2]
    ϵh′ = errors[:,3]

    hs = 1.0 ./ Ns
    pg,ph,ph′ = error_slopes(hs, errors, ["δg", "δh", "δh′"], verbosity)

    hs,ϵg,ϵh,ϵh′,pg,ph,ph′
end

function norm_rot!(v)
    normalize!(v)
    vc = v.args[2]
    vc[:] *= sign(vc[1])
    v
end

function diagonalize_hamiltonian(H, R::B, nev, σ; method=:arnoldi_shift_invert) where {B<:AbstractQuasiMatrix}
    A,target = if method == :arnoldi_shift_invert
        ShiftAndInvert(H, R, σ), LR()
    else
        H, SR()
    end
    schurQR,history = partialschur(A, nev=nev, which=target)

    ϕ = [norm_rot!(R*schurQR.Q[:,j]) for j = 1:nev]

    θ = schurQR.eigenvalues
    λ = if method == :arnoldi_shift_invert
        σ .+ inv.(θ)
    else
        θ
    end

    λ,ϕ
end

function get_kinetic_operator(R::B, Z=0, ℓ=0) where {B<:AbstractQuasiMatrix}
    D̃ = Derivative(Base.axes(R,1))
    D = !iszero(Z) ? CoulombDerivative(D̃, Z, ℓ) : D̃
    apply(*, R', D', D, R) / -2
end

function test_particle_in_a_box(R, L, nev; kwargs...) where {B<:AbstractQuasiMatrix}
    Tm = get_kinetic_operator(R)
    λ,ϕ = diagonalize_hamiltonian(Tm, R, nev, 0.0; kwargs...)

    r̃ = CompactBases.locs(R)

    n = 1:nev
    δλ = abs.(λ - n.^2*π^2/(2L^2))

    # Could/should also test eigenvectors

    λ,ϕ,r̃,R,δλ
end

function test_singular_scheme(R, Z, ℓ, nev; kwargs...) where {B<:AbstractQuasiMatrix}
    r = axes(R,1)

    n = 1:nev
    λₐ = -inv.(2(n .+ ℓ).^2)

    coulomb(r) = -inv(r)
    centrifugal(r) = ℓ*(ℓ+1)/2r^2

    Tm = get_kinetic_operator(R, Z, ℓ) + R'*QuasiDiagonal(centrifugal.(r))*R
    V = R'*QuasiDiagonal(coulomb.(r))*R
    H = Tm + V
    λ,ϕ = diagonalize_hamiltonian(H, R, nev, 1.1λₐ[1]; kwargs...)

    r̃ = CompactBases.locs(R)

    δλ = abs.(λ - λₐ)

    # Could/should also test eigenvectors

    λ,ϕ,r̃,R,δλ
end

function compute_diagonalization_errors(fun::Function, problem::Function, Ns, args...; verbosity=0, kwargs...) where B
    p = Progress(length(Ns))
    errors = map(Ns) do N
        t = time()
        λ,ϕ,r̃,R,δλ = problem(fun(N), args...; kwargs...)
        elapsed = time()-t
        ProgressMeter.next!(p)
        vcat(δλ,elapsed)'
    end |> e -> vcat(e...)

    errors,elapsed = errors[:,1:end-1],errors[:,end]

    hs = 1.0 ./ Ns
    slopes = error_slopes(hs, errors, ["$j" for j in 1:size(errors,2)], verbosity)

    errors,slopes,elapsed
end
