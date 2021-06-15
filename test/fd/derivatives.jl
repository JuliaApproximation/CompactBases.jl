@testset "Derivatives" begin
    @testset "Normal finite-differences" begin
        B = FiniteDifferences(5,1.0)
        D = Derivative(axes(B,1))

        ∇ = applied(*, B', D, B)
        ∇² = applied(*, B', D', D, B)

        ∂ = materialize(∇)
        T = materialize(∇²)

        @test ∂ isa Tridiagonal
        @test T isa SymTridiagonal

        @test all(diag(∂) .== 0)
        @test all(diag(∂,1) .== 0.5)
        @test all(diag(∂,-1) .== -0.5)

        @test all(diag(T) .== -2)
        @test all(diag(T,1) .== 1)
        @test all(diag(T,-1) .== 1)
    end

    @testset "Restricted bases" begin
        N = 10
        ρ = 1.0

        R = StaggeredFiniteDifferences(N, ρ)
        r = axes(R, 1)
        D = Derivative(r)

        ∇ = R'D*R
        ∇² = R'D'D*R

        for (sela,selb) in Iterators.product([1:10, 3:6, 8:10, 5:10, 4:5], [1:10, 3:6, 8:10, 5:10, 4:5])
            # In case we have partially overlapping restrictions, the
            # bandwidth will not be the full three-point stencil, and if we
            # have no overlap at all, the bandwidth will be zero.
            overlap = sela ∩ selb
            expected_bandwidth = isempty(overlap) ? 0 : min(3, length(overlap)+1)

            Ra = R[:,sela]
            Rb = R[:,selb]

            ∂ = Ra'D*Rb
            @test ∂ == Matrix(∇)[sela,selb]
            if sela == selb
                @test bandwidths(∂) == (1,1)
            end
            @test ∂ isa BandedMatrix
            @test length(bandrange(∂)) == expected_bandwidth

            ∂² = Ra'D'*D*Rb
            @test ∂² == Matrix(∇²)[sela,selb]
            if sela == selb
                @test bandwidths(∂) == (1,1)
            end
            @test ∂² isa BandedMatrix
            @test length(bandrange(∂²)) == expected_bandwidth
        end
    end

    @testset "Implicit finite differences" begin
        N = 10
        ρ = 0.3
        Z = 1.0
        R = ImplicitFiniteDifferences(N, ρ)
        r = axes(R, 1)
        @testset "Singular origin=$(singular_origin)" for singular_origin = [true, false]
            D̃ = Derivative(r)
            D = singular_origin ? CoulombDerivative(D̃, Z, 0) : D̃

            λ,δβ₁ = if singular_origin
                √3 - 2, -Z*ρ*inv(12 - 10*Z*ρ)
            else
                0,0
            end

            ∂ = R'D*R
            @test ∂ isa CompactBases.ImplicitDerivative

            @test ∂.Δ isa Tridiagonal
            @test ∂.Δ[1,1] ≈ λ/2ρ
            @test all(iszero, ∂.Δ.d[2:end])
            @test all(∂.Δ.dl .≈ -1/2ρ)
            @test all(∂.Δ.du .≈ 1/2ρ)

            @test ∂.M isa SymTridiagonal
            @test ∂.M[1,1] ≈ (4+λ)/6
            @test all(∂.M.dv[2:end] .≈ 4/6)
            @test all(∂.M.ev .≈ 1/6)


            ∂² = R'D'D*R
            @test ∂ isa CompactBases.ImplicitDerivative

            @test ∂².Δ isa SymTridiagonal
            @test ∂².Δ[1,1] ≈ -2*(1+δβ₁)/ρ^2
            @test all(∂².Δ.dv[2:end] .≈ -2/ρ^2)
            @test all(∂².Δ.ev .≈ 1/ρ^2)

            @test ∂².M isa SymTridiagonal
            @test ∂².M[1,1] ≈ (-1/2)*(-(10-2δβ₁)/6)
            @test all(∂².M.dv[2:end] .≈ (-1/2)*(-10/6))
            @test all(∂².M.ev .≈ (-1/2)*(-1/6))
        end

        @testset "Restricted bases" begin
            R̃ = R[:, 1:5]
            D = Derivative(axes(R̃, 1))

            @test_throws DimensionMismatch R̃'D*R
            @test_throws DimensionMismatch R'D*R̃
            @test_throws DimensionMismatch R̃'D'D*R
            @test_throws DimensionMismatch R'D'D*R̃

            ∂ = R̃'D*R̃
            @test ∂ isa CompactBases.ImplicitDerivative
            @test size(∂) == (5,5)

            ∂² = R̃'D'D*R̃
            @test ∂² isa CompactBases.ImplicitDerivative
            @test size(∂²) == (5,5)

            @warn "Need to test correctness of ImplicitDerivative in restricted bases"
        end
    end
end

@testset "Derivative accuracy" begin
    Ns = 2 .^ (5:20)

    d = 1.0
    f,g,h,a,b = derivative_test_functions(d)

    @testset "kind = $B" for (order,B) in [(2,FiniteDifferences),
                                           (4,ImplicitFiniteDifferences)]
        @info "$B derivative accuracy"
        hs,ϵg,ϵh,ϵh′,pg,ph,ph′ = compute_derivative_errors(Ns, f, g, h, 1) do N
            B(range(a, stop=b, length=N+2)[2:end-1])
        end

        @test isapprox(pg, order, atol=0.03) || pg > order
        @test isapprox(ph, order, atol=0.03) || ph > order
    end
    @testset "kind = StaggeredFiniteDifferences" begin
        dd = b-a
        @testset "$(label)" for (order,label,fun) = [(2.0, "Uniform", N -> (N,dd/(N+1))),
                                                     (2.0, "Non-uniform", N -> (dd*(1 .- cos.(π*range(0,stop=1,length=N+2)[2:end-1]))/2,))]
            @info "$(label) StaggeredFiniteDifferences derivative accuracy"
            hs,ϵg,ϵh,ϵh′,pg,ph,ph′ = compute_derivative_errors(Ns, x -> f(x-dd/2), x -> g(x-dd/2), x -> h(x-dd/2), 1) do N
                StaggeredFiniteDifferences(fun(N)...)
            end

            @test isapprox(pg, order, atol=0.03) || pg > order
            @test isapprox(ph, order, atol=0.03) || ph > order
        end
    end
end

@testset "Particle-in-a-box" begin
    @testset "Eigenvalues convergence rate" begin
        Ns = 2 .^ (7:14)
        nev = 3
        L = 1.0
        @testset "kind = $B" for (order,B) in [(2,FiniteDifferences),
                                               (4,ImplicitFiniteDifferences)]
            @info "$B particle-in-a-box eigenvalues convergence rate"
            ϵλ,slopes,elapsed =
                compute_diagonalization_errors(N -> B(N, L/(N+1)), test_particle_in_a_box,
                                               Ns, L, nev,
                                               verbosity=1)
            for (p,o) in zip(slopes, [order,order-1,order-1])
                @test isapprox(p, o, atol=0.1) || p > o
            end
        end
    end
end

# This tests the discretization of the Laplacian, especially near
# the origin where the potential term is singular.
@testset "Hydrogen bound states" begin
    @testset "Eigenvalues and eigenvectors" begin
        rₘₐₓ = 300
        ρ = 0.25
        N = ceil(Int, rₘₐₓ/ρ + 1/2)

        R = StaggeredFiniteDifferences(N, ρ)
        r = axes(R, 1)

        D = CoulombDerivative(Derivative(Base.axes(R,1)), 1, 0)
        Tm = (R'D'D*R) / -2
        @test Tm isa SymTridiagonal

        V = R'*QuasiDiagonal((r -> -inv(r)).(r))*R

        H = Tm + V

        ee = eigen(H)

        n = 1:10
        λₐ = -inv.(2n.^2)

        abs_error = ee.values[n] - λₐ
        rel_error = abs_error ./ abs.(1e-10 .+ abs.(λₐ))

        @test abs_error[1] < 3e-5
        @test all(rel_error .< 1e-3)

        r̃ = locs(R)
        # Table 2.2 Foot (2005)
        Rₐ₀ = [2exp.(-r̃),
               -(1/2)^1.5 * 2 * (1 .- r̃/2).*exp.(-r̃/2), # Why the minus?
               (1/3)^1.5 * 2 * (1 .- 2r̃/3 .+ 2/3*(r̃/3).^2).*exp.(-r̃/3)]
        expected_errors = [1e-3,1e-3,2e-3]

        for i = 1:3
            v = ee.vectors[:,i]
            N = norm(v)*√ρ
            # The sign from the diagonalization is arbitrary; make max lobe positive
            N *= sign(v[argmax(abs.(v))])
            abs_error = v/N .- r̃.*Rₐ₀[i]
            @test norm(abs_error)/abs(1e-10+norm(r̃.*Rₐ₀[i])) < expected_errors[i]
        end
    end

    @testset "Eigenvalues convergence rate" begin
        Ns = 2 .^ (7:14)
        nev = 3
        rₘₐₓ = 100.0
        Z = 1.0
        ℓ = 0

        @testset "kind = $B" for (order,B) in [(2,StaggeredFiniteDifferences),
                                               (4,ImplicitFiniteDifferences)]
            @info "$B hydrogen eigenvalues convergence rate"
            ϵλ,slopes,elapsed = compute_diagonalization_errors(test_singular_scheme, Ns, Z, ℓ, nev, verbosity=1) do N
                ρ = B == StaggeredFiniteDifferences ? rₘₐₓ/(N+0.5) : rₘₐₓ/(N+1)
                B(N, ρ)
            end
            for p in slopes
                @test isapprox(p, order, atol=0.04) || p > order
            end
        end
    end
end
