@testset "Interpolation" begin
    @testset "T = $(T)" for (T,o,γ) in [(Float64,x -> one(Float64),1.0),(ComplexF64,x -> one(ComplexF64),exp(im*π/4))]
        L = 10.0
        f = x -> γ * sin(2π*x/L)

        N = 100

        rr = range(0, stop=L, length=1000)

        @testset "$(label)" for (label,R,tol) in [
            ("Finite-differences", FiniteDifferences(N, L/(N+1)), 0.01),
            ("Uniform, staggered finite-differences", StaggeredFiniteDifferences(N, L/(N-1/2)), 0.04),
            ("Non-uniform, staggered finite-differences", StaggeredFiniteDifferences(0.1, 0.2, 0.01, L), 0.04),
            ("Implicit finite-differences", ImplicitFiniteDifferences(N, L/(N+1)), 0.01),
            ("FE-DVR", FEDVR(range(0, stop=L, length=N), 4), 5e-6),
            ("B-splines", BSpline(LinearKnotSet(5, 0, L, N)), 5e-9)
        ]
            x = axes(R,1)
            fc = R \ f.(x)

            χ = R[rr,:]
            @test χ*fc ≈ f.(rr) atol=tol
            if !(R isa AbstractFiniteDifferences)
                # We don't bother testing the various
                # finite-differences for interpolations of constants,
                # since they do not yet support non-zero Dirichlet
                # boundary conditions.
                oc = R \ o.(x)
                @test norm(χ*oc .- 1) < 5e-14
            end
        end
    end
end
