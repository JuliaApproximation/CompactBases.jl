@testset "$(kind) operators" for (kind,Vfun) in [
    ("Linear", (R,fx) -> LinearOperator(R'*QuasiDiagonal(fx)*R, R)),
    ("Diagonal", (R,fx) -> DiagonalOperator(applied(*, R, R \ fx)))
]
    L = 10.0
    f = x -> sin.(2π*x/L)
    oo = x -> 1

    N = 100

    @testset "$(label)" for (label,R) in [
        ("Finite-differences", FiniteDifferences(N, L/(N+1))),
        ("Uniform, staggered finite-differences", StaggeredFiniteDifferences(N, L/(N-1/2))),
        ("Non-uniform, staggered finite-differences", StaggeredFiniteDifferences(0.1, 0.2, 0.01, L)),
        ("FE-DVR", FEDVR(range(0, stop=L, length=N), 4)),
        ("B-splines", BSpline(LinearKnotSet(5, 0, L, N)))
    ]
        x = axes(R,1)
        fx = f.(x)
        c = R \ fx
        o = R \ oo.(x)

        V = Vfun(R, fx)
        m = size(R,2)
        @test size(V) == (m,m)
        @test size(V,1) == size(V,2) == m
        @test eltype(V) == eltype(c)

        u = similar(o)
        mul!(u, V, o)
        @test u ≈ c
    end
end

@testset "ShiftAndInvert" begin
    L = 10.0
    N = 100

    @testset "$(label)" for (label,R,tol) in [
        ("Finite-differences", FiniteDifferences(N, L/(N+1)), 0.004),
        ("Uniform, staggered finite-differences", StaggeredFiniteDifferences(N, L/(N-1/2)), 0.05),
        ("Non-uniform, staggered finite-differences", StaggeredFiniteDifferences(0.1, 0.2, 0.01, L), 0.05),
        ("Implicity finite-differences", ImplicitFiniteDifferences(N, L/(N+1)), 5e-6),
        ("FE-DVR", FEDVR(range(0, stop=L, length=N), 4)[:,2:end-1], 8e-11),
        ("B-splines", BSpline(LinearKnotSet(5, 0, L, N))[:,2:end-1], 5e-13)
    ]
        λ,ϕ,r̃,R,δλ = test_particle_in_a_box(R, L, 5)
        @test norm(δλ) < tol
    end
end
