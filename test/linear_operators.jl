@testset "Linear operators" begin

end

@testset "Diagonal operators" begin
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
        c = R \ f.(x)
        o = R \ oo.(x)
        V = DiagonalOperator(applied(*, R, c))
        u = similar(o)
        mul!(u, V, o)
        @test u ≈ c
    end
end

@testset "ShiftAndInvert" begin

end
