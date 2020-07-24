@testset "Densities" begin
    f = x -> sin(2π*x)
    g = x -> x*exp(-x)
    h = x -> f(x)*g(x)

    rmin = 0.01
    rmax = 10.0
    ρ = 0.05
    N = ceil(Int, rmax/ρ)
    k = 7
    Nn = 71

    ρmin=rmin
    ρmax=0.5
    α=0.01

    @testset "$R" for (R,kind,rtol) in [
        (FiniteDifferences(N, ρ), :orthogonal_uniform, 1e-14),
        (StaggeredFiniteDifferences(ρmin, ρmax, α, rmax, 0.0), :orthogonal_non_uniform, 1e-14),
        (FEDVR(range(0, stop=rmax, length=Nn), k), :orthogonal_non_uniform, 1e-14),
        (BSpline(LinearKnotSet(k, 0, rmax, Nn)), :non_orthogonal, 1e-5),
    ]
        r = axes(R,1)

        cf = R \ f.(r)
        cg = R \ g.(r)
        ch = R \ h.(r)

        ρ = Density(applied(*,R,cf), applied(*,R,cg))
        ch2 = ρ.ρ

        @test eltype(ρ) == eltype(cf)

        if kind == :orthogonal_uniform
            @test ρ.LV == I
            @test ρ.RV == I
            @test ρ.C == I
        elseif kind == :orthogonal_non_uniform
            @test ρ.LV == I
            @test ρ.RV == I
        end

        @test ch2 ≈ ch atol=1e-14 rtol=rtol
    end

    @testset "Pretty-printing" begin
        R = FiniteDifferences(99, 0.1)
        u = R*rand(size(R,2))
        v = R*rand(size(R,2))
        ρ = Density(u, v)
        fp = FunctionProduct{false}(u, v)

        @test string(ρ) == "99 .* 99 -> 99 FunctionProduct Float64, conjugated (<=> Density)"
        @test string(fp) == "99 .* 99 -> 99 FunctionProduct Float64"

        buf = IOBuffer()
        show(buf, MIME"text/plain"(), ρ)
        @test String(take!(buf)) == """
            99 .* 99 -> 99 FunctionProduct Float64, conjugated (<=> Density); L .* R -> R, with
              L: $(string(R))
              R: $(string(R))"""

        buf = IOBuffer()
        show(buf, MIME"text/plain"(), fp)
        @test String(take!(buf)) == """
            99 .* 99 -> 99 FunctionProduct Float64; L .* R -> R, with
              L: $(string(R))
              R: $(string(R))"""
    end
end
