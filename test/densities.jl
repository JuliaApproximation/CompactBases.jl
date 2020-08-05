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

    @testset "Multiple functions" begin
        L = 6.0
        N = 1000
        ρ = L/(N+1)
        @testset "$R" for R in [
            FiniteDifferences(N, ρ),
            StaggeredFiniteDifferences(0.01, 0.2, 0.1, L),
            FEDVR(range(0, stop=L, length=20), 7)[:,2:end-1],
            BSpline(LinearKnotSet(7, 0, L, 200))
        ]
            r = axes(R,1)
            f = exp.(-(r .- L/2).^2/0.5)
            gs = [
                sin.(1π/L*r),
                sin.(2π/L*r),
                sin.(3π/L*r),
                sin.(4π/L*r),
                sin.(5π/L*r),
                sin.(6π/L*r)
            ]
            w = r -> 1 / r

            cf = R \ f
            cgs = reduce(hcat, [R \ g for g in gs])
            wr = w.(r)

            u = applied(*, R, cf)
            v = applied(*, R, cgs)

            # One times many functions
            ρuv = Density(u, v, w=w)
            uv_refs = reduce(hcat, [R \ (f .* g .* wr) for g in gs])
            @test ρuv.ρ ≈ uv_refs

            # Many times many functions
            ρvv = Density(v, v, w=w)
            vv_refs = reduce(hcat, [R \ (g .* g .* wr) for g in gs])
            @test ρvv.ρ ≈ vv_refs

            @test_throws DimensionMismatch Density(applied(*, R, rand(size(R, 2), 2)), applied(*, R, rand(size(R, 2), 3)))
        end
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

        u = R*rand(size(R,2), 1)
        v = R*rand(size(R,2), 5)
        ρ = Density(u, v)
        @test string(ρ) == "(99, :) .* (99, :) -> (99, 5) FunctionProduct Float64, conjugated (<=> Density)"

        u = R*rand(size(R,2), 5)
        v = R*rand(size(R,2), 5)
        ρ = Density(u, v)
        @test string(ρ) == "(99, :) .* (99, :) -> (99, 5) FunctionProduct Float64, conjugated (<=> Density)"

        R = BSpline(LinearKnotSet(7, 0, 6.0, 200))

        u = R*rand(size(R,2))
        v = R*rand(size(R,2), 5)
        ρ = Density(u, v)
        @test string(ρ) == "206 .* (206, 5) -> (206, 5) FunctionProduct Float64, conjugated (<=> Density)"

        u = R*rand(size(R,2), 5)
        v = R*rand(size(R,2))
        ρ = Density(u, v)
        @test string(ρ) == "(206, 5) .* 206 -> (206, 5) FunctionProduct Float64, conjugated (<=> Density)"

        u = R*rand(size(R,2), 1)
        v = R*rand(size(R,2), 5)
        ρ = Density(u, v)
        @test string(ρ) == "(206, 1) .* (206, 5) -> (206, 5) FunctionProduct Float64, conjugated (<=> Density)"

        u = R*rand(size(R,2), 5)
        v = R*rand(size(R,2), 1)
        ρ = Density(u, v)
        @test string(ρ) == "(206, 5) .* (206, 1) -> (206, 5) FunctionProduct Float64, conjugated (<=> Density)"

        u = R*rand(size(R,2), 5)
        v = R*rand(size(R,2), 5)
        ρ = Density(u, v)
        @test string(ρ) == "(206, 5) .* (206, 5) -> (206, 5) FunctionProduct Float64, conjugated (<=> Density)"
    end
end
