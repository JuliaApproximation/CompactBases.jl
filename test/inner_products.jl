@testset "Inner products" begin
    rmax = 1.0
    ρ = 0.1
    N = ceil(Int, rmax/ρ)
    k = 7

    @testset "$(grid_type)" for grid_type = [:fd, :sfd, :implicit_fd, :fedvr, :bsplines]
        R = if grid_type == :fd
            FiniteDifferences(N, ρ)
        elseif grid_type == :sfd
            StaggeredFiniteDifferences(N, ρ)
        elseif grid_type == :implicit_fd
            ImplicitFiniteDifferences(N, ρ)
        elseif grid_type == :fedvr
            t = range(0, stop=rmax, length=3)
            FEDVR(t, k)
        elseif grid_type == :bsplines
            t = LinearKnotSet(k, 0, rmax, 3)
            BSpline(t)
        end

        r = axes(R,1)

        c = R \ sin.(2π*r)
        f = R*c
        d = R \ cos.(2π*r)
        g = R*d

        F = R*hcat(c,d)
        FF = F'F
        @test size(FF) == (2,2)
        @test isapprox(FF, 0.5I, rtol=1e-3)

        @test f isa CompactBases.FuncVector
        @test isapprox(norm(f), 1/√2, rtol=1e-3)
        f′ = deepcopy(f)
        normalize!(f′)
        @test isapprox(norm(f′), 1, rtol=1e-3)

        S = R'R

        @test apply(*, c', R', R, c) isa Real
        @test isapprox(apply(*, c', R', R, c), 0.5, rtol=1e-3)
        @test dot(c, S, c) isa Real
        @test isapprox(dot(c, S, c), 0.5, rtol=1e-3)
        @test dot(c, S, d) isa Real
        @test isapprox(dot(c, S, d), 0.0, atol=1e-14)
        @test f'f isa Real
        @test isapprox(f'f, 0.5, rtol=1e-3)
        @test f'g isa Real
        @test isapprox(f'g, 0.0, atol=1e-14)
        @test g'g isa Real
        @test isapprox(g'g, 0.5, rtol=1e-3)
    end

    @warn "Need to test inner products with restricted bases as well"
end
