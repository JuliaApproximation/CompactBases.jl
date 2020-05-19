@testset "Densities" begin
    R = FiniteDifferences(20,1.0)
    R⁻¹ = pinv(R)
    r̃ = locs(R)
    χ = R[r̃,:]

    r = axes(R,1)

    fu = r -> r^2*exp(-r)
    u = R*(R\fu.(r))

    fv = r -> r^6*exp(-r)
    v = R*(R\fv.(r))

    w = u .* v
    fw = r -> fu(r)*fv(r)

    @test norm(χ * (R⁻¹*w) - fw.(r̃)) == 0

    y = R*rand(ComplexF64, size(R,2))
    y² = y .* y
    @test all(isreal.(R⁻¹*y²))
    @test all(R⁻¹*y² .== abs2.(R⁻¹*y))

    @testset "Lazy densities" begin
        uv = u .⋆ v
        @test uv isa CompactBases.FDDensity

        w′ = similar(u)
        copyto!(w′, uv)
        @test norm(χ * (R⁻¹*w′) - fw.(r̃)) == 0

        uu = R*repeat(R⁻¹*u,1,2)
        vv = R*repeat(R⁻¹*v,1,2)
        uuvv = uu .⋆ vv
        ww′ = similar(uu)
        copyto!(ww′, uuvv)

        @test norm(χ * (R⁻¹*ww′) .- fw.(r̃)) == 0

        yy = y .⋆ y
        @test yy isa CompactBases.FDDensity
        wy = similar(y)
        copyto!(wy, yy)
        @test all(isreal.(R⁻¹*wy))
        @test all(R⁻¹*wy .== abs2.(R⁻¹*y))
    end

    @testset "Restricted bases" begin
        N = 100
        n = 10
        dr = 0.1

        R = StaggeredFiniteDifferences(N, dr)
        R̃ = R[:,1:n]

        f = 3

        ϕ = applied(*, R, f*ones(N))
        ϕ̃ = applied(*, R̃, f*ones(n))

        ρ = similar(ϕ)
        ρ̃ = similar(ϕ̃)

        @test ϕ .⋆ ϕ isa CompactBases.FDDensity
        @test ϕ̃ .⋆ ϕ̃ isa CompactBases.FDDensity

        copyto!(ρ, ϕ .⋆ ϕ)
        @test ρ.args[2] == f^2*ones(N)
        copyto!(ρ̃, ϕ̃ .⋆ ϕ̃)
        @test ρ̃.args[2] == f^2*ones(n)
    end
end
