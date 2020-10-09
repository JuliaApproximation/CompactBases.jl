@testset "Inner products" begin
    rₘₐₓ = 300
    ρ = rand()
    N = ceil(Int, rₘₐₓ/ρ + 1/2)

    for B in [StaggeredFiniteDifferences(N, ρ), FiniteDifferences(N, 1.0)]
        S = B'B
        for T in [Float64,ComplexF64]
            vv = rand(T, size(B,2))
            v = B*vv
            lv = B ⋆ vv
            normalize!(v)

            @test norm(v) ≈ 1.0
            @test applied(*, v'.args..., v.args...) isa CompactBases.FDInnerProduct # {T,Float64,StaggeredFiniteDifferences{Float64,Int}}
            @test v'v ≈ 1.0
            @test dot(vv,S,vv) ≈ 1.0

            lazyip = lv' ⋆ lv

            @test lazyip isa CompactBases.LazyFDInnerProduct
            @test materialize(lazyip) ≈ 1.0
        end
    end

    @testset "Restricted bases" begin
        N = 100
        n = 10
        dr = 0.1

        R = StaggeredFiniteDifferences(N, dr)
        R̃ = R[:,1:n]

        f = 3

        ϕ = R*(f*ones(N))
        ϕ̃ = R̃*(f*ones(n))

        @test norm(ϕ) ≈ f*√(N*dr)
        @test norm(ϕ̃) ≈ f*√(n*dr)

        @test ϕ'ϕ == (f^2*N*dr)
        @test ϕ̃'ϕ̃ == (f^2*n*dr)

        @test apply(*, ϕ', ϕ) == (f^2*N*dr)
        @test apply(*, ϕ̃', ϕ̃) == (f^2*n*dr)

        normalize!(ϕ)
        @test norm(ϕ) ≈ 1.0
        normalize!(ϕ̃)
        @test norm(ϕ̃) ≈ 1.0
    end
end
