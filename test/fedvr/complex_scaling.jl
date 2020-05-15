@testset "Complex scaling" begin
    B = FEDVR(range(0,stop=20,length=71), 10)
    C = FEDVR(range(0,stop=20,length=71), 10, t₀=10.0, ϕ=π/3)
    @test C.t₀ ≥ 10.0
    @test_throws ArgumentError FEDVR(range(0,stop=20,length=71), 10, t₀=30.0, ϕ=π/3)
    @test complex_rotate(5, B) == 5
    @test complex_rotate(15, B) == 15
    @test complex_rotate(5, C) == 5
    @test complex_rotate(15, C) ≈ 10 + 5*exp(im*π/3)
end

@testset "Real locations" begin
    rₘₐₓ = 20
    R = FEDVR(range(0,stop=rₘₐₓ,length=11), 10, t₀=rₘₐₓ/2, ϕ=π/3)
    R′ = FEDVR(range(0,stop=rₘₐₓ,length=11), 10)
    @test norm(real_locs(R)-real_locs(R′)) == 0
end
