@testset "Simple tests" begin
    B = FEDVR(range(0,stop=20,length=71), 10)
    C = FEDVR(range(0,stop=20,length=71), 10)

    @test B == C
    @test axes(B) == (Inclusion(0..20), 1:631)
    @test size(B) == (ℵ₁, 631)

    @test nel(B) == 70

    @test B[0.0,1] ≠ 0.0
    @test B[20.0,end] ≠ 0.0

    @testset "Restricted bases" begin
        B̃ = B[:, 2:end-1]
        @test axes(B̃) == (Inclusion(0..20), 1:629)
        @test size(B̃) == (ℵ₁, 629)

        @test leftendpoint(B̃) == 0
        @test rightendpoint(B̃) == 20

        @test B̃[0.0,1] == 0.0
        @test B̃[20.0,end] == 0.0

        @test B isa FEDVR
        @test B isa CompactBases.FEDVROrRestricted
        @test B' isa QuasiAdjoint{<:Any,<:FEDVR}
        @test B' isa CompactBases.AdjointFEDVROrRestricted
        @test !(B̃ isa FEDVR)
        @test B̃ isa CompactBases.RestrictedFEDVR
        @test B̃ isa CompactBases.FEDVROrRestricted
        @test B̃' isa CompactBases.AdjointRestrictedFEDVR
        @test B̃' isa CompactBases.AdjointFEDVROrRestricted

        @test unrestricted_basis(B) == B
        @test unrestricted_basis(B̃) == B

        @test restriction_extents(B) == (0,0)
        @test restriction_extents(B̃) == (1,1)

        @test locs(B̃) == locs(B)[2:end-1]
    end
end

@testset "Pretty printing" begin
    B = FEDVR(range(0,stop=20,length=71), 10)
    C = FEDVR(range(0,stop=20,length=71), 10, t₀=10.0, ϕ=π/3)
    @test occursin("FEDVR{Float64} basis with 70 elements on 0.0..20.0", string(B))
    @test occursin("FEDVR{Complex{Float64}} basis with 70 elements on 0.0..20.0 with ECS @ 60.00° starting at 10", string(C))
end

@testset "Element access" begin
    B = FEDVR(range(0,stop=20,length=71), 10)
    @test first(@elem(B, x, 1)) == 0
    @test last(@elem(B, x, 1)) == first(@elem(B, x, 2))
    @test last(@elem(B, x, 69)) == first(@elem(B, x, 70))
    @test_throws BoundsError @elem(B, x, 71)
end
