@testset "Basis functions" begin
    rₘₐₓ = 10
    ρ = rand()
    N = ceil(Int, rₘₐₓ/ρ + 1/2)

    for B in [StaggeredFiniteDifferences(N, ρ), FiniteDifferences(N, 1.0)]
        x = locs(B)
        χ = B[x, :]
        @test χ isa AbstractSparseMatrix
        @test χ ≈ Diagonal(ones(N))
    end

    @testset "Basis functions in restricted bases" begin
        N = 100
        L = 20.0
        δx = L/(N+1)

        B = FiniteDifferences((1:N) .- div(N,2), δx)
        xx = axes(B,1)
        a = xx.domain
        x = range(a.left, stop=a.right, length=max(N,101))
        χ = B[x,:]
        @test χ isa AbstractSparseMatrix

        sell = 1:round(Int, 0.6N)
        selv = round(Int, 0.4N):N

        Bl = B[:,sell]
        Bv = B[:,selv]

        χl = Bl[x, :]
        @test χl isa AbstractSparseMatrix
        @test χl == χ[:,sell]
        χv = Bv[x, :]
        @test χv == χ[:,selv]

        @test Bv[x, 4:8] == χ[:,selv[4:8]]
    end
end
