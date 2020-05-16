@testset "Scalar operators" begin
    B = FEDVR(1.0:7, 4)
    B̃ = B[:,2:end-1]
    r = axes(B,1)

    V = B'*QuasiDiagonal(identity.(r))*B
    @test V == Diagonal(B.x)

    Ṽ = apply(*, B̃', QuasiDiagonal(identity.(r)), B̃)
    @test Ṽ == Diagonal(B.x[2:end-1])
end

@testset "Mass matrices and inverses" begin
    t = range(0,stop=20,length=5)
    @testset "Order $k" for k ∈ [2,5]
        R = FEDVR(t, k)
        R̃ = R[:,2:end-1]

        m = size(R,2)
        n = size(R̃,2)

        @testset "Mass matrices" begin
            d = R'R
            @test d == Diagonal(Ones(m))

            d̃ = R̃'R̃
            @test d̃ == Diagonal(Ones(n))

            @test R'R̃ == (R̃'R)' == BandedMatrix((-1 => ones(Int,n),), (m,n), (1,-1))
            @test R̃'R[:,1:end-1] == BandedMatrix((1 => ones(Int,n),), (n,m-1), (-1,1))
        end

        @testset "Inverses" begin
            R⁻¹ = pinv(R)
            # R̃⁻¹ = pinv(R̃)

            @test R⁻¹*R === I
            @test R*R⁻¹ === I

            cu = rand(size(R,2))
            cv = rand(size(R,2))
            cuv = [cu cv]

            u = R*cu
            v = R*cv
            uv = R*cuv

            @test R⁻¹*u === cu
            @test R⁻¹*v === cv
            @test R⁻¹*uv === cuv

            ut = u'
            # Fails with: ERROR: MethodError: no method matching axes(::UniformScaling{Float64}, ::Int64)
            # @test ut*R⁻¹' === ut.args[1]

            @warn "Need to implement/test basis inverses for restricted bases"
        end
    end
end
