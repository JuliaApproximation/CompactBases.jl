@testset "Inner products" begin
    @testset "Vectors" begin
        t = range(0,stop=1,length=11)
        B = FEDVR(t, 4)
        B̃ = B[:,1:5]

        uv = ones(size(B,2))
        vv = ones(size(B̃,2))

        @testset "Direct inner products" begin
            u = B*uv
            v = B̃*vv

            @test u'v == 5
            @test (uv'*(B'B̃)*vv)[1] == 5
            @test norm(v) ≈ √5
        end
    end

    @testset "Matrices" begin
        t = range(-1.0,stop=1.0,length=5)
        R = FEDVR(t, 4)[:,2:end-1]

        n = 3
        Φ = rand(ComplexF64, size(R,2), n)
        for i = 1:n
            # Orthogonalize against all previous vectors
            for j = 1:i-1
                c = Φ[:,j]'Φ[:,i]/(Φ[:,j]'Φ[:,j])
                Φ[:,i] -= c*Φ[:,j]
            end
            Φ[:,i] /= norm(R*Φ[:,i])
        end

        @testset "Direct inner products" begin
            Φ̂ = R*Φ
            Φ̂'Φ̂
            @test Φ̂'Φ̂ ≈ I
        end
    end
end
