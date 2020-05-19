@testset "Splines" begin
    @testset "Discontinuous knot set" begin
        t = ArbitraryKnotSet(3, [0.0, 1, 1, 3, 4, 6], 1, 3)
        B = BSpline(t)

        @test all(B[0.0, :] .== 0)
        @test all(B[0.0:0.5:1.0, :][1, :] .== 0)

        @test B.B == B[B.x,:]
        # Reference generated using AnalyticBSplines.jl
        @test B.S ≈ [ 3//5    2//9     2//45    0//1    0//1
                      2//9    7//15   83//270   1//270  0//1
                      2//45  83//270  26//27   83//270  2//45
                      0//1    1//270  83//270   7//15   2//9
                      0//1    0//1     2//45    2//9    2//5]
    end

    @testset "Cardinal splines" begin
        for k = 1:5
            t = ArbitraryKnotSet(k, 1.0:6.0, 1, 1)
            R = BSpline(t)
            @test size(R.B,2) == 6-k
            @test size(R.S) == (6-k,6-k)
            @test (bandwidth(R.S,1),bandwidth(R.S,2)) == (k-1,k-1)
        end
    end

    @testset "Evaluate B-splines" begin
        k = 3
        t = LinearKnotSet(k, 0, 1, 2)
        B = BSpline(t)
        @testset "Eval on subintervals" begin
            x₁ = range(-1,stop=-0.5,length=10)
            χ₁ = B[x₁, :]
            @test norm(χ₁) == 0

            function testbasis(x)
                χ = B[x, :]
                B̃ = spzeros(Float64, length(x), 4)
                B̃[:,1] = (x .>= 0) .* (x .< 0.5) .* ((2x .- 1).^2)
                B̃[:,2] = (x .>= 0) .* (x .< 0.5) .* (2/3*(1 .- (3x .- 1).^2)) +
                    (x .>= 0.5) .* (x .< 1)  .* (2*(x .- 1).^2)
                B̃[:,3] = (x .>= 0) .* (x .< 0.5) .* (2*x.^2) +
                    (x .>= 0.5) .* (x .< 1)  .* (2/3*(1 .- (3x .- 2).^2))
                B̃[:,4] = (x .>= 0.5) .* (x .<= 1) .* ((2x .- 1).^2)
                for j = 1:4
                    δ,δr = vecdist(χ[:,j],B̃[:,j])
                    @test δ < 1e-15
                    @test δr < 1e-15
                end
            end
            testbasis(range(0,stop=1,length=50))
            testbasis(range(-0.5,stop=0.6,length=40))
            testbasis(range(0.5,stop=1.6,length=40))
        end
    end

    @testset "Restricted basis" begin
        t = LinearKnotSet(3, 0.0, 1.0, 6)
        B = BSpline(t,3)
        B̃ = B[:,2:end-1]
        @test size(B̃) == (ℵ₁, 6)
        χ̃ = B̃[CompactBases.locs(B), :]
        @test size(χ̃,2) == 6
        @test χ̃ == B.B[:,2:end-1]
    end

    @testset "Evaluate spline" begin
        @testset "Linear knot set" begin
            k = 4
            t = LinearKnotSet(k, 0, 1, 5)
            B = BSpline(t)
            x = range(0, stop=1, length=301)
            χ = B[x, :]

            i = 1:size(B,2)
            c = sin.(i)

            # Whether we evaluate the spline or the basis functions should
            # not make a difference.
            s = (B*c)[x]
            s′ = χ*c
            @test all(isfinite, s)
            @test all(isfinite, s′)
            @test all(!isnan, s)
            @test all(!isnan, s′)
            @test s == s′
        end
        @testset "Discontinuous knot set" begin
            t = ArbitraryKnotSet(3, [0.0, 1, 1, 3, 4, 6], 1, 3)
            B = BSpline(t,3)

            χ = B[B.x, :]
            @test all(isfinite, χ)
            @test all(!isnan, χ)
        end
    end

    @testset "Mass matrices" begin
        k = 4
        t = LinearKnotSet(k, 0, 1, 5)
        B = BSpline(t)
        @test B'B == B.S
        @test_broken B.S isa Symmetric

        B̃ = B[:,2:end-1]
        @test B̃'B̃ == B.S[2:end-1,2:end-1]

        @test B̃'B == B.S[2:end-1,:]

        # TODO: Test mass matrices between different B-splines bases
        # but same knot set/quadrature points.
    end

    @testset "Function interpolation" begin
        B = BSpline(LinearKnotSet(7, 0, 1, 10))
        x = axes(B,1)
        c = B \ sin.(2π*x)
        x̃ = range(0, stop=1, length=301)
        χ = B[x̃,:]
        @test all(vecdist(χ*c, sin.(2π*x̃)) .< (3e-6, 2e-7))

        B̃ = B[:,2:end-1]
        c̃ = B̃ \ sin.(2π*x)
        χ̃ = B̃[x̃,:]
        @test all(vecdist(χ̃*c̃, sin.(2π*x̃)) .< (3e-6, 2e-7))
    end
end
