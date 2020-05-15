@testset "Function interpolation" begin
    rₘₐₓ = 20
    @testset "t₀ = $(t₀), ϕ = $(ϕ)" for (t₀,ϕ) in [(0.0,0.0)# ,
                                                   # (rₘₐₓ/2, π/3)
                                                   ]
        R = FEDVR(range(0,stop=rₘₐₓ,length=11), 10, t₀=t₀, ϕ=ϕ)
        R⁻¹ = pinv(R)
        r = axes(R,1)
        r̃ = range(0,stop=rₘₐₓ,length=1001)
        χ = R[r̃,:]

        fu = r -> r^2*exp(-r)
        fv = r -> r^6*exp(-8r)

        u = R*(R \ fu.(r))
        @test norm(χ * (R⁻¹*u) - fu.(r̃)) < 1e-6
        v = R*(R \ fv.(r))
        @test norm(χ * (R⁻¹*v) - fv.(r̃)) < 5e-5

        @testset "Restricted basis" begin
            R̃ = R[:,2:end-1]
            χ̃ = R̃[r̃,:]

            # # Actually, we want
            # ũ = R̃*(R̃ \ fu.(r))
            # # but R̃*... expands the coefficient vector by two elements
            # # (which are "hidden" by the restriction matrix), which
            # # causes a dimension mismatch below.
            ũ = R̃ \ fu.(r)
            # @test norm(χ̃ * ũ.args[2] - fu.(r̃)) < 5e-6
            @test norm(χ̃ * ũ - fu.(r̃)) < 2e-6
            # ṽ = R̃*(R̃ \ fv.(r))
            ṽ = R̃ \ fv.(r)
            # @test norm(χ̃ * ṽ.args[2] - fv.(r̃)) < 1e-4
            @test norm(χ̃ * ṽ - fv.(r̃)) < 5e-5

            h = r -> (r-rₘₐₓ/2)^2
            c = R̃ \ h.(r)
            # Vandermonde matrix to compare with; its generation is
            # costly, which is why we prefer to compute overlaps with
            # basis functions instead.
            rl = locs(R̃)
            V = R̃[rl,:]
            c̃ = V \ h.(rl)

            @test c ≈ c̃ atol=1e-13
        end
    end
    @testset "Restricted interpolation" begin
        L = 20.0
        n = 10
        t = range(-L/2, stop=L/2, length=n)
        B = FEDVR(t, [4,3,4,2,5,4,2,4,8])
        N = size(B,2)
        xx = axes(B,1)

        Bl = B[:,1:round(Int, 0.6N)+1]
        Bv = B[:,round(Int, 0.4N):end]

        f = sin
        c = B \ f.(xx)
        cl = Bl \ f.(xx)
        cv = Bv \ f.(xx)

        @test length(cl) == size(Bl,2)
        @test c[1:length(cl)-1] == cl[1:end-1]

        @test length(cv) == size(Bv,2)
        @test c[end-length(cv)+1:end] == cv
    end
end
