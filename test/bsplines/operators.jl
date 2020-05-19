@testset "Operators" begin
    @testset "Coulomb operator" begin
        k = 7
        N = 31

        a,b = 0,70

        coulomb(r) = -1/r
        @testset "$name knot set" for (name, t, tol) in [("Linear", LinearKnotSet(k, a, b, N), 9e-3),
                                                         ("Exponential", ExpKnotSet(k, -1.0, log10(b), N), 2e-7)]
            B = BSpline(t)[:,2:end-1]
            S = B'B

            r = axes(B,1)

            f = B \ (r -> r^2 * exp(-r)).(r)
            g = B \ (r -> -r*exp(-r)).(r)

            V = B'*QuasiDiagonal(coulomb.(r))*B
            g̃ = S \ V*f

            @test g ≈ g̃ atol=tol
        end
    end
end
