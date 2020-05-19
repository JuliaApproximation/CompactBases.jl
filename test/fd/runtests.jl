@testset "Finite-differences" begin
    include("basis_functions.jl")
    include("scalar_operators.jl")
    include("inner_products.jl")

    @testset "Function interpolation" begin
        R = FiniteDifferences(20,1.0)
        R⁻¹ = pinv(R)
        r̃ = locs(R)
        χ = R[r̃,:]

        r = axes(R,1)

        y = Inclusion(0.0..5.0)
        @test_throws DimensionMismatch R \ y.^2

        fu = r -> r^2*exp(-r)
        u = R*(R\fu.(r))
        @test norm(χ * (R⁻¹*u) - fu.(r̃)) == 0

        fv = r -> r^6*exp(-r)
        v = R*(R\fv.(r))
        @test norm(χ * (R⁻¹*v) - fv.(r̃)) == 0
    end
    
    include("derivatives.jl")
    include("densities.jl")
end
