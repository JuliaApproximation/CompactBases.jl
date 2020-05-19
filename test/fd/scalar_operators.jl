@testset "Mass matrices and inverses" begin
    R = FiniteDifferences(20,0.2)
    @testset "Mass matrix" begin
        @test R'R == step(R)*I
    end
    @testset "Inverses" begin
        R⁻¹ = pinv(R)
        @test R⁻¹*R == I
        @test R*R⁻¹ == I

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
    end
end

@testset "Scalar operators" begin
    B = FiniteDifferences(5,1.0)
    x = axes(B,1)
    x² = B'QuasiDiagonal(x.^2)*B
    @test x² isa Diagonal

    v = ones(5)
    @test x²*v == locs(B).^2

    y = Inclusion(0.0..5.0)
    @test_throws DimensionMismatch B'QuasiDiagonal(y.^2)*B

    @testset "Restricted bases" begin
        N = 10
        ρ = 1.0

        R = FiniteDifferences(N, ρ)
        r = axes(R, 1)

        sel = 3:6
        sel2 = 8:10

        R̃ = R[:, sel]
        R′ = R[:, sel2]

        x = QuasiDiagonal(r)

        apply_obj = applied(*, R', x, R)
        @test LazyArrays.ApplyStyle(*, typeof(R'), typeof(x), typeof(R)) ==
            CompactBases.FiniteDifferencesStyle()

        A = apply(*, R', x, R)
        @test A isa Diagonal

        a = apply(*, R̃', x, R)
        @test a isa BandedMatrix
        @test a == Matrix(A)[sel,:]
        @test bandwidths(a) == (-2,2)

        a = apply(*, R', x, R̃)
        @test a isa BandedMatrix
        @test a == Matrix(A)[:,sel]
        @test bandwidths(a) == (2,-2)

        a = apply(*, R̃', x, R̃)
        @test a isa Diagonal
        @test a == Matrix(A)[sel,sel]

        # Non-overlapping restrictions
        a = apply(*, R̃', x, R′)
        @test a isa BandedMatrix
        @test size(a) == (length(sel), length(sel2))
        @test iszero(a)
        @test bandwidths(a) == (5,-5)
    end
end
