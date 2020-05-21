@testset "Lazy derivatives" begin
    for (t₀,ϕ) in [(1.0,0.0), (4.0,π/3)]
        B = FEDVR(1.0:7, 4, t₀=t₀, ϕ=ϕ)
        B̃ = B[:,2:end-1]
        D = Derivative(axes(B,1))

        BD = B'⋆D
        BDD = B'⋆D'⋆D

        @test BD*B == B'*D*B
        @test_broken BDD*B == B'*D'*D*B

        @test B'D*B == B'*D*B
        @test B'D'D*B == B'*D'*D*B
    end
end

@testset "Materialize derivatives" begin
    @testset "$style materialization" for (style,first_derivative,second_derivative) in [
        ("Infix", (B,D) -> B'*D*B, (B,D) -> B'*D'*D*B),
        ("applied", (B,D) -> materialize(applied(*, B', D, B)),
         (B,D) -> materialize(applied(*, B', D', D, B))),
        ("apply", (B,D) -> apply(*, B', D, B),
         (B,D) -> apply(*, B', D', D, B))
    ]
        B = FEDVR(1.0:7, 4)
        D = Derivative(axes(B,1))

        ∇ = first_derivative(B,D)
        ∇² = second_derivative(B,D)

        A′ = Matrix(undef, B)
        A′′ = Matrix(undef, B)

        CompactBases.derop!(A′, B, 1)
        CompactBases.derop!(A′′, B, 2)

        @test ∇ isa BlockSkylineMatrix
        @test ∇² isa BlockSkylineMatrix

        @test ∇ == A′
        @test ∇² == A′′
    end
end

@testset "Derivatives in restricted bases" begin
    @testset "$style materialization" for (style,first_derivative,second_derivative) in [
        ("Infix", (B,D) -> B'*D*B, (B,D) -> B'*D'*D*B),
        ("applied", (B,D) -> materialize(applied(*, B', D, B)),
         (B,D) -> materialize(applied(*, B', D', D, B))),
        ("apply", (B,D) -> apply(*, B', D, B),
         (B,D) -> apply(*, B', D', D, B))
    ]
        @testset "order = $(order), sel = $(sel), N = $N" for (order,sel,N) in [
            (4,2:12,5),(4,1:13,5),([5,4,2,2],5:9,5),([5,4,3,2],5:10,5),
            (2,2:3,5),(2,2:4,5),(2,2:5,5),(2,1:5,5),
            (5,2:4,2)
        ]
            t = range(0,stop=20,length=N)

            R = FEDVR(t, order)
            R̃ = R[:,sel]

            D = Derivative(axes(R̃,1))

            expT = (any(order .> 2) ? BlockSkylineMatrix : Tridiagonal)

            ∇ = first_derivative(R, D)
            ∇̃ = first_derivative(R̃, D)

            @test ∇ isa expT
            @test ∇̃ isa expT

            @test Matrix(∇)[sel,sel] == Matrix(∇̃)

            ∇² = second_derivative(R, D)
            ∇̃² = second_derivative(R̃, D)

            @test ∇² isa expT
            @test ∇̃² isa expT

            @test Matrix(∇²)[sel,sel] == Matrix(∇̃²)
        end
    end

    @testset "Gradually decreasing size" begin
        L = 20.0
        n = 10
        t = range(-L/2, stop=L/2, length=n)
        B = FEDVR(t, [4,3,4,2,5,4,2,4,8])
        x = axes(B,1)

        D = Derivative(x)
        Dm = apply(*, B', D, B)
        @test Dm isa BlockSkylineMatrix

        @testset "From the left, a = $a" for a = 1:size(B,2)
            B̃ = B[:,a:end]
            D̃m = apply(*, B̃', D, B̃)
            @test D̃m isa BlockSkylineMatrix
            @test Dm[a:end,a:end] == D̃m
        end

        @testset "From the right, b = $b" for b = size(B,2):-1:1
            B̃ = B[:,1:b]
            D̃m = apply(*, B̃', D, B̃)
            @test D̃m isa BlockSkylineMatrix
            @test Dm[1:b,1:b] == D̃m
        end
    end
end

function compute_fedvr_derivative_errors(a, b, Ns, order::Integer, s::Integer, e::Integer,
                                         f::Function, g::Function, h::Function,
                                         verbosity=0)
    compute_derivative_errors(Ns, f, g, h, verbosity) do N
        t = range(a, stop=b, length=N)
        R = FEDVR(t, order)
        if (s,e) != (0,0)
            n = size(R,2)
            R = R[:,(1+s):(n-e)]
        end
        R,1
    end
end

@testset "Derivative accuracy" begin
    (f,g,h,a,b),s,e = derivative_test_functions(1.0), 1, 1

    orders = [(2,1.5,3),
              (3,1.2,2.5),
              (4,1.3,2.5),
              (5,1.3,2.5),
              (6,1.3,2.5),
              (7,1,2.5),
              (8,1,2.5),
              (9,1,2.5),
              (10,1,2.5)]
    slopes = zeros(length(orders),3)

    @info "FE-DVR derivative accuracy"
    for (i,(order,oa,ob)) in enumerate(orders)
        @show order
        Ns = ceil.(Int, 10 .^ range(oa,stop=ob,length=10))
        @show extrema(Ns)
        hs,ϵg,ϵh,ϵh′,pg,ph,ph′ = compute_fedvr_derivative_errors(a, b, Ns, order, s, e, f, g, h, 1)
        slopes[i,:] = [pg ph ph′]
    end

    oo = first.(orders)

    println("Derivative convergence rates:")
    pretty_table([oo slopes], ["Order", "pg", "ph", "ph′"])
    println()

    plt = lineplot(oo, slopes[:,1], name="pg", xlabel="FE-DVR order", ylabel="Error slope")
    lineplot!(plt, oo, slopes[:,2], name="ph")
    lineplot!(plt, oo, slopes[:,3], name="ph′")
    display(plt)
    println()
    println()
    
    # The convergence rate should be order - 1 (polynomial order =
    # polynomial degree + 1), but since the error slope fitting is a
    # bit error prone, we require that it is greater than order - 2.5.
    @test_broken all(slopes[:,1] .> oo .- 2.5)
    @test all(slopes[:,1] .> oo .- 3.5)
    # Since the approximation to h is calculated by computing ∇∇f, we
    # lose one order extra, compared to ∇²f.
    @test_broken all(slopes[:,2] .> oo .- 3.5)
    @test all(slopes[:,2] .> oo .- 4.5)
    @test_broken all(slopes[:,3] .> oo .- 2.5)
    @test all(slopes[:,3] .> oo .- 4.0)
end
