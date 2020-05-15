function test_block_structure(t, o, l, u)
    B = FEDVR(t, o)
    M = Matrix(undef, B)

    @test M isa BlockSkylineMatrix
    bs = M.block_sizes
    @test length(bs.l) == length(l)+1
    @test length(bs.u) == length(u)+1

    # TODO Test block sizes

    @test all(bs.l[1:end-2] .== l[1:end-1])
    length(l) > 0 && @test bs.l[end-1] ≥ l[end]
    @test all(bs.u[3:end] .== u[2:end])
    length(u) > 0 && @test bs.u[2] ≥ u[1]
end

@testset "Block structure" begin
    @testset "Simple" begin
    test_block_structure(1.0:2, [4], [], [])
    test_block_structure(1.0:7, [2,2,3,4,2,4], [1,1,2,1,2,1,1,1], [1,1,1,2,1,2,1,1])
    test_block_structure(1.0:8, [2,2,3,4,2,4,2], [1,1,2,1,2,1,1,2,1,1], [1,1,1,2,1,2,1,1,2,1])
    test_block_structure(range(0,1,length=7), 1 .+ (1:6), [1,2,1,2,1,2,1,2,1,1], [1,1,2,1,2,1,2,1,2,1])
    test_block_structure(range(0,1,length=7), reverse(1 .+ (1:6)), [1,2,1,2,1,2,1,2,1,1], [1,1,2,1,2,1,2,1,2,1])
    test_block_structure(range(0,1,length=7), 2 .+ (1:6), [1,2,1,2,1,2,1,2,1,1], [1,1,2,1,2,1,2,1,2,1])
        @test Matrix(undef, FEDVR(1.0:7, 2)) isa Tridiagonal
    end
    @testset "Restrictions" begin
        L = 20.0
        n = 10
        t = range(-L/2, stop=L/2, length=n)
        B = FEDVR(t, [4,3,4,2,5,4,2,4,8])
        @test CompactBases.elements(B) == 1:9
        @test CompactBases.block_structure(B) == [3, 1, 1, 1, 2, 1, 1, 3, 1, 2, 1, 1, 2, 1, 7]
        # @test CompactBases.block_bandwidths(B) == ([1,2,1,2,1,1,2,1,2,1,1,2,1,1],
        #                                          [1,2,1,2,1,1,2,1,2,1,1,2,1,1])

        N = size(B,2)
        a = round(Int, 0.4N)
        b = round(Int, 0.6N)+1

        @test CompactBases.elements(B[:,1:b]) == 1:8
        @test CompactBases.block_structure(B[:,1:b]) == [3, 1, 1, 1, 2, 1, 1, 3, 1, 2, 1, 1]
        @test CompactBases.elements(B[:,2:b]) == 1:8
        @test CompactBases.block_structure(B[:,2:b]) == [2, 1, 1, 1, 2, 1, 1, 3, 1, 2, 1, 1]
        @test CompactBases.elements(B[:,3:b]) == 1:8
        @test CompactBases.block_structure(B[:,3:b]) == [1, 1, 1, 1, 2, 1, 1, 3, 1, 2, 1, 1]
        @test CompactBases.elements(B[:,4:b]) == 1:8
        @test CompactBases.block_structure(B[:,4:b]) == [2, 1, 2, 1, 1, 3, 1, 2, 1, 1]

        @test CompactBases.elements(B[:,1:b-1]) == 1:7
        @test CompactBases.block_structure(B[:,1:b-1]) == [3, 1, 1, 1, 2, 1, 1, 3, 1, 3]
        @test CompactBases.elements(B[:,1:b-2]) == 1:6
        @test CompactBases.block_structure(B[:,1:b-2]) == [3, 1, 1, 1, 2, 1, 1, 3, 1, 2]
        @test CompactBases.elements(B[:,1:b-3]) == 1:6
        @test CompactBases.block_structure(B[:,1:b-3]) == [3, 1, 1, 1, 2, 1, 1, 3, 1, 1]
        @test CompactBases.elements(B[:,1:b-4]) == 1:6
        @test CompactBases.block_structure(B[:,1:b-4]) == [3, 1, 1, 1, 2, 1, 1, 4]

        @test CompactBases.elements(B[:,a:end]) == 5:9
        @test CompactBases.block_structure(B[:,a:end]) == [3, 1, 2, 1, 1, 2, 1, 7]
        @test CompactBases.elements(B[:,a:end-1]) == 5:9
        @test CompactBases.block_structure(B[:,a:end-1]) == [3, 1, 2, 1, 1, 2, 1, 6]
        @test CompactBases.elements(B[:,a:end-7]) == 5:9
        @test CompactBases.block_structure(B[:,a:end-7]) == [3, 1, 2, 1, 1, 3]
        @test CompactBases.elements(B[:,a:end-8]) == 5:8
        @test CompactBases.block_structure(B[:,a:end-8]) == [3, 1, 2, 1, 1, 2]

        # @test CompactBases.block_bandwidths(B[:,a:end-1]) == ([1,2,1,1,2,1,1],[1,1,2,1,1,2,1])
        # @test CompactBases.block_bandwidths(B[:,a:end-7]) == ([1,2,1,1,1],[1,1,2,1,1])

        @testset "Single element" begin
            B = FEDVR(range(-L/2, stop=L/2, length=2), 8)
            @test CompactBases.block_structure(B) == [8]
            @test CompactBases.block_structure(B[:,1:end-1]) == [7]
            @test CompactBases.block_structure(B[:,2:end]) == [7]
            @test CompactBases.block_structure(B[:,2:end-1]) == [6]
            @test CompactBases.block_bandwidths(B) == ([0],[0])
            @test CompactBases.block_bandwidths(B[:,2:end-1]) == ([0],[0])
        end
    end
end

function test_blocks(f::Function, t, o)
    B = FEDVR(t, o)
    A = Matrix(undef, B)

    blocks = map(f, 1:CompactBases.nel(B))
    coords = vcat(1,1 .+ cumsum(B.order[1:end-1] .- 1))
    CompactBases.set_elements!(f, A, B)

    for (i,(b,c)) in enumerate(zip(blocks,coords))
        sel = c .+ (0:size(b,1)-1)
        b′ = copy(A[sel,sel])
        i > 1 && (b′[1,1] -= blocks[i-1][end,end])
        i < length(blocks) && (b′[end,end] -= blocks[i+1][1,1])
        @test all(b′ .== b)
    end
end

@testset "Set blocks" begin
    test_blocks(1.0:3,[2,3]) do i
        i*ones(i+1,i+1)
    end
    begin
        o = [2,2,3,4,2,4]
        test_blocks(1.0:7,o) do i
            i*ones(o[i],o[i])
        end
    end
    begin
        o = 4 .+ (1:6)
        test_blocks(range(0,1,length=7), o) do i
            i*ones(o[i],o[i])
        end
    end
end

@testset "Set blocks 2" begin
    L = 20.0
    n = 10
    t = range(-L/2, stop=L/2, length=n)
    B = FEDVR(t, [4,3,4,2,5,4,2,4,8])

    N = size(B,2)

    Bl = B[:,1:round(Int, 0.6N)+1]
    M = Matrix(undef, Bl, Int)
    @test M isa BlockSkylineMatrix
    @test eltype(M) == Int
    @test size(M) == (18,18)

    CompactBases.set_elements!(M, Bl) do i
        o = CompactBases.order(Bl, i)
        i*ones(Int,o,o)
    end

    for (v,ijs) in [1 => [(1,1), (1,2), (2,1)],
                    3 => [(2,2)],
                    2 => [(2,3), (2,4), (3,2), (4,2), (3,3), (3,4), (4,3)],
                    5 => [(4,4)],
                    3 => [(5,4), (6,4), (4,5), (4,6), (5,5), (6,5), (5,6)],
                    7 => [(6,6)],
                    4 => [(7,6), (6,7)],
                    9 => [(7,7)],
                    5 => [(8,7), (9,7), (7,8), (7,9), (8,8), (9,8), (8,9)],
                    11 => [(9,9)],
                    6 => [(10,9), (11,9), (9,10), (9,11), (10,10), (11,10), (10,11)],
                    13 => [(11,11)],
                    7 => [(12,11), (11,12)],
                    15 => [(12,12)]]
        for (i,j) in ijs
            @test all(M[Block(i,j)] .== v)
        end
    end

    Bl2 = B[:,4:round(Int, 0.6N)]
    M2 = Matrix(undef, Bl2)
    @time CompactBases.set_elements!(M2, Bl2) do i
        o = CompactBases.order(Bl2, i)
        i*ones(Int,o,o)
    end
    @test M2 == M[4:end-1,4:end-1]
end
