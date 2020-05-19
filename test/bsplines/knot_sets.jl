using IntervalSets
import CompactBases: RightContinuous, find_interval, within_interval, within_support

"""
    within_interval_linear(x, interval)

Return the range of indices of elements of the sorted vector/range `x`
that fall within `interval`. This loops through all elements of `x`
(linear complexity), so it should not be used for large number of
elements.
"""
function within_interval_linear(x, interval)
    sel = findall(e -> e ∈ interval, x)
    isempty(sel) ? (1:0) : (minimum(sel):maximum(sel))
end

function test_within_interval(x, interval, expected=nothing; reversed=true)
    if reversed
        x = reverse(x)
    end
    N = length(x)

    linear_expected = within_interval_linear(x, interval)
    if isnothing(expected)
        expected = linear_expected
    else
        if reversed
            expected = (N-expected[end]+1):(N-expected[1]+1)
        end
        # This is mostly a sanity check to test if the expectations
        # are actually correct.
        expected ≠ linear_expected &&
            throw(ArgumentError("Invalid expectations: $(expected), actual: $(linear_expected)"))
    end

    result = :(within_interval($x, $interval))
    expr = :($result == $expected || isempty($result) && isempty($expected))
    if !(@eval $expr)
        println("Looking for elements of $x ∈ $interval, got $(@eval $result), expected $expected")
        length(x) < 30 && println("    x = ", collect(enumerate(x)), "\n")
    end
    @eval @test $expr
end

@testset "Knot sets" begin
    @testset "Simple tests" begin
        @test ArbitraryKnotSet(3, [0.0, 1, 1, 3, 4, 6], 1, 3) == [0.0, 1, 1, 3, 4, 6, 6, 6]
        @test LinearKnotSet(3, 0, 1, 2) == [0,0,0,0.5,1,1,1]
        @test LinearKnotSet(2, 0, 1, 2, 1, 1) == [0,0.5,1]
        @test ExpKnotSet(2, -4, 2, 7) == [0,0,0.0001,0.001,0.01,0.1,1,10,100,100]
    end

    @testset "k = $k" for k ∈ 1:6
        @testset "Full multiplicity" begin
            t = LinearKnotSet(k, 0, 1, 1)
            @test length(t) == 2k
            @test collect(t) == vcat(fill(0, k), fill(1, k))
            @test length(nonempty_intervals(t)) == 1
        end
        @testset "Simple endpoints" begin
            t = LinearKnotSet(k, 0, 1, k, 1, 1)
            @test length(t) == k+1
            @test collect(t) == range(0, stop=1, length=k+1)
            @test length(nonempty_intervals(t)) == k
        end
    end
    @test length(nonempty_intervals(ArbitraryKnotSet(3, [0.0, 1, 1, 3, 4, 6], 1, 3))) == 4

    @testset "Invalid order/multiplicities" begin
        @testset "Invalid order" begin
            @test_throws ArgumentError LinearKnotSet(0, 0, 1, 4)
            @test_throws ArgumentError LinearKnotSet(-1, 0, 1, 4)
        end
        @testset "Invalid multipicities" begin
            @test_throws ArgumentError LinearKnotSet(1, 0, 1, 4, 0, 1)
            @test_throws ArgumentError LinearKnotSet(1, 0, 1, 4, 1, 0)
            @test_throws ArgumentError LinearKnotSet(1, 0, 1, 4, 1, 2)
        end
        @testset "Not enough interior points" begin
            @test_throws ArgumentError LinearKnotSet(2, 0, 1, 1, 1, 1)
        end
    end

    @testset "Properties" begin
        @testset "#intervals = $N" for N = 1:4
            @testset "k = $k" for k = 1:N
                @testset "Simple multiplicities" begin
                    t = LinearKnotSet(k, 0.0, 1.0, N, 1, 1)
                    @test numintervals(t) == N
                    @test numfunctions(t) == N - k + 1
                end
            end
        end
    end

    @testset "Compact support" begin
        a,b = 0,1
        x = range(a, stop=b, length=21)
        @testset "Interval coverage" begin
            @testset "Two intervals" begin
                @testset "Reversed: $reversed" for reversed in [false, true]
                    test_within_interval(x, 0..0.5, 1:11, reversed=reversed)
                end
                @testset "L=$L" for L=[:closed,:open]
                    @testset "R=$R" for R=[:closed,:open]
                        @testset "Reversed: $reversed" for reversed in [false, true]
                            test_within_interval(x, Interval{L,R}(0,0.5), reversed=reversed)
                            test_within_interval(x, Interval{L,R}(0.25,0.5), reversed=reversed)
                        end
                    end
                end
            end
            @testset "Three intervals" begin
                @testset "Reversed: $reversed" for reversed in [false, true]
                    test_within_interval(x, RightContinuous(0,1/3), 1:7, reversed=reversed)
                    test_within_interval(x, RightContinuous(1/3,2/3), 8:14, reversed=reversed)
                    test_within_interval(x, 2/3..1, 15:21, reversed=reversed)
                end
            end
            @testset "Open interval" begin
                @testset "Reversed: $reversed" for reversed in [false, true]
                    test_within_interval(x, OpenInterval(0.2,0.4), 6:8, reversed=reversed)
                end
            end
            @testset "Random intervals" begin
                @testset "L=$L" for L=[:closed,:open]
                    @testset "R=$R" for R=[:closed,:open]
                        @testset "Reversed: $reversed" for reversed in [false, true]
                            for i = 1:20
                                interval = Interval{L,R}(minmax(rand(),rand())...)
                                test_within_interval(x, interval, reversed=reversed)
                            end
                        end
                    end
                end
            end
        end

        @testset "Partially covered intervals" begin
            k = 3
            t = LinearKnotSet(k, 0, 1, 2)
            @testset "$name, x = $x" for (name,x) in [
                ("Outside left",range(-1,stop=-0.5,length=10)),
                ("Touching left",range(-1,stop=0,length=10)),
                ("Touching left-ϵ",range(-1,stop=0-eps(),length=10)),
                ("Touching left+ϵ",range(-1,stop=0+eps(),length=10)),

                ("Outside right",range(1.5,stop=2,length=10)),
                ("Touching right",range(1,stop=2,length=10)),
                ("Touching right-ϵ",range(1-eps(),stop=2,length=10)),
                ("Touching right+ϵ",range(1+eps(),stop=2,length=10)),

                ("Other right",range(0.5,stop=1,length=10)),
                ("Other right-ϵ",range(0.5-eps(0.5),stop=1,length=10)),
                ("Other right+ϵ",range(0.5+eps(0.5),stop=1,length=10)),

                ("Complete", range(0,stop=1,length=10)),
                ("Complete-ϵ", range(eps(),stop=1-eps(),length=10)),
                ("Complete+ϵ", range(-eps(),stop=1+eps(),length=10)),

                ("Left partial", range(-0.5,stop=0.6,length=10)),
                ("Left", range(-0.5,stop=1.0,length=10)),
                ("Right partial", range(0.5,stop=1.6,length=10)),
                ("Right", range(0,stop=1.6,length=10))]
                @testset "L=$L" for L=[:closed,:open]
                    @testset "R=$R" for R=[:closed,:open]
                        @testset "Reversed: $reversed" for reversed in [false, true]
                            for i in nonempty_intervals(t)
                                interval = Interval{L,R}(t[i], t[i+1])
                                test_within_interval(x, interval, reversed=reversed)
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "Finding intervals" begin
        t = ArbitraryKnotSet(3, [0.0, 1, 1, 3, 4, 6], 1, 3)
        R = BSpline(t)
        # Test without specifying initial interval
        @test find_interval(t, 0.5) == 1
        for (i,xs) in enumerate([[0, 0.5, 1.0-eps()],
                                 [],
                                 [1.0, 2.0, 3.0-eps(3.0)],
                                 [3, 3.5, 4.0-eps(4.0)],
                                 [4, 5, 6.0-eps(6.0),6]])
            for j in [1,i] # Test starting from the beginning, and from the interval itself
                for x in xs
                    result = :(find_interval($t, $x, $j))
                    res = @eval $result
                    expr = :($result == $i)
                    (@eval $expr) ||
                        println("Expected to find $x ∈ interval #$i, got $(isnothing(res) ? "nothing" : res)")
                    @eval @test $expr
                    # @test find_interval(t, x) ==
                    find_interval(t, x)
                    within_support(x:x, t, i)[1][2]
                end
            end
        end
        @test isnothing(find_interval(t, 0.5, 2))
    end

    @testset "Support of Heavyside splines" begin
        a,b = 0,1
        x = range(a, stop=b, length=21)
        t = LinearKnotSet(1, a, b, 4)
        supports = [within_support(x, t, j)
                    for j = 1:numfunctions(t)]
        @test length(supports) == 4
        # Each basis function should cover one interval only (since order = 1).
        @test all(length.(supports) .== 1)
        # Test that all elements of x are covered by the basis
        # functions, and that none of the basis functions overlap.
        @test first.(first.(supports)) == [1:5, 6:10, 11:15, 16:21]
    end

    @testset "Support of linear splines" begin
        a,b = 0,1
        x = range(a, stop=b, length=21)
        @testset "Simple multiplicity" begin
            t = LinearKnotSet(2, a, b, 2, 1, 1)
            supports = [within_support(x, t, j)
                        for j = 1:numfunctions(t)]
            @test length(supports) == 1
            @test supports[1] == [(1:10,1), (11:21,2)]
        end
        @testset "Full multiplicity" begin
            @testset "#intervals = 2" begin
                t = LinearKnotSet(2, a, b, 2)
                supports = [within_support(x, t, j)
                            for j = 1:numfunctions(t)]
                @test length(supports) == 3
                @test supports[1] == [(1:10,2)]
                @test supports[2] == [(1:10,2), (11:21,3)]
                @test supports[3] == [(11:21,3)]
            end
            @testset "#intervals = 3" begin
                t = LinearKnotSet(2, a, b, 3)
                supports = [within_support(x, t, j)
                            for j = 1:numfunctions(t)]
                @test length(supports) == 4
                @test supports[1] == [(1:7,2)]
                @test supports[2] == [(1:7,2), (8:14,3)]
                @test supports[3] == [(8:14,3), (15:21,4)]
                @test supports[4] == [(15:21,4)]
            end
        end
    end
end

@testset "Quadrature" begin
    t = LinearKnotSet(1, 0, 1, 2)
    x,w = CompactBases.lgwt(t, 2)
    @test all(w .== 1/4)
    @test x ≈ [-1,1,-1,1]/(4*√3) + [1,1,3,3]/4 rtol=1e-13

    @test_logs (:warn, "N = 1 quadrature point not enough to calculate overlaps between polynomials of order k = 3") CompactBases.lgwt(LinearKnotSet(3, 0, 1, 3), 1)

    @test_logs (:warn, "N = 2 quadrature points not enough to calculate overlaps between polynomials of order k = 3") CompactBases.lgwt(LinearKnotSet(3, 0, 1, 3), 2)
end
