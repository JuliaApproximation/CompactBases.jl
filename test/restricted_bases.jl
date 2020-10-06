@testset "Restrictions" begin
    @testset "$(label)" for (R,label) in [
        (FEDVR(range(0, stop=1.0, length=30), 6), "FE-DVR")
        (StaggeredFiniteDifferences(100, 0.1), "Finite-differences")
    ]
        for (A,B) in [(R,R), (R,R[:,6:20]), (R[:,6:20],R[:,3:80])]
            rA = CompactBases.restriction(A)
            rB = CompactBases.restriction(B)

            rAB = CompactBases.combined_restriction(A,B)
            rBA = CompactBases.combined_restriction(B,A)
            @test rAB == rA'rB
            @test rAB.data isa Ones
            @test rBA == rB'rA
            @test rBA.data isa Ones

            S = A'B
            @test S.data isa Fill
        end
    end
end
