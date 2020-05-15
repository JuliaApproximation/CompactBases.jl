import CompactBases: nel, complex_rotate, real_locs

@testset "FEDVR" begin
    include("basics.jl")
    include("complex_scaling.jl")
    include("block_structure.jl")
    include("scalar_operators.jl")
    include("inner_products.jl")
    include("function_interpolation.jl")
    include("derivatives.jl")
    include("densities.jl")
end
