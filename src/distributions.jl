abstract type NodeDistribution end
struct Uniform <: NodeDistribution end
struct NonUniform <: NodeDistribution end

distribution(::AbstractRange) = Uniform()
distribution(::AbstractVector) = NonUniform()
