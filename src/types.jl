const FuncArray{T,N,B<:BasisOrRestricted} = MulQuasiArray{T,N,<:Tuple{B,<:AbstractArray{T,N}}}
const FuncVector{T,B<:BasisOrRestricted} = FuncArray{T,1,B}
const FuncMatrix{T,B<:BasisOrRestricted} = FuncArray{T,2,B}
const FuncVecOrMat{T,B<:BasisOrRestricted} = Union{FuncVector{T,B},FuncMatrix{T,B}}

const AdjointFuncArray{T,N,B<:BasisOrRestricted} = MulQuasiArray{T,<:Any,<:Tuple{<:Adjoint{T,<:AbstractArray{T,N}},
                                                                    <:QuasiAdjoint{T,<:B}}}
const AdjointFuncVector{T,B<:BasisOrRestricted} = AdjointFuncArray{T,1,B}
const AdjointFuncMatrix{T,B<:BasisOrRestricted} = AdjointFuncArray{T,2,B}
const AdjointFuncVecOrMat{T,B<:BasisOrRestricted} = Union{AdjointFuncVector{T,B},AdjointFuncMatrix{T,B}}
