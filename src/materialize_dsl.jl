# Pending https://github.com/JuliaArrays/LazyArrays.jl/pull/91

# For unparametrized destination types
generate_copyto!_signature(dest, dest_type::Symbol, Msig) =
    :(Base.copyto!($(dest)::$(dest_type), applied_obj::$(Msig)))

# For parametrized destination types
function generate_copyto!_signature(dest, dest_type::Expr, Msig)
    dest_type.head == :curly ||
        throw(ArgumentError("Invalid destination specification $(dest)::$(dest_type)"))
    :(Base.copyto!($(dest)::$(dest_type), applied_obj::$(Msig)) where {$(dest_type.args[2:end]...)})
end

function generate_copyto!(body, factor_names, Msig)
    body.head == :(->) ||
        throw(ArgumentError("Invalid copyto! specification"))
    body.args[1].head == :(::) ||
        throw(ArgumentError("Invalid destination specification $(body.args[1])"))
    (dest,dest_type) = body.args[1].args
    copyto!_signature = generate_copyto!_signature(dest, dest_type, Msig)
    f_body = quote
        axes($dest) == axes(applied_obj) || throw(DimensionMismatch("axes must be same"))
        $(factor_names) = applied_obj.args
        $(body.args[2].args...)
        $(dest)
    end
    Expr(:function, copyto!_signature, f_body)
end

"""
    @materialize function op(args...)

This macro simplifies the setup of a few functions necessary for the
materialization of [`Applied`](@ref) objects:

- `ApplyStyle`, used to ensure dispatch of the applied object to the
  routines below

- `copyto!(dest::DestType, applied_obj::Applied{...,op})` performs the
  actual materialization of `applied_obj` into the destination object
  which has been generated by

- `similar` which usually returns a suitable matrix

- `materialize` which makes use of the above functions

# Example

```julia
@materialize function *(Ac::MyAdjointBasis,
                        O::MyOperator,
                        B::MyBasis)
    MyApplyStyle # An instance of this type will be returned by ApplyStyle
    T -> begin # generates similar
        A = parent(Ac)
        parent(A) == parent(B) ||
            throw(ArgumentError("Incompatible bases"))

        # There may be different matrices best representing different
        # situations:
        if ...
            Diagonal(Vector{T}(undef, size(B,1)))
        else
            Tridiagonal(Vector{T}(undef, size(B,1)-1),
                        Vector{T}(undef, size(B,1)),
                        Vector{T}(undef, size(B,1)-1))
        end
    end
    dest::Diagonal{T} -> begin # generate copyto!(dest::Diagonal{T}, ...) where T
        dest.diag .= 1
    end
    dest::Tridiagonal{T} -> begin # generate copyto!(dest::Tridiagonal{T}, ...) where T
        dest.dl .= -2
        dest.ev .= 1
        dest.du .= 3
    end
end
```
"""
macro materialize(expr)
    expr.head == :function || expr.head == :(=) || error("Must start with a function")
    @assert expr.args[1].head == :call
    op = expr.args[1].args[1]

    bodies = filter(e -> !(e isa LineNumberNode), expr.args[2].args)
    length(bodies) < 3 &&
        throw(ArgumentError("At least three blocks required (ApplyStyle, similar, and at least one copyto!)"))

    factor_types = :(<:Tuple{})
    factor_names = :(())
    apply_style = first(bodies)
    apply_style_fun = :(LazyArrays.ApplyStyle(::typeof($op)) = $(apply_style)())

    # Generate Applied signature
    for arg in expr.args[1].args[2:end]
        arg isa Expr && arg.head == :(::) ||
            throw(ArgumentError("Invalid argument specification $(arg)"))
        arg_name, arg_typ = arg.args
        push!(factor_types.args[1].args, :(<:$(arg_typ)))
        push!(factor_names.args, arg_name)
        push!(apply_style_fun.args[1].args, :(::Type{<:$(arg_typ)}))
    end
    Msig = :(LazyArrays.Applied{$(apply_style), typeof($op), $(factor_types)})

    sim_body = bodies[2]
    sim_body.head == :(->) ||
        throw(ArgumentError("Invalid similar specification"))
    T = first(sim_body.args)

    copytos! = map(body -> generate_copyto!(body, factor_names, Msig), bodies[3:end])

    f = quote
        $(apply_style_fun)

        function Base.similar(applied_obj::$Msig, ::Type{$T}=eltype(applied_obj)) where $T
            $(factor_names) = applied_obj.args
            $(sim_body.args[2])
        end

        $(copytos!...)

        LazyArrays.materialize(applied_obj::$Msig) =
            copyto!(similar(applied_obj, eltype(applied_obj)), applied_obj)
    end
    esc(f)
end
