using LinearAlgebra

using Test
using PrettyTables

function test_derivatives(R::AbstractQuasiMatrix,
                          f::Function, g::Function, h::Function)
    r = axes(R,1)

    D = Derivative(r)

    ∇ = R' * D * R
    ∇² = R' * D' * D * R

    fv = R \ f.(r)

    gv = ∇ * fv
    δgv = R \ g.(r) - gv

    hv = ∇ * gv
    δhv = R \ h.(r) - hv

    hv′ = ∇² * fv
    δhv′ = R \ h.(r) - hv′

    R,fv,gv,hv,hv′,δgv,δhv,δhv′
end

function error_slope(loghs,ϵ)
    # To avoid the effect of round-off errors on the order
    # estimation.
    i = argmin(abs.(ϵ))

    ([loghs[1:i] ones(i)] \ log10.(abs.(ϵ[1:i])))[1]
end

function derivative_test_functions(d)
    a,b = √d*[-1,1]

    # The functions vanish at the boundaries; Dirichlet0 boundary
    # conditions. They /cannot/ be evaluated at the boundaries, since
    # that results in NaNs.
    f = x -> exp(-1/(d-x^2))
    g = x -> -2*exp(-1/(d-x^2))*x/((d-x^2)^2)
    h = x -> -2*exp(-1/(d-x^2))*(d^2 + 2*(d-1)*x^2-3x^4)/((d-x^2)^4)

    f,g,h,a,b
end

function compute_derivative_errors(fun::Function,
                                   a, b, Ns,
                                   f::Function, g::Function, h::Function,
                                   verbosity=0)
    errors = map(Ns) do N
        R,μ = fun(N)
        R,fv,gv,hv,hv′,δgv,δhv,δhv′ = test_derivatives(R, f, g, h)
        [norm(δgv)*μ norm(δhv)*μ norm(δhv′)*μ]
    end |> e -> vcat(e...)

    ϵg = errors[:,1]
    ϵh = errors[:,2]
    ϵh′ = errors[:,3]

    hs = 1.0 ./ Ns

    loghs = log10.(hs)
    pg = error_slope(loghs, ϵg)
    ph = error_slope(loghs, ϵh)
    ph′ = error_slope(loghs, ϵh′)

    verbosity > 0 &&
        pretty_table([hs errors], ["h", "δg [$(pg)]", "δh [$(ph)]", "δh′ [$(ph′)]"])

    hs,ϵg,ϵh,ϵh′,pg,ph,ph′    
end
