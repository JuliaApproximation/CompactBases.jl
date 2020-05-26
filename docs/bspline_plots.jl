function Bspline_text(x, ϕ, j, k, color)
    xⱼ = mean_position(x, ϕ)
    txt = text(xⱼ,0.7maximum(ϕ),
               latexstring("\\mathrm{B}_{$j,$k}"),
               horizontalalignment="center",
               color=color)
    gray = 1.0 - mean_color(color)
    outline = RGB(lerp(gray, (gray<0.5 ? 0 : 1), 0.6)*[1,1,1]...)
    txt.set_path_effects([PathEffects.Stroke(linewidth=1,
                                             foreground="#$(hex(outline))"),
                          PathEffects.Normal()])
end

function bsplines_cardinal_splines()
    a,b = 1.0,6.0
    x = range(a, stop=b, length=1001)
    cfigure("cardinal splines",figsize=(7,9)) do
        t = 0
        for k = 1:5
            t = ArbitraryKnotSet(k, a:b, 1, 1)
            R = BSpline(t)
            χ = R[x,:]
            csubplot(5,1,k,nox=k<5) do
                for j = 1:size(χ,2)
                    ϕ = view(χ, :, j)
                    l=plot(x, ϕ)[1]
                    Bspline_text(x, ϕ, j, k, l.get_color())
                end
                ylabel(latexstring("k = $k"))
                margins(0.1,0.1)
            end
        end
        xlabel(L"x")
    end
    savedocfig("bsplines/cardinal-splines")
end

function bsplines_discontinuous_splines()
    t = ArbitraryKnotSet(3, [0.0, 1, 1, 3, 4, 6], 1, 3)
    r = range(first(t), stop=last(t), length=301)
    R = BSpline(t)
    χ = R[r, :]

    cfigure("basis functions",figsize=(7,6)) do
        csubplot(211,nox=true) do
            plot(r, χ, length(r) < 30 ? ".-" : "-")
            margins(0.1, 0.1)
        end
        csubplot(212) do
            rplot(t)
            margins(0.1, 0.1)
            xlabel(L"r")
            ylabel("Multiplicity")
        end
    end
    savedocfig("bsplines/discontinuous-splines")
end

function bsplines_full_multiplicity_splines()
    a,b = 1.0,6.0
    x = range(a, stop=b, length=1001)
    cfigure("splines",figsize=(7,9)) do
        t = 0
        for k = 1:5
            t = ArbitraryKnotSet(k, a:b)
            R = BSpline(t)
            χ = R[x,:]
            csubplot(5,1,k,nox=k<5) do
                for j = 1:size(χ,2)
                    ϕ = view(χ, :, j)
                    l=plot(x, ϕ)[1]
                    Bspline_text(x, ϕ, j, k, l.get_color())
                end
                ylabel(latexstring("k = $k"))
                margins(0.1,0.1)
            end
        end
        xlabel(L"x")
    end
    savedocfig("bsplines/full-multiplicity-splines")
end

function bsplines_spline1d()
    k = 4
    t = LinearKnotSet(k, 0, 1, 5)
    B = BSpline(t)
    x = range(0, stop=1, length=301)
    χ = B[x, :]

    i = 1:size(B,2)
    c = sin.(i)
    s = χ*c

    cfigure("spline 1d", figsize=(7,9)) do
        csubplot(311, nox=true) do
            plot(x, s)
            plot([mean_position(x, view(χ, :, j)) for j in i], c, "s-")
            ylabel(L"s(x)")
            margins(0.1,0.1)
        end
        csubplot(312, nox=true) do
            for j = 1:size(χ,2)
                ϕ = view(χ, :, j)
                l=plot(x, ϕ)[1]
                Bspline_text(x, ϕ, j, k, l.get_color())
            end
            margins(0.1,0.1)
        end
        csubplot(313) do
            rplot(t)
            margins(0.1,0.1)
            xlabel(L"x")
        end
    end

    savedocfig("bsplines/spline-1d")
end

function bsplines_spline2d()
    k = 4
    t = LinearKnotSet(k, 0, 1, 5)
    B = BSpline(t)
    x = range(0, stop=1, length=301)
    χ = B[x, :]

    i = 1:size(B,2)
    c = [sin.(i) tan.(i)]
    s = χ*c

    cfigure("spline 2d", figsize=(6,6)) do
        plot(s[:,1], s[:,2])
        plot(c[:,1], c[:,2], "s-")
        xlabel(L"x")
        ylabel(L"y")
        margins(0.1,0.1)
    end

    savedocfig("bsplines/spline-2d")
end

function bsplines_quadrature_points()
    t = ArbitraryKnotSet(3, [0.0, 1, 1, 3, 4, 6], 1, 3)
    x = range(first(t), stop=last(t), length=301)
    B = BSpline(t)
    χ = B[x, :]

    xl = nothing
    cfigure("quadrature points", figsize=(7,9)) do
        csubplot(411, nox=true) do
            plot(x, χ, length(x) < 30 ? ".-" : "-")
            margins(0.1, 0.1)
        end
        csubplot(412, nox=true) do
            rplot(t)
            margins(0.1, 0.1)
            xl = xlim()
            ylabel("Multiplicity")
        end
        csubplot(413, nox=true) do
            plot(B.x, 1:length(B.x), ".-")
            margins(0.1,0.1)
            xlim(xl)
            ylabel(L"i")
        end
        csubplot(414) do
            plot(B.x, B.w, ".-")
            margins(0.1,0.1)
            xlim(xl)
            xlabel(L"x")
            ylabel(L"w_i")
        end
    end
    savedocfig("bsplines/quadrature-points")
end

function bsplines_function_interpolation()
    t = LinearKnotSet(7, 0, 7, 10)
    B = BSpline(t)
    x = range(0, stop=7, length=301)
    χ = B[x, :]
    c = B \ sin.(axes(B,1))
    i = 1:size(B,2)

    cfigure("function bsplines_interpolation",figsize=(7,9)) do
        csubplot(311, nox=true) do
            plot(x, χ*c)
            plot([mean_position(x, view(χ, :, j)) for j in i], c, "s-")
        end
        csubplot(312, nox=true) do
            plot(x, χ*c-sin.(x))
            ylabel("Error")
        end
        csubplot(313) do
            plot(x, χ)
            xlabel(L"x")
        end
    end
    savedocfig("bsplines/function-interpolation")
end

function bsplines_restricted_basis_interpolation()
    t = LinearKnotSet(7, 0.0, 1.0, 6)
    x = range(first(t), stop=last(t), length=300)[2:end-1]

    B = BSpline(t)
    B̃ = B[:,2:end-1]

    f1 = x -> sin(2π*x)
    f2 = x -> cos(2π*x)

    xx = axes(B, 1)
    c1 = B \ f1.(xx)
    c2 = B \ f2.(xx)

    c̃1 = B̃ \ f1.(xx)
    c̃2 = B̃ \ f2.(xx)

    χ = B[x, :]
    χ̃ = B̃[x, :]
    x_avg = [mean_position(x, view(χ, :, j)) for j in 1:size(B,2)]
    x̃_avg = [mean_position(x, view(χ̃, :, j)) for j in 1:size(B̃,2)]

    yl = nothing
    cfigure("restricted basis interpolation", figsize=(7,9)) do
        csubplot(321, nox=true) do
            l=plot(x, χ*c1)[1]
            plot(x_avg, c1, "s:", color=l.get_color())
            l=plot(x, χ*c2)[1]
            plot(x_avg, c2, "s:", color=l.get_color())
        end
        csubplot(322, nox=true) do
            l=plot(x, χ̃*c̃1)[1]
            plot(x̃_avg, c̃1, "s:", color=l.get_color())
            l=plot(x, χ̃*c̃2)[1]
            plot(x̃_avg, c̃2, "s:", color=l.get_color())
            axes_labels_opposite(:y)
        end
        csubplot(324, nox=true) do
            semilogy(x, abs.(χ̃*c̃1 - f1.(x)))
            semilogy(x, abs.(χ̃*c̃2 - f2.(x)))
            axes_labels_opposite(:y)
            yl = ylim()
            ylabel("Error")
        end
        csubplot(323, nox=true) do
            semilogy(x, abs.(χ*c1 - f1.(x)))
            semilogy(x, abs.(χ*c2 - f2.(x)))
            ylabel("Error")
            ylim(yl)
        end
        csubplot(325) do
            plot(x, χ)
            xlabel(L"x")
            yl = ylim()
        end
        csubplot(326) do
            plot(x, χ̃)
            axes_labels_opposite(:y)
            ylim(yl)
            xlabel(L"x")
        end
    end
    savedocfig("bsplines/restricted-basis-interpolation")
end

function bsplines_smooth_interpolation()
    f = x -> sin(2π*x)

    rng = MersenneTwister(123);

    N = 10
    x = clamp.(sort(range(0, stop=1, length=N) + 0.1(2rand(rng,N) .- 1)), 0, 1);
    y = f.(x) + 0.1(2rand(rng,N) .- 1);

    t3 = LinearKnotSet(3, 0.0, 1.0, 6);
    t4 = LinearKnotSet(4, 0.0, 1.0, 6);
    B3 = BSpline(t3,k′=1)
    B4 = BSpline(t4,k′=1)

    c3 = B3[x,:] \ y
    c4 = B4[x,:] \ y

    r = range(first(t3), stop=last(t3), length=300)
    χ3 = B3[r,:]
    χ4 = B4[r,:]

    cfigure("smooth interpolation") do
        plot(x, y, "s", label="Samples")
        plot(r, χ3*c3, label="3rd order spline")
        plot(r, χ4*c4, label="4th order spline")
        plot(r, f.(r), "--", label=L"\sin(2\pi x)")
        legend()
    end
    savedocfig("bsplines/smooth-interpolation")
end

function bsplines_diagonal_operators()
    k = 7
    N = 31

    a,b = 0,70

    coulomb(r) = -1/r

    cfigure("V(x)", figsize=(7,9)) do
        for (j,(t,x)) in enumerate([(LinearKnotSet(k, a, b, N),
                                     range(a, stop=b, length=500)[2:end-1]),
                                    (ExpKnotSet(k, -1.0, log10(b), N),
                                     10 .^ range(-1.0, stop=log10(b), length=500)[2:end-1])])
            B = BSpline(t)[:,2:end-1]
            S = B'B

            χ = B[x,:]

            xx = axes(B,1)

            f = B \ (x -> x^2*exp(-x)).(xx)
            g = B \ (x -> -x*exp(-x)).(xx)

            V = B'*QuasiDiagonal(coulomb.(xx))*B
            g̃ = S \ V*f

            csubplot(3,2,(j-1)+1, nox=true) do
                plot(x, χ*f, label=L"f(x)")
                plot(x, χ*g̃, label=L"\tilde{g}(x)")
                plot(x, χ*g, "--", label=L"g(x)")
                yl=ylim()
                plot(x, coulomb.(x), label=L"V(x)")
                ylim(yl)
                legend(framealpha=0.75)
                xscale("log")
                iseven(j) && axes_labels_opposite(:y)
            end
            csubplot(3,2,(j-1)+3, nox=true) do
                plot(x, χ*(g-g̃), label=L"g(x)-\tilde{g}(x)")
                legend(framealpha=0.75)
                xscale("log")
                ylabel("Error")
                iseven(j) && axes_labels_opposite(:y)
            end
            csubplot(3,2,(j-1)+5) do
                plot(x, χ)
                xscale("log")
                xlabel(L"x")
                iseven(j) && axes_labels_opposite(:y)
            end
        end
    end
    savedocfig("bsplines/diagonal-operators")
end

function find_second_derivative(B, f::Function)
    S = B'*B
    x = axes(B,1)
    D = Derivative(x)
    ∇² = B'*D'*D*B

    # Project function onto B-spline basis
    cf = B \ f.(x)
    # Find derivative
    cg = S \ ∇²*cf

    cf,cg
end

function bsplines_sine_derivative()
    t = LinearKnotSet(10, 0, 10, 30);
    B = BSpline(t)[:,2:end-1]

    f = x -> sin(2π*x)
    g = x -> -4π^2*sin(2π*x)

    cf,cg = find_second_derivative(B, f)

    x = range(first(t), stop=last(t), length=1001)
    χ = B[x, :]

    cfigure("derivatives", figsize=(7,9)) do
        csubplot(411,nox=true) do
            l=plot(x, χ*cf, label=L"\tilde{f}(x)")[1]
            plot(x, f.(x), ":", label=L"f(x)")
            legend(framealpha=0.75)
        end
        csubplot(412,nox=true) do
            l=plot(x, (χ*cg), label=L"\tilde{g}(x)")[1]
            plot(x, g.(x), ":", label=L"g(x)")
            legend(framealpha=0.75)
        end
        csubplot(413,nox=true) do
            semilogy(x, abs.(χ*cg-g.(x)), label=L"|\tilde{g}(x)-g(x)|")
            legend(framealpha=0.75)
            ylim(1e-5,1)
            ylabel("Error")
        end
        csubplot(414) do
            plot(x, χ)
            xlabel(L"x")
        end
    end
    savedocfig("bsplines/sine-derivative")
end

function bsplines_ode_hookes_law(xₘₐₓ, kspring, k, N)
    t = LinearKnotSet(k, 0, xₘₐₓ, N)
    # By omitting the first basis function, we enforce V(0) = 0
    B = BSpline(t,k′=1)[:,2:end]
    xx = axes(B, 1)
    S = B'B

    D = Derivative(xx)
    ∇ = B'*D*B

    # Hooke's law
    F = x -> -kspring*x
    # Exact potential
    V = x -> kspring*x^2/2

    # Expand Hooke's law on B-splines
    cF = B \ F.(xx)
    # Solve for expansion coefficients of potential
    cV = -∇ \ S*cF

    x = range(first(t), stop=last(t), length=500)
    χ = B[x,:]
    x_avg = [mean_position(x, view(χ, :, j)) for j in 1:size(B,2)]

    cfigure("Hooke's law",figsize=(7,9)) do
        csubplot(411,nox=true) do
            l=plot(x, χ*cF, label=L"\tilde{F}(x)")[1]
            plot(-x, -χ*cF, "--", color=l.get_color())
            plot(x_avg, cF, ".:", color=l.get_color(), label=L"c_F")
            plot(x, F.(x), ":", label=L"F(x)")
            legend(framealpha=0.75)
        end
        csubplot(412,nox=true) do
            l=plot(x, χ*cV, label=L"\tilde{V}(x)")[1]
            plot(-x, χ*cV, "--", color=l.get_color())
            plot(x_avg, cV, ".:", color=l.get_color(), label=L"c_V")
            plot(x, V.(x), ":", label=L"V(x)")
            legend(framealpha=0.75)
        end
        csubplot(413,nox=true) do
            plot(x, χ*cV - V.(x))
            xl = xlim()
            xlim(-xl[2],xl[2])
            ylabel("Error")
        end
        csubplot(414) do
            plot(x, χ)
            xl = xlim()
            xlim(-xl[2],xl[2])
            xlabel(L"x")
        end
    end
    savedocfig("bsplines/hookes-law-$(k)-$(N)")
end

struct ShiftAndInvert{TA,TB,TT}
    A⁻¹::TA
    B::TB
    temp::TT
end

Base.size(S::ShiftAndInvert, args...) = size(S.A⁻¹, args...)
Base.eltype(S::ShiftAndInvert) = eltype(S.A⁻¹)

function LinearAlgebra.mul!(y,M::ShiftAndInvert,x)
    mul!(M.temp, M.B, x)
    ldiv!(y, M.A⁻¹, M.temp)
end

construct_linear_map(A,B,σ=0) =
    ShiftAndInvert(factorize(A-σ*B),B,Vector{eltype(A)}(undef, size(A,1)))

function bsplines_hydrogen_eigenstates()

    k = 7
    N = 31

    a,b = 0,70

    coulomb(r) = -1/r

    nev = 5
    σ = -0.5

    n = 1:nev

    cfigure("Hydrogen", figsize=(7,9)) do
        for (j,(t,x,tol)) in enumerate([(LinearKnotSet(k, a, b, N),
                                         range(a, stop=b, length=500)[2:end-1],
                                         9e-3),
                                        (ExpKnotSet(k, -1.0, log10(b), N),
                                         10 .^ range(-1.0, stop=log10(b), length=500)[2:end-1],
                                         2e-7)])
            B = BSpline(t)[:,2:end-1]
            xx = axes(B, 1)
            S = B'B

            χ = B[x,:]

            D = Derivative(xx)

            ∇² = B'*D'*D*B

            T = -∇²/2

            V = B'*QuasiDiagonal(coulomb.(xx))*B

            H = T + V

            schurQR,history = partialschur(construct_linear_map(H, S, σ), nev=nev)
            println(history)

            θ = schurQR.eigenvalues
            E = real(σ .+ inv.(θ))

            csubplot(3,2,(j-1)+1, nox=true) do
                for i = 1:nev
                    ϕ = schurQR.Q[:,i]
                    plot(x, real(E[i]) .+ abs2.(χ*ϕ)/(ϕ'S*ϕ)[1])
                end
                # legend(framealpha=0.75)
                xscale("log")
                iseven(j) && axes_labels_opposite(:y)
            end
            csubplot(3,2,(j-1)+3, nox=true) do
                plot(x, coulomb.(x))
                xscale("log")
                iseven(j) && axes_labels_opposite(:y)
            end
            csubplot(3,2,(j-1)+5) do
                plot(x, χ)
                plot(x, χ)
                xscale("log")
                xlabel(L"x")
                iseven(j) && axes_labels_opposite(:y)
            end
        end
    end
    savedocfig("bsplines/hydrogen-eigenstates")
end

mkpath(joinpath(fig_dir, "bsplines"))
@echo bsplines_cardinal_splines()
@echo bsplines_discontinuous_splines()
@echo bsplines_full_multiplicity_splines()
@echo bsplines_spline1d()
@echo bsplines_spline2d()
@echo bsplines_quadrature_points()
@echo bsplines_function_interpolation()
@echo bsplines_restricted_basis_interpolation()
@echo bsplines_smooth_interpolation()
@echo bsplines_diagonal_operators()
@echo bsplines_sine_derivative()
@echo bsplines_ode_hookes_law(3, 0.1, 7, 30)
@echo bsplines_ode_hookes_law(3, 0.1, 3, 1)
@echo bsplines_hydrogen_eigenstates()
