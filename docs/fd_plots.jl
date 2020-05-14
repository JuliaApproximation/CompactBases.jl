import CompactBases: locs, weights

function compare_staggered_non_uniform()
    N = 50
    rmax = 30.0

    B1 = StaggeredFiniteDifferences(rmax, N)

    ρmin = 0.001
    ρmax = 1.0
    α = 0.1
    B2 = StaggeredFiniteDifferences(ρmin, ρmax, α, rmax)

    S1 = B1'B1
    S2 = B2'B2

    x1 = locs(B1)
    x2 = locs(B2)
    w1 = weights(B1)
    w2 = weights(B2)

    cfigure("staggered non_uniform grid") do
        csubplot(211,nox=true) do
            plot(eachindex(x1), x1, ".-", label="Uniform")
            plot(eachindex(x2), x2, ".-", label="Non-uniform")
            ylabel("Node location")
            legend()
        end
        csubplot(212) do
            plot(eachindex(x1), w1, ".-")
            plot(eachindex(x2), w2, ".-")
            ylabel("Weights")
            xlabel("Node index")
        end
    end
    savefig("docs/src/figures/fd/compare_staggered_non_uniform_grid.svg")

    ξ = 10.0 .^ range(-3, stop=log10(25), length=1000)
    χ1 = B1[ξ, :]
    χ2 = B2[ξ, :]

    cfigure("staggered non_uniform basis") do
        csubplot(211,nox=true) do
            semilogx(ξ, χ1)
            ylabel("Uniform")
        end
        csubplot(212) do
            semilogx(ξ, χ2)
            ylabel("Non-uniform")
        end
    end
    savefig("docs/src/figures/fd/compare_staggered_non_uniform_basis.svg")

    f = x -> sin(2π*x/rmax)*exp(-4x/rmax)

    xx1 = axes(B1,1)
    c1 = B1 \ f.(xx1)
    xx2 = axes(B2,1)
    c2 = B2 \ f.(xx2)

    f1 = χ1*c1
    f2 = χ2*c2
    fe = f.(ξ)

    cfigure("staggered non_uniform reconstruction") do
        csubplot(211,nox=true) do
            l = plot(x1, c1, ".:", label="Uniform coeffs")[1]
            plot(ξ, f1, color=l.get_color(), label="Uniform reconstruction")
            l = plot(x2, c2, ".:", label="Non-uniform coeffs")[1]
            plot(ξ, f2, "--", color=l.get_color(), label="Non-uniform reconstruction")
            plot(ξ, fe, ":", label="Exact")
            xscale("log")
            legend(loc=2)
        end
        csubplot(212) do
            loglog(ξ, abs.(f1-fe), label="Uniform error")
            loglog(ξ, abs.(f2-fe), label="Non-uniform error")
            legend(loc=2)
            xlabel(L"x")
        end
    end
    savefig("docs/src/figures/fd/compare_staggered_non_uniform_reconstruction.svg")
end

mkpath("docs/src/figures/fd")
@echo compare_staggered_non_uniform()
