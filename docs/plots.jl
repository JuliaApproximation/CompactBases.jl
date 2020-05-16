using PyPlot
using Jagot.plotting
plot_style("ggplot")
using PyPlotRecipes
using PyCall
PathEffects = pyimport("matplotlib.patheffects")
using Statistics
using Random
using Colors

using LinearAlgebra
using ArnoldiMethod

using CompactBases

function mean_color(color::String)
    c = parse(Colorant, color)
    mean([c.r,c.g,c.b])
end

lerp(a,b,t) = (1-t)*a + t*b

mean_position(x, ϕ) = ϕ'*Diagonal(x)*ϕ/(ϕ'ϕ)

function savedocfig(name,dir="figures")
    fig = gcf()
    savefig(joinpath("docs/src",dir,"$(name).svg"),
            transparent=true,
            facecolor=fig.get_facecolor())
    close(fig)
end

function logo()
    t = ArbitraryKnotSet(3, [0.0, 1, 1, 3, 4, 6], 1, 3)
    r = range(first(t), stop=last(t), length=301)
    R = BSpline(t)
    χ = R[r, :]

    cfigure("logo",figsize=(7,3)) do
        plot(r, χ, length(r) < 30 ? ".-" : "-", linewidth=4)
        margins(0.1, 0.1)
        axis("off")
    end
    savefig("docs/src/assets/logo.svg", transparent=true)
end

function simple_example()
    a,b = 0.0,1.0 # Extents
    N = 13 # Number of nodes
    k = 5 # Order of FE-DVR/B-splines
    x = range(a, stop=b, length=1000)

    Δx = (b-a)/(N+1) # Grid spacing
    # Standard, uniform finite-differences
    fd = FiniteDifferences(N, Δx)
    # Staggered, uniform finite-differences
    sfd_uni = StaggeredFiniteDifferences(N, Δx)
    # Staggered, non-uniform finite-differences
    sfd_nonuni = StaggeredFiniteDifferences(0.01, 0.5, 0.1, b)

    # Finite-element boundaries
    tf = range(a, stop=b, length=N+2)
    # We can vary the polynomial order in each element
    forder = vcat(7, fill(4,length(tf)-2))
    # By indexing the second dimension, we can implement Dirichlet0
    # boundary conditions.
    fem = FEDVR(tf, forder)[:,2:end-1]

    tb = ExpKnotSet(k, -2.0, log10(b), N+1)
    splines = BSpline(tb)[:,2:end-1]

    bases = [fd, sfd_uni, sfd_nonuni, fem, splines]
    m = length(bases)

    cfigure("simple example", figsize=(6,9)) do
        for (i,basis) in enumerate(bases)
            csubplot(m,1,i,nox=i<m) do
                plot(x, basis[x,:])
            end
        end
    end
    savedocfig("simple_example")
end

function restricted_bases()
    B = FEDVR(range(0, stop=10.0, length=5), 7)
    B̃ = B[:, 3:17]

    x = range(0, stop=10.0, length=1000)
    χ = B[x, :]
    χ̃ = B̃[x, :]

    cfigure("restricted bases") do
        csubplot(211,nox=true) do
            plot(x, χ)
            ylabel("Original basis")
        end
        csubplot(212) do
            plot(x, χ̃)
            ylabel("Restricted basis")
        end
    end
    savedocfig("restricted_bases")
end

macro echo(expr)
    println(expr)
    :(@time $expr)
end

@info "Documentation plots"
mkpath("docs/src/figures")
@echo logo()
@echo simple_example()
@echo restricted_bases()
include("bspline_plots.jl")
include("fd_plots.jl")
