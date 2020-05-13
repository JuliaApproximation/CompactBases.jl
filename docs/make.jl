using Documenter, CompactBases

isdefined(Main, :NOPLOTS) && NOPLOTS || include("plots.jl")

DocMeta.setdocmeta!(CompactBases, :DocTestSetup, :(using CompactBases); recursive=true)
makedocs(;
    modules=[CompactBases],
    format=Documenter.HTML(assets = ["assets/latex.js"]),
    pages = [
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Usage" => [
            "Basis creation" => "usage.md",
            "Knot sets" => "knot_sets.md",
            "Splines" => [
                "Spline creation & evaluation" => "splines.md",
                "Function approximation" => "function_approximation.md",
            ],
            "Approximating operators" => "operators.md",
            "Examples" => [
                "Differentiating functions" => "differentiation.md",
                "Ordinary differential equations" => "odes.md",
                "Eigenproblems" => "eigenproblems.md"
            ]
        ],
    ],
    repo="https://github.com/JuliaApproximation/CompactBases.jl/blob/{commit}{path}#L{line}",
    sitename="CompactBases.jl",
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com>",
)

deploydocs(;
    repo="github.com/JuliaApproximation/CompactBases.jl",
)
