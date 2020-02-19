using Documenter, CompactBases

makedocs(;
    modules=[CompactBases],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/JuliaApproximation/CompactBases.jl/blob/{commit}{path}#L{line}",
    sitename="CompactBases.jl",
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/JuliaApproximation/CompactBases.jl",
)
