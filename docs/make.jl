using MultivariateInterpolation
using Documenter

DocMeta.setdocmeta!(MultivariateInterpolation, :DocTestSetup, :(using MultivariateInterpolation); recursive=true)

makedocs(;
    modules=[MultivariateInterpolation],
    authors="Viktor Korotynskiy <korotynskiy.viktor@gmail.com> and contributors",
    sitename="MultivariateInterpolation.jl",
    format=Documenter.HTML(;
        canonical="https://azoviktor.github.io/MultivariateInterpolation.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/azoviktor/MultivariateInterpolation.jl",
    devbranch="main",
)
