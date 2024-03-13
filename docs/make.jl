push!(LOAD_PATH,"../src/")

using MultivariateInterpolation
using Documenter

DocMeta.setdocmeta!(MultivariateInterpolation, :DocTestSetup, :(using MultivariateInterpolation); recursive=true)

makedocs(;
    modules=[MultivariateInterpolation],
    authors="Viktor Korotynskiy <korotynskiy.viktor@gmail.com> and contributors",
    repo="https://github.com/MultivariatePolynomialSystems/MultivariateInterpolation.jl/blob/{commit}{path}#{line}",
    sitename="MultivariateInterpolation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MultivariatePolynomialSystems.github.io/MultivariateInterpolation.jl",
        edit_link="main",
        assets=["assets/custom.css"],
        collapselevel=2
    ),
    pages=[
        "Introduction" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MultivariatePolynomialSystems/MultivariateInterpolation.jl.git",
    devbranch="main",
)
