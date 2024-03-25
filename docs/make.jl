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
        canonical="https://multivariatepolynomialsystems.github.io/MultivariateInterpolation.jl",
        edit_link="main",
        assets=["assets/custom.css"],
        collapselevel=2
    ),
    pages=[
        "Introduction" => "index.md",
        "Common interpolation techniques" => [
            "Interpolation basis" => "inter_basis.md",
            "Fixed-free interpolation" => "fixed_free.md"
        ],
        "Functions on Algebraic Varieties" => [
            "Algebraic varieties" => "on_varieties/varieties.md"
            "ExpressionMap" => "on_varieties/maps.md"
            "Sampling" => "on_varieties/sampling.md"
            "Interpolation" => "on_varieties/interpolation.md"
            # "Scaling symmetries" => "variety/scalings.md"
        ],
        "Functions on Affine Spaces" => [

        ]
    ],
)

deploydocs(;
    repo="github.com/MultivariatePolynomialSystems/MultivariateInterpolation.jl.git",
    devbranch="main",
)
