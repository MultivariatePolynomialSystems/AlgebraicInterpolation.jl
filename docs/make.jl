push!(LOAD_PATH,"../src/")

using AlgebraicInterpolation
using Documenter

DocMeta.setdocmeta!(AlgebraicInterpolation, :DocTestSetup, :(using AlgebraicInterpolation); recursive=true)

makedocs(;
    modules=[AlgebraicInterpolation],
    authors="Viktor Korotynskiy <korotynskiy.viktor@gmail.com> and contributors",
    repo="https://github.com/MultivariatePolynomialSystems/AlgebraicInterpolation.jl/blob/{commit}{path}#{line}",
    sitename="AlgebraicInterpolation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://multivariatepolynomialsystems.github.io/AlgebraicInterpolation.jl",
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
            "Lie symmetries" => "on_varieties/lie-symmetries.md"
            "Sampling" => "on_varieties/sampling.md"
            "Interpolation" => "on_varieties/interpolation.md"
        ],
        "Functions on Affine Spaces" => [

        ]
    ],
)

deploydocs(;
    repo="github.com/MultivariatePolynomialSystems/AlgebraicInterpolation.jl.git",
    devbranch="main",
)
