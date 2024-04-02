export AffineSpace


"""
    AffineSpace <: AbstractAlgebraicVariety

An [`AbstractAlgebraicVariety`](@ref) that represents an affine space ``\\mathbb{C}^n``.

# Constructors
```julia
AffineSpace(vars::AbstractVector{Variable})
AffineSpace(vars::AbstractArray)
```

# Examples
```jldoctest
julia> @var R[1:3,1:3] t[1:3]
(Variable[R₁₋₁ R₁₋₂ R₁₋₃; R₂₋₁ R₂₋₂ R₂₋₃; R₃₋₁ R₃₋₂ R₃₋₃], Variable[t₁, t₂, t₃])

julia> AffineSpace(t)
AffineSpace of dimension 3
 variables: t₁, t₂, t₃

julia> AffineSpace(R)
AffineSpace of dimension 9
 variables: R₁₋₁, R₂₋₁, R₃₋₁, R₁₋₂, R₂₋₂, R₃₋₂, R₁₋₃, R₂₋₃, R₃₋₃

julia> AffineSpace([R, t])
AffineSpace of dimension 12
 variables: R₁₋₁, R₂₋₁, R₃₋₁, R₁₋₂, R₂₋₂, R₃₋₂, R₁₋₃, R₂₋₃, R₃₋₃, t₁, t₂, t₃
```
"""
struct AffineSpace <: AbstractAlgebraicVariety
    vars::Vector{Variable}
    samples::Dict{FixedFreeVariables, FixedFreeSamples}

    function AffineSpace(vars, samples)
        if isempty(vars)
            error("Vector of variables must be nonempty")
        end
        return new(vars, samples)
    end
end

AffineSpace(
    vars::AbstractArray
) = AffineSpace(Variable.(collect(flatten(vars))), Dict())

variables(𝔸::AffineSpace) = 𝔸.vars
dimension(𝔸::AffineSpace; kwargs...) = nvariables(𝔸)

function Base.show(io::IO, 𝔸::AffineSpace)
    println(io, "AffineSpace of dimension $(dimension(𝔸))")
    print(io, " variables: ", join(variables(𝔸), ", "))
end

expressions(::AffineSpace) = [Expression(0)]
generate_sample(𝔸::AffineSpace) = rand_unit(ComplexF64, nvariables(𝔸))

tangent_space(
    𝔸::AffineSpace,
    args...;
    kwargs...
) = rand_unit(ComplexF64, nvariables(𝔸), nvariables(𝔸))

finite_dominant_projection(𝔸::AffineSpace; kwargs...) = ExpressionMap(𝔸, variables(𝔸))

sample(
    𝔸::AffineSpace,
    vars::FixedFreeVariables;
    nsamples::Int=1,
    kwargs...
) = FixedFreeSamples(
        rand_unit(ComplexF64, nfixed(vars)),
        rand_unit(ComplexF64, nfree(vars), nsamples)
    )

function sample!(
    𝔸::AffineSpace,
    vars::FixedFreeVariables;
    nsamples::Int=1,
    kwargs...
)
    s = sample(𝔸, vars; nsamples=nsamples)
    𝔸.samples[vars] = s
end

samples(
    𝔸::AffineSpace,
    vars::FixedFreeVariables
) = get(𝔸.samples, vars, nothing)