export AffineSpace


"""
    AffineSpace <: AbstractAlgebraicVariety
"""
struct AffineSpace <: AbstractAlgebraicVariety
    vars::Vector{Variable}
    samples::Dict{FixedFreeVariables, FixedFreeSamples}
end

variables(𝔸::AffineSpace) = 𝔸.vars
expressions(::AffineSpace) = [Expression(0)]
generate_sample(𝔸::AffineSpace) = rand_unit(ComplexF64, nvariables(𝔸))

tangent_space(
    𝔸::AffineSpace,
    args...;
    kwargs...
) = rand_unit(ComplexF64, nvariables(𝔸), nvariables(𝔸))

dimension(𝔸::AffineSpace; kwargs...) = nvariables(𝔸)
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