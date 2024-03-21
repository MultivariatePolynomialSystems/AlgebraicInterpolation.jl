export SampledVariety


"""
    SampledVariety <: AbstractSampledVariety
"""
struct SampledVariety <: AbstractSampledVariety
    system::DifferentiatedVariety
    samples::Dict{FixedFreeVariables, FixedFreeSamples}
    sample_generator::Union{Function, Nothing}
end

SampledVariety(
    F::System,
    sample_generator::Union{Function, Nothing}
) = SampledVariety(DifferentiatedVariety(F), Dict(), sample_generator)
SampledVariety(F::System) = SampledVariety(F, nothing)

HC.System(F::SampledVariety) = System(F.system)

unknowns(F::Union{DifferentiatedVariety, SampledVariety}) = unknowns(System(F))
HC.parameters(F::Union{DifferentiatedVariety, SampledVariety}) = parameters(System(F))
variables(F::Union{DifferentiatedVariety, SampledVariety}) = variables(System(F))
nunknowns(F::Union{DifferentiatedVariety, SampledVariety}) = nunknowns(System(F))
HC.nparameters(F::Union{DifferentiatedVariety, SampledVariety}) = nparameters(System(F))
nvariables(F::Union{DifferentiatedVariety, SampledVariety}) = nvariables(System(F))
HC.expressions(F::Union{DifferentiatedVariety, SampledVariety}) = expressions(System(F))

HC.subs(
    F::System,
    substitutions::Pair...
) = System(
        subs(expressions(F), substitutions...);
        variables=setdiff(variables(F), vcat(first.(substitutions)...))
    )

(F::SampledVariety)(x::AbstractVector) = System(F)(x)

samples(
    F::SampledVariety,
    vars::FixedFreeVariables
) = get(F.samples, vars, nothing)

function nsamples(F::SampledVariety, vars::FixedFreeVariables)
    s = samples(F, vars)
    !isnothing(s) && return size(free(s), 2)
    return 0
end

add_samples!(
    F::SampledVariety,
    vars::FixedFreeVariables,
    s::FixedFreeSamples;
    tol::Real=1e-10
) = F.samples[vars] = isnothing(samples(F, vars)) ? s : hcat(samples(F, vars), s; tol=tol)

function find_sample(F::SampledVariety; kwargs...)
    vars = FixedFreeVariables(variables(F))
    s = samples(F, vars)
    if isnothing(s) && !isnothing(F.sample_generator)
        s = F.sample_generator(vars=vars, nsamples=1)
        return isnothing(s) ? find_sample(System(F); kwargs...) : free(s)[:]
    end
    return isnothing(s) ? find_sample(System(F); kwargs...) : free(s)[:,rand(1:nsamples(s))]
end

full_jacobian(F::SampledVariety) = full_jacobian(F.system)

# TODO: add keyword arg for same fixed values?
"""
    sample!(F::Union{SampledSystem, MapGraph}, vars::FixedFreeVariables; kwargs...)

Generates samples using [`sample`](@ref) method and updates `F`. The computed samples can be
then obtained by running `samples(F, vars)`.
"""
function sample!(
    F::SampledVariety,
    vars::FixedFreeVariables;
    nsamples::Int=1,
    start_point::Union{AbstractVector, Nothing}=nothing,
    tols::Tolerances=Tolerances(),
    rand_method::Symbol=:rand_unit
)
    s = sample(
        F,
        vars;
        nsamples=nsamples,
        start_point=start_point,
        tols=tols,
        rand_method=rand_method
    )
    F.samples[vars] = s
end
