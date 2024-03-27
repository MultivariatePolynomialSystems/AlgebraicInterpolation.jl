export AlgebraicVariety


"""
    AlgebraicVariety <: AbstractAlgebraicVariety
"""
struct AlgebraicVariety <: AbstractAlgebraicVariety
    system::System
    jacobian::Matrix{Expression}
    samples::Dict{FixedFreeVariables, FixedFreeSamples}
    sample_generator::Union{Function, Nothing}
end

AlgebraicVariety(
    exprs::AbstractVector{Expression};
    variables::AbstractVector{Variable}=variables(exprs)
) = AlgebraicVariety(System(exprs; variables=variables), differentiate(exprs, variables), Dict(), nothing)

AlgebraicVariety(
    exprs::AbstractVector{Expression},
    sample_generator::Function;
    variables::AbstractVector{Variable}
) = AlgebraicVariety(System(exprs; variables=variables), differentiate(exprs, variables), Dict(), sample_generator)

expressions(X::AlgebraicVariety) = expressions(X.system)
variables(F::System) = vcat(HC.variables(F), HC.parameters(F)) # WARNING: redefine variables from HC
variables(X::AlgebraicVariety) = variables(X.system)
HC.System(X::AlgebraicVariety) = X.system

jacobian(X::AlgebraicVariety) = X.jacobian

samples(
    X::AlgebraicVariety,
    vars::FixedFreeVariables
) = get(X.samples, vars, nothing)

function nsamples(X::AlgebraicVariety, vars::FixedFreeVariables)
    s = samples(X, vars)
    !isnothing(s) && return size(free(s), 2)
    return 0
end

add_samples!(
    X::AlgebraicVariety,
    vars::FixedFreeVariables,
    s::FixedFreeSamples;
    tol::Real=1e-10
) = X.samples[vars] = isnothing(samples(X, vars)) ? s : hcat(samples(X, vars), s; tol=tol)

function generate_sample(X::AlgebraicVariety)
    vars = FixedFreeVariables(variables(X))
    s = samples(X, vars)
    if isnothing(s) && !isnothing(X.sample_generator)
        s = X.sample_generator(vars=vars, nsamples=1)
        return isnothing(s) ? generate_sample(System(X)) : free(s)[:]
    end
    return isnothing(s) ? generate_sample(System(X)) : free(s)[:,rand(1:nsamples(s))]
end

# HC.subs(
#     X::AlgebraicVariety,
#     substitutions::Pair...
# ) = AlgebraicVariety(
#         subs(expressions(F), substitutions...);
#         variables=setdiff(variables(F), vcat(first.(substitutions)...))
#     )