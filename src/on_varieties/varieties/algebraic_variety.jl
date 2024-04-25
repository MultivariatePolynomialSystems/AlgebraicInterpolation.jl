export AlgebraicVariety


"""
    AlgebraicVariety <: AbstractAlgebraicVariety

An [`AbstractAlgebraicVariety`](@ref) that represents basic algebraic variety defined
by polynomial equations.

# Constructors
```julia
AlgebraicVariety(
    exprs::AbstractVector{Expression},
    sample_generator::Function;
    variables::AbstractVector{Variable}
)
AlgebraicVariety(
    exprs::AbstractVector{Expression};
    variables::AbstractVector{Variable}
)
AlgebraicVariety(exprs::AbstractArray, sample_generator::Function; variables::AbstractArray)
AlgebraicVariety(exprs::AbstractArray; variables::AbstractArray)
```

# Examples
```jldoctest
julia> @var R[1:2,1:2] t[1:2]
(Variable[R₁₋₁ R₁₋₂; R₂₋₁ R₂₋₂], Variable[t₁, t₂])

julia> AlgebraicVariety([R'*R-I, det(R)-1]; variables=[R, t])
AlgebraicVariety X ⊂ ℂ⁶
 6 variables: R₁₋₁, R₂₋₁, R₁₋₂, R₂₋₂, t₁, t₂
 5 expressions: 
  -1 + R₁₋₁^2 + R₂₋₁^2
  R₁₋₁*R₁₋₂ + R₂₋₁*R₂₋₂
  R₁₋₁*R₁₋₂ + R₂₋₁*R₂₋₂
  -1 + R₁₋₂^2 + R₂₋₂^2
  -1 + R₁₋₁*R₂₋₂ - R₁₋₂*R₂₋₁
```
"""
struct AlgebraicVariety <: AbstractAlgebraicVariety
    system::System
    variables::Vector{Variable}
    jacobian::Matrix{Expression}
    samples::Dict{FixedFreeVariables, FixedFreeSamples}
    sample_generator::Union{Function, Nothing}
end

function AlgebraicVariety(
    exprs::AbstractArray,
    sample_generator::Union{Function, Nothing};
    variables::AbstractArray
) 
    exprs = Expression.(collect(flatten(exprs)))
    vars = Variable.(collect(flatten(variables)))
    F = System(exprs; variables=vars)
    jac = differentiate(exprs, vars)
    return AlgebraicVariety(F, vars, jac, Dict(), sample_generator)
end

AlgebraicVariety(
    exprs::AbstractArray;
    variables::AbstractArray
) = AlgebraicVariety(exprs, nothing; variables=variables)

expressions(X::AlgebraicVariety) = HC.expressions(X.system)
variables(F::System) = vcat(HC.variables(F), HC.parameters(F)) # WARNING: redefine variables from HC
variables(X::AlgebraicVariety) = X.variables
nvariables(X::AlgebraicVariety) = length(X.variables)

function Base.show(io::IO, X::AlgebraicVariety, offset::String)
    println(io, "$(offset)AlgebraicVariety X ⊂ ℂ$(superscript(nvariables(X)))")
    println(io, "$(offset) $(phrase(nvariables(X), "variable")): ", join(variables(X), ", "))
    println(io, "$(offset) $(phrase(nexpressions(X), "expression")): ")
    for (j, expr) in enumerate(expressions(X))
        print(io, "$(offset)  ", expr)
        j < nexpressions(X) && print(io, "\n")
    end
end

Base.show(io::IO, X::AlgebraicVariety) = show(io, X, "")

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

HC.subs(
    X::AlgebraicVariety,
    substitutions::Pair...
) = AlgebraicVariety(
        subs(expressions(F), substitutions...);
        variables=setdiff(variables(F), vcat(first.(substitutions)...))
    )