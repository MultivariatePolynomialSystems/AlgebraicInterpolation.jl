export AbstractAlgebraicVariety,
    variables,
    nvariables,
    expressions,
    nexpressions,
    generate_sample,
    jacobian,
    tangent_space,
    dimension,
    finite_dominant_projection,
    sample,
    sample!,
    add_samples!,
    samples


"""
    abstract type AbstractAlgebraicVariety end
"""
abstract type AbstractAlgebraicVariety end

"""
    variables(X::AbstractAlgebraicVariety) -> Vector{Variable}

Returns the variables of `X`.
"""
function variables(X::AbstractAlgebraicVariety)
    error("Not implemented")
end

"""
    nvariables(X::AbstractAlgebraicVariety) -> Int

Returns the number of variables of `X`.
"""
nvariables(X::AbstractAlgebraicVariety) = length(variables(X))

"""
    expressions(X::AbstractAlgebraicVariety) -> Vector{Expression}

Returns the expressions of `X`.
"""
function expressions(X::AbstractAlgebraicVariety)
    error("Not implemented")
end

"""
    nexpressions(X::AbstractAlgebraicVariety) -> Int

Returns the number of expressions of `X`.
"""
nexpressions(X::AbstractAlgebraicVariety) = length(expressions(X))

HC.System(X::AbstractAlgebraicVariety) = System(expressions(X); variables=variables(X))

"""
    generate_sample(X::AbstractAlgebraicVariety; <keyword arguments>) -> Vector{ComplexF64}

Generates a sample from `X`.
"""
function generate_sample(X::AbstractAlgebraicVariety; kwargs...)
    F = System(X)
    return generate_sample(F; kwargs...)
end

"""
    jacobian(X::AbstractAlgebraicVariety) -> Matrix{Expression}

Returns the symbolic jacobian of `X`.
"""
jacobian(X::AbstractAlgebraicVariety) = differentiate(expressions(X), variables(X))

"""
    jacobian(X::AbstractAlgebraicVariety, x::AbstractVector{<:Number}) -> Matrix{<:Number}

Returns the jacobian of `X` evaluated at `x`.
"""
jacobian(
    X::AbstractAlgebraicVariety,
    x::AbstractVector{<:Number}
) = evaluate(jacobian(X), variables(X) => x)

"""
    tangent_space(X::AbstractAlgebraicVariety, x::AbstractVector{<:Number}; <keyword arguments>) -> Matrix{<:Number}

Returns the tangent space of `X` at `x`.

# Keyword arguments
- `tols::Tolerances=Tolerances()`: tolerances for numerical computations.
"""
tangent_space(
    X::AbstractAlgebraicVariety,
    x::AbstractVector{<:Number};
    tols::Tolerances=Tolerances(),
    var_ids::Union{Vector{Int}, Colon}=:
) = nullspace(jacobian(X, x); atol=tols.nullspace_atol)[var_ids, :]

"""
    dimension(X::AbstractAlgebraicVariety; <keyword arguments>)

Computes the dimension of `X`.

# Keyword arguments
- `sample::Union{AbstractVector{<:Number}, Nothing}=nothing`: point that belongs to `X`.
- `tols::Tolerances=Tolerances()`: tolerances for numerical computations.
"""
function dimension(
    X::AbstractAlgebraicVariety;
    sample::Union{AbstractVector{<:Number}, Nothing}=nothing,
    tols::Tolerances=Tolerances()
)
    x = isnothing(x) ? generate_sample(X) : x
    if !isnothing(x)
        @assert norm(F(x)) < tols.common_tol
        J = full_jacobian(F, x)
        return nvariables(F) - rank(J; atol=tols.rank_atol)
    else
        error("Cannot generate a random sample of the algebraic variety, provide one!")
    end
end

"""
    finite_dominant_projection(X::AbstractAlgebraicVariety; <keyword arguments>) -> ExpressionMap

Returns a finite dominant projection from `X` to an affine space.

# Keyword arguments
- `sample::Union{AbstractVector, Nothing}=nothing`: point that belongs to `X`.
- `tols::Tolerances=Tolerances()`: tolerances for numerical computations.
"""
function finite_dominant_projection(
    X::AbstractAlgebraicVariety;
    sample::Union{AbstractVector{<:Number}, Nothing}=nothing,
    tols::Tolerances=Tolerances()
)
    x = isnothing(x) ? find_sample(F) : x
    proj_vars = Int[]
    dom_dim = dimension(F, x; tols=tols)
    im_dim = 0
    for i in 1:nvariables(F)
        push!(proj_vars, i)
        φ = ExpressionMap(F, proj_vars)
        im_dim_curr = image_dimension(φ, x; tols=tols)
        im_dim_curr == dom_dim && return φ
        if im_dim_curr > im_dim
            im_dim += 1
        else
            pop!(proj_vars)
        end
    end
    error("Couldn't find a finite dominant projection")
end

"""
    sample(X::AbstractAlgebraicVariety, vars::FixedFreeVariables; <keyword arguments>) -> FixedFreeSamples

Returns the samples of `X` in the given `vars`. 
Results in an error, if it is impossible or unreasonable to sample the given `vars` from `X`.

# Keyword arguments
- `nsamples::Integer=1`: number of samples.
- `start_point::Union{AbstractVector, Nothing}=nothing`: starting point for homotopy continuation.
- `tols::Tolerances=Tolerances()`: tolerances for numerical computations.
- `rand_method::Symbol=:rand_unit`: method for generating random samples.
"""
function sample(X::AbstractAlgebraicVariety, vars::FixedFreeVariables; kwargs...)
    error("Not implemented")
end

"""
    add_samples!(X::AbstractAlgebraicVariety, s::FixedFreeSamples)

Adds new samples `s` to `X`.
"""
function add_samples!(X::AbstractAlgebraicVariety, s::FixedFreeSamples)
    error("Not implemented")
end

"""
    sample!(X::AbstractAlgebraicVariety, vars::FixedFreeSamples; <keyword arguments>)

Samples `X` in the given `vars` and updates `X` with these samples.
"""
function sample!(X::AbstractAlgebraicVariety, vars::FixedFreeVariables; kwargs...)
    s = sample(X, vars; kwargs...)
    add_samples!(X, s)
end

"""
    samples(X::AbstractAlgebraicVariety, vars::FixedFreeVariables) -> FixedFreeSamples

Returns the saved samples of `X` in the given `vars`.
"""
function samples(X::AbstractAlgebraicVariety, vars::FixedFreeVariables)
    error("Not implemented")
end