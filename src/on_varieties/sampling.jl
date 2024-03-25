export AbstractSampledVariety,
    SampledSystem,
    unknowns,
    parameters,
    variables,
    nunknowns,
    nparameters,
    nvariables,
    expressions,
    samples,
    nsamples,
    full_jacobian,
    find_sample,
    dimension,
    finite_dominant_projection,
    sample,
    sample!,
    possible_to_sample


struct VariableEvaluation
    eval::Dict{Variable, ComplexF64}
end

VariableEvaluation(
    vars::Vector{Variable}, 
    vals::Vector{ComplexF64}
) = VariableEvaluation(Dict(zip(vars, vals)))

Base.getindex(
    eval::VariableEvaluation,
    vars::Vector{Variable}
) = [eval.eval[v] for v in vars]

# WARNING: redefine 'variables' from HC; this definition is IMHO semantically more correct
unknowns(F::System) = HC.variables(F)
variables(F::System) = vcat(unknowns(F), parameters(F))
nunknowns(F::System) = length(unknowns(F))
nvariables(F::System) = length(variables(F))

function find_sample(F::System; kwargs...)
    pair = HC.find_start_pair(F; kwargs...)
    isnothing(pair) && error("Couldn't generate a random sample of the system")
    x, p = pair
    p = isnothing(p) ? ComplexF64[] : p
    return vcat(x, p)
end

full_jacobian(F::System) = differentiate(expressions(F), variables(F))







dimension(
    F::System,
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
) = dimension(DifferentiatedVariety(F), x; tols=tols)

function finite_dominant_projection(
    F::AbstractDifferentiatedVariety,
    x::Union{AbstractVector, Nothing}=nothing;
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

function finite_dominant_projection(
    F::System,
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
)
    φ = finite_dominant_projection(DifferentiatedVariety(F), x; tols=tols)
    params = image_vars(φ)
    unknws = setdiff(variables(F), params)
    return System(expressions(F); variables=unknws, parameters=params)
end

function reasonable_to_sample(
    F::AbstractDifferentiatedVariety,
    vars::FixedFreeVariables,
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
)
    dim_fixed = 0
    if !isempty(fixed(vars))
        φ = ExpressionMap(F, fixed(vars))
        dim_fixed = image_dimension(φ, x; tols=tols)
    end

    φ = ExpressionMap(F, variables(vars))
    dim_all = image_dimension(φ, x; tols=tols)
    return dim_all - dim_fixed ≠ nfree(vars)
end

"""
    possible_to_sample(X::AbstractDifferentiatedVariety, vars::FixedFreeVariables; kwargs...)

Returns `false` if there exists a free variable in `vars` that is determined in 
finite many ways by all fixed variables in `vars`. 

Throws a WARNING (for keyword argument `logging=true`) if there are 
no constraints in free variables after fixing all the fixed ones.

*Keyword arguments*:
* 
* `tols::Tolerances=Tolerances()`: tolerances structure used for computations. Tolerances used 
in this method: same as in [`image_dimension`](@ref).
* `logging::Bool=true`
"""
function possible_to_sample(
    F::AbstractDifferentiatedVariety,
    vars::FixedFreeVariables,
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances(),
    logging::Bool=true
)
    @assert variables(vars) ⊆ variables(F)
    x = isnothing(x) ? find_sample(F) : x

    dim_fixed = 0
    if !isempty(fixed(vars))
        φ = ExpressionMap(F, fixed(vars))
        dim_fixed = image_dimension(φ, x; tols=tols)
    end

    if logging
        φ = ExpressionMap(F, variables(vars))
        dim_all = image_dimension(φ, x; tols=tols)
        if dim_all - dim_fixed == nfree(vars)
            if isempty(fixed(vars))
                @warn "There are no constraints in $(join(free(vars), ", "))"
            else
                @warn "There are no constraints in $(join(free(vars), ", ")) after fixing $(join(fixed(vars), ", "))"
            end
        end
    end

    isempty(fixed(vars)) && return true

    for var in free(vars)
        φ = ExpressionMap(F, vcat(fixed(vars), var))
        dim_var = image_dimension(φ, x; tols=tols)
        if dim_var == dim_fixed
            logging && @info "Variable $(var) is determined in finite many ways after fixing $(join(fixed(vars), ", "))"
            return false
        end
    end

    return true
end

# justified by the fact that the smallest degree is achieved in maximal fixed case
# follows from the fact that by fixing additional vars we don't increase the degree
function max_fixed_to_sample(
    F::AbstractDifferentiatedVariety,
    vars::AbstractArray{Variable},
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
)
    @assert !isempty(vars)
    @assert vars ⊆ variables(F)
    for var in vars

    end
end

# Supposes that each Result contains only 1 (if SUCCESS) or 0 (if FAIL) solutions
function extract_samples(
    results::Vector{Tuple{Result, Vector{ComplexF64}}},
    F::System;
    resample::Bool=false
)
    n_instances = length(results)
    samples = zeros(ComplexF64, nvariables(F), n_instances)
    k = 1
    for (res, p) in results
        sols = HC.solutions(res)
        if length(sols) == 1
            samples[1:nunknowns(F), k] = first(sols)
            samples[nunknowns(F)+1:end, k] = p
            k += 1
        elseif !resample
            error("Number of solutions in the $(k)-th result is $(length(sols)), expected $(n_tracked)")
        end
    end
    for i in k:n_instances
        while true
            instance_id = rand(1:i-1)  # TODO: what if i == 1?
            x₀ = samples[1:nunknowns(F), instance_id]
            p₀ = last(results[instance_id])
            p₁ = randn(ComplexF64, nparameters(F))
            res = HC.solve(
                F,
                x₀,
                start_parameters = p₀,
                target_parameters = p₁
            )
            sols = HC.solutions(res)
            if length(sols) == 1
                samples[1:nunknowns(F), i] = first(sols)
                samples[nunknowns(F)+1:end, i] = p₁
                break
            end
        end
    end
    return samples
end

# TODO: add keyword arg for same fixed values?
"""
    sample(F::Union{SampledSystem, MapGraph}, vars::FixedFreeVariables; kwargs...)

Returns samples for the given polynomial system `F` and variables `vars`. The return type
is `FixedFreeSamples`.

*Keyword arguments*:
* `nsamples::Integer=1`: defines number of samples to compute
* `start_point::Union{AbstractVector, Nothing}=nothing`: specifies the starting point of the variety
  defined by `F` from which to start collecting other samples.
* `tols::Tolerances=Tolerances()`: tolerances structure used for computations. Tolerances used 
  in this method: `rank_atol`, 
* `rand_method::Symbol=:rand_unit`: method for generating random samples.
"""
function sample(
    F::AbstractDifferentiatedVariety,
    vars::FixedFreeVariables;
    nsamples::Integer=1,
    start_point::Union{AbstractVector, Nothing}=nothing,
    tols::Tolerances=Tolerances(),
    rand_method::Symbol=:rand_unit
)
    x = isnothing(x) ? find_sample(F) : x
    if !possible_to_sample(F, vars, x; tols=tols, logging=false)
        error("Impossible to sample")
    end
    if !reasonable_to_sample(F, vars, x; tols=tols)
        error("Sampling is unnecessary")
    end

    if F isa SampledSystem
        if !isnothing(F.sample_generator)
            s = F.sample_generator(nsamples=nsamples, vars=vars) # TODO: add rand_method?
            !isnothing(s) && return s
        end
    end

    x = VariableEvaluation(variables(F), x)
    x_fixed = x[fixed(vars)]
    G = subs(System(F), fixed(vars) => x_fixed)

    vars_not_fixed = setdiff(variables(F), fixed(vars)) # TODO: variables(G)?
    G = finite_dominant_projection(G, x[vars_not_fixed]; tols=tols)

    ps = [eval(rand_method)(ComplexF64, nparameters(G)) for _ in 1:nsamples]
    res = HC.solve(
        G,
        [x[unknowns(G)]],
        start_parameters = x[parameters(G)],
        target_parameters = ps
    )
    s = extract_samples(res, G; resample=true)

    ids_free = indexin(free(vars), variables(G)) # TODO: is there problem?
    free_samples = s[ids_free, :]
    return FixedFreeSamples(x_fixed, free_samples)
end
