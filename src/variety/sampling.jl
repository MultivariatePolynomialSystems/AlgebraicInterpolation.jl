export FixedFreeVariables,
    FixedFreeSamples,
    fixed,
    free,
    nfixed,
    nfree,
    AbstractSampledSystem,
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

struct FixedFreeVariables
    fixed::Vector{Variable}
    free::Vector{Variable}

    function FixedFreeVariables(fixed::Vector{Variable}, free::Vector{Variable})
        if !isempty(fixed ∩ free)
            error("Nontrivial intersection of fixed and free variables")
        end
        if isempty(free)
            error("Array of free variables must be nonempty")
        end
        return new(fixed, free)
    end
end

function FixedFreeVariables(
    fixed::Union{Variable, AbstractArray},
    free::Union{Variable, AbstractArray}
)
    return FixedFreeVariables(collect(flatten(fixed)), collect(flatten(free)))
end
FixedFreeVariables(free::Union{Variable, AbstractArray}) = FixedFreeVariables(Variable[], free)

function Base.:(==)(v₁::FixedFreeVariables, v₂::FixedFreeVariables)
    return v₁.fixed == v₂.fixed && v₁.free == v₂.free
end

Base.hash(
    vars::FixedFreeVariables,
    u::UInt64
) = hash(vars.fixed, hash(vars.free, u))

fixed(vars::FixedFreeVariables) = vars.fixed
free(vars::FixedFreeVariables) = vars.free
nfixed(vars::FixedFreeVariables) = length(vars.fixed)
nfree(vars::FixedFreeVariables) = length(vars.free)
variables(vars::FixedFreeVariables) = vcat(vars.fixed, vars.free)

function Base.show(io::IO, vars::FixedFreeVariables)
    println(
        io,
        "FixedFreeVariables: ",
        "$(nfixed(vars)) fixed, ",
        "$(nfree(vars)) free"
    )
    println(io, " fixed: ", join(fixed(vars), ", "))
    print(io, " free: ", join(free(vars), ", "))
end

# TODO: remove?
struct FixedFreeIndices
    fixed::Vector{Int}
    free::Vector{Int}
end

fixed(ids::FixedFreeIndices) = ids.fixed
free(ids::FixedFreeIndices) = ids.free
nfixed(ids::FixedFreeIndices) = length(ids.fixed)
nfree(ids::FixedFreeIndices) = length(ids.free)

FixedFreeIndices(
    vars::FixedFreeVariables,
    all_vars::Vector{Variable}
) = FixedFreeVariables(indexin(fixed(vars), all_vars), indexin(free(vars), all_vars))

struct FixedFreeSamples
    fixed::Vector{ComplexF64}
    free::Matrix{ComplexF64}
end

FixedFreeSamples(
    nfixed::Int,
    nfree::Int,
    nsamples::Int
) = FixedFreeSamples(zeros(ComplexF64, nfixed), zeros(ComplexF64, nfree, nsamples))

fixed(s::FixedFreeSamples) = s.fixed
free(s::FixedFreeSamples) = s.free
nsamples(s::FixedFreeSamples) = size(s.free, 2)

function Base.hcat(s₁::FixedFreeSamples, s₂::FixedFreeSamples; tol::Real=1e-10) # TODO: make tol optional?
    @assert norm(fixed(s₁)-fixed(s₂)) < tol
    return FixedFreeSamples(fixed(s₁), hcat(free(s₁), free(s₂)))
end

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

abstract type AbstractDifferentiatedSystem end

struct DifferentiatedSystem <: AbstractDifferentiatedSystem
    system::System
    full_jacobian::Matrix{Expression}
end

DifferentiatedSystem(F::System) = DifferentiatedSystem(F, full_jacobian(F))
HC.System(F::DifferentiatedSystem) = F.system

(F::DifferentiatedSystem)(x::AbstractVector) = System(F)(x)
find_sample(F::DifferentiatedSystem; kwargs...) = find_sample(System(F); kwargs...)

full_jacobian(F::DifferentiatedSystem) = F.full_jacobian
full_jacobian(
    F::AbstractDifferentiatedSystem,
    x::AbstractVector
) = evaluate(full_jacobian(F), variables(F) => x)

tangent_space(
    F::AbstractDifferentiatedSystem,
    x::AbstractVector;
    tols::Tolerances=Tolerances(),
    var_ids::Union{Vector{Int}, Colon}=:
) = nullspace(full_jacobian(F, x); atol=tols.nullspace_atol)[var_ids, :]

function dimension(
    F::AbstractDifferentiatedSystem,
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
)
    x = isnothing(x) ? find_sample(F) : x
    if !isnothing(x)
        @assert norm(F(x)) < tols.common_tol
        J = full_jacobian(F, x)
        return nvariables(F) - rank(J; atol=tols.rank_atol)
    else
        error("Cannot generate a random sample of the system, provide one!")
    end
end

dimension(
    F::System,
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
) = dimension(DifferentiatedSystem(F), x; tols=tols)

function finite_dominant_projection(
    F::AbstractDifferentiatedSystem,
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
    φ = finite_dominant_projection(DifferentiatedSystem(F), x; tols=tols)
    params = image_vars(φ)
    unknws = setdiff(variables(F), params)
    return System(expressions(F); variables=unknws, parameters=params)
end

function reasonable_to_sample(
    F::AbstractDifferentiatedSystem,
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
    possible_to_sample(
        F::AbstractDifferentiatedSystem,
        vars::FixedFreeVariables,
        x::Union{AbstractVector, Nothing}=nothing;
        tols::Tolerances=Tolerances(),
        logging::Bool=true
    )

Returns false if there exists a free variable that is determined in finite many ways by all fixed variables. 
Throws a WARNING if there are no constraints in free variables after fixing all the fixed ones.
"""
function possible_to_sample(
    F::AbstractDifferentiatedSystem,
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
    F::AbstractDifferentiatedSystem,
    vars::AbstractArray{Variable},
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
)
    @assert !isempty(vars)
    @assert vars ⊆ variables(F)
    for var in vars

    end
end

abstract type AbstractSampledSystem <: AbstractDifferentiatedSystem end

struct SampledSystem <: AbstractSampledSystem
    system::DifferentiatedSystem
    samples::Dict{FixedFreeVariables, FixedFreeSamples}
    sample_generator::Union{Function, Nothing}
end

SampledSystem(
    F::System,
    sample_generator::Union{Function, Nothing}
) = SampledSystem(DifferentiatedSystem(F), Dict(), sample_generator)
SampledSystem(F::System) = SampledSystem(F, nothing)

HC.System(F::SampledSystem) = System(F.system)

unknowns(F::Union{DifferentiatedSystem, SampledSystem}) = unknowns(System(F))
HC.parameters(F::Union{DifferentiatedSystem, SampledSystem}) = parameters(System(F))
variables(F::Union{DifferentiatedSystem, SampledSystem}) = variables(System(F))
nunknowns(F::Union{DifferentiatedSystem, SampledSystem}) = nunknowns(System(F))
HC.nparameters(F::Union{DifferentiatedSystem, SampledSystem}) = nparameters(System(F))
nvariables(F::Union{DifferentiatedSystem, SampledSystem}) = nvariables(System(F))
HC.expressions(F::Union{DifferentiatedSystem, SampledSystem}) = expressions(System(F))

HC.subs(
    F::System,
    substitutions::Pair...
) = System(
        subs(expressions(F), substitutions...);
        variables=setdiff(variables(F), vcat(first.(substitutions)...))
    )

(F::SampledSystem)(x::AbstractVector) = System(F)(x)

samples(
    F::SampledSystem,
    vars::FixedFreeVariables
) = get(F.samples, vars, nothing)

function nsamples(F::SampledSystem, vars::FixedFreeVariables)
    s = samples(F, vars)
    !isnothing(s) && return size(free(s), 2)
    return 0
end

add_samples!(
    F::SampledSystem,
    vars::FixedFreeVariables,
    s::FixedFreeSamples;
    tol::Real=1e-10
) = F.samples[vars] = isnothing(samples(F, vars)) ? s : hcat(samples(F, vars), s; tol=tol)

function find_sample(F::SampledSystem; kwargs...)
    vars = FixedFreeVariables(variables(F))
    s = samples(F, vars)
    if isnothing(s) && !isnothing(F.sample_generator)
        s = F.sample_generator(vars=vars, nsamples=1)
        return isnothing(s) ? find_sample(System(F); kwargs...) : free(s)[:]
    end
    return isnothing(s) ? find_sample(System(F); kwargs...) : free(s)[:,rand(1:nsamples(s))]
end

full_jacobian(F::SampledSystem) = full_jacobian(F.system)

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
function sample(
    F::AbstractDifferentiatedSystem,
    vars::FixedFreeVariables=FixedFreeVariables(variables(F)),
    x::Union{AbstractVector, Nothing}=nothing;
    nsamples::Int=1,
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

# TODO: add keyword arg for same fixed values?
function sample!(
    F::SampledSystem,
    vars::FixedFreeVariables=FixedFreeVariables(variables(F)),
    x::Union{AbstractVector, Nothing}=nothing;
    nsamples::Int=1,
    tols::Tolerances=Tolerances(),
    rand_method::Symbol=:rand_unit
)
    s = sample(F, vars, x; nsamples=nsamples, tols=tols, rand_method=rand_method)
    F.samples[vars] = s
end