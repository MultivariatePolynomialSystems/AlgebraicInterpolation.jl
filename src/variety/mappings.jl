export ExpressionMap,
    domain,
    domain_vars,
    image_vars,
    domain_dimension,
    image_dimension,
    is_dominant,
    MapGraph,
    bottom_domain,
    extended_bottom_domain


MapDomain = Union{Vector{Variable}, AbstractDifferentiatedSystem}

struct ExpressionMap{T<:MapDomain}
    domain::T
    expr_vars::Vector{Variable}
    exprs::Vector{Expression}
    domain_image_vars::Vector{Int}
    domain_nonimage_vars::Vector{Int}
    exprs_jacobian::Union{Matrix{Expression}, Nothing}
end

ExpressionMap(
    domain::AbstractDifferentiatedSystem,
    expr_vars::Vector{Variable},
    exprs::Vector{Expression},
    domain_image_vars::Vector{Int}
) = ExpressionMap(
        domain,
        expr_vars,
        exprs,
        domain_image_vars,
        setdiff(1:nvariables(domain), domain_image_vars),
        isempty(exprs) ? nothing : differentiate(exprs, variables(domain))
    )

ExpressionMap(
    domain::Vector{Variable},
    expr_vars::Vector{Variable},
    exprs::Vector{Expression},
    domain_image_vars::Vector{Int}
) = ExpressionMap(
        domain,
        expr_vars,
        exprs,
        domain_image_vars,
        setdiff(1:length(domain), domain_image_vars),
        isempty(exprs) ? nothing : differentiate(exprs, domain)
    )

ExpressionMap(
    domain::System,
    expr_vars::Vector{Variable},
    exprs::Vector{Expression},
    domain_image_vars::Vector{Int}
) = ExpressionMap(DifferentiatedSystem(domain), expr_vars, exprs, domain_image_vars)

ExpressionMap(
    domain::Union{MapDomain, System},
    expr_vars::AbstractArray{Variable},
    exprs::AbstractArray{Expression}
) = ExpressionMap(domain, expr_vars[:], exprs[:], Int[])

ExpressionMap(
    domain::Union{MapDomain, System},
    domain_image_vars::Vector{Int}
) = ExpressionMap(domain, Variable[], Expression[], domain_image_vars)

ExpressionMap(
    domain::Union{AbstractDifferentiatedSystem, System},
    domain_image_vars::AbstractArray{Variable}
) = ExpressionMap(domain, [findfirst(var->var==v, variables(domain)) for v in domain_image_vars[:]])

ExpressionMap(
    domain::Vector{Variable},
    domain_image_vars::AbstractArray{Variable}
) = ExpressionMap(domain, [findfirst(var->var==v, domain) for v in domain_image_vars[:]])

domain(φ::ExpressionMap) = φ.domain
domain_vars(φ::ExpressionMap{Vector{Variable}}) = φ.domain
domain_vars(φ::ExpressionMap{<:AbstractDifferentiatedSystem}) = variables(φ.domain)
ndomain_vars(φ::ExpressionMap{Vector{Variable}}) = length(φ.domain)
ndomain_vars(φ::ExpressionMap{<:AbstractDifferentiatedSystem}) = nvariables(φ.domain)

expr_vars(φ::ExpressionMap) = φ.expr_vars
nexpr_vars(φ::ExpressionMap) = length(φ.expr_vars)
domain_image_vars(φ::ExpressionMap) = domain_vars(φ)[φ.domain_image_vars]
image_vars(φ::ExpressionMap) = vcat(φ.expr_vars, domain_image_vars(φ))
nimage_vars(φ::ExpressionMap) = length(φ.expr_vars) + length(φ.domain_image_vars)

# doesn't include repetitions
variables(φ::ExpressionMap) = vcat(domain_vars(φ), φ.expr_vars)
nvariables(φ::ExpressionMap) = ndomain_vars(φ) + length(φ.expr_vars)

HC.expressions(φ::ExpressionMap) = φ.exprs # TODO: define all_expressions?
(φ::ExpressionMap)(x::AbstractVector) = evaluate(vcat(φ.exprs, domain_image_vars(φ)), domain_vars(φ)=>x)

expr_dict(φ::ExpressionMap) = Dict(zip(φ.expr_vars, φ.exprs))

function Base.show(io::IO, φ::ExpressionMap)
    println(
        io,
        "ExpressionMap: ",
        "ℂ$(superscript(ndomain_vars(φ))) ⊇ X",
        " - - > ",
        "ℂ$(superscript(nimage_vars(φ)))"
    )
    println(io, " domain:") # TODO
    println(io, " action:")
    if isempty(φ.exprs)
        print(io, "  projection to ", join(domain_image_vars(φ), ", "))
    else
        # TODO
    end
end

domain_dimension(
    φ::ExpressionMap{Vector{Variable}},
    ::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
) = ndomain_vars(φ)

domain_dimension(
    φ::ExpressionMap{<:AbstractDifferentiatedSystem},
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
) = dimension(domain(φ), x; tols=tols)

find_sample(domain::Vector{Variable}) = rand_unit(ComplexF64, length(domain))

function image_dimension(
    φ::ExpressionMap{T},
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
) where {T<:MapDomain}
    x = isnothing(x) ? find_sample(domain(φ)) :  x
    if !isnothing(x)
        vars = domain_vars(φ)
        dφₓ = isnothing(φ.exprs_jacobian) ? nothing : evaluate(φ.exprs_jacobian[:,φ.domain_nonimage_vars], vars => x)
        if T isa Vector{Variable}
            if isnothing(dφₓ)
                return length(φ.domain_image_vars)
            else
                return rank(dφₓ; atol=tols.rank_atol) + length(φ.domain_image_vars) # TODO: consider rtol
            end
        else
            if isnothing(dφₓ)
                M = tangent_space(domain(φ), x; tols=tols, var_ids=φ.domain_image_vars)
            else
                TₓX = tangent_space(domain(φ), x; tols=tols)
                M = [dφₓ * TₓX[φ.domain_nonimage_vars,:]; TₓX[φ.domain_image_vars,:]]
            end
            return rank(M; atol=tols.rank_atol) # TODO: consider rtol
        end
    else
        error("Cannot generate a random sample in the domain, provide one!")
    end
end

is_dominant(
    φ::ExpressionMap,
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
) = image_dimension(φ, x; tols=tols) == nimage_vars(φ)

struct MapGraph{T<:ExpressionMap} <: AbstractSampledSystem
    map::T
    samples::Dict{FixedFreeVariables, FixedFreeSamples}
end

# TODO: check constructors
MapGraph(φ::ExpressionMap) = MapGraph(φ, Dict{FixedFreeVariables, FixedFreeSamples}())

domain(Γ::MapGraph) = domain(Γ.map)
domain_vars(Γ::MapGraph) = domain_vars(Γ.map)
ndomain_vars(Γ::MapGraph) = ndomain_vars(Γ.map)
image_vars(Γ::MapGraph) = image_vars(Γ.map)
nimage_vars(Γ::MapGraph) = nimage_vars(Γ.map)
expr_vars(Γ::MapGraph) = expr_vars(Γ.map)
nexpr_vars(Γ::MapGraph) = nexpr_vars(Γ.map)
variables(Γ::MapGraph) = variables(Γ.map)
nvariables(Γ::MapGraph) = nvariables(Γ.map)
HC.expressions(Γ::MapGraph) = vcat(expressions(domain(Γ)), expressions(Γ.map) .- expr_vars(Γ))
expr_dict(Γ::MapGraph) = expr_dict(Γ.map)

(Γ::MapGraph)(x::AbstractVector) = evaluate(expressions(Γ), variables(Γ) => x)

samples(
    Γ::MapGraph,
    vars::FixedFreeVariables
) = get(Γ.samples, vars, nothing)

function nsamples(Γ::MapGraph, vars::FixedFreeVariables)
    s = samples(Γ, vars)
    !isnothing(s) && return size(free(s), 2)
    return 0
end

function Base.show(io::IO, Γ::MapGraph)
    println(
        io,
        "MapGraph: Γ ⊂ ",
        "ℂ$(superscript(ndomain_vars(Γ)))",
        " × ",
        "ℂ$(superscript(nimage_vars(Γ)))"
    )
    println(io, " domain:")
    print(io, " map:")
    # TODO
end

function find_sample(Γ::MapGraph; kwargs...)
    s = samples(Γ, FixedFreeVariables(variables(Γ)))
    if !isnothing(s)
        return free(s)[:,rand(1:nsamples(s))]
    else
        s = find_sample(domain(Γ); kwargs...)
        return isnothing(s) ? nothing : vcat(s, Γ.map(s))
    end
end

function tangent_space(
    Γ::MapGraph,
    x::AbstractVector;
    tols::Tolerances=Tolerances(),
    var_ids::AbstractVector{Int}=1:nvariables(Γ)
)
    @assert norm(Γ(x)) < tols.common_tol
    if (length(var_ids) == nvariables(Γ))
        dom_ids = 1:ndomain_vars(Γ)
        expr_ids = 1:nexpr_vars(Γ)
    else
        tfs = var_ids .<= ndomain_vars(Γ)
        dom_ids = var_ids[tfs]
        expr_ids = var_ids[.!tfs] .- ndomain_vars(Γ)
    end
    x = x[1:ndomain_vars(Γ)]
    T = tangent_space(domain(Γ), x; tols=tols)
    isempty(expr_ids) && return T[dom_ids, :]
    B = evaluate(Γ.map.exprs_jacobian[expr_ids, :], domain_vars(Γ) => x)
    return vcat(T[dom_ids, :], B*T)
end

dimension(
    Γ::MapGraph,
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
) = dimension(domain(Γ), x[1:ndomain_vars(Γ)]; tols=tols)

# TODO: add keyword arg for same fixed values?
sample(
    domain::Vector{Variable},
    vars::FixedFreeVariables=FixedFreeVariables(domain),
    ::Union{AbstractVector, Nothing}=nothing;
    nsamples::Int=1,
    tols::Tolerances=Tolerances(),
    rand_method::Symbol=:rand_unit
) = FixedFreeSamples(
        eval(rand_method)(ComplexF64, nfixed(vars)),
        eval(rand_method)(ComplexF64, nfree(vars), nsamples)
    )

function bottom_domain(Γ::MapGraph)
    D = domain(Γ)
    D isa MapGraph && return bottom_domain(D)
    return D
end

function construct_expression(Γ::MapGraph, var::Variable, vars::Vector{Variable})
    expr = get(expr_dict(Γ), var, nothing)
    if isnothing(expr)
        return construct_expression(domain(Γ), var, vars)
    else
        expr_vars = HC.variables(expr)
        map_expr_vars = setdiff(expr_vars, vars)
        isempty(map_expr_vars) && return expr
        subs_exprs = [construct_expression(Γ, v, vars) for v in map_expr_vars]
        return subs(expr, map_expr_vars => subs_exprs)
    end
end

function var_expression(Γ::MapGraph, var::Variable, vars::Vector{Variable})
    if var ∉ variables(Γ) # TODO: remove?
        error("Variable doesn't exist")
    end
    var ∈ vars && return var
    return construct_expression(Γ, var, vars)
end

function extended_bottom_domain(Γ::MapGraph, vars::Vector{Variable})
    B = bottom_domain(Γ)
    bottom_vars = variables(B)
    map_vars = setdiff(vars, bottom_vars)
    var_exprs = [var - var_expression(Γ, var, bottom_vars) for var in map_vars]
    return System(vcat(expressions(B), var_exprs); variables=vcat(bottom_vars, map_vars))
end

function reconstruct_values(
    Γ::MapGraph,
    var::Variable,
    sampled_vars::FixedFreeVariables,
    sampled_values::FixedFreeSamples
)
    expr = construct_expression(Γ, var, variables(sampled_vars))
    expr_fixed = subs(expr, fixed(sampled_vars) => fixed(sampled_values))
    evals = zeros(ComplexF64, nsamples(sampled_values))
    for (j, vals) in enumerate(eachcol(free(sampled_values)))
        evals[j] = evaluate(expr_fixed, free(sampled_vars) => vals)
    end
    return evals
end

function reconstruct_values(
    Γ::MapGraph,
    vars::Vector{Variable},
    sampled_vars::FixedFreeVariables,
    sampled_values::FixedFreeSamples
)
    vals = zeros(ComplexF64, length(vars), nsamples(sampled_values))
    for (j, var) in enumerate(vars)
        vals[j, :] = reconstruct_values(Γ, var, sampled_vars, sampled_values)
    end
    return vals
end

# TODO: add keyword arg for same fixed values?
function sample(
    Γ::MapGraph,
    vars::FixedFreeVariables=FixedFreeVariables(variables(F)),
    x::Union{AbstractVector, Nothing}=nothing; # TODO: make keyword?
    nsamples::Int=1,
    tols::Tolerances=Tolerances(),
    rand_method::Symbol=:rand_unit
)
    x = isnothing(x) ? find_sample(Γ) : x
    if !possible_to_sample(Γ, vars, x; tols=tols, logging=false)
        error("Impossible to sample")
    end
    if !reasonable_to_sample(Γ, vars, x; tols=tols)
        error("Sampling is unnecessary")
    end

    # 1. add recursively reconstructed constraints for fixed_map_vars into bottom_domain (obtain F)
    F = extended_bottom_domain(Γ, fixed(vars))

    # 2. subs values for fixed(vars) into F and find finite dominant projection (obtain G)
    x = VariableEvaluation(variables(Γ), x)
    x_fixed = x[fixed(vars)]
    G = subs(F, fixed(vars) => x_fixed)
    G = finite_dominant_projection(G, x[variables(G)]; tols=tols)

    # 3. sample G; TODO: create separate method?
    ps = [eval(rand_method)(ComplexF64, nparameters(G)) for _ in 1:nsamples]
    res = HC.solve(
        G,
        [x[unknowns(G)]],
        start_parameters = x[parameters(G)],
        target_parameters = ps
    )
    samples_G = extract_samples(res, G; resample=true)

    # 4. extract values for bottom_vars in free(vars)
    free_bottom_vars = free(vars) ∩ variables(G)
    bottom_s = samples_G[indexin(free_bottom_vars, variables(G)), :]

    # 5. recursively compute values for map_vars in free(vars)
    free_map_vars = setdiff(free(vars), variables(G))
    map_s = reconstruct_values(
        Γ,
        free_map_vars,
        FixedFreeVariables(fixed(vars), variables(G)),
        FixedFreeSamples(x_fixed, samples_G)
    )

    # 6. create samples
    free_samples = zeros(ComplexF64, nfree(vars), nsamples)
    free_samples[indexin(free_bottom_vars, free(vars)), :] = bottom_s
    free_samples[indexin(free_map_vars, free(vars)), :] = map_s
    return FixedFreeSamples(x_fixed, free_samples)
end

# TODO: add keyword arg for same fixed values?
function sample!(
    Γ::MapGraph,
    vars::FixedFreeVariables=FixedFreeVariables(variables(F)),
    x::Union{AbstractVector, Nothing}=nothing; # TODO: make keyword?
    nsamples::Int=1,
    tols::Tolerances=Tolerances(),
    rand_method::Symbol=:rand_unit
)
    s = sample(Γ, vars, x; nsamples=nsamples, tols=tols, rand_method=rand_method)
    Γ.samples[vars] = s
end