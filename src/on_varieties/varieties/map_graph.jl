export MapGraph


"""
    MapGraph{T<:ExpressionMap} <: AbstractAlgebraicVariety

An [`AbstractAlgebraicVariety`](@ref) that represents a graph 
``\\Gamma = \\{(x, \\varphi(x)) \\;|\\; x \\in X\\}`` of an 
[`ExpressionMap`](@ref) ``\\varphi \\colon X \\dashrightarrow \\mathbb{C}^m``.

# Constructors
```julia
MapGraph(φ::ExpressionMap)
```

# Examples
```jldoctest
julia> @var R[1:2,1:2] t[1:2] s[1:2]
(Variable[R₁₋₁ R₁₋₂; R₂₋₁ R₂₋₂], Variable[t₁, t₂], Variable[s₁, s₂])

julia> X = AlgebraicVariety([R'*R-I, det(R)-1]; variables=[R, t]);

julia> Γ = MapGraph(ExpressionMap(X, s, R*t))
MapGraph Γ ⊂ ℂ⁶ × ℂ²
 domain part:
  AlgebraicVariety X ⊂ ℂ⁶
   6 variables: R₁₋₁, R₂₋₁, R₁₋₂, R₂₋₂, t₁, t₂
   5 expressions: 
    -1 + R₁₋₁^2 + R₂₋₁^2
    R₁₋₁*R₁₋₂ + R₂₋₁*R₂₋₂
    R₁₋₁*R₁₋₂ + R₂₋₁*R₂₋₂
    -1 + R₁₋₂^2 + R₂₋₂^2
    -1 + R₁₋₁*R₂₋₂ - R₁₋₂*R₂₋₁
 image part:
  s₁ = t₁*R₁₋₁ + t₂*R₁₋₂
  s₂ = t₁*R₂₋₁ + t₂*R₂₋₂
```
"""
struct MapGraph{T<:ExpressionMap} <: AbstractAlgebraicVariety
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
expressions(Γ::MapGraph) = vcat(expressions(domain(Γ)), expr_vars(Γ) .- expressions(Γ.map))
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

function Base.show(io::IO, Γ::MapGraph, offset::String)
    println(
        io,
        "$(offset)MapGraph Γ ⊂ ",
        "ℂ$(superscript(ndomain_vars(Γ)))",
        " × ",
        "ℂ$(superscript(nimage_vars(Γ)))"
    )
    println(io, "$(offset) domain part:")
    show(io, domain(Γ), "$(offset)  ")
    print(io, "\n")
    println(io, "$(offset) image part:")
    show_map_action(io, Γ.map, "$(offset)  ")
end

Base.show(io::IO, Γ::MapGraph) = show(io, Γ, "")

function generate_sample(Γ::MapGraph; kwargs...)
    s = samples(Γ, FixedFreeVariables(variables(Γ)))
    if !isnothing(s)
        return free(s)[:,rand(1:nsamples(s))]
    else
        s = generate_sample(domain(Γ); kwargs...)
        return isnothing(s) ? nothing : vcat(s, Γ.map(s))
    end
end

function tangent_space(
    Γ::MapGraph,
    x::AbstractVector{<:Number};
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

function dimension(
    Γ::MapGraph;
    sample::Union{AbstractVector{<:Number}, Nothing}=nothing,
    tols::Tolerances=Tolerances()
)
    isnothing(sample) && return dimension(domain(Γ); tols=tols)
    return dimension(domain(Γ); sample=sample[1:ndomain_vars(Γ)], tols=tols)
end

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
    x = isnothing(x) ? generate_sample(Γ) : x
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