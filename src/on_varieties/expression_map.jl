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


"""
    ExpressionMap{T <: AbstractAlgebraicVariety}
"""
struct ExpressionMap{T<:AbstractAlgebraicVariety}
    domain::T
    expr_vars::Vector{Variable}
    exprs::Vector{Expression}
    domain_image_vars::Vector{Int}
    domain_nonimage_vars::Vector{Int}
    exprs_jacobian::Union{Matrix{Expression}, Nothing}
end

ExpressionMap(
    domain::AbstractDifferentiatedVariety,
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
    domain::AbstractAlgebraicVariety,
    expr_vars::AbstractArray{Variable},
    exprs::AbstractArray{Expression}
) = ExpressionMap(domain, expr_vars[:], exprs[:], Int[])

ExpressionMap(
    domain::AbstractAlgebraicVariety,
    domain_image_vars::Vector{Int}
) = ExpressionMap(domain, Variable[], Expression[], domain_image_vars)

ExpressionMap(
    domain::AbstractDifferentiatedVariety,
    domain_image_vars::AbstractArray{Variable}
) = ExpressionMap(domain, [findfirst(var->var==v, variables(domain)) for v in domain_image_vars[:]])

ExpressionMap(
    domain::Vector{Variable},
    domain_image_vars::AbstractArray{Variable}
) = ExpressionMap(domain, [findfirst(var->var==v, domain) for v in domain_image_vars[:]])

domain(φ::ExpressionMap) = φ.domain
domain_vars(φ::ExpressionMap) = variables(φ.domain)
ndomain_vars(φ::ExpressionMap) = nvariables(φ.domain)

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

"""
    domain_dimension(φ::ExpressionMap; kwargs...)

Computes the dimension of the domain of ``\\varphi``.
"""
domain_dimension(
    φ::ExpressionMap,
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
) = dimension(domain(φ), x; tols=tols)

find_sample(domain::Vector{Variable}) = rand_unit(ComplexF64, length(domain))

"""
    image_dimension(φ::ExpressionMap; kwargs...)

Computes the dimension of the image of ``\\varphi``.
"""
function image_dimension(
    φ::ExpressionMap{T},
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
) where {T<:AbstractAlgebraicVariety}
    x = isnothing(x) ? find_sample(domain(φ)) :  x
    if !isnothing(x)
        vars = domain_vars(φ)
        dφₓ = isnothing(φ.exprs_jacobian) ? nothing : evaluate(φ.exprs_jacobian[:,φ.domain_nonimage_vars], vars => x)
        if T isa AffineSpace
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


"""
    is_dominant(φ::ExpressionMap; kwargs...)

Returns `true` if ``\\varphi \\colon X \\dashrightarrow \\mathbb{C}^m`` is dominant, i.e. if

```math
\\mathrm{dim}(\\mathrm{im}(\\varphi)) = m.
```
"""
is_dominant(
    φ::ExpressionMap,
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
) = image_dimension(φ, x; tols=tols) == nimage_vars(φ)

