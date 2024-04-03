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

A data type that represents a rational map 
``\\varphi \\colon X \\dashrightarrow \\mathbb{C}^m``.

# Constructors
```julia
ExpressionMap(
    domain::AbstractAlgebraicVariety,
    vars::AbstractArray,
    exprs::AbstractArray
)
ExpressionMap(
    domain::AbstractAlgebraicVariety,
    vars::AbstractArray
)
```

# Examples
```jldoctest
julia> @var R[1:3,1:3] t[1:3] E[1:3,1:3]
(Variable[R₁₋₁ R₁₋₂ R₁₋₃; R₂₋₁ R₂₋₂ R₂₋₃; R₃₋₁ R₃₋₂ R₃₋₃], Variable[t₁, t₂, t₃], Variable[E₁₋₁ E₁₋₂ E₁₋₃; E₂₋₁ E₂₋₂ E₂₋₃; E₃₋₁ E₃₋₂ E₃₋₃])

julia> X = AlgebraicVariety([R'*R-I, det(R)-1]; variables=[R, t]);

julia> tₓ = [0 -t[3] t[2]; t[3] 0 -t[1]; -t[2] t[1] 0]
3×3 Matrix{Expression}:
   0  -t₃   t₂
  t₃    0  -t₁
 -t₂   t₁    0

julia> φ = ExpressionMap(X, E, tₓ*R)
ExpressionMap: ℂ¹² ⊇ X - - > ℂ⁹
 domain:
  AlgebraicVariety X ⊂ ℂ¹²
   12 variables: R₁₋₁, R₂₋₁, R₃₋₁, R₁₋₂, R₂₋₂, R₃₋₂, R₁₋₃, R₂₋₃, R₃₋₃, t₁, t₂, t₃
   10 expressions:
    -1 + R₁₋₁^2 + R₂₋₁^2 + R₃₋₁^2
    R₁₋₁*R₁₋₂ + R₂₋₁*R₂₋₂ + R₃₋₂*R₃₋₁
    R₁₋₁*R₁₋₃ + R₂₋₁*R₂₋₃ + R₃₋₃*R₃₋₁
    R₁₋₁*R₁₋₂ + R₂₋₁*R₂₋₂ + R₃₋₂*R₃₋₁
    -1 + R₁₋₂^2 + R₂₋₂^2 + R₃₋₂^2
    R₁₋₂*R₁₋₃ + R₂₋₂*R₂₋₃ + R₃₋₃*R₃₋₂
    R₁₋₁*R₁₋₃ + R₂₋₁*R₂₋₃ + R₃₋₃*R₃₋₁
    R₁₋₂*R₁₋₃ + R₂₋₂*R₂₋₃ + R₃₋₃*R₃₋₂
    -1 + R₁₋₃^2 + R₂₋₃^2 + R₃₋₃^2
    -1 + (R₁₋₂*R₂₋₃ - R₁₋₃*R₂₋₂)*R₃₋₁ - (R₁₋₂*R₃₋₃ - R₁₋₃*R₃₋₂)*R₂₋₁ + (-R₃₋₂*R₂₋₃ + R₃₋₃*R₂₋₂)*R₁₋₁
 action:
  E₁₋₁ = t₂*R₃₋₁ - t₃*R₂₋₁
  E₂₋₁ = -t₁*R₃₋₁ + t₃*R₁₋₁
  E₃₋₁ = t₁*R₂₋₁ - t₂*R₁₋₁
  E₁₋₂ = t₂*R₃₋₂ - t₃*R₂₋₂
  E₂₋₂ = -t₁*R₃₋₂ + t₃*R₁₋₂
  E₃₋₂ = t₁*R₂₋₂ - t₂*R₁₋₂
  E₁₋₃ = t₂*R₃₋₃ - t₃*R₂₋₃
  E₂₋₃ = -t₁*R₃₋₃ + t₃*R₁₋₃
  E₃₋₃ = t₁*R₂₋₃ - t₂*R₁₋₃
```
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
    domain::AbstractAlgebraicVariety,
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

function ExpressionMap(
    domain::AbstractAlgebraicVariety,
    expr_vars::AbstractArray,
    exprs::AbstractArray
)
    expr_vars = Variable.(collect(flatten(expr_vars)))
    exprs = Expression.(collect(flatten(exprs)))
    return ExpressionMap(domain, expr_vars, exprs, Int[])
end

ExpressionMap(
    domain::AbstractAlgebraicVariety,
    domain_image_vars::Vector{Int}
) = ExpressionMap(domain, Variable[], Expression[], domain_image_vars)

function ExpressionMap(
    domain::AbstractAlgebraicVariety,
    domain_image_vars::AbstractArray
)
    domain_image_vars = Variable.(collect(flatten(domain_image_vars)))
    domain_image_ids = [findfirst(var->var==v, variables(domain)) for v in domain_image_vars]
    return ExpressionMap(domain, domain_image_ids)
end

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

expressions(φ::ExpressionMap) = φ.exprs # TODO: define all_expressions?
(φ::ExpressionMap)(x::AbstractVector) = evaluate(vcat(φ.exprs, domain_image_vars(φ)), domain_vars(φ)=>x)

expr_dict(φ::ExpressionMap) = Dict(zip(φ.expr_vars, φ.exprs))

function show_map_action(io::IO, φ::ExpressionMap, offset::String)
    if isempty(expressions(φ))
        print(io, "projection to ", join(domain_image_vars(φ), ", "))
    else
        for (j, (var, expr)) in enumerate(zip(expr_vars(φ), expressions(φ)))
            print(io, "$(offset)", var, " = ", expr)
            j < nexpr_vars(φ) && print(io, "\n")
        end
        if !isempty(domain_image_vars(φ))
            print(io, "\n")
            print(io, "$(offset)projection to ", join(domain_image_vars(φ), ", "))
        end
    end
end

function Base.show(io::IO, φ::ExpressionMap, offset::String)
    println(
        io,
        "$(offset)ExpressionMap: ",
        "ℂ$(superscript(ndomain_vars(φ))) ⊇ X",
        " - - > ",
        "ℂ$(superscript(nimage_vars(φ)))"
    )
    println(io, "$(offset) domain:")
    show(io, domain(φ), "$(offset)  ")
    print(io, "\n")
    println(io, "$(offset) action:")
    show_map_action(io, φ, "$(offset)  ")
end

Base.show(io::IO, φ::ExpressionMap) = show(io, φ, "")

"""
    domain_dimension(φ::ExpressionMap; kwargs...)

Computes the dimension of the domain of ``\\varphi``.
"""
domain_dimension(
    φ::ExpressionMap,
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
) = dimension(domain(φ), x; tols=tols)

generate_sample(domain::Vector{Variable}) = rand_unit(ComplexF64, length(domain))

"""
    image_dimension(φ::ExpressionMap; kwargs...)

Computes the dimension of the image of ``\\varphi``.
"""
function image_dimension(
    φ::ExpressionMap{T},
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
) where {T<:AbstractAlgebraicVariety}
    x = isnothing(x) ? generate_sample(domain(φ)) :  x
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

