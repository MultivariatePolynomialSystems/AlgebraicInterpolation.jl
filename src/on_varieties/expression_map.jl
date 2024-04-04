export ExpressionMap,
    domain,
    domain_vars,
    ndomain_vars,
    expression_vars,
    nexpression_vars,
    projection_vars,
    nprojection_vars,
    image_vars,
    nimage_vars,
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

# Constructor
```julia
ExpressionMap(
    domain::AbstractAlgebraicVariety;
    expressions::Pair{<:AbstractArray, <:AbstractArray}=Pair(Variable[], Expression[]),
    projection::AbstractArray=Variable[]
)
```

# Examples
```jldoctest
julia> @var R[1:2,1:2] t[1:2] s[1:2]
(Variable[R₁₋₁ R₁₋₂; R₂₋₁ R₂₋₂], Variable[t₁, t₂], Variable[s₁, s₂])

julia> X = AlgebraicVariety([R'*R-I, det(R)-1]; variables=[R, t]);

julia> φ = ExpressionMap(X; expressions=Pair(s, R*t), projection=t)
ExpressionMap: ℂ⁶ ⊇ X - - > ℂ⁴
 domain:
  AlgebraicVariety X ⊂ ℂ⁶
   6 variables: R₁₋₁, R₂₋₁, R₁₋₂, R₂₋₂, t₁, t₂
   5 expressions: 
    -1 + R₁₋₁^2 + R₂₋₁^2
    R₁₋₁*R₁₋₂ + R₂₋₁*R₂₋₂
    R₁₋₁*R₁₋₂ + R₂₋₁*R₂₋₂
    -1 + R₁₋₂^2 + R₂₋₂^2
    -1 + R₁₋₁*R₂₋₂ - R₁₋₂*R₂₋₁
 action:
  s₁ = t₁*R₁₋₁ + t₂*R₁₋₂
  s₂ = t₁*R₂₋₁ + t₂*R₂₋₂
  projection to t₁, t₂
```
"""
struct ExpressionMap{T<:AbstractAlgebraicVariety}
    domain::T
    expr_vars::Vector{Variable}
    exprs::Vector{Expression}
    proj_vars::Vector{Int}
    nonproj_vars::Vector{Int}
    exprs_jacobian::Union{Matrix{Expression}, Nothing}
end

ExpressionMap(
    domain::AbstractAlgebraicVariety,
    expr_vars::Vector{Variable},
    exprs::Vector{Expression},
    proj_vars::Vector{Int}
) = ExpressionMap(
        domain,
        expr_vars,
        exprs,
        proj_vars,
        setdiff(1:nvariables(domain), proj_vars),
        isempty(exprs) ? nothing : differentiate(exprs, variables(domain))
    )

function ExpressionMap(
    domain::AbstractAlgebraicVariety;
    expressions::Pair{<:AbstractArray, <:AbstractArray}=Pair(Variable[], Expression[]),
    projection::AbstractArray=Variable[]
)
    expr_vars = Variable.(collect(flatten(first(expressions))))
    exprs = Expression.(collect(flatten(last(expressions))))
    dom_im_vars = Variable.(collect(flatten(projection)))
    dom_im_ids = Int.([findfirst(var->var==v, variables(domain)) for v in dom_im_vars])
    return ExpressionMap(domain, expr_vars, exprs, dom_im_ids)
end

"""
    domain(φ::ExpressionMap{T<:AbstractAlgebraicVariety}) -> T

Return the domain of `φ`.

# Examples
```julia-repl
julia> domain(φ)
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
domain(φ::ExpressionMap) = φ.domain

"""
    domain_vars(φ::ExpressionMap) -> Vector{Variable}

Return the variables of the domain of `φ`.

# Examples
```julia-repl
julia> domain_vars(φ)
6-element Vector{Variable}:
 R₁₋₁
 R₂₋₁
 R₁₋₂
 R₂₋₂
   t₁
   t₂
```
"""
domain_vars(φ::ExpressionMap) = variables(φ.domain)

"""
    ndomain_vars(φ::ExpressionMap) -> Int

Return the number of variables in the domain of `φ`.

# Examples
```julia-repl
julia> ndomain_vars(φ)
6
```
"""
ndomain_vars(φ::ExpressionMap) = nvariables(φ.domain)

"""
    expression_vars(φ::ExpressionMap) -> Vector{Variable}

Return the expression variables of `φ`.

# Examples
```julia-repl
julia> expression_vars(φ)
2-element Vector{Variable}:
 s₁
 s₂
```
"""
expression_vars(φ::ExpressionMap) = φ.expr_vars

"""
    nexpression_vars(φ::ExpressionMap) -> Int

Return the number of expression variables of `φ`.

# Examples
```julia-repl
julia> nexpression_vars(φ)
2
```
"""
nexpression_vars(φ::ExpressionMap) = length(φ.expr_vars)

"""
    projection_vars(φ::ExpressionMap) -> Vector{Variable}

Return the projection variables of `φ`.

# Examples
```julia-repl
julia> projection_vars(φ)
2-element Vector{Variable}:
 t₁
 t₂
```
"""
projection_vars(φ::ExpressionMap) = domain_vars(φ)[φ.proj_vars]

"""
    nprojection_vars(φ::ExpressionMap) -> Int

Return the number of projection variables of `φ`.

# Examples
```julia-repl
julia> nprojection_vars(φ)
2
```
"""
nprojection_vars(φ::ExpressionMap) = length(φ.proj_vars)

"""
    image_vars(φ::ExpressionMap) -> Vector{Variable}

Return the concatenated vector of expression and projection variables of `φ`.

# Examples
```julia-repl
julia> image_vars(φ)
4-element Vector{Variable}:
 s₁
 s₂
 t₁
 t₂
```
"""
image_vars(φ::ExpressionMap) = vcat(expression_vars(φ), projection_vars(φ))

"""
    nimage_vars(φ::ExpressionMap) -> Int

Return the number of image variables of `φ`.

# Examples
```julia-repl
julia> nimage_vars(φ)
4
```
"""
nimage_vars(φ::ExpressionMap) = nexpression_vars(φ) + nprojection_vars(φ)

"""
    variables(φ::ExpressionMap) -> Vector{Variable}

Return the concatenated vector of domain and expression variables of `φ`.

# Examples
```julia-repl
julia> variables(φ)
8-element Vector{Variable}:
 R₁₋₁
 R₂₋₁
 R₁₋₂
 R₂₋₂
   t₁
   t₂
   s₁
   s₂
```
"""
variables(φ::ExpressionMap) = vcat(domain_vars(φ), expression_vars(φ))

"""
    nvariables(φ::ExpressionMap) -> Int

Return the number of variables of `φ`.

# Examples
```julia-repl
julia> nvariables(φ)
8
```
"""
nvariables(φ::ExpressionMap) = ndomain_vars(φ) + nexpression_vars(φ)

"""
    expressions(φ::ExpressionMap) -> Vector{Expression}

Return the expressions that define `φ`. Doesn't include the projection variables.

# Examples
```julia-repl
julia> expressions(φ)
2-element Vector{Expression}:
 t₁*R₁₋₁ + t₂*R₁₋₂
 t₁*R₂₋₁ + t₂*R₂₋₂
```
"""
expressions(φ::ExpressionMap) = φ.exprs # TODO: define all_expressions?
(φ::ExpressionMap)(x::AbstractVector) = evaluate(vcat(φ.exprs, projection_vars(φ)), domain_vars(φ)=>x)

expr_dict(φ::ExpressionMap) = Dict(zip(φ.expr_vars, φ.exprs))

function show_map_action(io::IO, φ::ExpressionMap, offset::String)
    if isempty(expressions(φ))
        print(io, "$(offset)projection to ", join(projection_vars(φ), ", "))
    else
        for (j, (var, expr)) in enumerate(zip(expression_vars(φ), expressions(φ)))
            print(io, "$(offset)", var, " = ", expr)
            j < nexpression_vars(φ) && print(io, "\n")
        end
        if !isempty(projection_vars(φ))
            print(io, "\n")
            print(io, "$(offset)projection to ", join(projection_vars(φ), ", "))
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
    domain_dimension(φ::ExpressionMap; <keyword_arguments>) -> Int

Compute the dimension of the domain of ``\\varphi``.

# Examples
```julia-repl
julia> domain_dimension(φ)
3
```
"""
domain_dimension(
    φ::ExpressionMap;
    domain_sample::Union{AbstractVector{<:Number}, Nothing}=nothing,
    tols::Tolerances=Tolerances()
) = dimension(domain(φ); sample=domain_sample, tols=tols)

"""
    image_dimension(φ::ExpressionMap; <keyword_arguments>) -> Int

Compute the dimension of the image of ``\\varphi``.

# Examples
```julia-repl
julia> image_dimension(φ)
3
```
"""
function image_dimension(
    φ::ExpressionMap{T};
    domain_sample::Union{AbstractVector{<:Number}, Nothing}=nothing,
    tols::Tolerances=Tolerances()
) where {T<:AbstractAlgebraicVariety}
    x = isnothing(domain_sample) ? generate_sample(domain(φ)) :  domain_sample
    if !isnothing(x)
        vars = domain_vars(φ)
        dφₓ = isnothing(φ.exprs_jacobian) ? nothing : evaluate(φ.exprs_jacobian[:,φ.nonproj_vars], vars => x)
        if T isa AffineSpace
            if isnothing(dφₓ)
                return nprojection_vars(φ)
            else
                return rank(dφₓ; atol=tols.rank_atol) + nprojection_vars(φ) # TODO: consider rtol
            end
        else
            if isnothing(dφₓ)
                M = tangent_space(domain(φ), x; tols=tols, var_ids=φ.proj_vars)
            else
                TₓX = tangent_space(domain(φ), x; tols=tols)
                M = [dφₓ * TₓX[φ.nonproj_vars,:]; TₓX[φ.proj_vars,:]]
            end
            return rank(M; atol=tols.rank_atol) # TODO: consider rtol
        end
    else
        error("Cannot generate a random sample in the domain, provide one!")
    end
end


"""
    is_dominant(φ::ExpressionMap; <keyword_arguments>)

Return `true` if ``\\varphi \\colon X \\dashrightarrow \\mathbb{C}^m`` is dominant, i.e. if

```math
\\mathrm{dim}(\\mathrm{im}(\\varphi)) = m.
```

# Examples
```julia-repl
julia> is_dominant(φ)
false
```
"""
is_dominant(
    φ::ExpressionMap;
    domain_sample::Union{AbstractVector{<:Number}, Nothing}=nothing,
    tols::Tolerances=Tolerances()
) = image_dimension(φ; domain_sample=domain_sample, tols=tols) == nimage_vars(φ)

