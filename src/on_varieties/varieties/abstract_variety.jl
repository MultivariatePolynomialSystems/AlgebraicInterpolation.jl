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
    samples


"""
    abstract type AbstractAlgebraicVariety end
"""
abstract type AbstractAlgebraicVariety end

"""
    variables(X::AbstractAlgebraicVariety) -> Vector{Variable}

Return the variables of `X`.

# Examples
```julia-repl
julia> variables(X)
6-element Vector{Variable}:
 R₁₋₁
 R₂₋₁
 R₁₋₂
 R₂₋₂
   t₁
   t₂

julia> variables(Γ)
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
function variables(X::AbstractAlgebraicVariety)
    error("Not implemented")
end

"""
    nvariables(X::AbstractAlgebraicVariety) -> Int

Return the number of variables of `X`.

# Examples
```julia-repl
julia> nvariables(X)
6

julia> nvariables(Γ)
8
```
"""
nvariables(X::AbstractAlgebraicVariety) = length(variables(X))

"""
    expressions(X::AbstractAlgebraicVariety) -> Vector{Expression}

Return the expressions of `X`.

# Examples
```julia-repl
julia> expressions(X)
5-element Vector{Expression}:
       -1 + R₁₋₁^2 + R₂₋₁^2
      R₁₋₁*R₁₋₂ + R₂₋₁*R₂₋₂
      R₁₋₁*R₁₋₂ + R₂₋₁*R₂₋₂
       -1 + R₁₋₂^2 + R₂₋₂^2
 -1 + R₁₋₁*R₂₋₂ - R₁₋₂*R₂₋₁

julia> expressions(Γ)
7-element Vector{Expression}:
       -1 + R₁₋₁^2 + R₂₋₁^2
      R₁₋₁*R₁₋₂ + R₂₋₁*R₂₋₂
      R₁₋₁*R₁₋₂ + R₂₋₁*R₂₋₂
       -1 + R₁₋₂^2 + R₂₋₂^2
 -1 + R₁₋₁*R₂₋₂ - R₁₋₂*R₂₋₁
   s₁ - (t₁*R₁₋₁ + t₂*R₁₋₂)
   s₂ - (t₁*R₂₋₁ + t₂*R₂₋₂)
```
"""
function expressions(X::AbstractAlgebraicVariety)
    error("Not implemented")
end

"""
    nexpressions(X::AbstractAlgebraicVariety) -> Int

Return the number of expressions of `X`.

# Examples
```julia-repl
julia> nexpressions(X)
5

julia> nexpressions(Γ)
7
```
"""
nexpressions(X::AbstractAlgebraicVariety) = length(expressions(X))

"""
    generate_sample(X::AbstractAlgebraicVariety) -> Vector{ComplexF64}

Generate a sample from `X`.

# Examples
```julia-repl
julia> generate_sample(X)
6-element Vector{ComplexF64}:
 -0.8044679533846167 - 0.37355606563379196im
  0.7966605813908512 - 0.3772169611682977im
 -0.7966605813908512 + 0.3772169611682977im
 -0.8044679533846167 - 0.37355606563379196im
 -0.1653103992303188 - 0.014776408348405386im
  0.4869739603874435 - 0.8113780972782965im

julia> generate_sample(Γ)
8-element Vector{ComplexF64}:
  -0.9484588204432428 + 0.27190339385225026im
  -0.5995036004937357 - 0.4301711816163053im
   0.5995036004937357 + 0.43017118161630535im
  -0.9484588204432428 + 0.27190339385225026im
  -0.5638302344774876 - 0.3707771509690897im
  0.38660464552523993 + 0.4470517026250209im
   0.6750474426323179 + 0.6326747874590103im
 -0.30971285075613664 + 0.14593473983144623im
```
"""
function generate_sample(X::AbstractAlgebraicVariety)
    error("Not implemented")
end

"""
    jacobian(X::AbstractAlgebraicVariety) -> Matrix{Expression}

Return the symbolic jacobian of `X`.

# Examples
```julia-repl
julia> jacobian(X)
5×6 Matrix{Expression}:
 2*R₁₋₁  2*R₂₋₁       0       0  0  0
   R₁₋₂    R₂₋₂    R₁₋₁    R₂₋₁  0  0
   R₁₋₂    R₂₋₂    R₁₋₁    R₂₋₁  0  0
      0       0  2*R₁₋₂  2*R₂₋₂  0  0
   R₂₋₂   -R₁₋₂   -R₂₋₁    R₁₋₁  0  0

julia> jacobian(Γ)
7×8 Matrix{Expression}:
 2*R₁₋₁  2*R₂₋₁       0       0      0      0  0  0
   R₁₋₂    R₂₋₂    R₁₋₁    R₂₋₁      0      0  0  0
   R₁₋₂    R₂₋₂    R₁₋₁    R₂₋₁      0      0  0  0
      0       0  2*R₁₋₂  2*R₂₋₂      0      0  0  0
   R₂₋₂   -R₁₋₂   -R₂₋₁    R₁₋₁      0      0  0  0
    -t₁       0     -t₂       0  -R₁₋₁  -R₁₋₂  1  0
      0     -t₁       0     -t₂  -R₂₋₁  -R₂₋₂  0  1
```
"""
jacobian(X::AbstractAlgebraicVariety) = differentiate(expressions(X), variables(X))

"""
    jacobian(X::AbstractAlgebraicVariety, x::AbstractVector{<:Number}) -> Matrix{<:Number}

Return the jacobian of `X` evaluated at `x`.

# Examples
```julia-repl
julia> x = generate_sample(X);

julia> jacobian(X, x)
5×6 Matrix{ComplexF64}:
  -1.53374+0.72824im    -1.62756-0.686261im        0.0+0.0im             0.0+0.0im       0.0+0.0im  0.0+0.0im
  0.813778+0.343131im  -0.766869+0.36412im   -0.766869+0.36412im   -0.813778-0.343131im  0.0+0.0im  0.0+0.0im
  0.813778+0.343131im  -0.766869+0.36412im   -0.766869+0.36412im   -0.813778-0.343131im  0.0+0.0im  0.0+0.0im
       0.0+0.0im             0.0+0.0im         1.62756+0.686261im   -1.53374+0.72824im   0.0+0.0im  0.0+0.0im
 -0.766869+0.36412im   -0.813778-0.343131im   0.813778+0.343131im  -0.766869+0.36412im   0.0+0.0im  0.0+0.0im

julia> x = generate_sample(Γ);

julia> jacobian(Γ, x)
7×8 Matrix{ComplexF64}:
   1.98004+0.441998im   1.01158-0.865159im        0.0+0.0im       …        0.0+0.0im             0.0+0.0im       0.0+0.0im  0.0+0.0im
 -0.505789+0.43258im   0.990022+0.220999im   0.990022+0.220999im           0.0+0.0im             0.0+0.0im       0.0+0.0im  0.0+0.0im
 -0.505789+0.43258im   0.990022+0.220999im   0.990022+0.220999im           0.0+0.0im             0.0+0.0im       0.0+0.0im  0.0+0.0im
       0.0+0.0im            0.0+0.0im        -1.01158+0.865159im           0.0+0.0im             0.0+0.0im       0.0+0.0im  0.0+0.0im
  0.990022+0.220999im  0.505789-0.43258im   -0.505789+0.43258im            0.0+0.0im             0.0+0.0im       0.0+0.0im  0.0+0.0im
  -0.85827+0.362656im       0.0+0.0im       -0.625043-0.817264im  …  -0.990022-0.220999im   0.505789-0.43258im   1.0+0.0im  0.0+0.0im
       0.0+0.0im       -0.85827+0.362656im        0.0+0.0im          -0.505789+0.43258im   -0.990022-0.220999im  0.0+0.0im  1.0+0.0im
```
"""
jacobian(
    X::AbstractAlgebraicVariety,
    x::AbstractVector{<:Number}
) = evaluate(jacobian(X), variables(X) => x)

"""
    tangent_space(X::AbstractAlgebraicVariety, x::AbstractVector{<:Number}; <keyword arguments>) -> Matrix{<:Number}

Return the tangent space of `X` at `x`.

# Keyword arguments
- `tols::Tolerances=Tolerances()`: tolerances for numerical computations.

# Examples
```julia-repl
julia> x = generate_sample(X);

julia> tangent_space(X, x)
6×3 Matrix{ComplexF64}:
  0.326172+0.12926im   0.0-0.0im  0.0-0.0im
  0.597181-0.142401im  0.0-0.0im  0.0-0.0im
 -0.597181+0.142401im  0.0-0.0im  0.0-0.0im
  0.326172+0.12926im   0.0-0.0im  0.0-0.0im
      -0.0-0.0im       1.0-0.0im  0.0-0.0im
       0.0-0.0im       0.0-0.0im  1.0-0.0im

julia> x = generate_sample(Γ);

julia> tangent_space(Γ, x)
8×3 Matrix{ComplexF64}:
 0.0363119-0.652215im         0.0-0.0im              0.0-0.0im
 -0.242325-0.120728im         0.0-0.0im              0.0-0.0im
  0.242325+0.120728im         0.0-0.0im              0.0-0.0im
 0.0363119-0.652215im         0.0-0.0im              0.0-0.0im
      -0.0-0.0im              1.0-0.0im              0.0-0.0im
       0.0-0.0im              0.0-0.0im              1.0-0.0im
  0.421251+0.0796069im  -0.201566+0.390146im    -1.05694-0.0744038im
  0.458737-0.15848im      1.05694+0.0744038im  -0.201566+0.390146im
```
"""
tangent_space(
    X::AbstractAlgebraicVariety,
    x::AbstractVector{<:Number};
    tols::Tolerances=Tolerances(),
    var_ids::Union{Vector{Int}, Colon}=:
) = nullspace(jacobian(X, x); atol=tols.nullspace_atol)[var_ids, :]

"""
    dimension(X::AbstractAlgebraicVariety; <keyword arguments>) -> Int

Computes the dimension of `X`.

# Keyword arguments
- `sample::Union{AbstractVector{<:Number}, Nothing}=nothing`: point that belongs to `X`.
- `tols::Tolerances=Tolerances()`: tolerances for numerical computations.

# Examples
```julia-repl
julia> dimension(X)
3

julia> dimension(Γ)
3
```
"""
function dimension(
    X::AbstractAlgebraicVariety;
    sample::Union{AbstractVector{<:Number}, Nothing}=nothing,
    tols::Tolerances=Tolerances()
)
    x = isnothing(sample) ? generate_sample(X) : sample
    if !isnothing(x)
        # TODO: @assert norm(F(x)) < tols.common_tol
        J = jacobian(X, x)
        return nvariables(X) - rank(J; atol=tols.rank_atol)
    else
        error("Cannot generate a random sample of the algebraic variety, provide one!")
    end
end

"""
    finite_dominant_projection(X::AbstractAlgebraicVariety; <keyword arguments>) -> ExpressionMap

Returns a finite dominant projection from `X` to an affine space.

# Keyword arguments
- `sample::Union{AbstractVector{<:Number}, Nothing}=nothing`: point that belongs to `X`.
- `tols::Tolerances=Tolerances()`: tolerances for numerical computations.

# Examples
```julia-repl
julia> finite_dominant_projection(X)
ExpressionMap: ℂ⁶ ⊇ X - - > ℂ³
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
  projection to R₁₋₁, t₁, t₂

julia> finite_dominant_projection(Γ)
ExpressionMap: ℂ⁸ ⊇ X - - > ℂ³
 domain:
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
 action:
  projection to R₁₋₁, t₁, t₂
```
"""
function finite_dominant_projection( # TODO: check efficiency for MapGraph
    X::AbstractAlgebraicVariety;
    sample::Union{AbstractVector{<:Number}, Nothing}=nothing,
    tols::Tolerances=Tolerances()
)
    x = isnothing(sample) ? generate_sample(X) : sample
    proj_vars = Int[]
    dom_dim = dimension(X; sample=x, tols=tols)
    if dom_dim == 0
        error("Zero-dimensional variety doesn't have a dominant projection")
    end
    im_dim = 0
    for i in 1:nvariables(X)
        push!(proj_vars, i)
        φ = ExpressionMap(X, proj_vars)
        im_dim_curr = image_dimension(φ; domain_sample=x, tols=tols)
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
"""
function sample(X::AbstractAlgebraicVariety, vars::FixedFreeVariables; kwargs...)
    error("Not implemented")
end

"""
    sample!(X::AbstractAlgebraicVariety, vars::FixedFreeSamples; <keyword arguments>) -> FixedFreeSamples

Samples `X` in the given `vars` and updates `X` with these samples.
"""
function sample!(X::AbstractAlgebraicVariety, vars::FixedFreeVariables; kwargs...)
    error("Not implemented")
end

"""
    samples(X::AbstractAlgebraicVariety, vars::FixedFreeVariables) -> FixedFreeSamples

Returns the saved samples of `X` in the given `vars`.
"""
function samples(X::AbstractAlgebraicVariety, vars::FixedFreeVariables)
    error("Not implemented")
end