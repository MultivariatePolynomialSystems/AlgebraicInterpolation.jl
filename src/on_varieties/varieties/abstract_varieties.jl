export AbstractAlgebraicVariety,
    AbstractDifferentiatedVariety,
    AbstractSampledVariety


"""
    AbstractAlgebraicVariety

*Interface*:

```julia
variables(X::AbstractAlgebraicVariety) -> AbstractVector{Variable}
```
```julia
nvariables(X::AbstractAlgebraicVariety) -> Integer
```
```julia
expressions(X::AbstractAlgebraicVariety) -> AbstractVector{Expression}
```
```julia
nexpressions(X::AbstractAlgebraicVariety) -> Integer
```
```julia
generate_sample(X::AbstractAlgebraicVariety) -> AbstractVector{ComplexF64}
```

"""
abstract type AbstractAlgebraicVariety end

"""
    AbstractDifferentiatedVariety <: AbstractAlgebraicVariety

*Interface*:

```julia
full_jacobian(X::AbstractDifferentiatedVariety) -> AbstractMatrix{Expression}
```
```julia
full_jacobian(
    X::AbstractDifferentiatedVariety,
    x::AbstractVector{ComplexF64}
) -> AbstractMatrix{ComplexF64}
```
```julia
tangent_space(
    X::AbstractDifferentiatedVariety,
    x::AbstractVector{ComplexF64}
) -> AbstractMatrix{ComplexF64}
```
```julia
dimension(X::AbstractDifferentiatedVariety) -> Integer
```
"""
abstract type AbstractDifferentiatedVariety <: AbstractAlgebraicVariety end

"""
    AbstractSampledVariety <: AbstractDifferentiatedVariety

Interface:
"""
abstract type AbstractSampledVariety <: AbstractDifferentiatedVariety end