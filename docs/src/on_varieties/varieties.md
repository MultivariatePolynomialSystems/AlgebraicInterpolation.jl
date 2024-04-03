# Algebraic varieties

Abstract type for algebraic varieties is given by

```@docs
AbstractAlgebraicVariety
```

## AffineSpace

```@docs
AffineSpace
```

## AlgebraicVariety

```@docs
AlgebraicVariety
```

## MapGraph

```@docs
MapGraph
```

## Methods

We will show the basic functionality of an [`AbstractAlgebraicVariety`](@ref) on 2 concrete examples
defined in [`AlgebraicVariety`](@ref) and [`MapGraph`](@ref).

```@docs
variables(::AbstractAlgebraicVariety)
nvariables(::AbstractAlgebraicVariety)
expressions
nexpressions
generate_sample
jacobian
tangent_space
dimension
finite_dominant_projection
```