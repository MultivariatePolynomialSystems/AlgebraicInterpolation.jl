# Interpolation basis

Abstract type for interpolation basis is given by

```@docs
AbstractInterpolationBasis
```

## Interface

In order to achieve the functionalities of this package that work on general [`AbstractInterpolationBasis`](@ref), one should implement the following interface:

```@docs
nelements
to_expressions
evaluate
```

## MonomialBasis

```@docs
MonomialBasis
```
