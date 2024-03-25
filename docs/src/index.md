# Introduction
MultivariateInterpolation.jl is a Julia package that provides methods for interpolating sparse and dense multivariate polynomial and rational functions from samples taken either from an algebraic variety or an affine space.

## Quick start

```@repl
using MultivariateInterpolation
@var R[1:3,1:3] t[1:3] E[1:3,1:3]
```

## Contents