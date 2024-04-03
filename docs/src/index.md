# Introduction
MultivariateInterpolation.jl is a Julia package that provides methods for interpolating sparse and dense multivariate polynomial and rational functions from samples taken either from an algebraic variety or an affine space.

## Quick start

```@repl
using MultivariateInterpolation
@var R[1:3,1:3] t[1:3] E[1:3,1:3]
X = AlgebraicVariety([R'*R-I, det(R)-1]; variables=[R, t]);
tₓ = [0 -t[3] t[2]; t[3] 0 -t[1]; -t[2] t[1] 0] # skew-symmetric matrix
φ = ExpressionMap(X, expressions=Pair(E, tₓ*R)); # map to essential matrices
image_dimension(φ)
Γ = MapGraph(φ)
```

## Contents