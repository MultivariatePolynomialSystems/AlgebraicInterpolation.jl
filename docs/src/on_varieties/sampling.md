# Sampling

In this Julia package we partially deal with polynomial and rational function defined on algebraic varieties. These varieties are defined by multivariate polynomial systems. We use [`HomotopyContinuation.jl`](https://www.juliahomotopycontinuation.org/) to sample these algebraic varieties. The samples are later used for the interpolation. 

An alternative would be that a user provides a method that will generate samples without using HomotopyContinuation. In practice it is usually possible to provide such a method, since we very often deal with structured polynomial systems.

```@docs
possible_to_sample
sample
sample!
```