module MultivariateInterpolation

import HomotopyContinuation
const HC = HomotopyContinuation
using HomotopyContinuation.ModelKit
export @var, Variable, Expression, System

using HomotopyContinuation: Result
using SparseArrays: SparseVector, SparseMatrixCSC, spzeros, AbstractSparseVector, findnz, sparse
using Combinatorics: partitions, multiset_permutations, combinations
using LinearAlgebra: nullspace, norm, rank, svdvals, svd
using Base.Iterators: flatten


include("utils.jl")
include("monomials.jl")
include("variety/sampling.jl")
include("variety/mappings.jl")
include("variety/scalings.jl")
include("variety/interpolation.jl")

end
