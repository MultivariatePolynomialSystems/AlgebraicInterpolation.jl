module MultivariateInterpolation

import HomotopyContinuation
const HC = HomotopyContinuation
using HomotopyContinuation.ModelKit
export @var, Variable, Expression, System

using HomotopyContinuation: Result
using SparseArrays: SparseVector, SparseMatrixCSC, spzeros, AbstractSparseVector, findnz, sparse
using Combinatorics: partitions, multiset_permutations, combinations

using LinearAlgebra: nullspace, norm, rank, svdvals, svd, det, I
export det, I

using Base.Iterators: flatten


include("utils.jl")

include("interpolation_bases/abstract_basis.jl")
include("interpolation_bases/monomials.jl")
include("fixed_free.jl")

include("on_varieties/varieties/abstract_variety.jl")
include("on_varieties/varieties/affine_space.jl")
include("on_varieties/varieties/algebraic_variety.jl")

include("on_varieties/sampling.jl")

include("on_varieties/expression_map.jl")
include("on_varieties/varieties/map_graph.jl")

# include("on_varieties/scalings.jl")
include("on_varieties/interpolation.jl")

end
