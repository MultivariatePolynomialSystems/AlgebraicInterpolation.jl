module MultivariateInterpolation

import HomotopyContinuation
const HC = HomotopyContinuation
using HomotopyContinuation.ModelKit
export @var, Variable, Expression, System, subs
using HomotopyContinuation: Result

import SparseArrays
using SparseArrays: SparseVector, SparseMatrixCSC, spzeros, AbstractSparseVector, findnz, sparse
using Combinatorics: partitions, multiset_permutations, combinations, with_replacement_combinations

using LinearAlgebra: nullspace, norm, rank, svdvals, svd, det, dot, I
export det, I

using Base.Iterators: flatten


include("interpolation_bases/abstract_basis.jl")
include("interpolation_bases/monomials.jl")
include("utils.jl")
include("fixed_free.jl")

include("on_varieties/varieties/abstract_variety.jl")
include("on_varieties/varieties/affine_space.jl")
include("on_varieties/varieties/algebraic_variety.jl")

include("on_varieties/sampling.jl")

include("on_varieties/expression_map.jl")
include("on_varieties/varieties/map_graph.jl")

include("on_varieties/lie-symmetries.jl")
include("on_varieties/scalings.jl")
include("on_varieties/interpolation.jl")

end
