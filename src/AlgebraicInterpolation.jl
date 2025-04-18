module AlgebraicInterpolation

import HomotopyContinuation
const HC = HomotopyContinuation
using HomotopyContinuation.ModelKit
export @var, Variable, Expression, System, subs, expand
using HomotopyContinuation: Result

import SparseArrays
using SparseArrays: SparseVector, SparseMatrixCSC, spzeros, AbstractSparseVector, findnz, sparse
using Combinatorics: partitions, multiset_permutations, combinations, with_replacement_combinations

import LinearAlgebra
using LinearAlgebra: nullspace, norm, rank, svdvals, svd, det, dot, I, Diagonal
export det, I, nullspace

using Base.Iterators: flatten
using Base: reduce


# 1. Basics of Interpolation
include("interpolation_bases/abstract_basis.jl")
include("interpolation_bases/monomials.jl")
include("utils.jl")
include("fixed_free.jl")

# 2. Algebraic Varieties
include("on_varieties/varieties/abstract_variety.jl")
include("on_varieties/varieties/affine_space.jl")
include("on_varieties/varieties/algebraic_variety.jl")

# 3. Sampling
include("on_varieties/sampling.jl")

# 4. Expression Map
include("on_varieties/expression_map.jl")
include("on_varieties/varieties/map_graph.jl")

# 5. Scaling symmetries
include("on_varieties/scalings.jl")

# 6. Interpolation
include("on_varieties/interpolation.jl")

end
