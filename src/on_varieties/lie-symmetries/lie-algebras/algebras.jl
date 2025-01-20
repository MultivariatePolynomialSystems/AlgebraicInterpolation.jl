export AbstractLieAlgebra,
    ScalingLieAlgebra,
    name,
    dim,
    exponents,
    rank,
    LieAlgebra,
    basis,
    cartan_subalgebra,
    positive_root_elements,
    negative_root_elements,
    weight_structure,
    weights,
    nweights,
    weight_spaces,
    SumLieAlgebra,
    ⊕


abstract type AbstractLieAlgebra end


struct ScalingLieAlgebra <: AbstractLieAlgebra
    name::String
    exps::SparseMatrixCSC # every row is a vector u = [u₁,...,uₖ] which acts on vars by λᵘ
end

function ScalingLieAlgebra(exps::Matrix{Int})
    name = size(exps, 1) == 1 ? "ℂ" : "ℂ$(superscript(size(exps, 1)))"
    return ScalingLieAlgebra(name, sparse(exps))
end

ScalingLieAlgebra(size::Int) = ScalingLieAlgebra("ℂ", sparse(ones(Int, 1, size)))

name(alg::ScalingLieAlgebra) = alg.name
dim(alg::ScalingLieAlgebra) = size(alg.exps, 1)
exponents(alg::ScalingLieAlgebra) = alg.exps
MultivariateInterpolation.rank(alg::ScalingLieAlgebra) = dim(alg)
Base.size(alg::ScalingLieAlgebra) = size(alg.exps, 2)

function basis(alg::ScalingLieAlgebra; as_matrices::Bool=false)
    as_matrices && return [Diagonal(e) for e in eachrow(alg.exps)]
    coeffs = eye(ComplexF64, dim(alg))
    return [LieAlgebraElem(alg, c) for c in eachcol(coeffs)]
end

cartan_subalgebra(alg::ScalingLieAlgebra) = basis(alg)
positive_root_elements(::ScalingLieAlgebra) = LieAlgebraElem{ScalingLieAlgebra}[]
negative_root_elements(::ScalingLieAlgebra) = LieAlgebraElem{ScalingLieAlgebra}[]

struct ChevalleyBasis
    std_basis::Vector{Matrix{ComplexF64}}
    cartan::Vector{Vector{ComplexF64}} # given by coefficients in std_basis
    positive::Vector{Vector{ComplexF64}} # given by coefficients in std_basis
    negative::Vector{Vector{ComplexF64}} # given by coefficients in std_basis
    positive_roots::Vector{Vector{Int}}
    negative_roots::Vector{Vector{Int}}
end

# Implements basic matrix Lie algebra
mutable struct LieAlgebra <: AbstractLieAlgebra
    name::String
    basis::ChevalleyBasis
    weight_structure::WeightStructure
    hw_spaces::Vector{WeightSpace}
end

function Base.show(io::IO, alg::LieAlgebra)
    println(io, "LieAlgebra $(name(alg))")
    println(io, " dimension: $(dim(alg))")
    print(io, " rank (dimension of Cartan subalgebra): $(rank(alg))")
end

function so3_lie_algebra()
    X₁ = [0 0 0; 0 0 -1; 0 1 0]
    X₂ = [0 0 1; 0 0 0; -1 0 0]
    X₃ = [0 -1 0; 1 0 0; 0 0 0]
    cartan = [[0, 0, im]]
    positive = [[im, -1, 0]]
    negative = [[im, 1, 0]]
    pos_roots = [[1]]
    neg_roots = [[-1]]
    ch_basis = ChevalleyBasis([X₁, X₂, X₃], cartan, positive, negative, pos_roots, neg_roots)
    ws = WeightStructure([-1, 0, 1], [[1, -im, 0], [0, 0, 1], [1, im, 0]])
    hw_spaces = [WeightSpace([1], [1, im, 0])]
    return LieAlgebra("𝖘𝖔(3,ℂ)", ch_basis, ws, hw_spaces)
end

# TODO
function LieAlgebra(name::String, size::Int)
    if name == "so" && size == 3
        return so3_lie_algebra()
    else
        error("Not implemented")
    end
end

name(alg::LieAlgebra) = alg.name
dim(alg::LieAlgebra) = length(alg.basis.std_basis)
MultivariateInterpolation.rank(alg::LieAlgebra) = length(alg.basis.cartan)
Base.size(alg::LieAlgebra) = size(alg.basis.std_basis[1], 1)

function basis(alg::LieAlgebra; as_matrices::Bool=false)
    if as_matrices
        return alg.basis.std_basis
    end
    coeffs = eye(ComplexF64, dim(alg))
    return [LieAlgebraElem(alg, c) for c in eachcol(coeffs)]
end

cartan_subalgebra(alg::LieAlgebra) = [LieAlgebraElem(alg, coeffs) for coeffs in alg.basis.cartan]
positive_root_elements(
    alg::LieAlgebra
) = [LieAlgebraElem(alg, coeffs, root) for (coeffs, root) in zip(alg.basis.positive, alg.basis.positive_roots)]
negative_root_elements(
    alg::LieAlgebra
) = [LieAlgebraElem(alg, coeffs, root) for (coeffs, root) in zip(alg.basis.negative, alg.basis.negative_roots)]
weight_structure(alg::LieAlgebra) = alg.weight_structure
weights(alg::LieAlgebra) = alg.weight_structure.weights
weights(alg::LieAlgebra, inds...) = getindex(alg.weight_structure.weights, inds...)
nweights(alg::LieAlgebra) = nweights(alg.weight_structure)
weight_spaces(alg::LieAlgebra) = alg.weight_structure.weight_spaces
weight_spaces(alg::LieAlgebra, inds...) = getindex(alg.weight_structure.weight_spaces, inds...)
hw_spaces(alg::LieAlgebra) = alg.hw_spaces


struct SumLieAlgebra <: AbstractLieAlgebra
    name::String
    algs::Vector{AbstractLieAlgebra}
end

name(g::SumLieAlgebra) = g.name
algebras(g::SumLieAlgebra) = g.algs
dim(g::SumLieAlgebra) = sum([dim(alg) for alg in algebras(g)])
MultivariateInterpolation.rank(g::SumLieAlgebra) = sum([rank(alg) for alg in algebras(g)])

⊕(
    alg₁::AbstractLieAlgebra,
    alg₂::AbstractLieAlgebra
) = SumLieAlgebra("$(name(alg₁)) ⊕ $(name(alg₂))", [alg₁, alg₂])

⊕(
    alg::SumLieAlgebra,
    alg₂::AbstractLieAlgebra
) = SumLieAlgebra("$(name(alg)) ⊕ $(name(alg₂))", [algebras(alg)..., alg₂])

function cartan_subalgebra(alg::SumLieAlgebra)
    cartan_elems = SumLieAlgebraElem[]
    for (i, a) in enumerate(algebras(alg))
        a_elems = cartan_subalgebra(a)
        alg_elems = [zero(alg) for _ in a_elems]
        for (j, elem) in enumerate(a_elems)
            alg_elems[j][i] = elem
        end
        append!(cartan_elems, alg_elems)
    end
    return cartan_elems
end

function positive_root_elements(alg::SumLieAlgebra)
    pos_root_elems = SumLieAlgebraElem[]
    for (i, a) in enumerate(algebras(alg))
        a_elems = positive_root_elements(a)
        alg_elems = [zero(alg) for _ in a_elems]
        for (j, elem) in enumerate(a_elems)
            alg_elems[j][i] = elem
            root = zeros(Int, rank(alg))
            root[i] = 1
            alg_elems[j].root = root
        end
        append!(pos_root_elems, alg_elems)
    end
    return pos_root_elems
end

function negative_root_elements(alg::SumLieAlgebra)
    neg_root_elems = SumLieAlgebraElem[]
    for (i, a) in enumerate(algebras(alg))
        a_elems = negative_root_elements(a)
        alg_elems = [zero(alg) for _ in a_elems]
        for (j, elem) in enumerate(a_elems)
            alg_elems[j][i] = elem
            root = zeros(Int, rank(alg))
            root[i] = -1
            alg_elems[j].root = root
        end
        append!(neg_root_elems, alg_elems)
    end
    return neg_root_elems
end