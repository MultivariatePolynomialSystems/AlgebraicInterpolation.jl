export AbstractLieAlgebra,
    ScalingLieAlgebra,
    name,
    dim,
    LieAlgebra,
    basis,
    cartan_subalgebra,
    positive_roots,
    negative_roots,
    set_chevalley_basis!,
    set_cartan_subalgebra!,
    set_positive_roots!,
    set_negative_roots!,
    weight_structure,
    weights,
    nweights,
    weight_spaces,
    SumLieAlgebra,
    ⊕,
    AbstractLieAlgebraElem,
    LieAlgebraElem,
    as_matrix,
    act,
    SumLieAlgebraElem,
    AbstractLieAlgebraAction,
    LieAlgebraAction,
    ScalingLieAction,
    SumLieAlgebraAction


abstract type AbstractLieAlgebra end

struct ScalingLieAlgebra <: AbstractLieAlgebra
    name::String
    dim::Int
end

name(alg::ScalingLieAlgebra) = alg.name
dim(alg::ScalingLieAlgebra) = alg.dim

# Implements basic matrix Lie algebra
mutable struct LieAlgebra <: AbstractLieAlgebra
    name::String
    basis::Vector{Matrix{ComplexF64}}
    chevalley_basis::Vector{Vector{Vector{ComplexF64}}} # given by coefficients in basis; [cartan, positive, negative]
    weight_structure::WeightStructure
end

function Base.show(io::IO, alg::LieAlgebra)
    println(io, "LieAlgebra $(name(alg))")
    print(io, " dimension: $(dim(alg))")
end

LieAlgebra(
    basis::Vector{Matrix{T}} where T <: Number
) = LieAlgebra("", basis, [[], [], []], WeightStructure())

LieAlgebra(
    name::String,
    basis::Vector{Matrix{T}} where T <: Number
) = LieAlgebra(name, basis, [[], [], []], WeightStructure())

function so3_lie_algebra()
    X₁ = [0 0 0; 0 0 -1; 0 1 0]
    X₂ = [0 0 1; 0 0 0; -1 0 0]
    X₃ = [0 -1 0; 1 0 0; 0 0 0]
    so3 = LieAlgebra("𝖘𝖔(3,ℂ)", [X₁, X₂, X₃])
    set_cartan_subalgebra!(so3, [[0, 0, im]])
    set_positive_roots!(so3, [[im, -1, 0]])
    set_negative_roots!(so3, [[im, 1, 0]])
    so3.weight_structure = WeightStructure([-1, 0, 1], [[1, -im, 0], [0, 0, 1], [1, im, 0]])
    return so3
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
dim(alg::LieAlgebra) = length(alg.basis)
Base.size(alg::LieAlgebra) = size(alg.basis[1], 1)

function basis(alg::LieAlgebra)
    coeffs = eye(ComplexF64, dim(alg))
    return [LieAlgebraElem(c, alg) for c in eachcol(coeffs)]
end

cartan_subalgebra(alg::LieAlgebra) = [LieAlgebraElem(coeffs, alg) for coeffs in alg.chevalley_basis[1]]
positive_roots(alg::LieAlgebra) = [LieAlgebraElem(coeffs, alg) for coeffs in alg.chevalley_basis[2]]
negative_roots(alg::LieAlgebra) = [LieAlgebraElem(coeffs, alg) for coeffs in alg.chevalley_basis[3]]

set_chevalley_basis!(
    alg::LieAlgebra,
    B::NTuple{3, Vector{Vector{T}}} where T <: Number
) = alg.chevalley_basis = B

set_cartan_subalgebra!(
    alg::LieAlgebra,
    B::Vector{Vector{T}} where T <: Number
) = alg.chevalley_basis[1] = B

set_positive_roots!(
    alg::LieAlgebra,
    B::Vector{Vector{T}} where T <: Number
) = alg.chevalley_basis[2] = B

set_negative_roots!(
    alg::LieAlgebra,
    B::Vector{Vector{T}} where T <: Number
) = alg.chevalley_basis[3] = B

weight_structure(alg::LieAlgebra) = alg.weight_structure
weights(alg::LieAlgebra) = alg.weight_structure.weights
weights(alg::LieAlgebra, inds...) = getindex(alg.weight_structure.weights, inds...)
nweights(alg::LieAlgebra) = nweights(alg.weight_structure)
weight_spaces(alg::LieAlgebra) = alg.weight_structure.weight_spaces
weight_spaces(alg::LieAlgebra, inds...) = getindex(alg.weight_structure.weight_spaces, inds...)


struct SumLieAlgebra <: AbstractLieAlgebra
    algs::Vector{AbstractLieAlgebra}
end

algebras(g::SumLieAlgebra) = g.algs
name(g::SumLieAlgebra) = join([name(alg) for alg in algebras(g)], " ⊕ ")

⊕(
    alg₁::LieAlgebra,
    alg₂::LieAlgebra
) = SumLieAlgebra([alg₁, alg₂])

abstract type AbstractLieAlgebraElem end

struct LieAlgebraElem <: AbstractLieAlgebraElem
    alg::LieAlgebra
    coeffs::Vector{ComplexF64}
end

Base.size(elem::LieAlgebraElem) = size(elem.algebra)

# called by Shift+Enter
function Base.show(io::IO, mime::MIME"text/plain", elem::LieAlgebraElem)
    println(io, "LieAlgebraElem of $(name(elem.algebra)):")
    println(io, " matrix representation:")
    show(io, mime, as_matrix(elem)) # TODO: add offset
end

# called by print and inside vectors/matrices
function Base.show(io::IO, elem::LieAlgebraElem)
    print(io, "LieAlgebraElem from $(name(elem.algebra)) with coefficients: ")
    show(io, elem.coeffs)
end

Base.randn(alg::LieAlgebra) = LieAlgebraElem(randn(ComplexF64, dim(alg)), alg)

Base.:*(a::Number, elem::LieAlgebraElem) = LieAlgebraElem(a*elem.coeffs, elem.algebra)
Base.:*(elem::LieAlgebraElem, a::Number) = LieAlgebraElem(a*elem.coeffs, elem.algebra)

function Base.:+(X::LieAlgebraElem, Y::LieAlgebraElem)
    @assert X.algebra == Y.algebra
    return LieAlgebraElem(X.coeffs+Y.coeffs, X.algebra)
end

function as_matrix(elem::LieAlgebraElem)
    X = zeros(ComplexF64, size(elem), size(elem))
    for (i, B) in enumerate(elem.algebra.basis)
        X += elem.coeffs[i] * B
    end
    return X
end

function act(elem::LieAlgebraElem, f::Expression, vars::Vector{Variable})
    @assert size(elem) == length(vars)
    X = as_matrix(elem)
    return dot(differentiate(f, vars), -X*vars)
end

function act(elem::LieAlgebraElem, f::Expression; var_groups::Vector{Vector{Variable}})
    X = as_matrix(elem)
    return sum([dot(differentiate(f, vars), -X*vars) for vars in var_groups])
end

function act(elem::LieAlgebraElem, mon::Monomial)
    @assert size(elem) == nvariables(mon)
    X = as_matrix(elem)
    return dot(differentiate(mon), -X*variables(mon))
end

function act(elem::LieAlgebraElem, mon::Monomial, var_groups::Vector{Vector{Variable}})
    X = as_matrix(elem)
    return sum([dot(differentiate(mon, vars), -X*vars) for vars in var_groups])
end

function as_matrix(elem::LieAlgebraElem, B::MonomialBasis)
    @assert size(elem) == nvariables(B)
    M = zeros(ComplexF64, length(B), length(B))
    for (i, mon) in enumerate(B)
        gMon = act(elem, mon)
        M[:, i] = coefficients(gMon, B)
    end
    return M
end

function as_matrix(elem::LieAlgebraElem, B::MonomialBasis, var_groups::Vector{Vector{Variable}})
    M = zeros(ComplexF64, length(B), length(B))
    for (i, mon) in enumerate(B)
        gMon = act(elem, mon, var_groups)
        M[:, i] = coefficients(gMon, B)
    end
    return M
end

struct SumLieAlgebraElem <: AbstractLieAlgebraElem
    alg::SumLieAlgebra
    elems::Vector{AbstractLieAlgebraElem}
end


abstract type AbstractLieAlgebraAction end

struct LieAlgebraAction <: AbstractLieAlgebraAction
    alg::LieAlgebra
    vars::Vector{Vector{Variable}}
end

LieAlgebraAction(
    alg::LieAlgebra,
    vars::AbstractVecOrMat{Variable}
) = LieAlgebraAction(alg, M2VV(hcat(vars)))

LieAlgebraAction(
    alg::LieAlgebra,
    vars::AbstractArray
) = LieAlgebraAction(alg, M2VV(hcat(vars...)))

algebra(g::LieAlgebraAction) = g.alg
action(g::LieAlgebraAction) = g.vars

function Base.show(io::IO, g::LieAlgebraAction)
    println(io, "LieAlgebraAction of $(name(algebra(g)))")
    print(io, " action: [", join([join(vars, ", ") for vars in action(g)], "], ["), "]")
end

# TODO
function are_commutative(g₁::LieAlgebraAction, g₂::LieAlgebraAction)
    return true
end

function ⊕(g₁::LieAlgebraAction, g₂::LieAlgebraAction)
    @assert are_commutative(g₁, g₂)
    alg = algebra(g₁) ⊕ algebra(g₂)
    return SumLieAlgebraAction(alg, [g₁, g₂])
end

struct ScalingLieAction{T<:AbstractMatrix{<:Integer}} <: AbstractLieAlgebraAction
    alg::ScalingLieAlgebra
    vars::Vector{Variable}
    exps::T # every row is a vector u = [u₁,...,uₖ] which acts on vars by λᵘ
end

struct SumLieAlgebraAction <: AbstractLieAlgebraAction
    alg::SumLieAlgebra
    actions::Vector{AbstractLieAlgebraAction}
end

algebra(g::SumLieAlgebraAction) = g.alg
actions(g::SumLieAlgebraAction) = g.actions
nsummands(g::SumLieAlgebraAction) = length(actions(g))

function Base.show(io::IO, g::SumLieAlgebraAction)
    println(io, "SumLieAlgebraAction of $(name(algebra(g)))")
    for (i, a) in enumerate(actions(g))
        print(io, " action of $(name(algebra(a))): [")
        print(io, join([join(vars, ", ") for vars in action(a)], "], ["), "]")
        i < nsummands(g) && print(io, "\n")
    end
end