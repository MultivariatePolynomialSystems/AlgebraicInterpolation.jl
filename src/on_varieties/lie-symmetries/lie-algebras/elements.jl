export AbstractLieAlgebraElem,
    LieAlgebraElem,
    algebra,
    root,
    as_matrix,
    act,
    SumLieAlgebraElem


abstract type AbstractLieAlgebraElem end

struct LieAlgebraElem{T<:Union{LieAlgebra, ScalingLieAlgebra}} <: AbstractLieAlgebraElem
    alg::T
    coeffs::Vector{ComplexF64}
    root::Union{Nothing, Vector{Int}}
end

function LieAlgebraElem(alg::T, coeffs::AbstractVector) where T<:Union{LieAlgebra, ScalingLieAlgebra}
    @assert length(coeffs) == dim(alg)
    return LieAlgebraElem{T}(alg, coeffs, nothing)
end

algebra(elem::LieAlgebraElem) = elem.alg
root(elem::LieAlgebraElem) = elem.root
Base.size(elem::LieAlgebraElem) = size(algebra(elem))
Base.zero(alg::Union{LieAlgebra, ScalingLieAlgebra}) = LieAlgebraElem(alg, zeros(ComplexF64, dim(alg)))
Base.randn(alg::Union{LieAlgebra, ScalingLieAlgebra}) = LieAlgebraElem(alg, randn(ComplexF64, dim(alg)))

# called by Shift+Enter
function Base.show(io::IO, mime::MIME"text/plain", elem::LieAlgebraElem)
    println(io, "LieAlgebraElem of $(name(algebra(elem))):")
    println(io, " matrix representation:")
    show(io, mime, as_matrix(elem)) # TODO: add offset
end

# called by print and inside vectors/matrices
function Base.show(io::IO, elem::LieAlgebraElem)
    print(io, "LieAlgebraElem of $(name(algebra(elem))) with coefficients: ")
    show(io, elem.coeffs)
end

Base.:*(a::Number, elem::LieAlgebraElem) = LieAlgebraElem(elem.alg, a*elem.coeffs)
Base.:*(elem::LieAlgebraElem, a::Number) = LieAlgebraElem(elem.alg, a*elem.coeffs)

function Base.:+(X::LieAlgebraElem, Y::LieAlgebraElem)
    @assert algebra(X) == algebra(Y)
    return LieAlgebraElem(algebra(X), X.coeffs+Y.coeffs)
end

function as_matrix(elem::LieAlgebraElem)
    return sum([elem.coeffs[i]*Bᵢ for (i, Bᵢ) in enumerate(basis(algebra(elem); as_matrices=true))])
end

function act(elem::LieAlgebraElem, mon::Monomial)
    @assert size(elem) == nvariables(mon)
    X = as_matrix(elem)
    return expand(dot(differentiate(mon), -X*variables(mon)))
end

function act(elem::LieAlgebraElem, f::Expression, vars::Vector{Variable})
    @assert size(elem) == length(vars)
    X = as_matrix(elem)
    return expand(dot(differentiate(f, vars), -X*vars))
end

mutable struct SumLieAlgebraElem <: AbstractLieAlgebraElem
    alg::SumLieAlgebra
    elems::Vector{AbstractLieAlgebraElem}
    root::Union{Nothing, Vector{Int}}
end

function SumLieAlgebraElem(alg::SumLieAlgebra, elems::Vector{<:AbstractLieAlgebraElem})
    @assert length(elems) == length(algebras(alg))
    return SumLieAlgebraElem(alg, elems, nothing)
end

algebra(elem::SumLieAlgebraElem) = elem.alg
elems(elem::SumLieAlgebraElem) = elem.elems
root(elem::SumLieAlgebraElem) = elem.root
Base.randn(alg::SumLieAlgebra) = SumLieAlgebraElem(alg, [randn(a) for a in algebras(alg)])
Base.zero(alg::SumLieAlgebra) = SumLieAlgebraElem(alg, [zero(a) for a in algebras(alg)])
Base.getindex(elem::SumLieAlgebraElem, i::Int) = elems(elem)[i]
Base.setindex!(elem::SumLieAlgebraElem, val, i) = (elems(elem)[i] = val)

function as_matrix(elem::AbstractLieAlgebraElem, B::MonomialBasis)
    @assert size(elem) == nvariables(B)
    M = zeros(ComplexF64, length(B), length(B))
    for (i, mon) in enumerate(B)
        gMon = act(elem, mon)
        M[:, i] = coefficients(gMon, B)
    end
    return M
end
