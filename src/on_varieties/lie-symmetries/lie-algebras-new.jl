abstract type AbstractLieAlgebra end

struct LieAlgebra <: AbstractLieAlgebra

end

struct SumLieAlgebra <: AbstractLieAlgebra
    algs::Vector{AbstractLieAlgebra}
end


abstract type AbstractLieAlgebraElem end

struct LieAlgebraElem <: AbstractLieAlgebraElem
    alg::LieAlgebra
    coeffs::Vector{ComplexF64}
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

struct SumLieAlgebraAction <: AbstractLieAlgebraAction
    alg::SumLieAlgebra
end

struct ScalingLieAction <: AbstractLieAlgebraAction

end