export Tolerances


abstract type AbstractMonomialVector end

# TODO: add general_tol for everything else?
@kwdef struct Tolerances
    common_tol::Float64=1e-10
    nullspace_atol::Float64=0
    nullspace_rtol::Float64=0
    rank_atol::Float64=0
    rank_rtol::Float64=0
    rref_tol::Float64=1e-5
    sparsify_tol::Float64=1e-5
end

function polynomial_function(
    coeffs::AbstractVector{<:Number},
    mons::AbstractMonomialVector;
    logging::Bool=false
)
    @assert length(coeffs) == length(mons)
    p = sum(to_expressions(mons).*coeffs)
    logging && println("polynomial = ", p)
    return p
end

function polynomial_functions(
    coeffs::AbstractMatrix{<:Number},
    mons::AbstractMonomialVector;
    logging::Bool=false
)
    return [polynomial_function(row, mons; logging=logging) for row in eachrow(coeffs)]
end