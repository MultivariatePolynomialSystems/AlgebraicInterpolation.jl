export interpolate_constraints


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

function interpolate_constraints(
    F::AbstractAlgebraicVariety,
    vars::FixedFreeVariables;
    degree::Integer,
    start_point::Union{AbstractVector, Nothing}=nothing,
    tols::Tolerances=Tolerances(),
    sample_rand_method::Symbol=:rand_unit
)
    mons = MonomialVector{Int8, Int16}(free(vars); degree=degree)
    s = sample(F, vars, start_point; nsamples=length(mons), rand_method=sample_rand_method)

    A = transpose(evaluate(mons, free(s))) # Vandermonde matrix
    if tols.nullspace_rtol == 0
        N = nullspace(A; atol=tols.nullspace_atol)
    else
        N = nullspace(A; atol=tols.nullspace_atol, rtol=tols.nullspace_rtol)
    end

    coeffs = transpose(N)
    size(coeffs, 1) == 0 && return Expression[]

    coeffs = rref(coeffs, tols.rref_tol)
    sparsify!(coeffs, tols.sparsify_tol; digits=1)
    return polynomial_functions(coeffs, mons)
end
