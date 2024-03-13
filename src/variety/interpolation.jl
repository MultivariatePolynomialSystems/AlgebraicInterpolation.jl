export interpolate_constraints


function interpolate_constraints(
    F::AbstractDifferentiatedSystem,
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
