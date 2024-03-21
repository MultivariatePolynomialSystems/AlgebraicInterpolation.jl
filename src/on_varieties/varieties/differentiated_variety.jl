export DifferentiatedVariety


"""
    DifferentiatedVariety <: AbstractDifferentiatedVariety
"""
struct DifferentiatedVariety <: AbstractDifferentiatedVariety
    system::System
    full_jacobian::Matrix{Expression}
end

DifferentiatedVariety(F::System) = DifferentiatedVariety(F, full_jacobian(F))
HC.System(F::DifferentiatedVariety) = F.system

(F::DifferentiatedVariety)(x::AbstractVector) = System(F)(x)
find_sample(F::DifferentiatedVariety; kwargs...) = find_sample(System(F); kwargs...)

full_jacobian(F::DifferentiatedVariety) = F.full_jacobian
full_jacobian(
    F::AbstractDifferentiatedVariety,
    x::AbstractVector
) = evaluate(full_jacobian(F), variables(F) => x)

tangent_space(
    F::AbstractDifferentiatedVariety,
    x::AbstractVector;
    tols::Tolerances=Tolerances(),
    var_ids::Union{Vector{Int}, Colon}=:
) = nullspace(full_jacobian(F, x); atol=tols.nullspace_atol)[var_ids, :]

function dimension(
    F::AbstractDifferentiatedVariety,
    x::Union{AbstractVector, Nothing}=nothing;
    tols::Tolerances=Tolerances()
)
    x = isnothing(x) ? find_sample(F) : x
    if !isnothing(x)
        @assert norm(F(x)) < tols.common_tol
        J = full_jacobian(F, x)
        return nvariables(F) - rank(J; atol=tols.rank_atol)
    else
        error("Cannot generate a random sample of the system, provide one!")
    end
end