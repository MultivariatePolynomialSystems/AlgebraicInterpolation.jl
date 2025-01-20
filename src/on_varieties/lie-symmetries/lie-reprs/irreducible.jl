export AbstractLieAlgebraRepresentation,
    IrreducibleRepresentation,
    action,
    space_basis,
    highest_weight_expression,
    nisotypic,
    to_expressions

abstract type AbstractLieAlgebraRepresentation end

struct IrreducibleRepresentation{T<:AbstractLieAlgebraAction} <: AbstractLieAlgebraRepresentation
    action::T
    hw_module::HighestWeightModule
end

action(ρ::IrreducibleRepresentation) = ρ.action
algebra(ρ::IrreducibleRepresentation) = algebra(ρ.action)
space_basis(ρ::IrreducibleRepresentation) = basis(ρ.hw_module)
highest_weight_module(ρ::IrreducibleRepresentation) = ρ.hw_module
highest_weight(ρ::IrreducibleRepresentation) = highest_weight(ρ.hw_module)
highest_weight_vector(ρ::IrreducibleRepresentation) = highest_weight_vector(ρ.hw_module)
highest_weight_expression(ρ::IrreducibleRepresentation) = dot(vector(highest_weight_vector(ρ)), to_expressions(space_basis(ρ)))
dim(ρ::IrreducibleRepresentation) = prod([2*j+1 for j in highest_weight(ρ)[1:2]]) # TODO

function Base.show(io::IO, ::MIME"text/plain", ρ::IrreducibleRepresentation)
    println(io, "IrreducibleRepresentation of dimension $(dim(ρ))")
    println(io, " Lie algebra: ", name(algebra(ρ)))
    print(io, " highest weight: ", highest_weight(ρ))
end

function Base.show(io::IO, ρ::IrreducibleRepresentation)
    print(io, "IrreducibleRepresentation of dimension $(dim(ρ))")
end

function orbit(
    wexpr::WeightExpression,
    action::AbstractLieAlgebraAction,
    processed_weights::Set{Vector{Int}};
    tol::Float64=1e-5
)
    c = coefficients(wexpr)
    norm(c) < tol && return WeightExpression[]
    push!(processed_weights, weight(wexpr))
    orbit_exprs = [wexpr]
    for neg_root_elem in negative_root_elements(algebra(action))
        new_weight = weight(wexpr) + root(neg_root_elem)
        if new_weight ∉ processed_weights
            new_expr = act(neg_root_elem, expression(wexpr), action)
            coeffs = div_by_lowest_magnitude(coefficients(new_expr, basis(wexpr)), tol)
            new_wexpr = WeightExpression(coeffs, basis(wexpr), new_weight)
            append!(orbit_exprs, orbit(new_wexpr, action, processed_weights))
        end
    end
    return orbit_exprs
end

orbit(
    wexpr::WeightExpression,
    action::AbstractLieAlgebraAction;
    tol::Float64=1e-5
) = orbit(wexpr, action, Set{Vector{Int}}(); tol=tol)

function to_expressions(ρ::IrreducibleRepresentation; tol::Float64=1e-5, in_rref::Bool=true)
    v = vector(highest_weight_vector(ρ))
    mons = space_basis(ρ)
    wexpr = WeightExpression(v, mons, highest_weight(ρ))
    orb = orbit(wexpr, action(ρ); tol=tol)
    if in_rref
        coeffs = hcat([coefficients(f) for f in orb]...)
        coeffs = rref(Matrix{ComplexF64}(transpose(coeffs)))
        sparsify!(coeffs, tol)
        return [dot(simplify_numbers(a), to_expressions(mons)) for a in eachrow(coeffs)]
    end
    return [expression(wexpr) for wexpr in orb]
end