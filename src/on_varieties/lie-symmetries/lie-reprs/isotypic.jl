export IsotypicComponent,
    irreducibles,
    nirreducible,
    highest_weight_expressions,
    highest_weight_vectors

struct IsotypicComponent{T<:AbstractLieAlgebraAction} <: AbstractLieAlgebraRepresentation
    action::T
    highest_weight::Vector{Int}
    irreds::Vector{IrreducibleRepresentation{T}}
end

action(ic::IsotypicComponent) = ic.action
algebra(ic::IsotypicComponent) = algebra(ic.action)
highest_weight(ic::IsotypicComponent) = ic.highest_weight
nirreducible(ic::IsotypicComponent) = length(ic.irreds)
dim(ic::IsotypicComponent) = nirreducible(ic)*dim(irreducibles(ic)[1])
irreducibles(ic::IsotypicComponent) = ic.irreds
highest_weight_expressions(ic::IsotypicComponent) = [highest_weight_expression(ρ) for ρ in irreducibles(ic)]
highest_weight_vectors(ic::IsotypicComponent) = [vector(highest_weight_vector(ρ)) for ρ in irreducibles(ic)]

function Base.show(io::IO, ::MIME"text/plain", ic::IsotypicComponent)
    println(io, "IsotypicComponent of dimension $(dim(ic))")
    println(io, " multiplicity of irreducible subrepresentation: ", nirreducible(ic))
    println(io, " dimension of irreducible subrepresentation: ", dim(irreducibles(ic)[1]))
    println(io, " Lie algebra: ", name(algebra(ic)))
    print(io, " highest weight: ", highest_weight(ic))
end

function Base.show(io::IO, ic::IsotypicComponent)
    print(io, "IsotypicComponent of dimension $(dim(ic))")
end

to_expressions(
    ic::IsotypicComponent
) = vcat([to_expressions(ρ) for ρ in irreducibles(ic)]...)