abstract type AbstractLieAlgebraRepresentation end

struct IrreducibleRepresentation{T<:AbstractLieAlgebraAction} <: AbstractLieAlgebraRepresentation
    action::T
    hw_module::HighestWeightModule
end

algebra(π::IrreducibleRepresentation) = algebra(π.action)
space_basis(π::IrreducibleRepresentation) = basis(π.hw_module)
highest_weight(π::IrreducibleRepresentation) = highest_weight(π.hw_module)
highest_weight_vector(π::IrreducibleRepresentation) = highest_weight_vector(π.hw_module)
dim(π::IrreducibleRepresentation) = prod([2*j+1 for j in highest_weight(π)]) # TODO: works only for so(3)

function to_expressions(π::IrreducibleRepresentation; tol::Float64=1e-5)
    exprs = Expression[]
    # TODO: extend the following to multiple negative roots
    πJ₋ = as_matrix(negative_roots(algebra(π))[1], space_basis(π), π.var_groups)
    expr_mons = to_expressions(space_basis(π))
    v = vector(highest_weight_vector(π))
    while norm(v) > tol
        v = div_by_lowest_magnitude(v, tol)
        sparsify!(v, tol)
        v = simplify_numbers(v)
        push!(exprs, dot(v, expr_mons))
        v = ComplexF64.(πJ₋*v)
    end
    return exprs
end


struct IsotypicComponent{T<:AbstractLieAlgebraAction} <: AbstractLieAlgebraRepresentation
    action::T
    highest_weight::Vector{Int}
    irreds::Vector{HighestWeightModule}
end

highest_weight(ic::IsotypicComponent) = ic.highest_weight
dim(ic::IsotypicComponent) = sum([dim(irr) for irr in ic.irreds])
mul(ic::IsotypicComponent) = length(ic.irreds)
irreducible_components(ic::IsotypicComponent) = ic.irreds

to_expressions(
    ic::IsotypicComponent;
    tol::Float64=1e-5
) = vcat([to_expressions(π; tol=tol) for π in ic.irreds]...)


struct PolynomialVectorSpace
    vars::Vector{Variable}
    degree::Int
    upto::Bool
end

PolynomialVectorSpace(;
    variables::Vector{Variable},
    degree::Int,
    upto::Bool=true
) = PolynomialVectorSpace(variables, degree, upto)

variables(V::PolynomialVectorSpace) = V.vars
nvariables(V::PolynomialVectorSpace) = length(V.vars)
degree(V::PolynomialVectorSpace) = V.degree
is_upto(V::PolynomialVectorSpace) = V.upto
dim(V::PolynomialVectorSpace) = num_mons(nvariables(V), degree(V); upto=is_upto(V))

Base.rand(
    V::PolynomialVectorSpace,
    n::Int
) = random_monomial_basis(length=n, nvars=nvariables(V), degree=degree(V), upto=upto(V))

Base.:(==)(
    V::PolynomialVectorSpace,
    W::PolynomialVectorSpace
) = (V.vars == W.vars) && (V.degree == W.degree) && (V.upto == W.upto)

function Base.show(io::IO, V::PolynomialVectorSpace)
    println(io, "PolynomialVectorSpace of dimension $(dim(V))")
    print(io, " variables: $(variables(V))")
end


struct LieAlgebraRepresentation{T<:AbstractLieAlgebraAction} <: AbstractLieAlgebraRepresentation
    action::T
    V::PolynomialVectorSpace
    isotypic::Vector{IsotypicComponent{T}}
end

function LieAlgebraRepresentation(
    alg::LieAlgebra,
    V::PolynomialVectorSpace,
    var_groups::Vector{Vector{Variable}},
    irreds::Vector{IrreducibleRepresentation}
)
    iso_dict = Dict{Vector{Int}, Vector{IrreducibleRepresentation}}()
    for irr in irreds
        irrs = get(iso_dict, highest_weight(irr), nothing)
        if isnothing(irrs)
            iso_dict[highest_weight(irr)] = [irr]
        else
            push!(irrs, irr)
        end
    end
    iso_comps = [IsotypicComponent(alg, hw, iso) for (hw, iso) in iso_dict]
    return LieAlgebraRepresentation(alg, V, var_groups, iso_comps)
end

algebra(π::LieAlgebraRepresentation) = π.alg
space(π::LieAlgebraRepresentation) = π.V
isotypic_components(π::LieAlgebraRepresentation) = π.isotypic
irreducible_components(π::LieAlgebraRepresentation) = vcat([irreducible_components(iso) for iso in π.isotypic]...)
dim(π::LieAlgebraRepresentation) = sum([dim(ic) for ic in isotypic_components(π)])

# called by Shift+Enter
function Base.show(io::IO, mime::MIME"text/plain", π::LieAlgebraRepresentation)
    println(
        io,
        "LieAlgebraRepresentation of $(name(π.alg)) ",
        "on the $(dim(π))-dimensional vector space:"
    )
    show(io, mime, π.var_groups)
    # print(io, " action on variables: $(π.var_groups)")
end

function nullspace_as_weight_vectors(
    M::Matrix{T},
    ws::WeightStructure
) where {T <: Number}
    wvs = WeightVector[]
    for (weight, weight_space) in zip(ws.weights, ws.weight_spaces)
        vs = M2VV(weight_space*nullspace(M*weight_space); copy=false)
        append!(wvs, [WeightVector(weight, v) for v in vs])
    end
    return wvs
end

function weight_module(
    alg::LieAlgebra,
    variables::Vector{Vector{Variable}},
    degrees::AbstractVector{<:Integer}
)
    v_mexps = [multiexponents(degree=Int8(d), nvars=Int16(length(vars)), upto=false) for (d, vars) in zip(degrees, variables)]
    tensor_basis = Base.Iterators.product([1:length(mexps) for mexps in v_mexps]...)
    ws = tensor_weight_structure(alg, degrees, v_mexps, tensor_basis)
    mon_bases = [MonomialBasis{Int8,Int16}(variables, i; degree=d, upto=false) for (i, d) in enumerate(degrees)]
    mons = ⊗(mon_bases, tensor_basis; equal_vars=true)
    return WeightModule(mons, ws)
end

function to_irreducible(
    alg::LieAlgebra,
    var_groups::Vector{Vector{Variable}},
    weight_module::WeightModule
)
    πXs = vcat([as_matrix(pos_root, basis(weight_module), var_groups) for pos_root in positive_roots(alg)]...)
    hw_vectors = nullspace_as_weight_vectors(πXs, weight_structure(weight_module))
    hw_modules = [HighestWeightModule(basis(weight_module), hwv) for hwv in hw_vectors]
    return [IrreducibleRepresentation(alg, var_groups, hwm) for hwm in hw_modules]
end

# TODO: supposes that all the vars in V occur in var_groups
function LieAlgebraRepresentation(
    alg::LieAlgebra,
    V::PolynomialVectorSpace;
    action::Vector{Vector{Variable}}
)
    @assert issetequal(vcat(action...), variables(V))
    groups_mexps = multiexponents(degree=degree(V), nvars=length(action), upto=is_upto(V))
    irreds = IrreducibleRepresentation[]
    for mexp in groups_mexps
        wm = weight_module(alg, action, mexp)
        append!(irreds, to_irreducible(alg, action, wm))
    end
    return LieAlgebraRepresentation(alg, V, action, irreds)
end