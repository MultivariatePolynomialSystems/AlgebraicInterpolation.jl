export LieAlgebraRepresentation,
    space,
    isotypic_components,
    highest_weights


struct LieAlgebraRepresentation{T<:AbstractLieAlgebraAction} <: AbstractLieAlgebraRepresentation
    action::T
    V::PolynomialVectorSpace
    isotypic::Vector{IsotypicComponent{T}}
end

function LieAlgebraRepresentation(
    action::T,
    V::PolynomialVectorSpace,
    irreds::Vector{IrreducibleRepresentation{T}}
) where T <: AbstractLieAlgebraAction
    iso_dict = Dict{Vector{Int}, Vector{IrreducibleRepresentation{T}}}()
    for irr in irreds
        irrs = get(iso_dict, highest_weight(irr), nothing)
        if isnothing(irrs)
            iso_dict[highest_weight(irr)] = [irr]
        else
            push!(irrs, irr)
        end
    end
    iso_comps = [IsotypicComponent(action, hw, iso) for (hw, iso) in iso_dict]
    return LieAlgebraRepresentation(action, V, iso_comps)
end

action(ρ::LieAlgebraRepresentation) = ρ.action
algebra(ρ::LieAlgebraRepresentation) = algebra(ρ.action)
space(ρ::LieAlgebraRepresentation) = ρ.V
nisotypic(ρ::LieAlgebraRepresentation) = length(ρ.isotypic)
isotypic_components(ρ::LieAlgebraRepresentation) = ρ.isotypic
irreducibles(ρ::LieAlgebraRepresentation) = vcat([irreducibles(iso) for iso in ρ.isotypic]...)
dim(ρ::LieAlgebraRepresentation) = dim(space(ρ))
highest_weights(ρ::LieAlgebraRepresentation) = [highest_weight(ic) for ic in isotypic_components(ρ)]

# called by Shift+Enter
function Base.show(io::IO, ::MIME"text/plain", ρ::LieAlgebraRepresentation)
    println(
        io,
        "LieAlgebraRepresentation of $(name(algebra(ρ))) ",
        "on the $(dim(ρ))-dimensional vector space"
    )
    println(io, " Lie algebra: ", name(algebra(ρ)))
    println(io, " number of isotypic components: ", nisotypic(ρ))
    print(io, " dimensions of isotypic components: ", join([dim(ic) for ic in isotypic_components(ρ)], ", "))
end

# called by print and inside vectors/matrices
function Base.show(io::IO, ρ::LieAlgebraRepresentation)
    print(
        io,
        "LieAlgebraRepresentation of $(name(algebra(ρ))) ",
        "on the $(dim(ρ))-dimensional vector space"
    )
end

function nullspace_as_weight_vectors(
    M::Matrix{T},
    ws::WeightStructure;
    tol::Real=1e-5
) where {T <: Number}
    wvs = WeightVector[]
    for weight_space in weight_spaces(ws)
        N = space(weight_space)*nullspace(M*space(weight_space))
        sparsify!(N, tol)
        vs = M2VV(N; copy=false)
        for v in vs
            v = div_by_lowest_magnitude(v, tol)
            push!(wvs, WeightVector(weight(weight_space), v))
        end
    end
    return wvs
end

sym_weight_structure(
    alg::LieAlgebra,
    d::Integer,
    mexps::Vector{<:SparseVector}
) = sym_weight_structure(weight_structure(alg), d, mexps)

sym_weight_structure(
    ws::WeightStructure,
    d::Integer,
    mons::MonomialBasis
) = sym_weight_structure(ws, d, multiexponents(mons))

tensor_weight_structure(
    ws::WeightStructure,
    ds::AbstractVector{<:Integer},
    v_mexps::Vector{<:Vector{<:SparseVector}},
    tensor_basis # TODO: add type
) = tensor_weight_structure([sym_weight_structure(ws, d, mexps) for (d, mexps) in zip(ds, v_mexps)], tensor_basis)

tensor_weight_structure(
    alg::LieAlgebra,
    ds::AbstractVector{<:Integer},
    v_mexps::Vector{<:Vector{<:SparseVector}},
    tensor_basis # TODO: add type
) = tensor_weight_structure([sym_weight_structure(alg, d, mexps) for (d, mexps) in zip(ds, v_mexps)], tensor_basis)

function weight_module(
    action::LieAlgebraAction,
    degrees::AbstractVector{<:Integer}
)
    v_mexps = [multiexponents(degree=Int8(d), nvars=Int16(length(vars)), upto=false) for (d, vars) in zip(degrees, var_groups(action))]
    tensor_basis = Base.Iterators.product([1:length(mexps) for mexps in v_mexps]...)
    ws = tensor_weight_structure(weight_structure(algebra(action)), degrees, v_mexps, tensor_basis)
    mon_bases = [MonomialBasis{Int8,Int16}(var_groups(action), i; degree=d, upto=false) for (i, d) in enumerate(degrees)]
    mons = ⊗(mon_bases, tensor_basis; equal_vars=true)
    return WeightModule(mons, ws)
end

function to_irreducible(
    action::AbstractLieAlgebraAction,
    wm::WeightModule
)
    ρXs = vcat([as_matrix(pos_root, basis(wm), action) for pos_root in positive_root_elements(algebra(action))]...)
    hw_vectors = nullspace_as_weight_vectors(ρXs, weight_structure(wm))
    hw_modules = [HighestWeightModule(basis(wm), hwv) for hwv in hw_vectors]
    return [IrreducibleRepresentation(action, hwm) for hwm in hw_modules]
end

# # TODO: supposes that all the vars in V occur in var_groups
# function LieAlgebraRepresentation(
#     action::LieAlgebraAction,
#     V::PolynomialVectorSpace
# )
#     var_grs = var_groups(action)
#     @assert issetequal(vcat(var_grs...), variables(V))
#     groups_mexps = multiexponents(degree=degree(V), nvars=length(var_grs), upto=is_upto(V))
#     irreds = IrreducibleRepresentation[]
#     for mexp in groups_mexps
#         wm = weight_module(action, mexp)
#         append!(irreds, to_irreducible(action, wm))
#     end
#     return LieAlgebraRepresentation(action, V, irreds)
# end

# TODO: compose into isotypic right away
function LieAlgebraRepresentation(
    action::AbstractLieAlgebraAction,
    V::PolynomialVectorSpace
)
    ws = weight_structure(action, variables(V))
    irreds = IrreducibleRepresentation{SumLieAlgebraAction}[]
    for d in degrees(V)
        mons = MonomialBasis{Int8,Int16}(variables(V); degree=d, upto=false)
        ws_d = sym_weight_structure(ws, d, mons)
        wm = WeightModule(mons, ws_d)
        append!(irreds, to_irreducible(action, wm))
    end
    return LieAlgebraRepresentation(action, V, irreds)
end
