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
    println(io, " number of isotypic components: ", nisotypic(ρ))
    print(io, " dimensions of isotypic components: ", join([dim(ic) for ic in isotypic_components(ρ)], ", "))
end

# 
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
    for (weight, weight_space) in zip(ws.weights, ws.weight_spaces)
        N = weight_space*nullspace(M*weight_space)
        sparsify!(N, tol)
        vs = M2VV(N; copy=false)
        for v in vs
            v = div_by_lowest_magnitude(v, tol)
            sparsify!(v, tol)
            push!(wvs, WeightVector(weight, v))
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
    action::LieAlgebraAction,
    weight_module::WeightModule
)
    ρXs = vcat([as_matrix(pos_root, basis(weight_module), var_groups(action)) for pos_root in positive_roots(algebra(action))]...)
    hw_vectors = nullspace_as_weight_vectors(ρXs, weight_structure(weight_module))
    hw_modules = [HighestWeightModule(basis(weight_module), hwv) for hwv in hw_vectors]
    return [IrreducibleRepresentation(action, hwm) for hwm in hw_modules]
end

# TODO: supposes that all the vars in V occur in var_groups
function LieAlgebraRepresentation(
    action::LieAlgebraAction,
    V::PolynomialVectorSpace
)
    var_grs = var_groups(action)
    @assert issetequal(vcat(var_grs...), variables(V))
    groups_mexps = multiexponents(degree=degree(V), nvars=length(var_grs), upto=is_upto(V))
    irreds = IrreducibleRepresentation[]
    for mexp in groups_mexps
        wm = weight_module(action, mexp)
        append!(irreds, to_irreducible(action, wm))
    end
    return LieAlgebraRepresentation(action, V, irreds)
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

# TODO:
# 1) Determine the highest weight spaces for each action:
# 1.1) Decompose the action on the variables into irreducibles, i.e. find the highest weight vectors
# 1.2) Compose isomorphic irreducibles together to form the highest weight spaces

# 2) Determine the highest weight spaces for the sum of 2 actions:
# 2.1) Compute all possible intersections of the highest weight spaces Vλ and Vτ of the 2 actions, Vλ ∩ Vτ corresponds
# to the weight [λ, τ]

# 3) Determine the highest weight spaces for the sum of more than 2 actions: apply the above procedure recursively

# 4) Determine the WeightStructure for the sum action from the highest weight spaces:
# 4.1) For each highest weight vector apply all possible combinations of lowering operators to get the weight vectors
# 4.2) Combine the weight vectors with equal weights into weight spaces and create the WeightStructure
function weight_structure(action::SumLieAlgebraAction)
    g₁, g₂, _ = actions(action)

    v = [1; im; 0; im; -1; 0; 0; 0; 0]
    J₋ = negative_root_elements(algebra(g₁))[1]
    
    @var E[1:3,1:3]
    mons = MonomialBasis{Int8,Int16}(E[:]; degree=1, upto=false)
    ρJ₋₁ = as_matrix(J₋, mons, g₁)
    ρJ₋₂ = as_matrix(J₋, mons, g₂)
    
    wv = [v, ρJ₋₂*v, ρJ₋₂^2*v, ρJ₋₁*v, ρJ₋₂*ρJ₋₁*v, ρJ₋₂^2*ρJ₋₁*v, ρJ₋₁^2*v, ρJ₋₂*ρJ₋₁^2*v, ρJ₋₂^2*ρJ₋₁^2*v]
    w = [[1,1,1], [1,0,1], [1,-1,1], [0,1,1], [0,0,1], [0,-1,1], [-1,1,1], [-1,0,1], [-1,-1,1]]
    
    return WeightStructure(w, wv)
end

function LieAlgebraRepresentation(
    action::SumLieAlgebraAction,
    V::PolynomialVectorSpace
)
    ws = weight_structure(action)
    irreds = IrreducibleRepresentation{SumLieAlgebraAction}[]
    for d in degrees(V)
        mons = MonomialBasis{Int8,Int16}(variables(V); degree=d, upto=false)
        ws_d = sym_weight_structure(ws, d, mons)
        wm = WeightModule(mons, ws_d)
        append!(irreds, to_irreducible(action, wm))
    end
    return LieAlgebraRepresentation(action, V, irreds)
end
