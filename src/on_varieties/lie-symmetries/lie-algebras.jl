export LieAlgebra,
    name,
    dim,
    size,
    basis,
    weight_structure,
    weights,
    weight_spaces,
    set_chevalley_basis!,
    set_cartan_subalgebra!,
    set_positive_roots!,
    set_negative_roots!,
    cartan_subalgebra,
    positive_roots,
    negative_roots,
    LieAlgebraElem,
    act,
    as_matrix,
    LieAlgebraRepresentation,
    IrreducibleLieAlgebraRepresentation,
    space_basis,
    to_irreducible,
    IsotypicComponent,
    mul,
    isotypic_components,
    irreducible_components,
    PolynomialVectorSpace,
    decompose,
    ×

struct WeightStructure
    weights::Vector{Vector{Int}}
    weight_spaces::Vector{Matrix{ComplexF64}}
end

WeightStructure() = WeightStructure([], [])

WeightStructure(
    weights::Vector{Int},
    weight_spaces::Vector{Matrix{T}} where T <: Number
) = WeightStructure([[w] for w in weights], weight_spaces)

WeightStructure(
    weights::Vector{Int},
    weight_vectors::Vector{Vector{T}} where T <: Number
) = WeightStructure(weights, [V2M(v) for v in weight_vectors])

weights(ws::WeightStructure) = ws.weights
weights(ws::WeightStructure, inds...) = getindex(ws.weights, inds...)
nweights(ws::WeightStructure) = length(ws.weights)
weight_spaces(ws::WeightStructure) = ws.weight_spaces
weight_spaces(ws::WeightStructure, inds...) = getindex(ws.weight_spaces, inds...)

# Base.getindex(ws::WeightStructure, weight::Vector{Int}) = ws.spaces[weight]

struct WeightVector
    weight::Vector{Int}
    vector::Vector{ComplexF64}
end

weight(wv::WeightVector) = wv.weight
vector(wv::WeightVector) = wv.vector

struct WeightModule
    basis::MonomialBasis
    weight_structure::WeightStructure
end

basis(wm::WeightModule) = wm.basis
weight_structure(wm::WeightModule) = wm.weight_structure

struct HighestWeightModule
    basis::MonomialBasis
    hw_vector::WeightVector
end

basis(hwm::HighestWeightModule) = hwm.basis
highest_weight(hwm::HighestWeightModule) = weight(hwm.hw_vector)
highest_weight_vector(hwm::HighestWeightModule) = hwm.hw_vector

# TODO: make parametric based on basis::Vector{Matrix{T}}?
mutable struct LieAlgebra
    name::String
    basis::Vector{Matrix{ComplexF64}}
    chevalley_basis::Vector{Vector{Vector{ComplexF64}}} # given by coefficients in basis; [cartan, positive, negative]
    weight_structure::WeightStructure # TODO: rename?
end

function Base.show(io::IO, alg::LieAlgebra)
    println(io, "LieAlgebra $(name(alg))")
    print(io, " dimension: $(dim(alg))")
end

LieAlgebra(
    basis::Vector{Matrix{T}} where T <: Number
) = LieAlgebra("", basis, [[], [], []], WeightStructure())

LieAlgebra(
    name::String,
    basis::Vector{Matrix{T}} where T <: Number
) = LieAlgebra(name, basis, [[], [], []], WeightStructure())

function so3_lie_algebra()
    X₁ = [0 0 0; 0 0 -1; 0 1 0]
    X₂ = [0 0 1; 0 0 0; -1 0 0]
    X₃ = [0 -1 0; 1 0 0; 0 0 0]
    so3 = LieAlgebra("𝖘𝖔(3,ℂ)", [X₁, X₂, X₃])
    set_cartan_subalgebra!(so3, [[0, 0, im]])
    set_positive_roots!(so3, [[im, -1, 0]])
    set_negative_roots!(so3, [[im, 1, 0]])
    so3.weight_structure = WeightStructure([-1, 0, 1], [[1, -im, 0], [0, 0, 1], [1, im, 0]])
    return so3
end

# TODO
function LieAlgebra(name::String, size::Int)
    if name == "so" && size == 3
        return so3_lie_algebra()
    else
        error("Not implemented")
    end
end

name(alg::LieAlgebra) = alg.name
dim(alg::LieAlgebra) = length(alg.basis)
Base.size(alg::LieAlgebra) = size(alg.basis[1], 1)

function basis(alg::LieAlgebra)
    coeffs = eye(ComplexF64, dim(alg))
    return [LieAlgebraElem(c, alg) for c in eachcol(coeffs)]
end

cartan_subalgebra(alg::LieAlgebra) = [LieAlgebraElem(coeffs, alg) for coeffs in alg.chevalley_basis[1]]
positive_roots(alg::LieAlgebra) = [LieAlgebraElem(coeffs, alg) for coeffs in alg.chevalley_basis[2]]
negative_roots(alg::LieAlgebra) = [LieAlgebraElem(coeffs, alg) for coeffs in alg.chevalley_basis[3]]

set_chevalley_basis!(
    alg::LieAlgebra,
    B::NTuple{3, Vector{Vector{T}}} where T <: Number
) = alg.chevalley_basis = B

set_cartan_subalgebra!(
    alg::LieAlgebra,
    B::Vector{Vector{T}} where T <: Number
) = alg.chevalley_basis[1] = B

set_positive_roots!(
    alg::LieAlgebra,
    B::Vector{Vector{T}} where T <: Number
) = alg.chevalley_basis[2] = B

set_negative_roots!(
    alg::LieAlgebra,
    B::Vector{Vector{T}} where T <: Number
) = alg.chevalley_basis[3] = B

weight_structure(alg::LieAlgebra) = alg.weight_structure
weights(alg::LieAlgebra) = alg.weight_structure.weights
weights(alg::LieAlgebra, inds...) = getindex(alg.weight_structure.weights, inds...)
nweights(alg::LieAlgebra) = nweights(alg.weight_structure)
weight_spaces(alg::LieAlgebra) = alg.weight_structure.weight_spaces
weight_spaces(alg::LieAlgebra, inds...) = getindex(alg.weight_structure.weight_spaces, inds...)

struct LieAlgebraElem
    coeffs::Vector{ComplexF64}
    algebra::LieAlgebra
end

Base.size(elem::LieAlgebraElem) = size(elem.algebra)

# called by Shift+Enter
function Base.show(io::IO, mime::MIME"text/plain", elem::LieAlgebraElem)
    println(io, "LieAlgebraElem of $(name(elem.algebra)):")
    println(io, " matrix representation:")
    show(io, mime, as_matrix(elem)) # TODO: add offset
end

# called by print and inside vectors/matrices
function Base.show(io::IO, elem::LieAlgebraElem)
    print(io, "LieAlgebraElem from $(name(elem.algebra)) with coefficients: ")
    show(io, elem.coeffs)
end

Base.:*(a::Number, elem::LieAlgebraElem) = LieAlgebraElem(a*elem.coeffs, elem.algebra)
Base.:*(elem::LieAlgebraElem, a::Number) = LieAlgebraElem(a*elem.coeffs, elem.algebra)

function Base.:+(X::LieAlgebraElem, Y::LieAlgebraElem)
    @assert X.algebra == Y.algebra
    return LieAlgebraElem(X.coeffs+Y.coeffs, X.algebra)
end

function as_matrix(elem::LieAlgebraElem)
    X = zeros(ComplexF64, size(elem), size(elem))
    for (i, B) in enumerate(elem.algebra.basis)
        X += elem.coeffs[i] * B
    end
    return X
end

function act(elem::LieAlgebraElem, f::Expression, vars::Vector{Variable})
    @assert size(elem) == length(vars)
    X = as_matrix(elem)
    return dot(differentiate(f, vars), -X*vars)
end

function act(elem::LieAlgebraElem, f::Expression; var_groups::Vector{Vector{Variable}})
    X = as_matrix(elem)
    return sum([dot(differentiate(f, vars), -X*vars) for vars in var_groups])
end

function act(elem::LieAlgebraElem, mon::Monomial)
    @assert size(elem) == nvariables(mon)
    X = as_matrix(elem)
    return dot(differentiate(mon), -X*variables(mon))
end

function act(elem::LieAlgebraElem, mon::Monomial, var_groups::Vector{Vector{Variable}})
    X = as_matrix(elem)
    return sum([dot(differentiate(mon, vars), -X*vars) for vars in var_groups])
end

function as_matrix(elem::LieAlgebraElem, B::MonomialBasis)
    @assert size(elem) == nvariables(B)
    M = zeros(ComplexF64, length(B), length(B))
    for (i, mon) in enumerate(B)
        gMon = act(elem, mon)
        M[:, i] = coefficients(gMon, B)
    end
    return M
end

function as_matrix(elem::LieAlgebraElem, B::MonomialBasis, var_groups::Vector{Vector{Variable}})
    M = zeros(ComplexF64, length(B), length(B))
    for (i, mon) in enumerate(B)
        gMon = act(elem, mon, var_groups)
        M[:, i] = coefficients(gMon, B)
    end
    return M
end

# symmetric power of vectors in vs
function ⊙(
    vs::Vector{<:AbstractVector{T}},
    tensor_basis::Vector{Vector{S}}
) where {T<:Number, S<:Integer}
    @assert !isempty(vs)
    @assert binomial(length(vs[1])+length(vs)-1, length(vs)) == length(tensor_basis)
    sp = zeros(T, length(tensor_basis))
    for (i, tnsr) in enumerate(tensor_basis)
        sym_basis = multiset_permutations(tnsr, length(tnsr))
        for inds in sym_basis
            sp[i] += prod([vs[j][ind] for (j, ind) in enumerate(inds)])
        end
    end
    return sp
end

function mexps_to_tensor_basis(mexps::Vector{<:SparseVector})
    return [vcat([fill(i, v) for (i, v) in zip(mexp.nzind, mexp.nzval)]...) for mexp in mexps]
end

⊙(
    vs::Vector{<:AbstractVector{T}} where T <: Number,
    mexps::Vector{<:SparseVector}
) = ⊙(vs, mexps_to_tensor_basis(mexps))

⊙(
    Ms::Vector{<:AbstractMatrix{T}} where T <: Number,
    mexps::Vector{<:SparseVector}
) = ⊙(vcat([M2VV(M; copy=false) for M in Ms]...), mexps)

function ⊙(
    Ms::Vector{Matrix{T}},
    multiplicities::Vector{Int}, # if [d₁, d₂], then take d₁ columns from M₁ and d₂ from M₂
    mexps::Vector{<:SparseVector}
) where {T<:Number}
    @assert length(Ms) == length(multiplicities)
    combs = [with_replacement_combinations(1:size(Ms[i], 2), multiplicities[i]) for i in 1:length(Ms)]
    prod_combs = Base.Iterators.product(combs...)
    sym_M = zeros(T, length(mexps), length(prod_combs))
    for (i, prod_comb) in enumerate(prod_combs)
        sym_M[:,i] = ⊙([view(Ms[j],:,col_ids) for (j, col_ids) in enumerate(prod_comb)], mexps)
    end
    return sym_M
end

function sym_weight_structure(alg::LieAlgebra, d::Integer, mexps::Vector{<:SparseVector})
    d = Int(d)
    d == 0 && return WeightStructure([0], [[1;;]])
    d == 1 && return weight_structure(alg)
    combs = multiexponents(; degree=d, nvars=nweights(alg))
    new_weights_dict = Dict{Vector{Int}, Vector{typeof(combs[1])}}()
    for comb in combs
        w = sum([comb.nzval[i]*weights(alg, comb.nzind[i]) for i in 1:length(comb.nzind)])
        val = get(new_weights_dict, w, nothing)
        if isnothing(val)
            new_weights_dict[w] = [comb]
        else
            push!(new_weights_dict[w], comb)
        end
    end
    new_weights = [zeros(Int, 0) for _ in 1:length(new_weights_dict)]
    new_weight_spaces = [zeros(ComplexF64, 0, 0) for _ in 1:length(new_weights_dict)]
    for (i, (weight, combs)) in enumerate(new_weights_dict)
        new_weights[i] = weight
        Ms = [⊙(weight_spaces(alg, comb.nzind), comb.nzval, mexps) for comb in combs]
        new_weight_spaces[i] = hcat(Ms...)
    end
    return WeightStructure(new_weights, new_weight_spaces)
end

# For vs = [v₁,...,vₖ], computes v₁ ⊗ ... ⊗ vₖ
function ⊗(
    vs::Vector{<:AbstractVector{T}},
    tensor_basis # TODO: add type
) where {T <: Number}
    v = zeros(T, length(tensor_basis))
    for (i, tensor) in enumerate(tensor_basis)
        v[i] = prod([vs[j][k] for (j,k) in enumerate(tensor)])
    end
    return v
end

# TODO: compare timings with kron
# For Ms = [M₁,...,Mₖ], computes tensor products v₁ ⊗ ... ⊗ vₖ of all possible combinations of columns vᵢ ∈ Mᵢ
function ⊗(
    Ms::Vector{<:AbstractMatrix{T}},
    tensor_basis # TODO: add type
) where {T <: Number}
    M = zeros(T, prod([size(M,1) for M in Ms]), prod([size(M,2) for M in Ms]))
    combs = Base.Iterators.product([axes(M,2) for M in Ms]...)
    for (i, comb) in enumerate(combs)
        M[:,i] = ⊗([view(Ms[j],:,col_id) for (j, col_id) in enumerate(comb)], tensor_basis)
    end
    return M
end

function tensor_weight_structure(
    v_ws::Vector{WeightStructure},
    tensor_basis # TODO: add type
)
    length(v_ws) == 1 && return v_ws[1]
    combs = Base.Iterators.product([1:nweights(ws) for ws in v_ws]...)
    new_weights_dict = Dict{Vector{Int}, Vector{NTuple{length(v_ws), Int}}}()
    for comb in combs
        w = sum([weights(v_ws[i], j) for (i, j) in enumerate(comb)])
        val = get(new_weights_dict, w, nothing)
        if isnothing(val)
            new_weights_dict[w] = [comb]
        else
            push!(new_weights_dict[w], comb)
        end
    end
    new_weights = [zeros(Int, 0) for _ in 1:length(new_weights_dict)]
    new_weight_spaces = [zeros(ComplexF64, 0, 0) for _ in 1:length(new_weights_dict)]
    for (i, (weight, combs)) in enumerate(new_weights_dict)
        new_weights[i] = weight
        Ms = [⊗([weight_spaces(v_ws[i], j) for (i, j) in enumerate(comb)], tensor_basis) for comb in combs]
        new_weight_spaces[i] = hcat(Ms...)
    end
    return WeightStructure(new_weights, new_weight_spaces)
end

tensor_weight_structure(
    alg::LieAlgebra,
    ds::Vector{<:Integer},
    v_mexps::Vector{<:Vector{<:SparseVector}},
    tensor_basis # TODO: add type
) = tensor_weight_structure([sym_weight_structure(alg, d, mexps) for (d, mexps) in zip(ds, v_mexps)], tensor_basis)

struct IrreducibleLieAlgebraRepresentation
    alg::LieAlgebra
    var_groups::Vector{Vector{Variable}}
    hw_module::HighestWeightModule
end

algebra(π::IrreducibleLieAlgebraRepresentation) = π.alg
space_basis(π::IrreducibleLieAlgebraRepresentation) = basis(π.hw_module)
highest_weight(π::IrreducibleLieAlgebraRepresentation) = highest_weight(π.hw_module)
highest_weight_vector(π::IrreducibleLieAlgebraRepresentation) = highest_weight_vector(π.hw_module)
dim(π::IrreducibleLieAlgebraRepresentation) = 2*highest_weight(π)[1]+1 # TODO: works only for so(3)

function to_expressions(π::IrreducibleLieAlgebraRepresentation; tol::Float64=1e-5)
    exprs = Expression[]
    # TODO: extend the following to multiple negative roots
    πJ₋ = as_matrix(negative_roots(π.alg)[1], space_basis(π), π.var_groups)
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

struct IsotypicComponent
    alg::LieAlgebra
    weight::Vector{Int}
    irreds::Vector{IrreducibleLieAlgebraRepresentation} # TODO: Vector{HighestWeightModule} ?
end

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

Base.:(==)(
    V::PolynomialVectorSpace,
    W::PolynomialVectorSpace
) = (V.vars == W.vars) && (V.degree == W.degree) && (V.upto == W.upto)

function Base.show(io::IO, V::PolynomialVectorSpace)
    println(io, "PolynomialVectorSpace of dimension $(dim(V))")
    print(io, " variables: $(variables(V))")
end

struct LieAlgebraRepresentation
    alg::LieAlgebra
    V::PolynomialVectorSpace
    var_groups::Vector{Vector{Variable}} # defines an action of a Lie algebra
    isotypic::Vector{IsotypicComponent}
end

function LieAlgebraRepresentation(
    alg::LieAlgebra,
    V::PolynomialVectorSpace,
    var_groups::Vector{Vector{Variable}},
    irreds::Vector{IrreducibleLieAlgebraRepresentation}
)
    iso_dict = Dict{Vector{Int}, Vector{IrreducibleLieAlgebraRepresentation}}()
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
    degrees::Vector{<:Integer}
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
    return [IrreducibleLieAlgebraRepresentation(alg, var_groups, hwm) for hwm in hw_modules]
end

# TODO: supposes that all the vars in V occur in var_groups
function LieAlgebraRepresentation(
    alg::LieAlgebra,
    V::PolynomialVectorSpace;
    action::Vector{Vector{Variable}}
)
    @assert issetequal(vcat(action...), variables(V))
    groups_mexps = multiexponents(degree=Int8(degree(V)), nvars=Int16(length(action)), upto=is_upto(V))
    irreds = IrreducibleLieAlgebraRepresentation[]
    for mexp in groups_mexps
        wm = weight_module(alg, action, Vector(mexp))
        append!(irreds, to_irreducible(alg, action, wm))
    end
    return LieAlgebraRepresentation(alg, action, irreds)
end

Base.randn(alg::LieAlgebra) = LieAlgebraElem(randn(ComplexF64, dim(alg)), alg)

# Having two representations π₁: g₁ → V and π₂: g₂ → V such that there
# actions commute, computes a representation π₁ ⊞ π₂: g₁ ⊕ g₂ → V
function ⊞(π₁::LieAlgebraRepresentation, π₂::LieAlgebraRepresentation)
    @assert space(π₁) == space(π₂)
    
end