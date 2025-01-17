export WeightStructure,
    weights,
    nweights,
    weight_spaces,
    WeightVector,
    weight,
    vector,
    WeightModule,
    basis,
    weight_structure,
    HighestWeightModule,
    highest_weight,
    highest_weight_vector


# TODO: change to Vector{WeightSpace}?
struct WeightStructure
    weights::Vector{Vector{Int}}
    weight_spaces::Vector{Matrix{ComplexF64}}
end

WeightStructure(
    weights::Vector{Vector{Int}},
    weight_vectors::Vector{Vector{T}} where T <: Number
) = WeightStructure(weights, [V2M(v) for v in weight_vectors])

WeightStructure() = WeightStructure([], [])

WeightStructure(
    weights::Vector{Int},
    weight_spaces::Vector{Matrix{T}} where T <: Number
) = WeightStructure([[w] for w in weights], weight_spaces)

WeightStructure(
    weights::Vector{Int},
    weight_vectors::Vector{Vector{T}} where T <: Number
) = WeightStructure(weights, [V2M(v) for v in weight_vectors])

function Base.show(io::IO, ws::WeightStructure)
    println(io, "WeightStructure of $(space_dim(ws))-dimensional vector space")
    println(io, " $(nweights(ws)) weights: ", join(weights(ws), ", "))
    println(io, " dimensions of $(nweights(ws)) weight spaces: ", join([size(M,2) for M in weight_spaces(ws)], ", "))
    print(io, " max weight space dimension: ", max([size(M,2) for M in weight_spaces(ws)]...))
end

weight_size(ws::WeightStructure) = length(ws.weights[1])
space_dim(ws::WeightStructure) = size(ws.weight_spaces[1], 1)
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

struct WeightSpace
    weight::Vector{Int}
    space::Matrix{ComplexF64}
end

weight(ws::WeightSpace) = ws.weight
space(ws::WeightSpace) = ws.space

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

function sym_weight_structure(ws::WeightStructure, d::Integer, mexps::Vector{<:SparseVector})
    d = Int(d)
    d == 0 && return WeightStructure([zeros(Int, weight_size(ws))], [[1;;]])
    d == 1 && return ws
    combs = multiexponents(; degree=d, nvars=nweights(ws))
    new_weights_dict = Dict{Vector{Int}, Vector{typeof(combs[1])}}()
    for comb in combs
        w = sum([comb.nzval[i]*weights(ws, comb.nzind[i]) for i in 1:length(comb.nzind)])
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
        Ms = [⊙(weight_spaces(ws, comb.nzind), comb.nzval, mexps) for comb in combs]
        new_weight_spaces[i] = hcat(Ms...)
    end
    return WeightStructure(new_weights, new_weight_spaces)
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
