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
