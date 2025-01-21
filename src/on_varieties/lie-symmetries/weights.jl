export WeightStructure,
    weights,
    nweights,
    weight_spaces,
    WeightVector,
    weight,
    vector,
    WeightSpace,
    WeightModule,
    basis,
    weight_structure,
    HighestWeightModule,
    highest_weight,
    highest_weight_vector,
    sum_weight_structure


struct WeightVector
    weight::Vector{Int}
    vector::Vector{ComplexF64}
end

weight(wv::WeightVector) = wv.weight
vector(wv::WeightVector) = wv.vector
Base.setindex!(wv::WeightVector, n::Number, i::Integer) = wv.vector[i] = n

mutable struct WeightSpace
    weight::Vector{Int}
    matrix::Matrix{ComplexF64}
end

WeightSpace(wv::WeightVector) = WeightSpace(weight(wv), V2M(vector(wv)))
WeightSpace(weight::AbstractVector, v::AbstractVector) = WeightSpace(weight, V2M(v))

weight(ws::WeightSpace) = ws.weight
weight_size(ws::WeightSpace) = length(weight(ws))
matrix(ws::WeightSpace) = ws.matrix
Base.size(ws::WeightSpace) = size(ws.matrix, 1)
dim(ws::WeightSpace) = size(ws.matrix, 2)
Base.getindex(ws::WeightSpace, inds...) = ws.matrix[inds...]
Base.setindex!(ws::WeightSpace, M::AbstractVecOrMat{T}, inds...) where {T <: Number} = ws.matrix[inds...] = M

function Base.hcat(ws1::WeightSpace, ws2::WeightSpace)
    if weight(ws1) != weight(ws2)
        throw(ArgumentError("Weights of the two spaces must be equal"))
    end
    if size(ws1) != size(ws2)
        throw(ArgumentError("Sizes of the two spaces must be equal"))
    end
    return WeightSpace(weight(ws1), hcat(matrix(ws1), matrix(ws2)))
end

function Base.show(io::IO, ::MIME"text/plain", ws::WeightSpace)
    println(io, "WeightSpace of dimension $(dim(ws))")
    print(io, " weight: $(weight(ws))")
end

function Base.show(io::IO, ws::WeightSpace)
    print(io, "WeightSpace of dimension $(dim(ws))")
end

function Base.iterate(ws::WeightSpace, state=1)
    state > dim(ws) && return nothing
    return (WeightVector(weight(ws), matrix(ws)[:, state]), state+1)
end

function Base.:∩(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T <: Number}
    N = nullspace(hcat(A, B))
    N = hcat([div_by_lowest_magnitude(N[:,i], 1e-8) for i in 1:size(N, 2)]...)
    sparsify!(N, 1e-8)
    return A*N[1:size(A,2), :]
end

function Base.:∩(ws1::WeightSpace, ws2::WeightSpace)
    new_weight = vcat(weight(ws1), weight(ws2))
    new_matrix = ∩(matrix(ws1), matrix(ws2))
    iszero(size(new_matrix, 2)) && return nothing
    return WeightSpace(new_weight, new_matrix)
end


struct WeightStructure
    weights::Vector{Vector{Int}}
    dict::Dict{Vector{Int}, WeightSpace} # weight => space
end

WeightStructure(
    weights::Vector{Vector{Int}},
    weight_spaces::Vector{Matrix{T}} where T <: Number
) = WeightStructure(
    weights,
    Dict(zip(weights, [WeightSpace(w, M) for (w, M) in zip(weights, weight_spaces)]))
)

WeightStructure(
    weights::Vector{Vector{Int}},
    weight_vectors::Vector{Vector{T}} where T <: Number
) = WeightStructure(weights, [V2M(v) for v in weight_vectors])

WeightStructure() = WeightStructure([], Dict())

WeightStructure(
    weights::Vector{Int},
    weight_spaces::Vector{Matrix{T}} where T <: Number
) = WeightStructure([[w] for w in weights], weight_spaces)

WeightStructure(
    weights::Vector{Int},
    weight_vectors::Vector{Vector{T}} where T <: Number
) = WeightStructure(weights, [V2M(v) for v in weight_vectors])

function WeightStructure(w_spaces::Vector{WeightSpace})
    ws = WeightStructure()
    for w_space in w_spaces
        push!(ws, w_space)
    end
    return ws
end

weight_spaces(
    ws::WeightStructure;
    as_matrices::Bool=false
) = as_matrices ? [matrix(ws[w]) for w in weights(ws)] : [ws[w] for w in weights(ws)]

weight_spaces(
    ws::WeightStructure,
    inds...;
    as_matrices::Bool=false
) = weight_spaces(ws; as_matrices=as_matrices)[inds...]

weights(ws::WeightStructure) = ws.weights
weight_size(ws::WeightStructure) = length(first(weights(ws)))
space_dim(ws::WeightStructure) = size(ws[first(weights(ws))])
weight(ws::WeightStructure, i::Integer) = weights(ws)[i]
nweights(ws::WeightStructure) = length(weights(ws))
dims(ws::WeightStructure) = [dim(M) for M in weight_spaces(ws)]

function Base.show(io::IO, ws::WeightStructure)
    println(io, "WeightStructure of $(space_dim(ws))-dimensional vector space")
    println(io, " $(nweights(ws)) weights: ", join(weights(ws), ", "))
    println(io, " dimensions of $(nweights(ws)) weight spaces: ", join(dims(ws), ", "))
    print(io, " max weight space dimension: ", maximum(dims(ws)))
end

Base.length(ws::WeightStructure) = nweights(ws)
Base.getindex(ws::WeightStructure, weight::Vector{Int}) = ws.dict[weight]
Base.getindex(ws::WeightStructure, i::Integer) = ws[weight(ws, i)]
Base.setindex!(ws::WeightStructure, ws_new::WeightSpace, weight::Vector{Int}) = ws.dict[weight] = ws_new
Base.haskey(ws::WeightStructure, weight::Vector{Int}) = haskey(ws.dict, weight)

function Base.iterate(ws::WeightStructure, state=1)
    state > nweights(ws) && return nothing
    return (ws[state], state+1)
end

function Base.push!(ws::WeightStructure, w_space::WeightSpace)
    if haskey(ws, weight(w_space))
        ws[weight(w_space)] = hcat(ws[weight(w_space)], w_space)
    else
        push!(weights(ws), weight(w_space))
        ws[weight(w_space)] = w_space
    end
    return ws
end

Base.push!(ws::WeightStructure, wv::WeightVector) = push!(ws, WeightSpace(wv))

function Base.:∩(ws1::WeightStructure, ws2::WeightStructure)
    new_spaces = [w_space₁ ∩ w_space₂ for w_space₁ in ws1, w_space₂ in ws2]
    return WeightStructure([ws for ws in new_spaces if !isnothing(ws)])
end

Base.:∩(ws::Vector{WeightStructure}) = reduce(∩, ws)

# function sum_weight_structure(ws::WeightStructure, d::Integer)
#     new_weight_spaces = [zeros(ComplexF64, space_dim(ws)*d, size(M, 2)*d) for M in weight_spaces(ws)]
#     for (M, new_M) in zip(weight_spaces(ws), new_weight_spaces)
#         for j in 1:d
#             new_M[(1:space_dim(ws)) .+ (j-1)*space_dim(ws), (1:size(M,2)) .+ (j-1)*size(M,2)] = M
#         end
#     end
#     return WeightStructure(weights(ws), new_weight_spaces)
# end

function sym_weight_structure(ws::WeightStructure, d::Integer, mexps::Vector{<:SparseVector})
    d = Int(d)
    d == 0 && return WeightStructure([zeros(Int, weight_size(ws))], [[1;;]])
    d == 1 && return ws
    combs = multiexponents(; degree=d, nvars=nweights(ws))
    new_weights_dict = Dict{Vector{Int}, Vector{typeof(combs[1])}}()
    for comb in combs
        w = sum([comb.nzval[i]*weight(ws, comb.nzind[i]) for i in 1:length(comb.nzind)])
        if haskey(new_weights_dict, w)
            push!(new_weights_dict[w], comb)
        else
            new_weights_dict[w] = [comb]
        end
    end
    new_weights = [zeros(Int, 0) for _ in 1:length(new_weights_dict)]
    new_weight_spaces = [zeros(ComplexF64, 0, 0) for _ in 1:length(new_weights_dict)]
    for (i, (weight, combs)) in enumerate(new_weights_dict)
        new_weights[i] = weight
        Ms = [⊙(weight_spaces(ws, comb.nzind; as_matrices=true), comb.nzval, mexps) for comb in combs]
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


# TODO: change to SparseVector
struct WeightExpression
    coeffs::Vector{ComplexF64}
    basis::MonomialBasis
    weight::Vector{Int}
end

WeightExpression(
    expr::Expression,
    basis::MonomialBasis,
    weight::Vector{Int}
) = WeightExpression(coefficients(expr, basis), basis, weight)

WeightExpression(
    v::WeightVector,
    basis::MonomialBasis
) = WeightExpression(vector(v), basis, weight(v))

coefficients(we::WeightExpression) = we.coeffs
basis(we::WeightExpression) = we.basis
weight(we::WeightExpression) = we.weight
expression(we::WeightExpression) = dot(coefficients(we), to_expressions(basis(we)))
WeightVector(we::WeightExpression) = WeightVector(weight(we), coefficients(we))