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