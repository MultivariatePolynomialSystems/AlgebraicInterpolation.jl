export sparsify!,
    simplify_numbers,
    eye, a2p, p2a, xx,
    num_mons, num_mons_upto,
    Tolerances,
    M2VV, V2M, M2V,
    coefficients

a2p(M::AbstractMatrix{<:Number}) = [M; ones(eltype(M), 1, size(M, 2))]
p2a(M::AbstractMatrix{<:Number}) = (M./M[end:end,:])[1:end-1,:]

M2VV(M::AbstractMatrix; copy::Bool=true) = copy ? [M[:,i] for i in axes(M, 2)] : [view(M,:,i) for i in axes(M,2)]

V2M(v::AbstractVector) = reshape(v, length(v), 1) # TODO: do I need this? I can access a Vector with [i,1] too...
M2V(M::AbstractMatrix) = M[:,1] # TODO: same here, I can access elements of 1-column matrix M by M[i] too...

function Base.copyto!(M::AbstractMatrix{T}, v::AbstractVector{AbstractVector{T}}; dim::Integer) where {T}
    for i in eachindex(v)
        copyto!(selectdim(M, dim, i), v[i])
    end
end

xx(v) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
xx2v(xx) = [-xx[2,3], xx[1,3], -xx[1,2]]
eye(T, n::Integer) = Matrix{T}(I(n))

prodpow(v::AbstractVector, e::AbstractSparseVector) = prod(v[e.nzind].^e.nzval)

function num_mons(n::Integer, d::Integer; upto::Bool=false)
    upto && return n > 0 ? binomial(Int(n + d), Int(d)) : 0
    return n > 0 ? binomial(Int(n - 1 + d), Int(d)) : 0
end

# TODO: test this
function sparsify!(v::AbstractVector{<:Number}, tol::Real; digits::Integer=0)
    for j in eachindex(v)
        if abs(imag(v[j])) < tol
            v[j] = real(v[j])
        elseif abs(round(imag(v[j]); digits=digits) - imag(v[j])) < tol
            v[j] = real(v[j]) + round(imag(v[j]); digits=digits)*im
        end
        if abs(real(v[j])) < tol
            v[j] = imag(v[j])*im
        elseif abs(round(real(v[j]); digits=digits) - real(v[j])) < tol
            v[j] = round(real(v[j]); digits=digits) + imag(v[j])*im
        end
    end
end

function sparsify!(M::AbstractMatrix{<:Number}, tol::Real; digits::Integer=0)
    for r in eachrow(M)
        sparsify!(r, tol; digits=digits)
    end
end

function simplify_numbers(v::AbstractVector{<:Number})
    v = Vector{Number}(v)
    for (i, vᵢ) in enumerate(v)
        try
            v[i] = Integer(vᵢ)
        catch
            try
                v[i] = Real(vᵢ)
            catch
                try
                    v[i] = Complex{Integer}(vᵢ)
                catch
                end
            end
        end
    end
    return v
end

function div_by_lowest_magnitude(v::AbstractVector{<:Number}, tol::Float64)
    j = 1
    while norm(v[j]) < tol
        j += 1
        if j > length(v)
            return v
        end
    end
    a = v[j]
    for vᵢ in v
        if tol < norm(vᵢ) < norm(a)
            a = vᵢ
        end
    end
    return v./a
end

function to_ordinal(n::Integer)::String
    if mod(n, 10) == 1
        mod(n, 100) == 11 && return "$(n)th"
        return "$(n)st"
    end
    if mod(n, 10) == 2
        mod(n, 100) == 12 && return "$(n)th"
        return "$(n)nd"
    end
    if mod(n, 10) == 3
        mod(n, 100) == 13 && return "$(n)th"
        return "$(n)rd"
    end
    return "$(n)th"
end

function subscript(n::Integer)::String
    c = n < 0 ? [Char(0x208B)] : []
    for d in reverse(digits(abs(n)))
        push!(c, Char(0x2080+d))
    end
    return join(c)
end

function superscript(n::Integer)::String
    c = n < 0 ? [Char(0x207B)] : []
    for d in reverse(digits(abs(n)))
        if d == 0 push!(c, Char(0x2070)) end
        if d == 1 push!(c, Char(0x00B9)) end
        if d == 2 push!(c, Char(0x00B2)) end
        if d == 3 push!(c, Char(0x00B3)) end
        if d > 3 push!(c, Char(0x2070+d)) end
    end
    return join(c)
end

# TODO: eachcol?
Base.findfirst(
    v::AbstractVector{<:Number},
    M::AbstractMatrix{<:Number};
    tol::Real=1e-5
) = findfirst(i->norm(M[:,i]-v)<tol, axes(M,2))

take_rows(
    f::Function,
    M::AbstractMatrix{T}
) where {T} = M[[f(r) for r in eachrow(M)], :]

function column_diffs(M::AbstractMatrix{T}) where {T<:Number}
    M = M - M[:,1]*ones(T, 1, size(M,2))
    return M[:,2:end]
end

phrase(i::Integer, word::String) = i == 1 ? "$(i) $(word)" : "$(i) $(word)s"

function rand_unit(T::Type, dims...)
    A = randn(T, dims...)
    return A./norm.(A)
end

# TODO: add general_tol for everything else?
@kwdef struct Tolerances
    common_tol::Float64=1e-10
    nullspace_atol::Float64=0
    nullspace_rtol::Float64=0
    rank_atol::Float64=0
    rank_rtol::Float64=0
    rref_tol::Float64=1e-5
    sparsify_tol::Float64=1e-5
end

function coefficients(f::Expression, mons::MonomialBasis)
    if iszero(f)
        return zeros(ComplexF64, length(mons))
    end
    dict = to_dict(expand(f), variables(mons))
    coeffs = zeros(ComplexF64, length(mons))
    for (mexp, coeff) in dict
        idx = findfirst(mexp, mons)
        if !isnothing(idx)
            coeffs[idx] = coeff
        else
            error("The given expression is not representable in the given basis")
        end
    end
    return coeffs
end

SparseArrays.sparse(
    Tv::Type{<:Integer},
    Ti::Type{<:Integer},
    v::AbstractVector{<:Integer}
) = SparseVector{Tv, Ti}(sparse(v))