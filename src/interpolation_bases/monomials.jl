export multiexponents,
    MonomialVector,
    iterate,
    evaluate,
    to_classes

"""
    MonomialVector{Tv<:Integer,Ti<:Integer} <: AbstractMonomialVector
"""
mutable struct MonomialVector{Tv<:Integer,Ti<:Integer} <: AbstractMonomialVector
    mexps::Vector{SparseVector{Tv,Ti}}
    vars::Vector{Variable}
end

# Generates list of multiexponents of degree @degree in @nvars variables
function multiexponents(; degree::Tv, nvars::Ti) where {Tv<:Integer,Ti<:Integer}
    mexps = [spzeros(Tv, Ti, nvars) for _ in 1:num_mons(nvars, degree)]
    k = 1
    for n in 1:nvars
        for part::Vector{Tv} in partitions(degree, n)
            for vals in multiset_permutations(part, n)
                for inds in combinations(Ti.(1:nvars), n)
                    mexps[k][inds] = vals
                    k += 1
                end
            end
        end
    end
    return mexps
end

function MonomialVector{Tv,Ti}(
    vars::Vector{Variable};
    degree::Integer=-1
) where {Tv<:Integer,Ti<:Integer}
    degree == -1 && return MonomialVector{Tv,Ti}([], vars)
    degree = convert(Tv, degree)
    mexps = [spzeros(Tv, Ti, length(vars))]
    for d::Tv in 1:degree
        append!(mexps, multiexponents(degree=d, nvars=Ti(length(vars))))
    end
    return MonomialVector{Tv,Ti}(mexps, vars)
end

variables(mons::MonomialVector) = mons.vars
nvariables(mons::MonomialVector) = length(variables(mons))
Base.length(mons::MonomialVector) = length(mons.mexps)

to_expression(
    mexp::AbstractSparseVector,
    vars::Vector{Variable}
) = prodpow(vars, mexp)

function to_expressions(mons::MonomialVector)
    return [to_expression(mexp, mons.vars) for mexp in mons]
end

function Base.show(io::IO, mons::MonomialVector)
    println(io, "$(length(mons))-element $(typeof(mons))")
    print(io, "[", join(to_expressions(mons), ", "), "]")
end

function Base.iterate(mons::MonomialVector, state=1)
    state > length(mons) && return nothing
    return (mons.mexps[state], state+1)
end

function Base.push!(
    mons::MonomialVector{Tv,Ti},
    mexp::SparseVector{Tv,Ti}
) where {Tv<:Integer,Ti<:Integer}
    push!(mons.mexps, mexp)
end

# TODO: test timings w/ and w/o view
prodpow(
    M::AbstractMatrix,
    e::AbstractSparseVector
) = prod(view(M,e.nzind,:).^e.nzval; dims=1)

function HC.evaluate(
    mons::MonomialVector,
    samples::Matrix{ComplexF64}
)
    evals = zeros(ComplexF64, length(mons), size(samples, 2))
    for (i, mexp) in enumerate(mons)
        evals[i, :] = prodpow(samples, mexp)
    end
    return evals
end
