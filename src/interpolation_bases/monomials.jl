export multiexponents,
    MonomialBasis,
    Monomial,
    multiexponents,
    iterate,
    evaluate,
    to_classes,
    ⊗


SPARSE_VECTOR_LINK = "https://docs.julialang.org/en/v1/stdlib/SparseArrays/#Sparse-Vector-Storage"

"""
    struct MonomialBasis{Tv<:Integer,Ti<:Integer} <: AbstractInterpolationBasis
        mexps::Vector{SparseVector{Tv,Ti}}
        vars::Vector{Variable}
    end

An [`AbstractInterpolationBasis`](@ref) that consists of monomials. Parametric
type `Tv` defines the type of exponents in multiexponents, `Ti` defines the type
of non-zero exponent indicies. See also [`SparseVector`]($(SPARSE_VECTOR_LINK)).

```julia
MonomialBasis{Tv<:Integer, Ti<:Integer}(; variables::Vector{Variable}, degree::Integer)
monomials(variables::Vector{Variable}, degree::Integer)
```

# Examples
```julia-repl
julia> @var x y z
(x, y, z)

julia> mons = MonomialBasis{Int8, Int16}(variables=[x,y,z], degree=2)
10-element MonomialBasis{Int8, Int16}
[1, x, y, z, x^2, y^2, z^2, x*y, x*z, y*z]

julia> samples = randn(ComplexF64, 3, 2)
3×2 Matrix{ComplexF64}:
  0.299344-0.238374im  -0.527805-0.360128im
 -0.114638+1.89994im    0.127791-0.846475im
  0.303708+1.24025im   0.0363844-0.264417im

julia> evaluate(mons, samples)
10×2 Matrix{ComplexF64}:
       1.0+0.0im              1.0+0.0im
  0.299344-0.238374im   -0.527805-0.360128im
 -0.114638+1.89994im     0.127791-0.846475im
  0.303708+1.24025im    0.0363844-0.264417im
 0.0327848-0.142712im    0.148886+0.380155im
  -3.59664-0.43561im    -0.700189-0.216343im
  -1.44598+0.753349im  -0.0685926-0.0192413im
  0.418581+0.596063im   -0.372288+0.400752im
  0.386557+0.298866im   -0.114428+0.126458im
  -2.39122+0.434849im   -0.219173-0.0645886im
```
"""
struct MonomialBasis{Tv<:Integer,Ti<:Integer} <: AbstractInterpolationBasis
    mexps::Vector{SparseVector{Tv,Ti}}
    vars::Vector{Variable}
end

# Generates a list of multiexponents of degree @degree in @nvars variables
function multiexponents(degree::Tv, nvars::Ti) where {Tv<:Integer,Ti<:Integer}
    mexps = [spzeros(Tv, Ti, nvars) for _ in 1:num_mons(nvars, degree)]
    iszero(degree) && return mexps
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

function multiexponents(; degree::Tv, nvars::Ti, upto::Bool=false) where {Tv<:Integer,Ti<:Integer}
    !upto && return multiexponents(degree, nvars)
    mexps = [spzeros(Tv, Ti, nvars)]
    for d::Tv in 1:degree
        append!(mexps, multiexponents(d, nvars))
    end
    return mexps
end

function MonomialBasis{Tv,Ti}(
    variables::Vector{Variable};
    degree::Integer=-1,
    upto::Bool=true
) where {Tv<:Integer,Ti<:Integer}
    degree == -1 && return MonomialBasis{Tv,Ti}([], variables)
    degree = convert(Tv, degree)
    mexps = multiexponents(degree=degree, nvars=Ti(length(variables)), upto=upto)
    return MonomialBasis{Tv,Ti}(mexps, variables)
end

function MonomialBasis{Tv,Ti}(
    variables::Vector{Vector{Variable}},
    i::Int;
    degree::Integer=-1,
    upto::Bool=true
) where {Tv<:Integer,Ti<:Integer}
    group_size = length(variables[1])
    vars = vcat(variables...)
    degree == -1 && return MonomialBasis{Tv,Ti}([], vars)
    degree = convert(Tv, degree)
    mexps = multiexponents(degree=degree, nvars=Ti(group_size), upto=upto)
    full_mexps = [spzeros(Tv, Ti, length(vars)) for _ in 1:length(mexps)]
    rng = ((i-1)*group_size+1):(i*group_size)
    for (j, mexp) in enumerate(mexps)
        full_mexps[j][rng] = mexp
    end
    return MonomialBasis{Tv,Ti}(full_mexps, vars)
end

multiexponents(mons::MonomialBasis) = mons.mexps
variables(mons::MonomialBasis) = mons.vars
nvariables(mons::MonomialBasis) = length(variables(mons))
Base.length(mons::MonomialBasis) = length(mons.mexps)
Base.:(==)(
    mons₁::MonomialBasis,
    mons₂::MonomialBasis
) = (mons₁.mexps == mons₂.mexps) && (mons₁.vars == mons₂.vars)

function Base.show(io::IO, mons::MonomialBasis)
    println(io, "$(length(mons))-element $(typeof(mons))")
    print(io, "[", join(to_expressions(mons), ", "), "]")
end

function ⊗(
    mon_bases::Vector{MonomialBasis{Tv,Ti}},
    tensor_basis; # TODO: add type
    equal_vars::Bool=false
) where {Tv<:Integer,Ti<:Integer}
    if equal_vars
        vars = variables(mon_bases[1])
        mexps = [spzeros(Tv, Ti, 0) for _ in 1:length(tensor_basis)]
        for (i, tensor) in enumerate(tensor_basis)
            mexps[i] = sum([multiexponents(mon_bases[j])[k] for (j, k) in enumerate(tensor)])
        end
        return MonomialBasis(mexps, vars)
    else
        return # TODO
    end
end

struct Monomial{Tv<:Integer,Ti<:Integer}
    mexp::SparseVector{Tv,Ti}
    vars::Vector{Variable}
    dict::Dict{Variable, Tv} # for a given var gives its exponent (only nonzero are saved)
end

function Monomial(
    mexp::SparseVector{Tv,Ti},
    vars::Vector{Variable}
) where {Tv<:Integer,Ti<:Integer}
    dict = Dict{Variable, Tv}(zip(vars[mexp.nzind], mexp.nzind))
    return Monomial{Tv, Ti}(mexp, vars, dict)
end

multiexponent(mon::Monomial) = mon.mexp
variables(mon::Monomial) = mon.vars
nvariables(mon::Monomial) = length(mon.vars)
to_expression(mon::Monomial) = prodpow(mon.vars, mon.mexp)

Base.:(==)(m₁::Monomial, m₂::Monomial) = to_expression(m₁) == to_expression(m₂)
Base.hash(m::Monomial, h::UInt) = hash(to_expression(m), h)

function Base.show(io::IO, mon::Monomial)
    println(io, "$(typeof(mon))")
    print(io, "$(to_expression(mon))")
end

to_expression(
    mexp::AbstractSparseVector,
    vars::Vector{Variable}
) = prodpow(vars, mexp)

function HC.differentiate(mon::Monomial, var::Variable)
    idx = get(mon.dict, var, nothing)
    if !isnothing(idx)
        mexp = copy(mon.mexp)
        mexp[idx] -= 1
        return (mon.mexp[idx])*to_expression(mexp, mon.vars)
    else
        return Expression(0)
    end
end

HC.differentiate(mon::Monomial) = [differentiate(mon, var) for var in mon.vars]
HC.differentiate(
    mon::Monomial,
    vars::Vector{Variable}
) = [differentiate(mon, var) for var in vars]


function to_expressions(mons::MonomialBasis)
    return [to_expression(mexp, mons.vars) for mexp in mons.mexps]
end

to_monomials(mons::MonomialBasis) = [mon for mon in mons]

function Base.iterate(mons::MonomialBasis, state=1)
    state > length(mons) && return nothing
    return (Monomial(mons.mexps[state], mons.vars), state+1)
end

function Base.push!(
    mons::MonomialBasis{Tv,Ti},
    mexp::SparseVector{Tv,Ti}
) where {Tv<:Integer,Ti<:Integer}
    push!(mons.mexps, mexp)
end

# TODO: test timings w/ and w/o view
prodpow(
    M::AbstractMatrix,
    e::AbstractSparseVector
) = prod(view(M,e.nzind,:).^e.nzval; dims=1)

HC.evaluate(
    mons::MonomialBasis,
    sample::AbstractVector
) = [prodpow(sample, mexp) for mexp in mons.mexps]

function HC.evaluate(
    mons::MonomialBasis,
    samples::AbstractMatrix{T}
) where {T <: Number}
    evals = zeros(T, length(mons), size(samples, 2))
    for (i, mexp) in enumerate(mons.mexps)
        evals[i, :] = prodpow(samples, mexp)
    end
    return evals
end

Base.getindex(
    mons::MonomialBasis{Tv,Ti},
    inds...
) where {Tv<:Integer,Ti<:Integer} = MonomialBasis{Tv,Ti}(getindex(mons.mexps, inds...), mons.vars)
Base.getindex(mons::MonomialBasis, ind::Int) = Monomial(mons.mexps[ind], mons.vars)

Base.findfirst(
    mexp::Vector{<:Integer},
    mons::MonomialBasis{Tv,Ti}
) where {Tv<:Integer,Ti<:Integer} = findfirst(x -> x==sparse(Tv, Ti, mexp), multiexponents(mons))