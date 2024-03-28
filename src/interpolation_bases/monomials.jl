export multiexponents,
    MonomialBasis,
    iterate,
    evaluate,
    to_classes


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

function MonomialBasis{Tv,Ti}(;
    variables::Vector{Variable},
    degree::Integer=-1
) where {Tv<:Integer,Ti<:Integer}
    degree == -1 && return MonomialBasis{Tv,Ti}([], variables)
    degree = convert(Tv, degree)
    mexps = [spzeros(Tv, Ti, length(variables))]
    for d::Tv in 1:degree
        append!(mexps, multiexponents(degree=d, nvars=Ti(length(variables))))
    end
    return MonomialBasis{Tv,Ti}(mexps, variables)
end

variables(mons::MonomialBasis) = mons.vars
nvariables(mons::MonomialBasis) = length(variables(mons))
Base.length(mons::MonomialBasis) = length(mons.mexps)

to_expression(
    mexp::AbstractSparseVector,
    vars::Vector{Variable}
) = prodpow(vars, mexp)

function to_expressions(mons::MonomialBasis)
    return [to_expression(mexp, mons.vars) for mexp in mons]
end

function Base.show(io::IO, mons::MonomialBasis)
    println(io, "$(length(mons))-element $(typeof(mons))")
    print(io, "[", join(to_expressions(mons), ", "), "]")
end

function Base.iterate(mons::MonomialBasis, state=1)
    state > length(mons) && return nothing
    return (mons.mexps[state], state+1)
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

function HC.evaluate(
    mons::MonomialBasis,
    samples::AbstractMatrix{ComplexF64}
)
    evals = zeros(ComplexF64, length(mons), size(samples, 2))
    for (i, mexp) in enumerate(mons)
        evals[i, :] = prodpow(samples, mexp)
    end
    return evals
end

# TODO: REMOVE LATER
Base.getindex(mons::MonomialBasis, inds...) = getindex(mons.mexps, inds...)