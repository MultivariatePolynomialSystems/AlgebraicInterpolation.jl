export FixedFreeVariables,
    FixedFreeSamples,
    fixed,
    free,
    nfixed,
    nfree


"""
    FixedFreeVariables
"""
struct FixedFreeVariables
    fixed::Vector{Variable}
    free::Vector{Variable}

    function FixedFreeVariables(fixed::Vector{Variable}, free::Vector{Variable})
        if !isempty(fixed ∩ free)
            error("Nontrivial intersection of fixed and free variables")
        end
        if isempty(free)
            error("Array of free variables must be nonempty")
        end
        return new(fixed, free)
    end
end

function FixedFreeVariables(
    fixed::Union{Variable, AbstractArray},
    free::Union{Variable, AbstractArray}
)
    return FixedFreeVariables(collect(flatten(fixed)), collect(flatten(free)))
end
FixedFreeVariables(free::Union{Variable, AbstractArray}) = FixedFreeVariables(Variable[], free)

function Base.:(==)(v₁::FixedFreeVariables, v₂::FixedFreeVariables)
    return v₁.fixed == v₂.fixed && v₁.free == v₂.free
end

Base.hash(
    vars::FixedFreeVariables,
    u::UInt64
) = hash(vars.fixed, hash(vars.free, u))

fixed(vars::FixedFreeVariables) = vars.fixed
free(vars::FixedFreeVariables) = vars.free
nfixed(vars::FixedFreeVariables) = length(vars.fixed)
nfree(vars::FixedFreeVariables) = length(vars.free)
variables(vars::FixedFreeVariables) = vcat(vars.fixed, vars.free)

function Base.show(io::IO, vars::FixedFreeVariables)
    println(
        io,
        "FixedFreeVariables: ",
        "$(nfixed(vars)) fixed, ",
        "$(nfree(vars)) free"
    )
    println(io, " fixed: ", join(fixed(vars), ", "))
    print(io, " free: ", join(free(vars), ", "))
end


"""
    FixedFreeSamples
"""
struct FixedFreeSamples
    fixed::Vector{ComplexF64}
    free::Matrix{ComplexF64}
end

FixedFreeSamples(
    nfixed::Int,
    nfree::Int,
    nsamples::Int
) = FixedFreeSamples(zeros(ComplexF64, nfixed), zeros(ComplexF64, nfree, nsamples))

fixed(s::FixedFreeSamples) = s.fixed
free(s::FixedFreeSamples) = s.free
nsamples(s::FixedFreeSamples) = size(s.free, 2)

function Base.hcat(s₁::FixedFreeSamples, s₂::FixedFreeSamples; tol::Real=1e-10) # TODO: make tol optional?
    @assert norm(fixed(s₁)-fixed(s₂)) < tol
    return FixedFreeSamples(fixed(s₁), hcat(free(s₁), free(s₂)))
end