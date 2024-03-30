export FixedFreeVariables,
    FixedFreeSamples,
    fixed,
    free,
    nfixed,
    nfree


"""
    FixedFreeVariables

Data type for creating fixed and free variables used for interpolation.

```julia
FixedFreeVariables(fixed::Vector{Variable}, free::Vector{Variable})
FixedFreeVariables(
    fixed::Union{Variable, AbstractArray},
    free::Union{Variable, AbstractArray}
)
FixedFreeVariables(free::Union{Variable, AbstractArray})
```

# Examples
```jldoctest
julia> @var x y[1:2] z[1:2,1:3]
(x, Variable[y₁, y₂], Variable[z₁₋₁ z₁₋₂ z₁₋₃; z₂₋₁ z₂₋₂ z₂₋₃])

julia> FixedFreeVariables(x)
FixedFreeVariables: 0 fixed, 1 free
 fixed: 
 free: x

julia> FixedFreeVariables([x, y])
FixedFreeVariables: 0 fixed, 3 free
 fixed: 
 free: x, y₁, y₂

julia> FixedFreeVariables([x, y], z)
FixedFreeVariables: 3 fixed, 6 free
 fixed: x, y₁, y₂
 free: z₁₋₁, z₂₋₁, z₁₋₂, z₂₋₂, z₁₋₃, z₂₋₃

julia> FixedFreeVariables([x, y], [y, z])
ERROR: Nontrivial intersection of fixed and free variables
[...]

julia> FixedFreeVariables([x, y, z], [])
ERROR: Array of free variables must be nonempty
[...]
```
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
    return FixedFreeVariables(Variable.(collect(flatten(fixed))), Variable.(collect(flatten(free))))
end
FixedFreeVariables(free::Union{Variable, AbstractArray}) = FixedFreeVariables(Variable[], free)

function Base.:(==)(v₁::FixedFreeVariables, v₂::FixedFreeVariables)
    return v₁.fixed == v₂.fixed && v₁.free == v₂.free
end

Base.hash(
    vars::FixedFreeVariables,
    u::UInt64
) = hash(vars.fixed, hash(vars.free, u))

"""
    fixed(vars::FixedFreeVariables)

Return the fixed variables in `vars`.
"""
fixed(vars::FixedFreeVariables) = vars.fixed

"""
    free(vars::FixedFreeVariables)

Return the free variables in `vars`.
"""
free(vars::FixedFreeVariables) = vars.free

"""
    nfixed(vars::FixedFreeVariables)

Return the number of fixed variables in `vars`.
"""
nfixed(vars::FixedFreeVariables) = length(vars.fixed)

"""
    nfree(vars::FixedFreeVariables)

Return the number of free variables in `vars`.
"""
nfree(vars::FixedFreeVariables) = length(vars.free)

"""
    variables(vars::FixedFreeVariables)

Return all the variables in `vars`.
"""
variables(vars::FixedFreeVariables) = vcat(vars.fixed, vars.free)

"""
    nvariables(vars::FixedFreeVariables)

Return the number of all the variables in `vars`.
"""
nvariables(vars::FixedFreeVariables) = nfixed(vars) + nfree(vars)

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

Data type contains samples from a variety of some [`FixedFreeVariables`](@ref).

```julia
FixedFreeSamples(fixed::Vector{ComplexF64}, free::Matrix{ComplexF64})
```

# Examples

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

"""
    fixed(s::FixedFreeSamples)

Return the fixed samples in `s`.
"""
fixed(s::FixedFreeSamples) = s.fixed

"""
    free(s::FixedFreeSamples)

Return the free samples in `s`.
"""
free(s::FixedFreeSamples) = s.free

"""
    nsamples(s::FixedFreeSamples) -> Int

Return the number of samples in `s`.
"""
nsamples(s::FixedFreeSamples) = size(s.free, 2)

function Base.hcat(s₁::FixedFreeSamples, s₂::FixedFreeSamples; tol::Real=1e-10) # TODO: make tol optional?
    @assert norm(fixed(s₁)-fixed(s₂)) < tol
    return FixedFreeSamples(fixed(s₁), hcat(free(s₁), free(s₂)))
end