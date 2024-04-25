export FixedFreeVariables,
    FixedFreeSamples,
    fixed,
    free,
    nfixed,
    nfree,
    nsamples


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
FixedFreeVariables(free::Union{Variable, AbstractArray}) = FixedFreeVariables([], free)

function Base.:(==)(v₁::FixedFreeVariables, v₂::FixedFreeVariables)
    return v₁.fixed == v₂.fixed && v₁.free == v₂.free
end

Base.hash(
    vars::FixedFreeVariables,
    u::UInt64
) = hash(vars.fixed, hash(vars.free, u))

"""
    fixed(vars::FixedFreeVariables) -> Vector{Variable}

Return the fixed variables in `vars`.
"""
fixed(vars::FixedFreeVariables) = vars.fixed

"""
    free(vars::FixedFreeVariables) -> Vector{Variable}

Return the free variables in `vars`.
"""
free(vars::FixedFreeVariables) = vars.free

"""
    nfixed(vars::FixedFreeVariables) -> Int

Return the number of fixed variables in `vars`.
"""
nfixed(vars::FixedFreeVariables) = length(vars.fixed)

"""
    nfree(vars::FixedFreeVariables) -> Int

Return the number of free variables in `vars`.
"""
nfree(vars::FixedFreeVariables) = length(vars.free)

"""
    variables(vars::FixedFreeVariables) -> Vector{Variable}

Return the concatenated vector of fixed and free variables in `vars`.
"""
variables(vars::FixedFreeVariables) = vcat(vars.fixed, vars.free)

"""
    nvariables(vars::FixedFreeVariables) -> Int

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

Data type contains samples from a variety for some [`FixedFreeVariables`](@ref).

# Constructor
```julia
FixedFreeSamples(fixed::Vector{ComplexF64}, free::Matrix{ComplexF64})
```

# Examples
```julia-repl
julia> v = randn(ComplexF64, 2);

julia> M = randn(ComplexF64, 3, 4);

julia> s = FixedFreeSamples(v, M)
FixedFreeSamples with 4 samples
 fixed:
2-element Vector{ComplexF64}:
  -0.7366318944925887 - 0.5245827233156782im
 -0.05897853939946192 + 0.4503548970705814im
 free:
3×4 Matrix{ComplexF64}:
 0.0193034-0.466551im   -0.277376+0.281461im  -0.0982807+1.60069im      1.0759+0.113134im
 -0.677777+0.110846im  -0.0734675-0.812262im   -0.802048+0.643131im  -0.322837+0.54686im
 -0.221659-0.897734im    -1.06199-0.43677im     0.669207+0.868208im   -0.02822-0.733581im
```
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
    fixed(s::FixedFreeSamples) -> Vector{ComplexF64}

Return the fixed samples in `s`.
"""
fixed(s::FixedFreeSamples) = s.fixed

"""
    free(s::FixedFreeSamples) -> Matrix{ComplexF64}

Return the free samples in `s`.
"""
free(s::FixedFreeSamples) = s.free

"""
    nsamples(s::FixedFreeSamples) -> Int

Return the number of samples in `s`.
"""
nsamples(s::FixedFreeSamples) = size(s.free, 2)

function Base.show(io::IO, s::FixedFreeSamples)
    println(io, "FixedFreeSamples with $(nsamples(s)) samples")
    println(io, " fixed:")
    display(fixed(s))
    println(io, " free:")
    display(free(s))
end

function Base.hcat(s₁::FixedFreeSamples, s₂::FixedFreeSamples; tol::Real=1e-10) # TODO: make tol optional?
    @assert norm(fixed(s₁)-fixed(s₂)) < tol
    return FixedFreeSamples(fixed(s₁), hcat(free(s₁), free(s₂)))
end