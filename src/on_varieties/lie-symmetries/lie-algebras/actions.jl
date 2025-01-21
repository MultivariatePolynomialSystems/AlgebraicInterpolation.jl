export AbstractLieAlgebraAction,
    LieAlgebraAction,
    ScalingLieAction,
    SumLieAlgebraAction,
    var_groups,
    hw_spaces


abstract type AbstractLieAlgebraAction end


struct ScalingLieAction <: AbstractLieAlgebraAction
    alg::ScalingLieAlgebra
    vars::Vector{Variable}
end

function ScalingLieAction(action_vars::Vector{Variable}, all_vars::Vector{Variable})
    ids = [findfirst(isequal(var), all_vars) for var in action_vars]
    exps = zeros(Int, 1, length(all_vars))
    exps[1, ids] .= 1
    return ScalingLieAction(ScalingLieAlgebra(exps), all_vars)
end

function ScalingLieAction(action_vars::AbstractArray; variables::AbstractArray=action_vars)
    act_vars = Variable.(collect(flatten(action_vars)))
    all_vars = Variable.(collect(flatten(variables)))
    return ScalingLieAction(act_vars, all_vars)
end

algebra(a::ScalingLieAction) = a.alg
variables(a::ScalingLieAction) = a.vars
exponents(a::ScalingLieAction) = a.alg.exps

function Base.show(io::IO, a::ScalingLieAction)
    println(io, "ScalingLieAction of $(name(algebra(a)))")
    print(io, " action:")
    U = exponents(a)
    if size(U, 1) == 1
        @var λ
        λ = [λ]
    else
        @var λ[1:size(U₀, 1)]
    end
    action = Vector{Vector{Tuple{Variable, Expression}}}([[] for _ in axes(U, 1)])
    vars = variables(a)
    for j in axes(U, 1)
        nzind, nzval = findnz(U[j, :])
        exprs = (λ[j].^nzval).*vars[nzind]
        action[j] = collect(zip(vars[nzind], exprs))
    end
    for free_action in action
        print(io, "\n  ")
        for (j, (var, expr)) in enumerate(free_action)
            print(io, var, " ↦ ", expr)
            j < length(free_action) && print(io, ", ")
        end
    end
end

function act(elem::LieAlgebraElem{ScalingLieAlgebra}, f::Union{Expression, Monomial}, action::ScalingLieAction)
    X = as_matrix(elem)
    return expand(dot(differentiate(f, variables(action)), -X*variables(action)))
end

function inv_weight_space(a::AbstractLieAlgebraAction, vars::Vector{Variable})
    inv_vars = setdiff(vars, variables(a))
    isempty(inv_vars) && return nothing
    ws = WeightSpace(zeros(Int, rank(algebra(a))), zeros(ComplexF64, length(vars), length(inv_vars)))
    for (i, var) in enumerate(inv_vars)
        ws[vars_dict[var], i] = 1
    end
    return ws
end

function weight_structure(a::ScalingLieAction, vars::Vector{Variable})
    @assert variables(a) ⊆ vars
    vars_dict = Dict(zip(vars, 1:length(vars)))
    hw_struct = WeightStructure()
    for (i, var) in enumerate(variables(a))
        weight = exponents(a)[:, i]
        hwv = WeightVector(weight, zeros(ComplexF64, length(vars)))
        hwv[vars_dict[var]] = 1
        push!(hw_struct, hwv)
    end
    inv_ws = inv_weight_space(a, vars)
    !isnothing(inv_ws) && push!(hw_struct, inv_ws)
    return hw_struct
end

hw_spaces(a::ScalingLieAction, vars::Vector{Variable}) = weight_structure(a, vars)


struct LieAlgebraAction <: AbstractLieAlgebraAction
    alg::LieAlgebra
    var_groups::Vector{Vector{Variable}}
end

LieAlgebraAction(
    alg::LieAlgebra,
    var_groups::AbstractVecOrMat{Variable}
) = LieAlgebraAction(alg, M2VV(hcat(var_groups)))

LieAlgebraAction(
    alg::LieAlgebra,
    var_groups::AbstractArray
) = LieAlgebraAction(alg, M2VV(hcat(var_groups...)))

algebra(a::LieAlgebraAction) = a.alg
var_groups(a::LieAlgebraAction) = a.var_groups
nvar_groups(a::LieAlgebraAction) = length(var_groups(a))
variables(a::LieAlgebraAction) = vcat(var_groups(a)...)
weights(a::LieAlgebraAction) = weights(a.alg)
weights(a::LieAlgebraAction, inds...) = weights(a.alg, inds...)
nweights(a::LieAlgebraAction) = nweights(a.alg)

function Base.show(io::IO, a::LieAlgebraAction)
    println(io, "LieAlgebraAction of $(name(algebra(a)))")
    print(io, " action: [", join([join(vars, ", ") for vars in var_groups(a)], "], ["), "]")
end

function act(elem::LieAlgebraElem{LieAlgebra}, f::Union{Expression, Monomial}, action::LieAlgebraAction)
    X = as_matrix(elem)
    return expand(sum([dot(differentiate(f, vars), -X*vars) for vars in var_groups(action)]))
end

function weight_structure(a::LieAlgebraAction, vars::Vector{Variable}; as_hw_spaces::Bool=false)
    @assert variables(a) ⊆ vars
    vars_dict = Dict(zip(vars, 1:length(vars)))
    ws = WeightStructure()
    alg_ws = as_hw_spaces ? hw_spaces(algebra(a)) : weight_structure(algebra(a))
    for w_space in alg_ws
        new_w_space = WeightSpace(weight(w_space), zeros(ComplexF64, length(vars), dim(w_space)*nvar_groups(a)))
        for (i, vars_group) in enumerate(var_groups(a))
            for (j, var) in enumerate(vars_group)
                new_w_space[vars_dict[var], (i-1)*dim(w_space)+1:i*dim(w_space)] = w_space[j, :]
            end
        end
        push!(ws, new_w_space)
    end
    inv_ws = inv_weight_space(a, vars)
    !isnothing(inv_ws) && push!(ws, inv_ws)
    return ws
end

hw_spaces(a::LieAlgebraAction, vars::Vector{Variable}) = weight_structure(a, vars, as_hw_spaces=true)


struct SumLieAlgebraAction <: AbstractLieAlgebraAction
    alg::SumLieAlgebra
    actions::Vector{AbstractLieAlgebraAction}
end

algebra(a::SumLieAlgebraAction) = a.alg
actions(a::SumLieAlgebraAction) = a.actions
nsummands(a::SumLieAlgebraAction) = length(actions(a))
Base.getindex(a::SumLieAlgebraAction, i::Int) = actions(a)[i]

function Base.show(io::IO, a::SumLieAlgebraAction)
    println(io, "SumLieAlgebraAction of $(name(algebra(a)))")
    # for (i, a) in enumerate(actions(g))
    #     print(io, " action of $(name(algebra(a))): [")
    #     print(io, join([join(vars, ", ") for vars in var_groups(a)], "], ["), "]")
    #     i < nsummands(g) && print(io, "\n")
    # end
    print(io, " action: ")
end

# TODO
function are_commutative(a₁::AbstractLieAlgebraAction, a₂::AbstractLieAlgebraAction)
    return true
end

function ⊕(a₁::AbstractLieAlgebraAction, a₂::AbstractLieAlgebraAction)
    @assert are_commutative(a₁, a₂)
    alg = algebra(a₁) ⊕ algebra(a₂)
    return SumLieAlgebraAction(alg, [a₁, a₂])
end

function ⊕(a₁::SumLieAlgebraAction, a₂::AbstractLieAlgebraAction)
    @assert are_commutative(a₁, a₂)
    alg = algebra(a₁) ⊕ algebra(a₂)
    return SumLieAlgebraAction(alg, [actions(a₁)..., a₂])
end

function act(elem::SumLieAlgebraElem, f::Union{Expression, Monomial}, action::SumLieAlgebraAction)
    return expand(sum([act(elem[i], f, action[i]) for i in 1:nsummands(action)]))
end

hw_spaces(
    a::SumLieAlgebraAction,
    vars::Vector{Variable}
) = ∩([hw_spaces(action, vars) for action in actions(a)])

function weight_structure(a::AbstractLieAlgebraAction, vars::Vector{Variable})
    hws = hw_spaces(a, vars)
    mons = MonomialBasis{Int8,Int16}(vars, degree=1, upto=false)
    ws = WeightStructure()
    for hw_space in hws
        for hwv in hw_space
            wexpr = WeightExpression(hwv, mons)
            orb = orbit(wexpr, a)
            [push!(ws, WeightVector(we)) for we in orb]
        end
    end
    return ws
end

function as_matrix(elem::AbstractLieAlgebraElem, B::MonomialBasis, action::AbstractLieAlgebraAction)
    @assert algebra(elem) == algebra(action)
    M = zeros(ComplexF64, length(B), length(B))
    for (i, mon) in enumerate(B)
        gMon = act(elem, mon, action)
        M[:, i] = coefficients(gMon, B)
    end
    return M
end