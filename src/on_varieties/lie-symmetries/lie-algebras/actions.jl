export AbstractLieAlgebraAction,
    LieAlgebraAction,
    ScalingLieAction,
    SumLieAlgebraAction,
    var_groups


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

algebra(g::ScalingLieAction) = g.alg
variables(g::ScalingLieAction) = g.vars
exponents(g::ScalingLieAction) = g.alg.exps

function Base.show(io::IO, g::ScalingLieAction)
    println(io, "ScalingLieAction of $(name(algebra(g)))")
    print(io, " action:")
    U = exponents(g)
    if size(U, 1) == 1
        @var λ
        λ = [λ]
    else
        @var λ[1:size(U₀, 1)]
    end
    action = Vector{Vector{Tuple{Variable, Expression}}}([[] for _ in axes(U, 1)])
    vars = variables(g)
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

algebra(g::LieAlgebraAction) = g.alg
var_groups(g::LieAlgebraAction) = g.var_groups
variables(g::LieAlgebraAction) = vcat(var_groups(g)...)
weights(g::LieAlgebraAction) = weights(g.alg)
weights(g::LieAlgebraAction, inds...) = weights(g.alg, inds...)
nweights(g::LieAlgebraAction) = nweights(g.alg)

function Base.show(io::IO, g::LieAlgebraAction)
    println(io, "LieAlgebraAction of $(name(algebra(g)))")
    print(io, " action: [", join([join(vars, ", ") for vars in var_groups(g)], "], ["), "]")
end


struct SumLieAlgebraAction <: AbstractLieAlgebraAction
    alg::SumLieAlgebra
    actions::Vector{AbstractLieAlgebraAction}
end

algebra(g::SumLieAlgebraAction) = g.alg
actions(g::SumLieAlgebraAction) = g.actions
nsummands(g::SumLieAlgebraAction) = length(actions(g))
Base.getindex(g::SumLieAlgebraAction, i::Int) = actions(g)[i]

function Base.show(io::IO, g::SumLieAlgebraAction)
    println(io, "SumLieAlgebraAction of $(name(algebra(g)))")
    # for (i, a) in enumerate(actions(g))
    #     print(io, " action of $(name(algebra(a))): [")
    #     print(io, join([join(vars, ", ") for vars in var_groups(a)], "], ["), "]")
    #     i < nsummands(g) && print(io, "\n")
    # end
    print(io, " action: ")
end

# TODO
function are_commutative(g₁::AbstractLieAlgebraAction, g₂::AbstractLieAlgebraAction)
    return true
end

function ⊕(g₁::AbstractLieAlgebraAction, g₂::AbstractLieAlgebraAction)
    @assert are_commutative(g₁, g₂)
    alg = algebra(g₁) ⊕ algebra(g₂)
    return SumLieAlgebraAction(alg, [g₁, g₂])
end

function ⊕(g₁::SumLieAlgebraAction, g₂::AbstractLieAlgebraAction)
    @assert are_commutative(g₁, g₂)
    alg = algebra(g₁) ⊕ algebra(g₂)
    return SumLieAlgebraAction(alg, [actions(g₁)..., g₂])
end

function act(elem::LieAlgebraElem{LieAlgebra}, f::Union{Expression, Monomial}, action::LieAlgebraAction)
    X = as_matrix(elem)
    return expand(sum([dot(differentiate(f, vars), -X*vars) for vars in var_groups(action)]))
end

function act(elem::LieAlgebraElem{ScalingLieAlgebra}, f::Union{Expression, Monomial}, action::ScalingLieAction)
    X = as_matrix(elem)
    return expand(dot(differentiate(f, variables(action)), -X*variables(action)))
end

function act(elem::SumLieAlgebraElem, f::Union{Expression, Monomial}, action::SumLieAlgebraAction)
    return expand(sum([act(elem[i], f, action[i]) for i in 1:nsummands(action)]))
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