export AbstractInterpolationBasis,
    nelements,
    to_expressions,
    evaluate


EXPRESSION_LINK = "https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable/model_kit/#HomotopyContinuation.ModelKit.Expression"

"""
    abstract type AbstractInterpolationBasis end
"""
abstract type AbstractInterpolationBasis end

"""
    nelements(B::AbstractInterpolationBasis) -> Integer

Returns the number of elements in `B`.
"""
function nelements(B::AbstractInterpolationBasis)
    error("Not implemented")
end

"""
    Base.iterate(B::AbstractInterpolationBasis, state=nothing)

TBW
"""
function Base.iterate(B::AbstractInterpolationBasis, state=nothing)
    error("Not implemented")
end

"""
    to_expressions(B::AbstractInterpolationBasis) -> Vector{Expression}

Returns the elements of `B` converted to [`Expression`]($(EXPRESSION_LINK))s.
"""
function to_expressions(B::AbstractInterpolationBasis)
    error("Not implemented")
end

"""
    HC.evaluate(B::AbstractInterpolationBasis, samples::AbstractMatrix{<:Number})

TBW
"""
function HC.evaluate(B::AbstractInterpolationBasis, samples::AbstractMatrix{<:Number})
    error("Not implemented")
end