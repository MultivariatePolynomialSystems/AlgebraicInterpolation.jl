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

Return the number of elements in `B`.
"""
function nelements(B::AbstractInterpolationBasis)
    error("Not implemented")
end

"""
    to_expressions(B::AbstractInterpolationBasis) -> AbstractVector{Expression}

Return the elements of `B` converted to [`Expression`]($(EXPRESSION_LINK))s.
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