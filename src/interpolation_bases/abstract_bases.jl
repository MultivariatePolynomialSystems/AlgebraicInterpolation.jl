export AbstractInterpolationBasis,
    AbstractMonomialVector


"""
    AbstractInterpolationBasis

Interface:
"""
abstract type AbstractInterpolationBasis end

"""
    AbstractMonomialVector <: AbstractInterpolationBasis

Interface:
"""
abstract type AbstractMonomialVector <: AbstractInterpolationBasis end