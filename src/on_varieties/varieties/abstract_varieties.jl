export AbstractAlgebraicVariety,
    AbstractDifferentiatedVariety,
    AbstractSampledVariety


"""
    AbstractAlgebraicVariety

Interface:
"""
abstract type AbstractAlgebraicVariety end

"""
    AbstractDifferentiatedVariety <: AbstractAlgebraicVariety

Interface:
"""
abstract type AbstractDifferentiatedVariety <: AbstractAlgebraicVariety end

"""
    AbstractSampledVariety <: AbstractDifferentiatedVariety

Interface:
"""
abstract type AbstractSampledVariety <: AbstractDifferentiatedVariety end