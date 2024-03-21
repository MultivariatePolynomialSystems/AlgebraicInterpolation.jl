export AffineSpace


"""
    AffineSpace <: AbstractAlgebraicVariety
"""
struct AffineSpace <: AbstractAlgebraicVariety
    vars::Vector{Variable}
end