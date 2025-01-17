export PolynomialVectorSpace,
    is_upto

struct PolynomialVectorSpace
    vars::Vector{Variable}
    degree::Int
    upto::Bool
end

PolynomialVectorSpace(;
    variables::AbstractArray,
    degree::Int,
    upto::Bool=true
) = PolynomialVectorSpace(Variable.(collect(flatten(variables))), degree, upto)

variables(V::PolynomialVectorSpace) = V.vars
nvariables(V::PolynomialVectorSpace) = length(V.vars)
degree(V::PolynomialVectorSpace) = V.degree
is_upto(V::PolynomialVectorSpace) = V.upto
degrees(V::PolynomialVectorSpace) = is_upto(V) ? range(0, degree(V)) : [degree(V)]
dim(V::PolynomialVectorSpace) = num_mons(nvariables(V), degree(V); upto=is_upto(V))

Base.rand(
    V::PolynomialVectorSpace,
    n::Int
) = random_monomial_basis(length=n, nvars=nvariables(V), degree=degree(V), upto=upto(V))

Base.:(==)(
    V::PolynomialVectorSpace,
    W::PolynomialVectorSpace
) = (V.vars == W.vars) && (V.degree == W.degree) && (V.upto == W.upto)

function Base.show(io::IO, V::PolynomialVectorSpace)
    println(io, "PolynomialVectorSpace of dimension $(dim(V))")
    println(io, " $(nvariables(V)) variables: ", join(variables(V), ", "))
    println(io, " degree of polynomials: $(degree(V))")
    print(io, " homogeneous: $(!is_upto(V))")
end