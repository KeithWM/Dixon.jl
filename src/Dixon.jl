module Dixon

using SpecialFunctions: beta
using OffsetArrays: OffsetArray
import Base.push!

export OffsetArray

pi3 = beta(1//3, 1//3)
export pi3

abstract type AbstractFunction end
struct DixonElliptic <: AbstractFunction end

abstract type AbstractCoefficients{S<:Number, T<:AbstractFunction} end
struct Coefficients{S, DixonElliptic} <: AbstractCoefficients{S, DixonElliptic}
    cm::Vector{S}
    sm::Vector{S}
end

function Coefficients{S, T}(
    n::Int
) where {S<:Number, T<:AbstractFunction}
    initial_vectors = initial_coefficients(T)
    coefficients = create_coefficients(S, n, initial_vectors, T)
    return Coefficients{S, T}(coefficients...)
end

function create_coefficients(
    number_type::Type{<:Number},
    n::Int,
    initial_vectors::NamedTuple{NT,NTuple{M,Vector{S}}},
    recursion::Type{<:AbstractFunction},
) where {NT,M,S<:Number}
    lengths = map(length, initial_vectors)
    all(l -> l == lengths[1], lengths) ||
        error("Not all initial vectors have the same length")
    n_initial = first(lengths) 
    vectors = NamedTuple{NT,NTuple{M,Vector{number_type}}}(initial_vectors)
     for k in n_initial:n-1
        push!(vectors, recursion)
    end
    return vectors
end

"""
The generic case when each vector is independent
"""
function push!(
    vectors::NamedTuple{NT,NTuple{M,V}}, recursion::Type{AbstractFunction}
) where {NT,M,V<:Vector}
    return map(k->push!(vectors[k], recursion, Val{k}()), NT)
end

function push!(
    vector::Vector, recursion::R, ::Val
) where {R<:Type{AbstractFunction}}
    error("Recurion type $R must implement the push! method.")
end

"""
Methods for Dixon's elliptic formula
"""
function push!(
    vectors::NamedTuple{NT,NTuple{M,V}}, ::Type{DixonElliptic}
) where {NT,M,V<:Vector}
    n = length(vectors[:cm]) + 1 - 1
    push!(
        vectors[:cm], 
        - 1//(3*n) * sum(vectors[:sm] .* vectors[:sm][end:-1:begin])
        )
    push!(
        vectors[:sm], 
        1//(3*n + 1) * sum(vectors[:cm] .* vectors[:cm][end:-1:begin])
    )
end

initial_coefficients(::Type{DixonElliptic}) = (cm=[1], sm=[1])

function compute(x::Number, coefficients::Coefficients{<:Number, DixonElliptic}, ::Val{:cm})
    sum(map((k, c) -> c^(3*(k-1)), enumerate(coefficients.cm)))
end

function compute(x::Number, coefficients::Coefficients{<:Number, DixonElliptic}, ::Val{:sm})
    sum(map((k, c) -> c^(3*(k-1) + 1), enumerate(coefficients.sm)))
end

export create_coefficients, compute
export Recursion, DixonRecursion

end # module Dixon