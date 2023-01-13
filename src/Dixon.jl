module Dixon

using SpecialFunctions: beta
using OffsetArrays: OffsetArray
import Base.push!

export OffsetArray

pi3 = beta(1//3, 1//3)
export pi3

abstract type Recursion end

function create_coefficients(
    number_type::Type,
    n::Int,
    recursion::Type{<:Recursion}
)
    initial_vectors = initial_coefficients(recursion)
    return create_coefficients(number_type, n, initial_vectors, recursion)
end

function create_coefficients(
    number_type::Type,
    n::Int,
    initial_vectors::NamedTuple{NT,NTuple{M,Vector{S}}},
    recursion::Type{<:Recursion},
) where {NT,M,S<:Number}
    # return Coefficients{M,N,T,NT}(initial_vector)
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
    vectors::NamedTuple{NT,NTuple{M,V}}, recursion::Type{<:Recursion}
) where {NT,M,V<:Vector}
    return map(k->push!(vectors[k], recursion, Val{k}()), NT)
end

function push!(
    vector::Vector, recursion::R, ::Val
) where {R<:Type{Recursion}}
    error("Recurion type $R must implement the push! method.")
end

"""
Methods for Dixon's elliptic formula
"""
struct DixonRecursion <: Recursion end

function push!(
    vectors::NamedTuple{NT,NTuple{M,V}}, recursion::Type{<:DixonRecursion}
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

initial_coefficients(::Type{DixonRecursion}) = (cm=[1], sm=[1])

export create_coefficients
export Recursion, DixonRecursion

end # module Dixon