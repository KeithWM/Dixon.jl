module Dixon

using SpecialFunctions: beta
using OffsetArrays: OffsetArray
using TaylorSeries: Taylor1, evaluate
import Base.push!

export OffsetArray

pi3 = beta(1//3, 1//3)
export pi3

abstract type AbstractFunction end
initial_coefficients(::T) where T<:AbstractFunction = error("AbstractFunction $T does not implement initial_coefficients method.")
power(::Int, ::T, ::Val) where T<:AbstractFunction = error("AbstractFunction $T does not implement power method.")

struct DixonElliptic <: AbstractFunction end

abstract type AbstractCoefficients{S<:Number, T<:AbstractFunction} end
struct Coefficients{S, DixonElliptic} <: AbstractCoefficients{S, DixonElliptic}
    cm::Taylor1{S}
    sm::Taylor1{S}
end

function Coefficients{S, T}(
    n::Int
) where {S<:Number, T<:AbstractFunction}
    initial_vectors = initial_coefficients(T)
    coefficients = create_coefficients(S, n, initial_vectors, T)
    taylors = create_taylors(coefficients, T)
    return Coefficients{S, T}(taylors...)
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

function create_taylors(
    coefficients::NamedTuple{NT,NTuple{M,Vector{S}}},
    recursion::Type{<:AbstractFunction},
) where {NT,M,S<:Number}
    map(nt->create_taylors(coefficients[nt], recursion, nt), NT)
end

function create_taylors(
    coefficients::Vector{S},
    recursion::Type{<:AbstractFunction},
    nt::Symbol
) where {S<:Number}
    powers = power.(1:length(coefficients), recursion, Val{nt}())
    taylor = zeros(S, maximum(powers) + 1)
    for (k, c) in enumerate(coefficients)
        taylor[powers[k] + 1] = c
    end
    return Taylor1(taylor, maximum(powers))
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

function compute(z::Number, coefficients::Coefficients{<:Number, DixonElliptic}, ::Val{:cm})
    return compute(z, coefficients, :cm, :sm)
end

function compute(z::Number, coefficients::Coefficients{<:Number, DixonElliptic}, ::Val{:sm})
    return compute(z, coefficients, :sm, :cm)
end

function compute(z::Number, coefficients::Coefficients{<:Number, DixonElliptic}, this::Symbol, other::Symbol)
    z = center(z, DixonElliptic)  # it might be inconvenient that this could get called a second time?
    if real(z) <= pi3/6
        # use _this_ expansion
        return getproperty(coefficients, this)(z)
    else
        # revert to using the other (sm for cm and vice versa)
        return compute(pi3/3 - z, coefficients, other, this)
    end
end

function center(z::Complex, ::Type{DixonElliptic})
    omega = exp(2im / 3 * pi) * pi3
    n_imag = round(imag(z) / imag(omega)) |> Int64
    z -= n_imag * omega
    n_real = round(real(z) / pi3) |> Int64
    z -= n_real * pi3
    return z
end

function center(z::Real, ::Type{DixonElliptic})
    n_real = round(z / pi3) |> Int64
    z -= n_real
    return z
end

power(k::Int, ::Type{DixonElliptic}, ::Val{:cm}) = 3 * (k - 1)
power(k::Int, ::Type{DixonElliptic}, ::Val{:sm}) = 3 * (k - 1) + 1

export create_coefficients, compute
export Recursion, DixonElliptic

end # module Dixon