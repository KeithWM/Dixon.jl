module Dixon

using SpecialFunctions: beta
using OffsetArrays: OffsetArray
using TaylorSeries: Taylor1
using Symbolics: @variables, @syms
using Symbolics: Differential
using Symbolics: substitute, value
import Base.push!

export OffsetArray

pi3 = beta(1//3, 1//3)
export pi3

abstract type AbstractFunction end
function initial_coefficients(::T) where {T<:AbstractFunction}
    return error("AbstractFunction $T does not implement initial_coefficients method.")
end
function power(::Int, ::T, ::Val) where {T<:AbstractFunction}
    return error("AbstractFunction $T does not implement power method.")
end

struct DixonElliptic <: AbstractFunction end

abstract type AbstractCoefficients{S<:Number,T<:AbstractFunction} end
struct Coefficients{S,DixonElliptic} <: AbstractCoefficients{S,DixonElliptic}
    cm::Taylor1{S}
    sm::Taylor1{S}
end

function Coefficients{S,T}(
    f, n::Int, z_0::U, y_0::NTuple{N,U}
) where {S<:Number,T<:AbstractFunction,U,N}
    taylors = create_taylors(T, f, n, z_0, y_0...)
    return Coefficients{S,T}(taylors...)
end

function create_taylors(
    ::Type{DixonElliptic}, f, n::Int, z_0::U, cm_0::U, sm_0::U
) where {U<:Number}
    @variables z::T
    @syms cm(z::T)::T sm(z::T)::T

    fp = derivatives{DixonElliptic}(f, n)
    evaluated = evaluate{DixonElliptic}(fp, n, z_0, cm_0, sm_0)
    return Taylor1(evaluated, n)
end

function derivatives(::Type{DixonElliptic}, f, n::Int)
    dz = Differential(z)
    subs = Dict(dz(cm(z)) => -sm(z)^2, dz(sm(z)) => cm(z)^2)
    fp = [substitute(expand_derivatives(f), subs)]
    for k in 1:n
        push!(fp, substitute(expand_derivatives(dz(fp[end])), subs) * (1 / k))
    end
    return fp
end

function evaluate(::Type{DixonElliptic}, fp, cm_0::U, sm_0::U) where {U<:Number}
    subbed = map(sub -> substitute(sub, Dict(cm(z) => cm_0, sm(z) => sm_0)), fp)
    evaluated = Symbolics.value.(subbed)
    return evaluated
end

export create_coefficients, compute
export DixonElliptic

end # module Dixon