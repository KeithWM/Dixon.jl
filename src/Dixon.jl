module Dixon

using SpecialFunctions: beta
using OffsetArrays: OffsetArray
using TaylorSeries: Taylor1
using Symbolics: @variables, @syms
using Symbolics: Differential
using Symbolics: expand_derivatives, substitute, value
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

abstract type AbstractCoefficients5{S<:Number,T<:AbstractFunction} end
struct Coefficients5{S,DixonElliptic} <: AbstractCoefficients5{S,DixonElliptic}
    z_0::S
    cm::Taylor1{S}
    sm::Taylor1{S}
end

function Coefficients5{W,T}(;
    n::Int, z_0::W, y_0::NTuple{N,V}
) where {W<:Number,T<:AbstractFunction,N,V<:Number}
    taylors = map(
        f -> create_taylors(T, f, n, z_0, y_0...), (cm=(c, s) -> c, sm=(c, s) -> s)
    )
    return Coefficients5{W,T}(z_0, taylors...)
end

function create_taylors(
    ::Type{DixonElliptic}, f, n::Int, z_0::U, cm_0::V, sm_0::V
) where {U<:W,V<:W,N} where {W<:Number}
    T = Real
    @variables z::T
    @syms cm(z::T)::T sm(z::T)::T

    fp = derivatives(DixonElliptic, f(cm(z), sm(z)), cm, sm, z, n)
    evaluated = evaluate(DixonElliptic, fp, cm, sm, z, cm_0, sm_0)
    return Taylor1(evaluated, n)
end

function derivatives(::Type{DixonElliptic}, f, cm, sm, z, n::Int)
    dz = Differential(z)
    subs = Dict(dz(cm(z)) => -sm(z)^2, dz(sm(z)) => cm(z)^2)
    fp = [substitute(expand_derivatives(f), subs)]
    for k in 1:n
        push!(fp, substitute(expand_derivatives(dz(fp[end])), subs) / BigInt(k))
    end
    return fp
end

function evaluate(::Type{DixonElliptic}, fp, cm, sm, z, cm_0::V, sm_0::V) where {V<:Number}
    subbed = map(sub -> substitute(sub, Dict(cm(z) => cm_0, sm(z) => sm_0)), fp)
    evaluated = value.(subbed)
    return evaluated
end

function compute(z::Number, coefficients::Coefficients5{<:Number,DixonElliptic}, ::Val{:cm})
    return compute(z, coefficients.z_0, coefficients, :cm, :sm)
end

function compute(z::Number, coefficients::Coefficients5{<:Number,DixonElliptic}, ::Val{:sm})
    return compute(z, coefficients.z_0, coefficients, :sm, :cm)
end

function compute(
    z::Number,
    z_0::Number,
    coefficients::Coefficients5{<:Number,DixonElliptic},
    this::Symbol,
    other::Symbol,
)
    coefficients.cm |> typeof
    z = center(z, DixonElliptic)  # it might be inconvenient that this could get called a second time?
    if real(z) <= pi3 / 6
        # use _this_ expansion
        return getproperty(coefficients, this)(z - z_0)
    else
        # revert to using the other (sm for cm and vice versa)
        return compute(pi3 / 3 - z, z_0, coefficients, other, this)
    end
end

function center(z::Complex, ::Type{DixonElliptic})
    omega = exp(2im / 3 * pi) * pi3
    n_imag = Int64(round(imag(z) / imag(omega)))
    z -= n_imag * omega
    n_real = Int64(round(real(z) / pi3))
    z -= n_real * pi3
    return z
end

function center(z::Real, ::Type{DixonElliptic})
    n_real = Int64(round(z / pi3 - 1 / 6))
    z -= n_real * pi3
    return z
end

export create_coefficients, compute
export DixonElliptic

end # module Dixon