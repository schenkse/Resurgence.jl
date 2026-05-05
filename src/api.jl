# Unified resummation API.
#
# Each subtype of `AbstractResummation` packages a method together with its
# parameters; `resum(method, a)` then dispatches to the canonical per-method
# function. Per-method functions remain the public, documented entry points;
# this layer is a convenience for code that needs to compare methods or treat
# them uniformly.

"""
    AbstractResummation

Abstract supertype for resummation method tags. See concrete subtypes
[`Shanks`](@ref), [`Pade`](@ref), [`BorelPade`](@ref), [`BorelLeRoyPade`](@ref),
and [`ConformalBorelPade`](@ref).
"""
abstract type AbstractResummation end

"""
    Shanks(n; depth = 1)

Shanks transformation tag. `resum(Shanks(n), a)` calls `shanks(a, n; depth)`.
"""
struct Shanks <: AbstractResummation
    n::Int
    depth::Int
    Shanks(n::Integer; depth::Integer = 1) = new(Int(n), Int(depth))
end

"""
    Richardson(n; depth = 1)

Richardson extrapolation tag. `resum(Richardson(n), a)` calls
`richardson(a, n; depth)`.
"""
struct Richardson <: AbstractResummation
    n::Int
    depth::Int
    Richardson(n::Integer; depth::Integer = 1) = new(Int(n), Int(depth))
end

"""
    Pade(n, m; x = 1)

Padé approximant tag. `resum(Pade(n, m; x), a)` calls `pade_value(a, n, m, x)`.
Argument order is numerator-first to match the Borel–Padé family.
"""
struct Pade{X} <: AbstractResummation
    n::Int
    m::Int
    x::X
    Pade(n::Integer, m::Integer; x = 1) = new{typeof(x)}(Int(n), Int(m), x)
end

"""
    BorelPade(n, m; x = 1, kwargs...)

Borel–Padé tag. `resum(BorelPade(n, m; x, kwargs...), a)` calls
`borel_pade(a; n, m, x, kwargs...)`.
"""
struct BorelPade{X,K} <: AbstractResummation
    n::Int
    m::Int
    x::X
    kwargs::K
    BorelPade(n::Integer, m::Integer; x = 1, kwargs...) =
        new{typeof(x),typeof(kwargs)}(Int(n), Int(m), x, kwargs)
end

"""
    BorelLeRoyPade(n, m; b = -1//2, x = 1, kwargs...)

Borel–Le Roy–Padé tag. `resum(BorelLeRoyPade(n, m; b, x, kwargs...), a)` calls
`borel_leroy_pade(a; n, m, b, x, kwargs...)`.
"""
struct BorelLeRoyPade{X,B,K} <: AbstractResummation
    n::Int
    m::Int
    b::B
    x::X
    kwargs::K
    BorelLeRoyPade(n::Integer, m::Integer; b::Real = -1//2, x = 1, kwargs...) =
        new{typeof(x),typeof(b),typeof(kwargs)}(Int(n), Int(m), b, x, kwargs)
end

"""
    ConformalBorelPade(n, m; x = 1, sing = 1, kwargs...)

Conformal-Borel–Padé tag. `resum(ConformalBorelPade(n, m; x, sing, kwargs...), a)`
calls `conformal_borel_pade(a; n, m, x, sing, kwargs...)`.
"""
struct ConformalBorelPade{X,S,K} <: AbstractResummation
    n::Int
    m::Int
    x::X
    sing::S
    kwargs::K
    ConformalBorelPade(n::Integer, m::Integer; x = 1, sing::Real = 1, kwargs...) =
        new{typeof(x),typeof(sing),typeof(kwargs)}(Int(n), Int(m), x, sing, kwargs)
end

"""
    resum(method::AbstractResummation, a)

Apply `method` to the formal power series with coefficients `a`.
"""
resum(s::Shanks, a) = shanks(a, s.n; depth = s.depth)
resum(r::Richardson, a) = richardson(a, r.n; depth = r.depth)
resum(p::Pade, a) = pade_value(a, p.n, p.m, p.x)
resum(bp::BorelPade, a) = borel_pade(a; n = bp.n, m = bp.m, x = bp.x, bp.kwargs...)
resum(bl::BorelLeRoyPade, a) =
    borel_leroy_pade(a; n = bl.n, m = bl.m, b = bl.b, x = bl.x, bl.kwargs...)
resum(cb::ConformalBorelPade, a) =
    conformal_borel_pade(a; n = cb.n, m = cb.m, x = cb.x, sing = cb.sing, cb.kwargs...)
