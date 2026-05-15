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
    WynnEps(n; depth = 1)

Wynn ε-algorithm tag. `resum(WynnEps(n; depth), a)` calls
`wynn_eps(a, n; depth)`.
"""
struct WynnEps <: AbstractResummation
    n::Int
    depth::Int
    WynnEps(n::Integer; depth::Integer = 1) = new(Int(n), Int(depth))
end

"""
    BrezinskiTheta(n; depth = 1)

Brezinski θ-algorithm tag. `resum(BrezinskiTheta(n; depth), a)` calls
`theta_brezinski(a, n; depth)`.
"""
struct BrezinskiTheta <: AbstractResummation
    n::Int
    depth::Int
    BrezinskiTheta(n::Integer; depth::Integer = 1) = new(Int(n), Int(depth))
end

"""
    BrezinskiRho(n; depth = 1)

Brezinski ρ-algorithm tag. `resum(BrezinskiRho(n; depth), a)` calls
`rho_brezinski(a, n; depth)`.
"""
struct BrezinskiRho <: AbstractResummation
    n::Int
    depth::Int
    BrezinskiRho(n::Integer; depth::Integer = 1) = new(Int(n), Int(depth))
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
    PadeCF(n, m; x = 1)

Continued-fraction Padé tag. `resum(PadeCF(n, m; x), a)` calls
`pade_cf_value(a, n, m, x)`. Restricted to `n == m` or `n + 1 == m`.
"""
struct PadeCF{X} <: AbstractResummation
    n::Int
    m::Int
    x::X
    PadeCF(n::Integer, m::Integer; x = 1) = new{typeof(x)}(Int(n), Int(m), x)
end

"""
    HermitePade(n, m, l; x = 1, branch = nothing)

Hermite/quadratic Padé tag. `resum(HermitePade(n, m, l; x, branch), a)`
calls `hermite_pade_value(a, n, m, l, x; branch)`.
"""
struct HermitePade{X,B} <: AbstractResummation
    n::Int
    m::Int
    l::Int
    x::X
    branch::B
    HermitePade(n::Integer, m::Integer, l::Integer; x = 1, branch = nothing) =
        new{typeof(x),typeof(branch)}(Int(n), Int(m), Int(l), x, branch)
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
    BorelPadeLateral(n, m; x = 1, side = +1, kwargs...)

Lateral Borel–Padé tag. `resum(BorelPadeLateral(n, m; x, side, kwargs...), a)`
calls `borel_pade_lateral(a; n, m, x, side, kwargs...)`.
"""
struct BorelPadeLateral{X,K} <: AbstractResummation
    n::Int
    m::Int
    x::X
    side::Int
    kwargs::K
    BorelPadeLateral(n::Integer, m::Integer; x = 1, side::Integer = +1, kwargs...) =
        new{typeof(x),typeof(kwargs)}(Int(n), Int(m), x, Int(side), kwargs)
end

"""
    BorelPadeMedian(n, m; x = 1, kwargs...)

Median Borel–Padé tag. `resum(BorelPadeMedian(n, m; x, kwargs...), a)` calls
`borel_pade_median(a; n, m, x, kwargs...)`.
"""
struct BorelPadeMedian{X,K} <: AbstractResummation
    n::Int
    m::Int
    x::X
    kwargs::K
    BorelPadeMedian(n::Integer, m::Integer; x = 1, kwargs...) =
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
    ConformalBorelPadePair(n, m; x = 1, sing = 1, kwargs...)

Conformal-Borel–Padé tag for a complex-conjugate Borel singularity pair at
`t = ±i·sing`. `resum(ConformalBorelPadePair(n, m; x, sing, kwargs...), a)`
calls `conformal_borel_pade_pair(a; n, m, x, sing, kwargs...)`.
"""
struct ConformalBorelPadePair{X,S,K} <: AbstractResummation
    n::Int
    m::Int
    x::X
    sing::S
    kwargs::K
    ConformalBorelPadePair(n::Integer, m::Integer; x = 1, sing::Real = 1, kwargs...) =
        new{typeof(x),typeof(sing),typeof(kwargs)}(Int(n), Int(m), x, sing, kwargs)
end

"""
    MeijerG(n; x = 1, kwargs...)

Meijer-G resummation tag. `resum(MeijerG(n; x, kwargs...), a)` calls
`borel_meijerg(a; n, x, kwargs...)`.
"""
struct MeijerG{X,K} <: AbstractResummation
    n::Int
    x::X
    kwargs::K
    MeijerG(n::Integer; x = 1, kwargs...) =
        new{typeof(x),typeof(kwargs)}(Int(n), x, kwargs)
end

"""
    Cesaro(n; depth = 1)

Cesàro-summation tag. `resum(Cesaro(n; depth), a)` calls `cesaro(a, n; depth)`.
"""
struct Cesaro <: AbstractResummation
    n::Int
    depth::Int
    Cesaro(n::Integer; depth::Integer = 1) = new(Int(n), Int(depth))
end

"""
    Abel(; x = 1)

Abel-summation tag. `resum(Abel(; x), a)` calls `abel(a; x)`.
"""
struct Abel{X} <: AbstractResummation
    x::X
    Abel(; x = 1) = new{typeof(x)}(x)
end

"""
    Levin(n; depth = nothing, variant = :u, β = 1)

Levin sequence-transformation tag. `resum(Levin(n; depth, variant, β), a)`
calls `levin(a, n; depth, variant, β)`. With `depth = nothing` (default), the
function-level default `length(a) - n - 1` is used at dispatch time.
"""
struct Levin{B,K} <: AbstractResummation
    n::Int
    depth::Union{Int,Nothing}
    variant::Symbol
    β::B
    kwargs::K
    Levin(n::Integer; depth::Union{Integer,Nothing} = nothing,
          variant::Symbol = :u, β::Real = 1, kwargs...) =
        new{typeof(β),typeof(kwargs)}(Int(n),
            depth === nothing ? nothing : Int(depth), variant, β, kwargs)
end

"""
    Weniger(n; depth = nothing, β = 1)

Weniger δ-transformation tag. `resum(Weniger(n; depth, β), a)` calls
`weniger(a, n; depth, β)`. With `depth = nothing` (default), the function-
level default `length(a) - n - 1` is used at dispatch time.
"""
struct Weniger{B,K} <: AbstractResummation
    n::Int
    depth::Union{Int,Nothing}
    β::B
    kwargs::K
    Weniger(n::Integer; depth::Union{Integer,Nothing} = nothing,
            β::Real = 1, kwargs...) =
        new{typeof(β),typeof(kwargs)}(Int(n),
            depth === nothing ? nothing : Int(depth), β, kwargs)
end

"""
    SidiS(n; depth = nothing, variant = :u, β = 1)

Sidi S-transformation tag. `resum(SidiS(n; depth, variant, β), a)` calls
`sidi_s(a, n; depth, variant, β)`. With `depth = nothing` (default), the
function-level default `length(a) - n - 1` is used at dispatch time.
"""
struct SidiS{B,K} <: AbstractResummation
    n::Int
    depth::Union{Int,Nothing}
    variant::Symbol
    β::B
    kwargs::K
    SidiS(n::Integer; depth::Union{Integer,Nothing} = nothing,
          variant::Symbol = :u, β::Real = 1, kwargs...) =
        new{typeof(β),typeof(kwargs)}(Int(n),
            depth === nothing ? nothing : Int(depth), variant, β, kwargs)
end

"""
    Hyperasymptotic(; x = 1, level = 1, kwargs...)

Hyperasymptotic-truncation tag. `resum(Hyperasymptotic(; x, level, kwargs...), a)`
calls `hyperasymptotic(a; x, level, kwargs...)`. `kwargs` may include
`action`, `β`, `A` to override the defaults derived from `stokes_fit`.
"""
struct Hyperasymptotic{X,K} <: AbstractResummation
    x::X
    level::Int
    kwargs::K
    Hyperasymptotic(; x = 1, level::Integer = 1, kwargs...) =
        new{typeof(x),typeof(kwargs)}(x, Int(level), kwargs)
end

"""
    resum(method::AbstractResummation, a)

Apply `method` to the formal power series with coefficients `a`.
"""
resum(s::Shanks, a) = shanks(a, s.n; depth = s.depth)
resum(r::Richardson, a) = richardson(a, r.n; depth = r.depth)
resum(w::WynnEps, a) = wynn_eps(a, w.n; depth = w.depth)
resum(b::BrezinskiTheta, a) = theta_brezinski(a, b.n; depth = b.depth)
resum(b::BrezinskiRho, a) = rho_brezinski(a, b.n; depth = b.depth)
resum(p::Pade, a) = pade_value(a, p.n, p.m, p.x)
resum(p::PadeCF, a) = pade_cf_value(a, p.n, p.m, p.x)
resum(h::HermitePade, a) = hermite_pade_value(a, h.n, h.m, h.l, h.x; branch = h.branch)
resum(bp::BorelPade, a) = borel_pade(a; n = bp.n, m = bp.m, x = bp.x, bp.kwargs...)
resum(bp::BorelPadeLateral, a) =
    borel_pade_lateral(a; n = bp.n, m = bp.m, x = bp.x, side = bp.side, bp.kwargs...)
resum(bp::BorelPadeMedian, a) =
    borel_pade_median(a; n = bp.n, m = bp.m, x = bp.x, bp.kwargs...)
resum(bl::BorelLeRoyPade, a) =
    borel_leroy_pade(a; n = bl.n, m = bl.m, b = bl.b, x = bl.x, bl.kwargs...)
resum(cb::ConformalBorelPade, a) =
    conformal_borel_pade(a; n = cb.n, m = cb.m, x = cb.x, sing = cb.sing, cb.kwargs...)
resum(cb::ConformalBorelPadePair, a) =
    conformal_borel_pade_pair(a; n = cb.n, m = cb.m, x = cb.x, sing = cb.sing, cb.kwargs...)
resum(mg::MeijerG, a) = borel_meijerg(a; n = mg.n, x = mg.x, mg.kwargs...)
resum(c::Cesaro, a) = cesaro(a, c.n; depth = c.depth)
resum(ab::Abel, a) = abel(a; x = ab.x)
resum(l::Levin, a) = l.depth === nothing ?
    levin(a, l.n; variant = l.variant, β = l.β, l.kwargs...) :
    levin(a, l.n; depth = l.depth, variant = l.variant, β = l.β, l.kwargs...)
resum(w::Weniger, a) = w.depth === nothing ?
    weniger(a, w.n; β = w.β, w.kwargs...) :
    weniger(a, w.n; depth = w.depth, β = w.β, w.kwargs...)
resum(s::SidiS, a) = s.depth === nothing ?
    sidi_s(a, s.n; variant = s.variant, β = s.β, s.kwargs...) :
    sidi_s(a, s.n; depth = s.depth, variant = s.variant, β = s.β, s.kwargs...)
resum(h::Hyperasymptotic, a) = hyperasymptotic(a; x = h.x, level = h.level, h.kwargs...)
