"""
    borel_pade(a; n, m, x = 1, regularize_poles = false, ε = nothing,
               side = +1, return_error = false, quad_kwargs...)

Borel–Padé resummation of the formal power series with coefficients `a`,
evaluated at argument `x`:

1. Borel-transform `a` to obtain `B[k] = a[k] / k!`.
2. Replace `B(t)` by its Padé approximant `[n/m]` `Bₚ(t) / B_q(t)`.
3. Laplace-transform `Bₚ(t) / B_q(t)` along the positive real axis:

        ∫₀^∞ (Bₚ(x t) / B_q(x t)) e^{-t} dt.

When `regularize_poles = true`, any poles of the Padé denominator on the
positive real axis are first shifted by `side · iε` (default
`ε = 100·eps(real(T))`, default `side = +1`) to push them off the contour.
This branch is only meaningful for `x > 0` because the integration runs along
positive `t`; for `x ≤ 0` an `ArgumentError` is raised. Selecting `side = +1`
or `side = −1` produces the two lateral Borel sums; see
[`borel_pade_lateral`](@ref).

`quad_kwargs` are forwarded to `QuadGK.quadgk` (e.g. `rtol`, `atol`, `order`).
With `return_error = true`, returns a `(value, abserr)` tuple; otherwise just
the value.

Requires `length(a) ≥ n + m + 1`.
"""
function borel_pade(a::AbstractVector{T};
                    n::Integer, m::Integer, x = 1,
                    regularize_poles::Bool = false,
                    ε::Union{Real,Nothing} = nothing,
                    side::Integer = +1,
                    return_error::Bool = false,
                    quad_kwargs...) where {T<:Number}
    length(a) ≥ n + m + 1 ||
        throw(ArgumentError("borel_pade needs length(a) ≥ n+m+1 = $(n+m+1)"))
    if regularize_poles && !(x isa Real && x > 0)
        throw(ArgumentError("regularize_poles=true requires positive real x (got $x)"))
    end
    Ba = borel_transform(a[1:n+m+1])
    Bp, Bq = pade(Ba, n, m)

    integrand = if regularize_poles
        ε_eff = ε === nothing ? 100 * eps(real(T)) : ε
        Bqr = poles_regularized(Polynomials.coeffs(Bq), ε_eff; side = side)
        # default to order=24 for the regularized contour unless caller overrides
        if !haskey(quad_kwargs, :order)
            quad_kwargs = (; quad_kwargs..., order = 24)
        end
        t -> Bp(x * t) / Bqr(x * t) * exp(-t)
    else
        t -> Bp(x * t) / Bq(x * t) * exp(-t)
    end

    val, err = quadgk(integrand, 0, Inf; quad_kwargs...)
    return return_error ? (val, err) : val
end

"""
    borel_leroy_pade(a; n, m, b = -1//2, x = 1, regularize_poles = false,
                     ε = nothing, side = +1, return_error = false, quad_kwargs...)

Borel–Le Roy–Padé resummation:

1. Borel–Le Roy transform `a` with parameter `b` to obtain `B[k] = a[k] / Γ(k+1+b)`.
2. Replace `B(t)` by its Padé approximant `[n/m]`.
3. Laplace-transform with the Le Roy weight:

        ∫₀^∞ t^b (Bₚ(x t) / B_q(x t)) e^{-t} dt.

`b` may be any real number (including non-integer values); the default `-1//2`
matches the original code. Pole regularization (`regularize_poles`, `ε`,
`side`) works as in [`borel_pade`](@ref); selecting opposite sides produces the
two lateral Le Roy–Borel sums. `quad_kwargs` are forwarded to `QuadGK.quadgk`.
"""
function borel_leroy_pade(a::AbstractVector{T};
                          n::Integer, m::Integer,
                          b::Real = -1//2, x = 1,
                          regularize_poles::Bool = false,
                          ε::Union{Real,Nothing} = nothing,
                          side::Integer = +1,
                          return_error::Bool = false,
                          quad_kwargs...) where {T<:Number}
    length(a) ≥ n + m + 1 ||
        throw(ArgumentError("borel_leroy_pade needs length(a) ≥ n+m+1 = $(n+m+1)"))
    if regularize_poles && !(x isa Real && x > 0)
        throw(ArgumentError("regularize_poles=true requires positive real x (got $x)"))
    end
    Ba = borel_leroy_transform(a[1:n+m+1], b)
    Bp, Bq = pade(Ba, n, m)

    integrand = if regularize_poles
        ε_eff = ε === nothing ? 100 * eps(real(T)) : ε
        Bqr = poles_regularized(Polynomials.coeffs(Bq), ε_eff; side = side)
        if !haskey(quad_kwargs, :order)
            quad_kwargs = (; quad_kwargs..., order = 24)
        end
        t -> t^b * Bp(x * t) / Bqr(x * t) * exp(-t)
    else
        t -> t^b * Bp(x * t) / Bq(x * t) * exp(-t)
    end

    val, err = quadgk(integrand, 0, Inf; quad_kwargs...)
    return return_error ? (val, err) : val
end

"""
    mittag_leffler_borel_pade(a; n, m, α, x = 1, return_error = false, quad_kwargs...)

Mittag-Leffler / generalised-Borel–Padé resummation of order `α > 0`:

1. Compute the Mittag-Leffler-Borel transform `B[k] = a[k] / Γ(α·k + 1)`
   (see [`mittag_leffler_borel_transform`](@ref)).
2. Replace `B(t)` by its Padé approximant `[n/m]`.
3. Integrate against the Mittag-Leffler kernel, in the change of variable
   `u = t^{1/α}` that linearises the exponential weight:

       f(x) ≈ ∫₀^∞ (Bₚ(x·u^α) / B_q(x·u^α)) · e^{−u} du.

For `α = 1` this reduces to [`borel_pade`](@ref). For `α > 1` it tolerates
super-factorial coefficient growth (e.g. `(2k)!`) that defeats the standard
Borel kernel; for `α < 1` it sharpens convergence on sub-factorial series.

No pole-regularisation / lateral variants for v1: the `u^α` change of variable
maps real-axis poles of the Padé denominator to fractional-power surfaces in
`u`, which calls for separate treatment.

`quad_kwargs` are forwarded to `QuadGK.quadgk` (e.g. `rtol`, `atol`, `order`).
Requires `length(a) ≥ n + m + 1`.
"""
function mittag_leffler_borel_pade(a::AbstractVector{T};
                                   n::Integer, m::Integer,
                                   α::Real, x = 1,
                                   return_error::Bool = false,
                                   quad_kwargs...) where {T<:Number}
    length(a) ≥ n + m + 1 ||
        throw(ArgumentError("mittag_leffler_borel_pade needs length(a) ≥ n+m+1 = $(n+m+1)"))
    α > 0 || throw(ArgumentError("mittag_leffler_borel_pade needs α > 0 (got $α)"))
    Ba = mittag_leffler_borel_transform(a[1:n+m+1], α)
    Bp, Bq = pade(Ba, n, m)

    αR = real(T)(α)
    integrand = u -> Bp(x * u^αR) / Bq(x * u^αR) * exp(-u)
    val, err = quadgk(integrand, 0, Inf; quad_kwargs...)
    return return_error ? (val, err) : val
end

"""
    conformal_borel_pade(a; n, m, x = 1, sing = 1, return_error = false, quad_kwargs...)

Conformal-Borel–Padé resummation:

1. Borel-transform `a`.
2. Re-expand the Borel transform under the conformal map
   `t = 4·sing·w / (1-w)²` (see [`conformal_reseries`](@ref)).
3. Padé-approximate the conformal series in `w`.
4. Map back via `w = (√(1+t/sing) - 1)/(√(1+t/sing) + 1)` and Laplace-integrate
   along the positive real `t` axis.

`sing > 0` is the assumed location of the nearest singularity of the Borel
transform on the negative real axis (i.e. the singularity sits at `t = -sing`).
"""
function conformal_borel_pade(a::AbstractVector{T};
                              n::Integer, m::Integer, x = 1,
                              sing::Real = 1,
                              return_error::Bool = false,
                              quad_kwargs...) where {T<:Number}
    length(a) ≥ n + m + 1 ||
        throw(ArgumentError("conformal_borel_pade needs length(a) ≥ n+m+1 = $(n+m+1)"))
    Ba = borel_transform(a[1:n+m+1])
    Bw = conformal_reseries(Ba, sing, n + m)
    Bp, Bq = pade(Bw, n, m)
    integrand = function (t)
        w = conformal_map(x * t; a = sing)
        return Bp(w) / Bq(w) * exp(-t)
    end
    val, err = quadgk(integrand, 0, Inf; quad_kwargs...)
    return return_error ? (val, err) : val
end

"""
    conformal_borel_pade_pair(a; n, m, x = 1, sing = 1, return_error = false, quad_kwargs...)

Conformal-Borel–Padé resummation for a complex-conjugate Borel-plane
singularity pair at `t = ±i·sing`:

1. Borel-transform `a`.
2. Re-expand the Borel transform under the singularity-pair conformal map
   `t = 2·sing·v / (1 - v²)` (see [`conformal_reseries_pair`](@ref)).
3. Padé-approximate the conformal series in `v`.
4. Map back via `v = (√(t² + sing²) − sing) / t` (numerically as
   `v = t / (√(t² + sing²) + sing)`) and Laplace-integrate along the
   positive real `t` axis.

`sing > 0` is the assumed magnitude of the imaginary singularity pair
(singularities at `t = ±i·sing`). PT-symmetric problems and the quartic
anharmonic oscillator are typical use cases.
"""
function conformal_borel_pade_pair(a::AbstractVector{T};
                                   n::Integer, m::Integer, x = 1,
                                   sing::Real = 1,
                                   return_error::Bool = false,
                                   quad_kwargs...) where {T<:Number}
    length(a) ≥ n + m + 1 ||
        throw(ArgumentError("conformal_borel_pade_pair needs length(a) ≥ n+m+1 = $(n+m+1)"))
    Ba = borel_transform(a[1:n+m+1])
    Bv = conformal_reseries_pair(Ba, sing, n + m)
    Bp, Bq = pade(Bv, n, m)
    integrand = function (t)
        v = conformal_map_pair(x * t; a = sing)
        return Bp(v) / Bq(v) * exp(-t)
    end
    val, err = quadgk(integrand, 0, Inf; quad_kwargs...)
    return return_error ? (val, err) : val
end

"""
    borel_pade_lateral(a; n, m, x = 1, side = +1, ε = nothing, kwargs...)

Lateral Borel–Padé sum: shorthand for
`borel_pade(a; n, m, x, regularize_poles = true, side, ε, kwargs...)`. The
positive-real-axis poles of the Padé denominator are pushed off the contour
by `side · iε`, giving one of the two lateral Borel sums. For real `a` and
real `x > 0` the result is generally complex.

Convention: `side = +1` shifts poles to `+iε` (upper half-plane), so the
real-axis integration contour passes below them — this matches the lateral
sum traditionally written `S_θ⁻` in the resurgence literature; `side = -1`
gives `S_θ⁺`. Requires `x > 0`.

Defaults `ε` to `√eps(real(eltype(a)))`. The lateral integrand has a
near-singularity of height `∼1/ε` at the shifted pole, so the very tight
default of `borel_pade` (`100·eps`) makes `quadgk` underflow; `√eps` keeps
both the pole shift small and the integrand resolvable.
"""
function borel_pade_lateral(a::AbstractVector{T};
                            n::Integer, m::Integer, x = 1,
                            side::Integer = +1,
                            ε::Union{Real,Nothing} = nothing,
                            kwargs...) where {T<:Number}
    ε_use = ε === nothing ? sqrt(eps(real(T))) : ε
    return borel_pade(a; n = n, m = m, x = x,
                      regularize_poles = true, side = side, ε = ε_use, kwargs...)
end

# (L⁺, L⁻) combinators shared by the median and discontinuity wrappers.
_median(Lp, Lm) = (Lp + Lm) / 2
_discontinuity(Lp, Lm) = (Lp - Lm) / (2 * im)

# Run `lateral_fn(a; side = +1, kwargs...)` and `lateral_fn(a; side = -1, kwargs...)`,
# then return `combine(Lp, Lm)`. Used to build median/discontinuity wrappers from
# their lateral building blocks without copy-pasting the signature.
function _lateral_combine(combine, lateral_fn, a; kwargs...)
    Lp = lateral_fn(a; side = +1, kwargs...)
    Lm = lateral_fn(a; side = -1, kwargs...)
    return combine(Lp, Lm)
end

"""
    borel_pade_median(a; n, m, x = 1, ε = nothing, kwargs...)

Median Borel–Padé sum `(L⁺ + L⁻)/2`, the average of the two lateral sums. For
real `a` the result is purely real (up to roundoff and `quadgk` tolerance);
this is the ambiguity-free real-valued resummation of a non-Borel-summable
series.
"""
borel_pade_median(a::AbstractVector;
                  n::Integer, m::Integer, x = 1,
                  ε::Union{Real,Nothing} = nothing,
                  kwargs...) =
    _lateral_combine(_median, borel_pade_lateral, a;
                     n = n, m = m, x = x, ε = ε, kwargs...)

"""
    borel_pade_discontinuity(a; n, m, x = 1, ε = nothing, kwargs...)

Stokes discontinuity `(L⁺ − L⁻) / (2i)` of the Borel–Padé sum. Encodes the
ambiguity of the lateral sums and is proportional to the leading instanton
contribution `Stokes_constant · e^{−|S|/x}`. For real `a` the result is real.
"""
borel_pade_discontinuity(a::AbstractVector;
                         n::Integer, m::Integer, x = 1,
                         ε::Union{Real,Nothing} = nothing,
                         kwargs...) =
    _lateral_combine(_discontinuity, borel_pade_lateral, a;
                     n = n, m = m, x = x, ε = ε, kwargs...)

"""
    borel_leroy_pade_lateral(a; n, m, b = -1//2, x = 1, side = +1, ε = nothing, kwargs...)

Lateral Borel–Le Roy–Padé sum. Same lateral conventions as
[`borel_pade_lateral`](@ref); the only difference is the Le Roy weight `t^b`
in the Laplace integrand. `ε` defaults to `√eps(real(eltype(a)))`, large
enough to keep the integrand resolvable across the shifted pole.
"""
function borel_leroy_pade_lateral(a::AbstractVector{T};
                                  n::Integer, m::Integer,
                                  b::Real = -1//2, x = 1,
                                  side::Integer = +1,
                                  ε::Union{Real,Nothing} = nothing,
                                  kwargs...) where {T<:Number}
    ε_use = ε === nothing ? sqrt(eps(real(T))) : ε
    return borel_leroy_pade(a; n = n, m = m, b = b, x = x,
                            regularize_poles = true, side = side,
                            ε = ε_use, kwargs...)
end

"""
    borel_leroy_pade_median(a; n, m, b = -1//2, x = 1, ε = nothing, kwargs...)

Median Borel–Le Roy–Padé sum `(L⁺ + L⁻) / 2`.
"""
borel_leroy_pade_median(a::AbstractVector;
                        n::Integer, m::Integer,
                        b::Real = -1//2, x = 1,
                        ε::Union{Real,Nothing} = nothing,
                        kwargs...) =
    _lateral_combine(_median, borel_leroy_pade_lateral, a;
                     n = n, m = m, b = b, x = x, ε = ε, kwargs...)

"""
    borel_leroy_pade_discontinuity(a; n, m, b = -1//2, x = 1, ε = nothing, kwargs...)

Stokes discontinuity `(L⁺ − L⁻) / (2i)` of the Borel–Le Roy–Padé sum.
"""
borel_leroy_pade_discontinuity(a::AbstractVector;
                               n::Integer, m::Integer,
                               b::Real = -1//2, x = 1,
                               ε::Union{Real,Nothing} = nothing,
                               kwargs...) =
    _lateral_combine(_discontinuity, borel_leroy_pade_lateral, a;
                     n = n, m = m, b = b, x = x, ε = ε, kwargs...)

"""
    borel_leroy_pade_odm(a; n, m, x = 1, b_grid = range(-0.4, 0.4; length = 17),
                         return_b = false, kwargs...)

Order-dependent-mapping (ODM) Borel–Le Roy–Padé sum: sweep the Le Roy
parameter `b` over `b_grid`, locate the stationary point of the
resummation result vs. `b`, and return the value there. Variational
selection of `b` à la Bender–Boettcher: `b` is treated as a tunable
parameter that the result should depend on weakly.

The "stationary point" is the interior grid index that minimises the
forward-difference quotient `|Δresult / Δb|`; ties are broken by
`|Δ²result / Δb²|` (flattest). With `return_b = true`, returns a
`(value, b_star)` tuple; otherwise just the value.

`b_grid` must have at least three points so that a meaningful difference
quotient can be formed. The default range stays well inside `(-1, 0)`,
where the Borel–Le Roy weight is regular and `Γ(k+1+b)` has no integer
poles.

`kwargs...` are forwarded to [`borel_leroy_pade`](@ref) (e.g. `rtol`,
`atol` for the integration).
"""
function borel_leroy_pade_odm(a::AbstractVector;
                              n::Integer, m::Integer, x = 1,
                              b_grid = range(-0.4, 0.4; length = 17),
                              return_b::Bool = false,
                              kwargs...)
    length(b_grid) ≥ 3 ||
        throw(ArgumentError("borel_leroy_pade_odm needs length(b_grid) ≥ 3 (got $(length(b_grid)))"))
    bs = collect(b_grid)
    # Each grid point is an independent Padé fit + quadgk; the b values share
    # nothing reusable since the Le Roy weight Γ(k+1+b) re-scales every Borel
    # coefficient. Spawn one task per point and gather; degenerates to
    # sequential at JULIA_NUM_THREADS=1.
    tasks = map(bs) do b
        Threads.@spawn borel_leroy_pade(a; n = n, m = m, b = b, x = x, kwargs...)
    end
    vals = map(fetch, tasks)
    bstar, vstar = _odm_stationary(bs, vals)
    return return_b ? (vstar, bstar) : vstar
end

# Pick the interior grid index minimising |Δvalue/Δparam|; ties broken by
# the smaller |Δ²value/Δparam²| (flattest).
function _odm_stationary(params::AbstractVector, vals::AbstractVector)
    length(params) == length(vals) ||
        throw(ArgumentError("params and vals must have the same length"))
    N = length(params)
    N ≥ 3 || throw(ArgumentError("need ≥ 3 points (got $N)"))
    # Forward first differences along the interior of the grid.
    bestidx = 0
    bestdv = Inf
    bestcurv = Inf
    @inbounds for i in 2:N-1
        dp_fwd = params[i+1] - params[i]
        dp_bwd = params[i] - params[i-1]
        # central first difference (matches forward at uniform spacing).
        dv = abs((vals[i+1] - vals[i-1]) / (params[i+1] - params[i-1]))
        # central second difference for tie-break.
        d2v = abs((vals[i+1] - 2 * vals[i] + vals[i-1]) / (dp_fwd * dp_bwd))
        if dv < bestdv || (dv == bestdv && d2v < bestcurv)
            bestdv = dv
            bestcurv = d2v
            bestidx = i
        end
    end
    return params[bestidx], vals[bestidx]
end
