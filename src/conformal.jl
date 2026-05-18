# Conformal mapping for the Borel plane.
#
# Given a Borel transform B(t) with its nearest singularity on the negative real
# axis at t = -a (a > 0) — the typical situation for asymptotic series in
# perturbative QFT — the map
#
#     t(w) = 4·a·w / (1 - w)²,    w(t) = (√(1 + t/a) - 1) / (√(1 + t/a) + 1)
#
# sends the unit disk |w| < 1 in the w-plane to the cut t-plane (cut along
# (-∞, -a]). Re-expanding B(t) in powers of w produces a series with a strictly
# wider domain of convergence than the original Borel series, suitable for
# Padé–Laplace processing.

"""
    conformal_map(t; a = 1)

Forward conformal map `w(t) = (√(1 + t/a) - 1) / (√(1 + t/a) + 1)`. The
parameter `a > 0` is the (assumed) location of the nearest Borel-plane
singularity on the negative real axis at `t = -a`.

For real `t` in the cut region `t < -a`, the result is complex (the map
sends the cut t-plane to the unit disk). Real `t ≥ -a` returns a real value;
complex `t` is handled directly.
"""
function conformal_map(t::Real; a::Real = 1)
    a > 0 || throw(ArgumentError("a must be positive (got $a)"))
    arg = 1 + t / a
    s = arg ≥ 0 ? sqrt(arg) : sqrt(complex(arg))
    return (s - 1) / (s + 1)
end

function conformal_map(t::Number; a::Real = 1)
    a > 0 || throw(ArgumentError("a must be positive (got $a)"))
    s = sqrt(1 + t / a)
    return (s - 1) / (s + 1)
end

"""
    inverse_conformal_map(w; a = 1)

Inverse conformal map `t(w) = 4·a·w / (1 - w)²`.
"""
function inverse_conformal_map(w::Number; a::Real = 1)
    a > 0 || throw(ArgumentError("a must be positive (got $a)"))
    return 4 * a * w / (1 - w)^2
end

"""
    conformal_reseries(b, a, N) -> Vector

Re-expand the Borel-plane series with coefficients `b` (constant term first)
under the conformal substitution `t = 4·a·w / (1 - w)²`, producing a vector of
`N + 1` coefficients in the variable `w`.

Internally this is power-series composition truncated at order `N`. No
intermediate symbolic engine is required.
"""
function conformal_reseries(b::AbstractVector{T}, a::Real, N::Integer) where {T<:Number}
    N ≥ 0 || throw(ArgumentError("N must be non-negative"))
    aT = T(a)
    # t(w) = 4a · w · (1-w)^(-2) = 4a · Σ_{k≥0} (k+1) w^{k+1}
    # Build the coefficient vector of t(w) truncated at order N.
    tw = zeros(T, N + 1)
    @inbounds for k in 0:N-1
        tw[k+2] = T(4) * aT * T(k + 1)
    end

    # Accumulate B(t(w)) = Σ_k b_k · t(w)^k, truncating each power at order N
    # via a direct convolution that never materialises the discarded high-order
    # coefficients. Two scratch buffers `pow` and `nextpow` are swapped each
    # iteration so no allocation happens inside the loop.
    out = zeros(T, N + 1)
    out[1] = b[1]                         # k = 0 contribution
    pow = zeros(T, N + 1); pow[1] = one(T)   # t(w)^0
    nextpow = zeros(T, N + 1)
    @inbounds for k in 1:length(b)-1
        _mul_trunc!(nextpow, pow, tw, N)
        pow, nextpow = nextpow, pow
        bk = b[k+1]
        for i in eachindex(pow)
            out[i] += bk * pow[i]
        end
    end
    return out
end

# In-place truncated multiplication: out ← p · q (coefficient convolution),
# truncated at total degree N. Caller owns `out`; it is zeroed on entry.
function _mul_trunc!(out::AbstractVector{T}, p::AbstractVector{T},
                     q::AbstractVector{T}, N::Integer) where {T}
    fill!(out, zero(T))
    @inbounds for i in eachindex(p)
        pi_ = p[i]
        iszero(pi_) && continue
        ki_max = min(N - (i - 1), length(q) - 1)
        for j in 1:ki_max+1
            out[(i - 1) + j] += pi_ * q[j]
        end
    end
    return out
end

# Non-mutating wrapper kept for callers that want the convenience form.
function _mul_trunc(p::AbstractVector{T}, q::AbstractVector{T}, N::Integer) where {T}
    return _mul_trunc!(zeros(T, N + 1), p, q, N)
end

# Conformal map for a complex-conjugate singularity pair at t = ±i·a
# (a > 0). The map
#
#     v(t) = t / (√(t² + a²) + a),    t(v) = 2 a v / (1 - v²)
#
# is the 1-to-1 conformal map of the cut t-plane (cuts emanating from
# ±i·a out to ±i·∞) to the open unit disk in v. The singularities at
# t = ±i·a map to v = ±i on the boundary; t = ±∞ maps to v = ±1. The
# forward form is rationalised so that v(0) = 0 is hit cleanly without
# a 0/0.
#
# Used for Borel transforms whose nearest singularities sit at ±i·a — a
# common situation in PT-symmetric problems and the quartic anharmonic
# oscillator. Re-expanding B(t) as B(t(v)) and Padé-approximating in v
# enlarges the domain of convergence past the imaginary singularities.

"""
    conformal_map_pair(t; a = 1)

Forward conformal map `v(t) = t / (√(t² + a²) + a)` for a complex-conjugate
Borel-plane singularity pair at `t = ±i·a` (`a > 0`). Maps the cut t-plane
to the open unit disk in `v`; the singularities at `t = ±i·a` map to
`v = ±i` on the unit circle, and `t = ±∞` map to `v = ±1`.

For real `t` the result is real (the map is regular on the real axis); for
complex `t`, the principal-branch square root is used.
"""
function conformal_map_pair(t::Real; a::Real = 1)
    a > 0 || throw(ArgumentError("a must be positive (got $a)"))
    return t / (sqrt(t^2 + a^2) + a)
end

function conformal_map_pair(t::Number; a::Real = 1)
    a > 0 || throw(ArgumentError("a must be positive (got $a)"))
    return t / (sqrt(t^2 + a^2) + a)
end

"""
    inverse_conformal_map_pair(v; a = 1)

Inverse conformal map `t(v) = 2 a v / (1 - v²)` for the singularity-pair
map. Analytic on `|v| < 1`.
"""
function inverse_conformal_map_pair(v::Number; a::Real = 1)
    a > 0 || throw(ArgumentError("a must be positive (got $a)"))
    return 2 * a * v / (1 - v^2)
end

"""
    conformal_reseries_pair(b, a, N) -> Vector

Re-expand the Borel-plane series with coefficients `b` (constant term
first) under the singularity-pair conformal substitution
`t = 2 a v / (1 - v²)`, producing a vector of `N + 1` coefficients in `v`.

Identical structure to [`conformal_reseries`](@ref); only the polynomial
of `t(v)` truncated to order `N` differs. Power-series composition,
truncated at every step, so the high-order discarded coefficients never
materialise.
"""
function conformal_reseries_pair(b::AbstractVector{T}, a::Real, N::Integer) where {T<:Number}
    N ≥ 0 || throw(ArgumentError("N must be non-negative"))
    aT = T(a)
    # t(v) = 2a · v · (1 + v² + v⁴ + …); truncated coefficients are
    # nonzero only at odd positions.
    tv = zeros(T, N + 1)
    @inbounds for k in 0:2:N-1
        tv[k+2] = T(2) * aT
    end

    out = zeros(T, N + 1)
    out[1] = b[1]
    pow = zeros(T, N + 1); pow[1] = one(T)
    nextpow = zeros(T, N + 1)
    @inbounds for k in 1:length(b)-1
        _mul_trunc!(nextpow, pow, tv, N)
        pow, nextpow = nextpow, pow
        bk = b[k+1]
        for i in eachindex(pow)
            out[i] += bk * pow[i]
        end
    end
    return out
end
