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
    # coefficients.
    out = zeros(T, N + 1)
    out[1] = b[1]                         # k = 0 contribution
    pow = zeros(T, N + 1); pow[1] = one(T)   # t(w)^0
    @inbounds for k in 1:length(b)-1
        pow = _mul_trunc(pow, tw, N)
        bk = b[k+1]
        for i in eachindex(pow)
            out[i] += bk * pow[i]
        end
    end
    return out
end

function _mul_trunc(p::AbstractVector{T}, q::AbstractVector{T}, N::Integer) where {T}
    out = zeros(T, N + 1)
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
