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
"""
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
    R = real(T)
    aT = T(a)
    # t(w) = 4a · w · (1-w)^(-2) = 4a · Σ_{k≥0} (k+1) w^{k+1}
    # Build the polynomial representation of t(w) truncated at order N.
    tw_coeffs = zeros(T, N + 1)
    @inbounds for k in 0:N-1
        tw_coeffs[k+2] = T(4) * aT * T(k + 1)
    end
    tw = Polynomial(tw_coeffs)

    # Accumulate B(t(w)) = Σ_k b_k · t(w)^k, truncating each power at order N.
    out = zeros(T, N + 1)
    out[1] = b[1]                         # k = 0 contribution
    pow = Polynomial([one(T)])            # t(w)^0
    @inbounds for k in 1:length(b)-1
        pow = _truncate(pow * tw, N)
        c = Polynomials.coeffs(pow)
        # add b[k+1] * c into out, padding if shorter than out
        for i in eachindex(c)
            out[i] += b[k+1] * c[i]
        end
    end
    return out
end

function _truncate(p::Polynomial{T}, N::Integer) where {T}
    c = Polynomials.coeffs(p)
    if length(c) ≤ N + 1
        return p
    end
    return Polynomial(c[1:N+1])
end
