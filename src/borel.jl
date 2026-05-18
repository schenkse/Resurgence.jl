"""
    borel_transform(a) -> Vector

Compute the Borel transform of the formal power series with coefficients `a`:

    B[k] = a[k] / k!

The result has the same length and element type as `a`. Inverse factorials are
generated iteratively (`f[k+1] = f[k] / k`), so this is overflow-free for
arbitrary-precision element types and avoids the ad-hoc `BigInt` conversion
needed when calling `factorial(::Int)` for `k > 20`.
"""
function borel_transform(a::AbstractVector{T}) where {T<:Number}
    invf = inv_factorials(T, length(a))
    return invf .* a
end

"""
    borel_leroy_transform(a, b) -> Vector

Borel–Le Roy transform with parameter `b`:

    B[k] = a[k] / Γ(k + 1 + b)

Reduces to the ordinary Borel transform when `b = 0`. The parameter `b` is
promoted to `eltype(a)` to keep the result type-stable.
"""
function borel_leroy_transform(a::AbstractVector{T}, b::Real) where {T<:Number}
    R = real(T)
    bT = R(b)
    out = similar(a)
    isempty(a) && return out
    # Γ(k+1+b) via the recurrence Γ(z+1) = z·Γ(z): seed with Γ(1+b), then
    # multiply by (k+b) at each step. One `gamma` call instead of length(a),
    # which is a big win for BigFloat where each call is expensive.
    g = gamma(one(R) + bT)
    @inbounds out[1] = a[1] / g
    @inbounds for k in 1:length(a)-1
        g *= R(k) + bT
        out[k+1] = a[k+1] / g
    end
    return out
end

"""
    mittag_leffler_borel_transform(a, α) -> Vector

Mittag-Leffler-Borel transform of order `α > 0`:

    B[k] = a[k] / Γ(α·k + 1).

The associated inverse transform (see [`mittag_leffler_borel_pade`](@ref)) uses
the kernel `ϕ_α(t) = (1/α)·t^{(1−α)/α}·e^{−t^{1/α}}` whose `k`-th moment is
`Γ(α·k + 1)`. With α = 1 this reduces to the standard Borel transform; α > 1
tolerates super-factorial growth such as `(2k)!`; α < 1 sharpens convergence
for sub-factorial growth.

Unlike [`borel_leroy_transform`](@ref) which can amortise `Γ` via the
recurrence `Γ(k+1+b) = (k+b)·Γ(k+b)`, the `Γ(α·k+1)` factors do not share a
common increment for non-integer `α`, so one `gamma` call is issued per
coefficient. For BigFloat element types this remains dominated by the
`gamma` cost.
"""
function mittag_leffler_borel_transform(a::AbstractVector{T}, α::Real) where {T<:Number}
    α > 0 || throw(ArgumentError("mittag_leffler_borel_transform needs α > 0 (got $α)"))
    R = real(T)
    αR = R(α)
    out = similar(a)
    @inbounds for k in 0:length(a)-1
        out[k+1] = a[k+1] / gamma(αR * R(k) + one(R))
    end
    return out
end

"""
    borel_ratios(b) -> Vector

Return the consecutive-coefficient ratios `b[k+1] / b[k]` of a Borel transform
(or any sequence). The length of the result is `length(b) - 1`. Element type is
preserved (no silent demotion to `Float64`).
"""
function borel_ratios(b::AbstractVector{T}) where {T<:Number}
    length(b) ≥ 2 || throw(ArgumentError("need length(b) ≥ 2"))
    any(iszero, @view b[1:end-1]) &&
        throw(ArgumentError("borel_ratios: divisor is zero (b[k] = 0 for some k < length(b))"))
    return b[2:end] ./ b[1:end-1]
end
