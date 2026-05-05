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
    borel_ratios(b) -> Vector

Return the consecutive-coefficient ratios `b[k+1] / b[k]` of a Borel transform
(or any sequence). The length of the result is `length(b) - 1`. Element type is
preserved (no silent demotion to `Float64`).
"""
function borel_ratios(b::AbstractVector{T}) where {T<:Number}
    length(b) ≥ 2 || throw(ArgumentError("need length(b) ≥ 2"))
    return b[2:end] ./ b[1:end-1]
end
