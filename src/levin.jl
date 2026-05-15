"""
    levin(a, n; depth = length(a) - n - 1, variant = :u, β = 1)

Apply the Levin sequence transformation to the partial sums of `a` at index
`n` with transform order `depth = k` and remainder estimate selected by
`variant`. Frequently outperforms [`shanks`](@ref) and [`wynn_eps`](@ref) on
monotone divergent or asymptotic series — the regime resurgence cares about.

`a[k]` is the k-th term of the series; with partial sums `Aⱼ = sum(a[1:j])`,

    Lₖ⁽ⁿ⁾ = Σⱼ (−1)ʲ C(k,j) (n+j+β)^{k−1} Aₙ₊ⱼ / ωₙ₊ⱼ
          ─────────────────────────────────────────────
            Σⱼ (−1)ʲ C(k,j) (n+j+β)^{k−1}      / ωₙ₊ⱼ

with `j = 0, …, k`. The remainder-estimate weight `ωₙ` depends on `variant`:

- `:t` —  `ωₙ = aₙ`                   (Smith–Ford t-transformation).
- `:u` —  `ωₙ = (n + β) aₙ`           (Smith–Ford u-transformation, default).
- `:v` —  `ωₙ = aₙ aₙ₊₁ / (aₙ₊₁ − aₙ)` (Smith–Ford v-transformation).

The default `depth` consumes all available data and is safe for every
variant; the `:v` variant additionally requires `a[n + depth + 1]` because
`ωₙ₊depth` reads one term beyond the partial-sum index.

When a term's weight is zero (a vanishes locally) that term is dropped from
both sums; if every weight is zero the partial sum `Aₙ₊depth` is returned.
"""
function levin(a::AbstractVector{T}, n::Integer;
               depth::Integer = length(a) - n - 1,
               variant::Symbol = :u, β::Real = 1) where {T<:Number}
    depth ≥ 1 || throw(ArgumentError("depth must be ≥ 1"))
    n ≥ 1 || throw(ArgumentError("n must be ≥ 1"))
    variant in (:u, :t, :v) ||
        throw(ArgumentError("variant must be :u, :t, or :v; got :$variant"))
    needed = variant === :v ? n + depth + 1 : n + depth
    needed ≤ length(a) || throw(ArgumentError(
        "need length(a) ≥ $needed for variant=:$variant, depth=$depth; got $(length(a))"))

    R = promote_type(T, typeof(one(T) * β))
    Np = n + depth
    S = Vector{R}(undef, Np)
    acc = zero(R)
    @inbounds for k in 1:Np
        acc += a[k]
        S[k] = acc
    end

    invω = _levin_inv_omega(R, a, n, depth, variant, β)
    return _levin_combine(S, invω, n, depth, R, β)
end

function _levin_inv_omega(::Type{R}, a, n, depth, variant, β) where {R}
    invω = Vector{R}(undef, depth + 1)
    if variant === :t
        @inbounds for j in 0:depth
            an = R(a[n+j])
            invω[j+1] = iszero(an) ? zero(R) : inv(an)
        end
    elseif variant === :u
        @inbounds for j in 0:depth
            an = R(a[n+j])
            invω[j+1] = iszero(an) ? zero(R) : inv(R(n + j + β) * an)
        end
    else # :v
        @inbounds for j in 0:depth
            an = R(a[n+j])
            ap = R(a[n+j+1])
            prod = an * ap
            invω[j+1] = iszero(prod) ? zero(R) : (ap - an) / prod
        end
    end
    return invω
end

function _levin_combine(S::AbstractVector{R}, invω::AbstractVector{R},
                        n::Integer, depth::Integer, ::Type{R}, β) where {R}
    k = depth
    num = zero(R)
    den = zero(R)
    C = one(R)            # C(k, 0)
    sign = one(R)         # (-1)^0
    @inbounds for j in 0:k
        x = R(n + j + β)
        w = x^(k - 1)     # k ≥ 1, so this is well-defined
        coef = sign * C * w
        num += coef * S[n+j] * invω[j+1]
        den += coef * invω[j+1]
        if j < k
            C = C * R(k - j) / R(j + 1)
            sign = -sign
        end
    end
    iszero(den) && return S[n+depth]
    return num / den
end

"""
    weniger(a, n; depth = length(a) - n - 1, β = 1)

Apply Weniger's δ-transformation to the partial sums of `a` at index `n`
with transform order `depth = k`. Often the single best black-box accelerator
for factorially divergent series — the textbook driver is the Stieltjes /
Euler series `Σ (-1)ᵏ k!`.

`a[k]` is the k-th term of the series; with partial sums `Aⱼ = sum(a[1:j])`,

    δₖ⁽ⁿ⁾ = Σⱼ (−1)ʲ C(k,j) (β+n+j)_{k−1} · Aₙ₊ⱼ / aₙ₊ⱼ₊₁
          ─────────────────────────────────────────────────
            Σⱼ (−1)ʲ C(k,j) (β+n+j)_{k−1} ·          1 / aₙ₊ⱼ₊₁

with `j = 0, …, k` and `(x)_p = x(x+1)…(x+p−1)` the rising factorial. The
remainder estimate `ωₙ = aₙ₊₁` makes this a Levin-style transform with
factorial rather than power weighting; `length(a) ≥ n + depth + 1` is
required because the j = depth term reads `aₙ₊depth₊₁`.

The constant overall factor `1/(β+n+k)_{k−1}` cancels between numerator and
denominator and is therefore dropped from the implementation.
"""
function weniger(a::AbstractVector{T}, n::Integer;
                 depth::Integer = length(a) - n - 1,
                 β::Real = 1) where {T<:Number}
    depth ≥ 1 || throw(ArgumentError("depth must be ≥ 1"))
    n ≥ 1 || throw(ArgumentError("n must be ≥ 1"))
    needed = n + depth + 1
    needed ≤ length(a) || throw(ArgumentError(
        "need length(a) ≥ $needed for depth=$depth; got $(length(a))"))

    R = promote_type(T, typeof(one(T) * β))
    Np = n + depth
    S = Vector{R}(undef, Np)
    acc = zero(R)
    @inbounds for k in 1:Np
        acc += a[k]
        S[k] = acc
    end

    k = depth
    num = zero(R)
    den = zero(R)
    C = one(R)              # C(k, 0)
    sign = one(R)           # (-1)^0
    # Incremental Pochhammer: P_j = (β+n+j)_{k-1} for k ≥ 1.
    P = one(R)
    @inbounds for i in 0:k-2
        P *= R(β + n + i)   # P_0
    end
    @inbounds for j in 0:k
        an1 = R(a[n+j+1])
        invω = iszero(an1) ? zero(R) : inv(an1)
        coef = sign * C * P
        num += coef * S[n+j] * invω
        den += coef * invω
        if j < k
            # advance binomial, Pochhammer, and sign for the next j
            C = C * R(k - j) / R(j + 1)
            sign = -sign
            if k ≥ 2
                # P_{j+1} = P_j * (β+n+j+k-1) / (β+n+j)
                P = P * R(β + n + j + k - 1) / R(β + n + j)
            end
        end
    end
    iszero(den) && return S[Np]
    return num / den
end
