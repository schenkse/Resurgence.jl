"""
    shanks(a, n; depth = 1)

Apply the Shanks transformation to the partial sums of `a`, evaluated at index
`n`. With `depth = k`, the transformation is iterated `k` times, which is the
standard practice for accelerating slowly convergent or alternating series.

`a[k]` is the k-th term of the series; the partial sum `Aₙ = sum(a[1:n])`. The
single-application formula is

    Sₙ = Aₙ₊₁ − (Aₙ₊₁ − Aₙ)² / ((Aₙ₊₁ − Aₙ) − (Aₙ − Aₙ₋₁))

For higher `depth`, the transformation is applied to the sequence of Shanks
values obtained from the previous round.

When consecutive partial-sum differences match exactly (the series has
effectively converged, or has a constant tail), the Shanks denominator
vanishes; the implementation returns the current partial sum in that case
rather than `NaN`/`Inf`.
"""
function shanks(a::AbstractVector{T}, n::Integer; depth::Integer = 1) where {T<:Number}
    depth ≥ 1 || throw(ArgumentError("depth must be ≥ 1"))
    n ≥ 2 || throw(ArgumentError("n must be ≥ 2"))
    n + 1 ≤ length(a) || throw(ArgumentError("need length(a) ≥ n+1, got $(length(a)) for n=$n"))
    if depth == 1
        return _shanks_once(a, n)
    end
    # Iterate: build the Shanks sequence around `n` to feed the next round.
    # We need 2·depth + 1 partial sums to take `depth` Shanks levels at index n.
    return _shanks_iter(a, n, depth)
end

function _shanks_once(a::AbstractVector{T}, n::Integer) where {T<:Number}
    Anm = sum(@view a[1:n-1])
    An  = Anm + a[n]
    Anp = An  + a[n+1]
    d1 = Anp - An
    d2 = An - Anm
    den = d1 - d2
    iszero(den) && return Anp
    return Anp - d1^2 / den
end

function _shanks_iter(a::AbstractVector{T}, n::Integer, depth::Integer) where {T<:Number}
    # Compute partial sums S[1..N] = A_1 .. A_N where N = n + depth.
    N = n + depth
    N ≤ length(a) || throw(ArgumentError("need length(a) ≥ n+depth = $N, got $(length(a))"))
    S = Vector{T}(undef, N)
    acc = zero(T)
    @inbounds for k in 1:N
        acc += a[k]
        S[k] = acc
    end
    # Apply the eps-algorithm-style Shanks recursion `depth` times to S.
    # After one application, S_new[k] = Shanks of (S[k-1], S[k], S[k+1]).
    cur = S
    for _ in 1:depth
        m = length(cur)
        m ≥ 3 || throw(ArgumentError("not enough partial sums to iterate Shanks"))
        nxt = Vector{T}(undef, m - 2)
        @inbounds for k in 2:m-1
            d1 = cur[k+1] - cur[k]
            d2 = cur[k] - cur[k-1]
            den = d1 - d2
            nxt[k-1] = iszero(den) ? cur[k+1] : cur[k+1] - d1^2 / den
        end
        cur = nxt
    end
    return last(cur)
end

"""
    richardson(a, n; depth = 1)

Apply Richardson extrapolation to the partial sums of `a`, evaluated at index
`n`. With `depth = d`, the transformation is iterated `d` times, eliminating
the leading `1/n, 1/n², …, 1/nᵈ` tails of the asymptotic expansion of the
partial-sum sequence.

`a[k]` is the k-th term of the series; the partial sum `Aₙ = sum(a[1:n])`.
The single-application formula is

    Rₙ = (n+1) Aₙ₊₁ − n Aₙ

For higher `depth`, the iterative recursion

    T⁽ᵏ⁺¹⁾ⱼ = ((j + k + 1) T⁽ᵏ⁾ⱼ₊₁ − j T⁽ᵏ⁾ⱼ) / (k + 1)

is applied with `T⁽⁰⁾ⱼ = Aⱼ`. The result is `T⁽ᵈ⁾ₙ`, identical to the
Salzer / Bender-Orszag (§8.1.30) closed form but free of its alternating-
sign cancellation. Richardson is the standard companion to [`shanks`](@ref)
for monotone, power-law-convergent series.
"""
function richardson(a::AbstractVector{T}, n::Integer; depth::Integer = 1) where {T<:Number}
    depth ≥ 1 || throw(ArgumentError("depth must be ≥ 1"))
    n ≥ 1 || throw(ArgumentError("n must be ≥ 1"))
    n + depth ≤ length(a) || throw(ArgumentError("need length(a) ≥ n+depth = $(n+depth), got $(length(a))"))
    if depth == 1
        return _richardson_once(a, n)
    end
    return _richardson_iter(a, n, depth)
end

function _richardson_once(a::AbstractVector{T}, n::Integer) where {T<:Number}
    An  = sum(@view a[1:n])
    Anp = sum(@view a[1:n+1])
    return (n + 1) * Anp - n * An
end

function _richardson_iter(a::AbstractVector{T}, n::Integer, depth::Integer) where {T<:Number}
    # Build partial sums S[1..N] = A_1 .. A_N where N = n + depth, then run
    # the Richardson recursion over a shrinking buffer indexed from `n`.
    N = n + depth
    N ≤ length(a) || throw(ArgumentError("need length(a) ≥ n+depth = $N, got $(length(a))"))
    cur = Vector{T}(undef, depth + 1)
    acc = zero(T)
    @inbounds for k in 1:N
        acc += a[k]
        if k ≥ n
            cur[k - n + 1] = acc
        end
    end
    # cur[i] currently holds T⁽⁰⁾_{n + i - 1} for i = 1 .. depth+1.
    # After level k, cur[i] holds T⁽ᵏ⁾_{n + i - 1} for i = 1 .. depth+1-k.
    for k in 0:depth-1
        @inbounds for i in 1:(depth - k)
            j = n + i - 1
            cur[i] = ((j + k + 1) * cur[i + 1] - j * cur[i]) / (k + 1)
        end
    end
    return cur[1]
end
