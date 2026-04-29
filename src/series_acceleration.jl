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
    An  = sum(@view a[1:n])
    Anp = sum(@view a[1:n+1])
    d1 = Anp - An
    d2 = An - Anm
    return Anp - d1^2 / (d1 - d2)
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
            nxt[k-1] = cur[k+1] - d1^2 / (d1 - d2)
        end
        cur = nxt
    end
    return last(cur)
end
