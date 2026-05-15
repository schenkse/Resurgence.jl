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

"""
    wynn_eps(a, n; depth = 1)

Apply the Wynn ε-algorithm to the partial sums of `a` and return the entry at
even column `2·depth` of the ε-tableau.

`a[k]` is the k-th term of the series; with partial sums `Aⱼ = sum(a[1:j])`,
the recursion is

    ε^{(n)}_{-1} = 0,    ε^{(n)}_0 = Aₙ,
    ε^{(n)}_{k+1} = ε^{(n+1)}_{k-1} + 1 / (ε^{(n+1)}_k − ε^{(n)}_k).

The even-column entries are the meaningful approximants; the odd columns are
auxiliary. Column `2·depth` realises the `depth`-th-order Shanks transform via
its determinantal form. At `depth = 1` this is one Aitken-Δ² step and matches
[`shanks`](@ref)`(a, n; depth = 1)` exactly; at higher `depth` it differs
entry-by-entry from iterated [`shanks`](@ref) (which iterates the Δ² operator)
but shares the same limit, and is typically better conditioned.

When two consecutive column-`k` entries coincide (the sequence has effectively
converged at that position), the recursion would divide by zero; the
implementation returns that entry instead of `NaN`/`Inf`.
"""
function wynn_eps(a::AbstractVector{T}, n::Integer; depth::Integer = 1) where {T<:Number}
    depth ≥ 1 || throw(ArgumentError("depth must be ≥ 1"))
    n ≥ 2 || throw(ArgumentError("n must be ≥ 2"))
    N = n + depth
    N ≤ length(a) || throw(ArgumentError("need length(a) ≥ n+depth = $N, got $(length(a))"))
    S = Vector{T}(undef, N)
    acc = zero(T)
    @inbounds for k in 1:N
        acc += a[k]
        S[k] = acc
    end
    return _wynn_eps_iter(S, depth)
end

function _wynn_eps_iter(S::AbstractVector{T}, depth::Integer) where {T<:Number}
    N = length(S)
    N ≥ 2 * depth + 1 || throw(ArgumentError("not enough partial sums to reach ε column $(2*depth)"))
    # Column -1 is identically zero; pad to length N so col_prev2[j+1] is valid.
    col_prev2 = zeros(T, N)
    col_prev  = collect(S)
    for _ in 1:(2 * depth)
        L = length(col_prev) - 1
        nxt = Vector{T}(undef, L)
        @inbounds for j in 1:L
            d = col_prev[j+1] - col_prev[j]
            iszero(d) && return col_prev[j+1]
            nxt[j] = col_prev2[j+1] + inv(d)
        end
        col_prev2 = col_prev
        col_prev  = nxt
    end
    return last(col_prev)
end

"""
    theta_brezinski(a, n; depth = 1)

Apply Brezinski's θ-algorithm to the partial sums of `a` and return the entry
at even column `2·depth` of the θ-tableau.

`a[k]` is the k-th term of the series; with partial sums `Aⱼ = sum(a[1:j])`,
the recursion is

    θ_{-1}^{(n)} = 0,    θ_0^{(n)} = Aₙ,
    θ_{2k+1}^{(n)} = θ_{2k-1}^{(n+1)} + 1 / (θ_{2k}^{(n+1)} − θ_{2k}^{(n)}),
    θ_{2k+2}^{(n)} = θ_{2k}^{(n+1)}
                     + (θ_{2k}^{(n+2)} − θ_{2k}^{(n+1)})
                       · (θ_{2k+1}^{(n+2)} − θ_{2k+1}^{(n+1)})
                     / (θ_{2k+1}^{(n+2)} − 2 θ_{2k+1}^{(n+1)} + θ_{2k+1}^{(n)}).

Even columns carry the accelerated estimates; odd columns are auxiliary. The
θ-algorithm complements [`wynn_eps`](@ref) on sequences where the ε-algorithm
breaks down (denominator vanishes), at the cost of three new data points per
acceleration level — hence the requirement `length(a) ≥ n + 3·depth`.

When a denominator vanishes the implementation returns the previous-column
entry at that position instead of `NaN`/`Inf`.
"""
function theta_brezinski(a::AbstractVector{T}, n::Integer;
                         depth::Integer = 1) where {T<:Number}
    depth ≥ 1 || throw(ArgumentError("depth must be ≥ 1"))
    n ≥ 1 || throw(ArgumentError("n must be ≥ 1"))
    N = n + 3 * depth
    N ≤ length(a) || throw(ArgumentError(
        "need length(a) ≥ n+3·depth = $N, got $(length(a))"))
    S = Vector{T}(undef, N)
    acc = zero(T)
    @inbounds for k in 1:N
        acc += a[k]
        S[k] = acc
    end
    return _theta_brezinski_iter(S, depth)
end

function _theta_brezinski_iter(S::AbstractVector{T}, depth::Integer) where {T<:Number}
    col_minus1 = zeros(T, length(S))
    col_prev = collect(S)
    for _ in 1:depth
        L = length(col_prev) - 1
        L ≥ 3 || throw(ArgumentError("not enough partial sums to iterate θ-algorithm"))
        col_odd = Vector{T}(undef, L)
        @inbounds for j in 1:L
            d = col_prev[j+1] - col_prev[j]
            iszero(d) && return col_prev[j+1]
            col_odd[j] = col_minus1[j+1] + inv(d)
        end
        L2 = L - 2
        col_even = Vector{T}(undef, L2)
        @inbounds for j in 1:L2
            Δprev = col_prev[j+2] - col_prev[j+1]   # Δθ_{2k}^{(n+1)}
            Δodd  = col_odd[j+2] - col_odd[j+1]      # Δθ_{2k+1}^{(n+1)}
            Δ2odd = col_odd[j+2] - 2 * col_odd[j+1] + col_odd[j]
            iszero(Δ2odd) && return col_prev[j+1]
            col_even[j] = col_prev[j+1] + Δprev * Δodd / Δ2odd
        end
        col_minus1 = col_odd
        col_prev = col_even
    end
    return last(col_prev)
end

"""
    rho_brezinski(a, n; depth = 1)

Apply Brezinski's ρ-algorithm to the partial sums of `a` and return the
entry at even column `2·depth` of the ρ-tableau.

`a[k]` is the k-th term of the series; with partial sums `Aⱼ = sum(a[1:j])`,

    ρ_{-1}^{(n)} = 0,    ρ_0^{(n)} = Aₙ,
    ρ_{k+1}^{(n)} = ρ_{k-1}^{(n+1)} + (k + 1) / (ρ_k^{(n+1)} − ρ_k^{(n)}).

Structurally identical to [`wynn_eps`](@ref) but with `(k + 1)` numerator,
which makes the tableau converge geometrically on logarithmically
convergent sequences `Aₙ ≈ A + C/nᵖ` with non-integer `p` — the regime
where [`richardson`](@ref) (designed for integer power-law tails) stalls.

When two consecutive column-`k` entries coincide the recursion would divide
by zero; the implementation returns that entry instead of `NaN`/`Inf`.
"""
function rho_brezinski(a::AbstractVector{T}, n::Integer;
                       depth::Integer = 1) where {T<:Number}
    depth ≥ 1 || throw(ArgumentError("depth must be ≥ 1"))
    n ≥ 1 || throw(ArgumentError("n must be ≥ 1"))
    N = n + 2 * depth
    N ≤ length(a) || throw(ArgumentError(
        "need length(a) ≥ n+2·depth = $N, got $(length(a))"))
    S = Vector{T}(undef, N)
    acc = zero(T)
    @inbounds for k in 1:N
        acc += a[k]
        S[k] = acc
    end
    return _rho_brezinski_iter(S, depth)
end

function _rho_brezinski_iter(S::AbstractVector{T}, depth::Integer) where {T<:Number}
    N = length(S)
    col_prev2 = zeros(T, N)
    col_prev = collect(S)
    for k in 1:(2 * depth)
        L = length(col_prev) - 1
        nxt = Vector{T}(undef, L)
        @inbounds for j in 1:L
            d = col_prev[j+1] - col_prev[j]
            iszero(d) && return col_prev[j+1]
            nxt[j] = col_prev2[j+1] + T(k) / d
        end
        col_prev2 = col_prev
        col_prev = nxt
    end
    return last(col_prev)
end

"""
    cesaro(a, n; depth = 1)

Cesàro mean of the partial sums of `a` at index `n`. With `depth = d`, the
arithmetic-mean averaging is iterated `d` times; `depth = 1` is the standard
(C, 1) Cesàro sum.

`a[k]` is the k-th term of the series; with partial sums `Aⱼ = sum(a[1:j])`,

    C⁽¹⁾ₙ = (1/n) Σⱼ₌₁ⁿ Aⱼ,    C⁽ᵏ⁺¹⁾ₙ = (1/n) Σⱼ₌₁ⁿ C⁽ᵏ⁾ⱼ.

Useful mainly as a baseline against more aggressive accelerators; recovers the
limit of `(C, k)`-summable sequences with `depth = k`.
"""
function cesaro(a::AbstractVector{T}, n::Integer; depth::Integer = 1) where {T<:Number}
    depth ≥ 1 || throw(ArgumentError("depth must be ≥ 1"))
    n ≥ 1 || throw(ArgumentError("n must be ≥ 1"))
    n ≤ length(a) || throw(ArgumentError("need length(a) ≥ n, got $(length(a))"))
    cur = Vector{T}(undef, n)
    acc = zero(T)
    @inbounds for k in 1:n
        acc += a[k]
        cur[k] = acc
    end
    for _ in 1:depth
        s = zero(T)
        @inbounds for k in 1:n
            s += cur[k]
            cur[k] = s / T(k)
        end
    end
    return last(cur)
end

"""
    abel(a; x = 1)

Truncated Abel-style evaluation `Σₖ a[k+1] xᵏ` via Horner's method. The
classical Abel sum is the radial limit `x → 1⁻`; this function evaluates the
finite power series at a user-supplied `x` (default `x = 1`, which reduces to
the plain partial sum). Useful as a baseline and as a building block when
sweeping `x` towards 1.
"""
function abel(a::AbstractVector{T}; x = one(real(T))) where {T<:Number}
    R = promote_type(T, typeof(x))
    isempty(a) && return zero(R)
    xR = convert(R, x)
    acc = zero(R)
    @inbounds for k in length(a):-1:1
        acc = acc * xR + a[k]
    end
    return acc
end
