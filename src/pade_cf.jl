"""
    pade_cf(a, n, m) -> (num::Polynomial, den::Polynomial)

Continued-fraction Padé approximant of the formal power series with
coefficients `a` (constant term first), built via the qd
(quotient–difference) algorithm and a 3-term polynomial recurrence on the
convergents of the resulting Stieltjes-type S-fraction. Returns a
numerator polynomial of degree `n` and a denominator polynomial of degree
`m` such that

    num(z) / den(z) = a₀ + a₁ z + a₂ z² + … + O(z^{n+m+1})

The qd algorithm naturally produces only the diagonal staircase Padé
sequence `[k/k]` and `[k/k+1]`; this implementation therefore restricts to
`n == m` or `n + 1 == m`. For general `(n, m)` use [`pade`](@ref).

Sometimes more numerically stable than the LU/pinv path of [`pade`](@ref)
on ill-conditioned moment problems and the natural object for J-fractions
and S-fractions in Stieltjes-series theory. Conversely, it has no graceful
rank-deficient fallback: division-by-zero anywhere in the qd table raises
`ArgumentError` directing users back to [`pade`](@ref).

Requires `length(a) ≥ n + m + 1` and `a[1], …, a[n+m]` non-zero.
"""
function pade_cf(a::AbstractVector{T}, n::Integer, m::Integer) where {T<:Number}
    n ≥ 0 || throw(ArgumentError("n must be ≥ 0"))
    m ≥ 0 || throw(ArgumentError("m must be ≥ 0"))
    (n == m || n + 1 == m) || throw(ArgumentError(
        "pade_cf supports only n == m or n + 1 == m; got n=$n, m=$m. Use `pade(a, n, m)` for general degrees."))
    length(a) ≥ n + m + 1 || throw(ArgumentError(
        "Padé [$n/$m] needs length(a) ≥ $(n+m+1), got $(length(a))"))

    K = n + m
    if K == 0
        return Polynomial(T[a[1]]), Polynomial(T[one(T)])
    end

    α = _qd_alphas(a, K)

    # 3-term recurrence on the convergents A_k(z), B_k(z) of
    #     U(z) = 1 - α₁z / (1 - α₂z / (1 - α₃z / …))
    # in standard CF form with b_k = 1 and a_k = -α_k z. Then
    # f(z) = a[1] / U(z) ≈ a[1] · B_K(z) / A_K(z).
    z = Polynomial(T[zero(T), one(T)])
    A_prev = Polynomial(T[one(T)])
    A_curr = Polynomial(T[one(T)])
    B_prev = Polynomial(T[zero(T)])
    B_curr = Polynomial(T[one(T)])
    for k in 1:K
        αz = α[k] * z
        A_new = A_curr - αz * A_prev
        B_new = B_curr - αz * B_prev
        A_prev, A_curr = A_curr, A_new
        B_prev, B_curr = B_curr, B_new
    end
    return a[1] * B_curr, A_curr
end

"""
    pade_cf_value(a, n, m, x)

Evaluate the continued-fraction Padé approximant of `a` at `x` directly via
a scalar 3-term recurrence on the CF convergents — no `Polynomial`
construction. Subject to the same `n == m` or `n + 1 == m` restriction as
[`pade_cf`](@ref).

Each call rebuilds the qd table from scratch; to sweep over many `x` for
the same series, call `(num, den) = pade_cf(a, n, m)` once and evaluate
`num(x) / den(x)` directly.
"""
function pade_cf_value(a::AbstractVector{T}, n::Integer, m::Integer, x) where {T<:Number}
    n ≥ 0 || throw(ArgumentError("n must be ≥ 0"))
    m ≥ 0 || throw(ArgumentError("m must be ≥ 0"))
    (n == m || n + 1 == m) || throw(ArgumentError(
        "pade_cf_value supports only n == m or n + 1 == m; got n=$n, m=$m. Use `pade_value(a, n, m, x)` for general degrees."))
    length(a) ≥ n + m + 1 || throw(ArgumentError(
        "Padé [$n/$m] needs length(a) ≥ $(n+m+1), got $(length(a))"))

    K = n + m
    S = promote_type(T, typeof(x))
    if K == 0
        return a[1] * one(S) / one(S)
    end

    α = _qd_alphas(a, K)
    A_prev = one(S); A_curr = one(S)
    B_prev = zero(S); B_curr = one(S)
    @inbounds for k in 1:K
        αx = α[k] * x
        A_new = A_curr - αx * A_prev
        B_new = B_curr - αx * B_prev
        A_prev, A_curr = A_curr, A_new
        B_prev, B_curr = B_curr, B_new
    end
    return a[1] * B_curr / A_curr
end

# Build the qd-table α-coefficients of the S-fraction
#     f(z) = a[1] / (1 - α₁ z / (1 - α₂ z / (1 - α₃ z / …)))
# Returns `α` of length `K = n + m`, with α_{2k-1} = q^{(0)}_k and
# α_{2k} = e^{(0)}_k pulled from the top row of the qd table. Throws
# ArgumentError on division-by-zero.
#
# Storage: only the previous and current columns of the qd table are kept,
# rolling forward. Column at level k has length K - k + 1 entries.
function _qd_alphas(a::AbstractVector{T}, K::Integer) where {T<:Number}
    @inbounds for j in 1:K
        iszero(a[j]) && throw(ArgumentError(
            "pade_cf: zero divisor in qd algorithm (a[$j] = 0); use `pade()` instead"))
    end

    α = Vector{T}(undef, K)

    # Column 1 (q's): q_col[j] = q^{(j-1)}_1 = a[j+1] / a[j], j = 1..K
    q_col = T[a[j+1] / a[j] for j in 1:K]
    α[1] = q_col[1]
    K == 1 && return α

    # Column 2 (e's): e_col[j] = e^{(j-1)}_1 = q^{(j)}_1 - q^{(j-1)}_1, j = 1..K-1
    # (uses e^{(j)}_0 = 0)
    e_col = T[q_col[j+1] - q_col[j] for j in 1:K-1]
    α[2] = e_col[1]

    for k in 3:K
        if isodd(k)
            # New q-column at level kq = (k+1)/2 from current (q_col, e_col).
            #   q^{(j-1)}_{kq} = e^{(j)}_{kq-1} / e^{(j-1)}_{kq-1} · q^{(j)}_{kq-1}
            new_len = K - k + 1
            q_col_new = Vector{T}(undef, new_len)
            @inbounds for j in 1:new_len
                e_jm1 = e_col[j]
                iszero(e_jm1) && throw(ArgumentError(
                    "pade_cf: zero divisor in qd table (e column at level $((k-1) >> 1), row $j); use `pade()` instead"))
                q_col_new[j] = e_col[j+1] / e_jm1 * q_col[j+1]
            end
            q_col = q_col_new
            α[k] = q_col[1]
        else
            # New e-column at level ke = k/2 from current (q_col, e_col).
            #   e^{(j-1)}_{ke} = q^{(j)}_{ke} - q^{(j-1)}_{ke} + e^{(j)}_{ke-1}
            new_len = K - k + 1
            e_col_new = Vector{T}(undef, new_len)
            @inbounds for j in 1:new_len
                e_col_new[j] = q_col[j+1] - q_col[j] + e_col[j+1]
            end
            e_col = e_col_new
            α[k] = e_col[1]
        end
    end

    return α
end
