# Meijer-G resummation.
#
# The classical Borel-summation route: fit the consecutive ratios of the Borel
# transform to a rational P(k)/Q(k), identify the resulting Borel transform as
# a generalized hypergeometric series, then Borel–Laplace-transform to the
# resummed value. The result is canonically a Meijer G-function; under the
# standard Slater identity it collapses to a single ₚFq evaluation, which
# `HypergeometricFunctions.pFq` already handles for divergent inputs via the
# Weniger sequence transformation.

"""
    _normalized_series_odd(a) -> Vector

Drop the leading entry of `a` and rescale by the new leading entry, so the
result starts with `1`. This is the "shifted Borel transform" used by the
even-`n` branch of [`borel_meijerg`](@ref).

Errors if `a[2] == 0`. Generic over the element type.
"""
function _normalized_series_odd(a::AbstractVector{T}) where {T<:Number}
    length(a) ≥ 2 || throw(ArgumentError("_normalized_series_odd needs length(a) ≥ 2"))
    iszero(a[2]) && throw(ArgumentError("_normalized_series_odd: a[2] is zero, cannot normalize"))
    return @views a[2:end] ./ a[2]
end

"""
    _construct_meijerg_matrix(r, N) -> Matrix

Build the `(2l+1)×(2l+1)` linear system used by [`_determine_pq`](@ref) to fit
the Borel ratios `r` to a rational `P(k)/Q(k)` of degree `≤ l`, where
`l = (N-1) ÷ 2`. Columns `1..l+1` carry the numerator basis `(k)^j`; columns
`l+2..2l+1` carry `-r[i]·(k)^j`.
"""
function _construct_meijerg_matrix(r::AbstractVector{T}, N::Integer) where {T<:Number}
    l = (N - 1) ÷ 2
    dim = 2l + 1
    A = Matrix{T}(undef, dim, dim)
    @inbounds for i in 1:dim, j in 1:dim
        if j ≤ l + 1
            A[i, j] = T(i - 1)^(j - 1)
        else
            A[i, j] = -r[i] * T(i - 1)^(j - 1 - l)
        end
    end
    return A
end

"""
    _determine_pq(r) -> (p, q)

Fit `r[i] = P(i-1) / Q(i-1)` where `P(k) = p[1] + p[2]·k + … + p[l+1]·k^l` and
`Q(k) = 1 + q[1]·k + … + q[l]·k^l`, with `l = (length(r) - 1) ÷ 2`. The length
of `r` must be odd. Returns `(p, q)` as separate vectors (`p` length `l+1`,
`q` length `l`).

Falls back to `pinv` for rank-deficient systems (the LU fast path throws
`SingularException`), so degenerate fits return a reduced-degree approximant
rather than crashing — same pattern as [`pade`](@ref).
"""
function _determine_pq(r::AbstractVector{T}) where {T<:Number}
    N = length(r)
    isodd(N) ||
        throw(ArgumentError("_determine_pq requires odd-length input (got length $N)"))
    l = (N - 1) ÷ 2
    A = _construct_meijerg_matrix(r, N)
    pqv = try
        A \ collect(r)
    catch e
        e isa LinearAlgebra.SingularException || rethrow()
        pinv(A) * collect(r)
    end
    chop!(pqv)
    return pqv[1:l+1], pqv[l+2:N]
end

"""
    borel_meijerg(a; n, x = 1, pFq_kwargs...)

Meijer-G resummation of the formal power series with coefficients `a`,
evaluated at argument `x`.

# Algorithm

1. Take the Borel transform `B[k] = a[k] / k!`.
2. Fit the consecutive Borel ratios `B[k+1]/B[k]` to a rational `P(k)/Q(k)`
   using a `(2l+1)×(2l+1)` linear system, with `l = (n-1) ÷ 2` for odd `n`
   and `l = (n-2) ÷ 2` for even `n` (the even branch first strips and
   normalises by `a[2]`).
3. Extract `α = -roots(P)` and `β = -roots(Q)`. The Borel transform is then
   identified with `_{l+1}F_l(α, 1; β; γt)` with `γ = leading(P)/leading(Q)`.
4. The Borel sum is the Laplace integral, which by the Slater identity is

       a[1] · _{l+2}F_l(α, 1, 1; β; γx)               (odd n)
       a[1] + a[2]·x · _{l+2}F_l(α, 1, 2; β; γx)     (even n)

   evaluated by [`HypergeometricFunctions.pFq`](@ref). For divergent
   arguments (which is the typical case here, since `p > q + 1`) `pFq`
   uses the Weniger sequence transformation internally.

# Arguments

- `n`: number of Borel-ratio terms used for the rational fit. Must be `≥ 3`.
  Odd `n` and `n+1` use the same number of ratios `2l+1` with `l = (n-1)÷2`.
- `x`: evaluation point. Must be nonzero.
- `pFq_kwargs`: forwarded to `HypergeometricFunctions.pFq` (e.g. for tuning
  the Weniger transformation).

No `return_error` keyword: unlike the `quadgk`-based methods, `pFq` does not
expose a residual estimate.

Requires `length(a) ≥ n + 1`.
"""
function borel_meijerg(a::AbstractVector{T};
                       n::Integer, x = 1,
                       pFq_kwargs...) where {T<:Number}
    n ≥ 3 ||
        throw(ArgumentError("borel_meijerg needs n ≥ 3 (got $n) — fewer ratios make the fit degenerate"))
    length(a) ≥ n + 1 ||
        throw(ArgumentError("borel_meijerg needs length(a) ≥ n+1 = $(n+1) (got $(length(a)))"))
    iszero(x) && throw(ArgumentError("borel_meijerg: x must be nonzero (got 0)"))

    Ba = borel_transform(a)

    if iseven(n)
        Ba_use = _normalized_series_odd(Ba)
        rBa = borel_ratios(Ba_use)[1:n-1]
    else
        rBa = borel_ratios(Ba)[1:n]
    end

    pv, qv = _determine_pq(rBa)
    pl = last(pv)
    ql = last(qv)
    iszero(pl) && throw(ArgumentError("borel_meijerg: P polynomial is degenerate (leading coefficient zero)"))
    iszero(ql) && throw(ArgumentError("borel_meijerg: Q polynomial is degenerate (leading coefficient zero)"))

    R = real(T)
    αs = -PolynomialRoots.roots(pv)
    qpoly = vcat(one(eltype(qv)), qv)
    βs = -PolynomialRoots.roots(qpoly)

    γ = pl / ql

    # Build pFq parameter vectors. Promote everything through Complex{R} since
    # PolynomialRoots returns complex roots regardless of input.
    CT = Complex{R}
    α_full = Vector{CT}(undef, length(αs) + 2)
    @inbounds for i in eachindex(αs)
        α_full[i] = CT(αs[i])
    end
    α_full[end-1] = one(CT)
    α_full[end] = iseven(n) ? CT(2) : one(CT)
    β_full = Vector{CT}(undef, length(βs))
    @inbounds for i in eachindex(βs)
        β_full[i] = CT(βs[i])
    end

    z_arg = CT(γ * x)
    pfq_val = HypergeometricFunctions.pFq(α_full, β_full, z_arg; pFq_kwargs...)

    raw = if iseven(n)
        a[1] + a[2] * x * pfq_val
    else
        a[1] * pfq_val
    end

    # Real input ⇒ commit to a real return type. The pFq path goes through
    # complex parameters (PolynomialRoots.roots is complex-valued), but for
    # a real series at a real x the answer is mathematically real and the
    # imaginary residue is at machine epsilon.
    if T <: Real && x isa Real
        return real(raw)
    else
        return raw
    end
end
