# Stokes / large-order coefficient extraction.
#
# Given a divergent perturbative series with coefficients `a[k+1]` (k = 0, 1, …)
# and the expected asymptotic
#
#     a[k+1] ∼ A · Γ(k + β) / S^{k + β} · (1 + c₁/k + c₂/k² + …)
#
# we extract the instanton action `S`, the exponent `β`, the prefactor `A`,
# and (optionally) subleading coefficients `c₁ … c_subleading`. Each of `S`,
# `β`, `A` is read off from the limit of a sequence built from consecutive
# coefficient ratios:
#
#     rₖ  = a[k+2] / (k · a[k+1])             →  1/S    + O(1/k)
#     βₖ  = S · a[k+2]/a[k+1] − k             →  β      + O(1/k)
#     Aₖ  = a[k+1] · S^{k+β} / Γ(k + β)       →  A      + O(1/k)
#
# Each sequence has an asymptotic expansion in powers of 1/k, so Richardson
# extrapolation eliminates the leading tails. We reuse the partial-sum-based
# `richardson` from `series_acceleration.jl` by converting each sequence
# `s` to its first differences `b` (with `b[1] = s[1]`) so that
# `sum(b[1:n]) = s[n]`.

"""
    stokes_action(a; depth = nothing, n = nothing)

Extract the instanton action `S` from the large-order behaviour of the
formal-power-series coefficients `a`. Builds the sequence

    rₖ = a[k+2] / (k · a[k+1]),    k = 1, 2, …, length(a) - 2,

which converges to `1/S` with `O(1/k)` corrections, and Richardson-extra-
polates the tail.

`depth` is the number of Richardson levels (default
`min(5, (length(a) - 2) ÷ 2 - 1)`); `n` is the tail anchor (default: the
largest index that admits the chosen depth). For real `a` and a
sign-alternating series the result is negative; for a non-Borel-summable
series with positive-real-axis singularity it is positive.
"""
function stokes_action(a::AbstractVector{T};
                       depth::Union{Integer,Nothing} = nothing,
                       n::Union{Integer,Nothing} = nothing) where {T<:Number}
    K = length(a) - 2
    K ≥ 2 || throw(ArgumentError("stokes_action: need length(a) ≥ 4"))
    r = _stokes_ratio_seq(a)
    d, idx = _resolve_depth_n(K, depth, n)
    inv_S = _richardson_seq(r, idx; depth = d)
    return one(T) / inv_S
end

"""
    stokes_exponent(a, S; depth = nothing, n = nothing)

Extract the exponent `β` once `S` is known. Builds

    βₖ = S · a[k+2]/a[k+1] − k,

which tends to `β` with `O(1/k)` corrections, and Richardson-extrapolates.
"""
function stokes_exponent(a::AbstractVector{T}, S::Number;
                         depth::Union{Integer,Nothing} = nothing,
                         n::Union{Integer,Nothing} = nothing) where {T<:Number}
    K = length(a) - 2
    K ≥ 2 || throw(ArgumentError("stokes_exponent: need length(a) ≥ 4"))
    U = promote_type(T, typeof(S))
    βs = Vector{U}(undef, K)
    @inbounds for k in 1:K
        iszero(a[k+1]) && throw(ArgumentError("stokes_exponent: a[$(k+1)] is zero"))
        βs[k] = U(S) * (U(a[k+2]) / U(a[k+1])) - U(k)
    end
    d, idx = _resolve_depth_n(K, depth, n)
    return _richardson_seq(βs, idx; depth = d)
end

"""
    stokes_constant(a, S, β; depth = nothing, n = nothing)

Extract the prefactor `A` once `S` and `β` are known. Builds

    Aₖ = a[k+1] · S^{k+β} / Γ(k + β),

which tends to `A` with `O(1/k)` corrections, and Richardson-extrapolates.
"""
function stokes_constant(a::AbstractVector{T}, S::Number, β::Number;
                         depth::Union{Integer,Nothing} = nothing,
                         n::Union{Integer,Nothing} = nothing) where {T<:Number}
    K = length(a) - 2
    K ≥ 2 || throw(ArgumentError("stokes_constant: need length(a) ≥ 4"))
    U = promote_type(T, typeof(S), typeof(β))
    As = Vector{U}(undef, K)
    @inbounds for k in 1:K
        kplusβ = U(k) + U(β)
        As[k] = U(a[k+1]) * U(S)^kplusβ / gamma(kplusβ)
    end
    d, idx = _resolve_depth_n(K, depth, n)
    return _richardson_seq(As, idx; depth = d)
end

"""
    stokes_fit(a; subleading = 0, depth = nothing, n = nothing)

One-call combined extraction. Returns a `NamedTuple` `(S, β, A, c)` where
`c` is a vector of subleading coefficients of length `subleading`
(possibly empty) fit by least squares to the residual

    a[k+1] · S^{k+β} / (A · Γ(k + β)) − 1  ≈  c₁/k + c₂/k² + … + c_subleading / k^subleading

over a tail window of the available coefficients.
"""
function stokes_fit(a::AbstractVector{T};
                    subleading::Integer = 0,
                    depth::Union{Integer,Nothing} = nothing,
                    n::Union{Integer,Nothing} = nothing) where {T<:Number}
    subleading ≥ 0 || throw(ArgumentError("subleading must be ≥ 0 (got $subleading)"))
    S = stokes_action(a; depth = depth, n = n)
    β = stokes_exponent(a, S; depth = depth, n = n)
    A = stokes_constant(a, S, β; depth = depth, n = n)
    U = promote_type(T, typeof(S), typeof(β), typeof(A))
    c = subleading == 0 ?
        Vector{U}(undef, 0) :
        _stokes_subleading(a, S, β, A, subleading)
    return (S = S, β = β, A = A, c = c)
end

# === helpers ============================================================

function _stokes_ratio_seq(a::AbstractVector{T}) where {T<:Number}
    K = length(a) - 2
    r = Vector{T}(undef, K)
    @inbounds for k in 1:K
        iszero(a[k+1]) && throw(ArgumentError("stokes ratio: a[$(k+1)] is zero"))
        r[k] = a[k+2] / (T(k) * a[k+1])
    end
    return r
end

function _resolve_depth_n(K::Int,
                          depth::Union{Integer,Nothing},
                          n::Union{Integer,Nothing})
    d = depth === nothing ? max(1, min(5, K ÷ 2 - 1)) : Int(depth)
    d ≥ 1 || throw(ArgumentError("depth must be ≥ 1 (got $d)"))
    idx = n === nothing ? K - d : Int(n)
    idx ≥ 1 || throw(ArgumentError("n must be ≥ 1 (got $idx)"))
    idx + d ≤ K ||
        throw(ArgumentError("n + depth = $(idx + d) exceeds available range K = $K"))
    return d, idx
end

# Sequence-level Richardson: convert s to first differences b so that
# sum(b[1:n]) = s[n], then call the partial-sum-based `richardson`.
function _richardson_seq(s::AbstractVector{T}, n::Integer;
                         depth::Integer = 1) where {T<:Number}
    b = similar(s)
    @inbounds begin
        b[1] = s[1]
        for i in 2:length(s)
            b[i] = s[i] - s[i-1]
        end
    end
    return richardson(b, n; depth = depth)
end

function _stokes_subleading(a::AbstractVector{T}, S::Number, β::Number,
                            A::Number, subleading::Integer) where {T<:Number}
    K = length(a) - 1
    Nfit = max(2 * subleading, K ÷ 2)
    Nfit > K && (Nfit = K)
    kstart = K - Nfit + 1
    U = promote_type(T, typeof(S), typeof(β), typeof(A))
    M = Matrix{U}(undef, Nfit, subleading)
    rhs = Vector{U}(undef, Nfit)
    @inbounds for (i, k) in enumerate(kstart:K)
        kplusβ = U(k) + U(β)
        ρ = U(a[k+1]) * U(S)^kplusβ / (U(A) * gamma(kplusβ)) - one(U)
        rhs[i] = ρ
        for j in 1:subleading
            M[i, j] = one(U) / U(k)^j
        end
    end
    return M \ rhs
end
