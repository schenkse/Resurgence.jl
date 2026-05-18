# High-level diagnostic / discovery helpers built on top of the per-method
# functions. These are aimed at lowering the entry barrier — picking a likely
# method for a given series, or running a small set of methods side-by-side —
# not at producing authoritative numbers (the underlying tools remain that).

"""
    diagnose(a) -> NamedTuple

Quick profile of a formal-power-series coefficient vector. Composes a tail
ratio inspection, [`stokes_fit`](@ref), and [`optimal_truncation`](@ref), and
tacks on a short ordered list of recommended resummation methods.

Returns a `NamedTuple` with fields:

- `growth ∈ (:convergent, :geometric, :factorial, :unknown)` —
  classification from `|a[k+1]/a[k]|` over the tail.
- `alternating::Bool` — every neighbouring pair in the tail window has
  opposite real-part sign.
- `S, β, A` — large-order data from `stokes_fit`, or `missing` if `a` is
  too short or `stokes_fit` throws.
- `Nstar, partial_opt, εN` — the `optimal_truncation` triple.
- `recommended::Vector{Symbol}` — ordered list of suggested per-method
  function names (e.g. `:borel_pade`, `:shanks`).

The classifier is heuristic; reach for the underlying tools
([`stokes_fit`](@ref), [`optimal_truncation`](@ref), [`borel_ratios`](@ref))
when you want the authoritative numbers.
"""
function diagnose(a::AbstractVector{T}) where {T<:Number}
    isempty(a) && throw(ArgumentError("diagnose: a must be non-empty"))
    growth      = _classify_growth(a)
    alternating = _is_alternating_tail(a)
    sk          = _safe_stokes_fit(a)
    Nstar, partial_opt, εN = optimal_truncation(a)
    recommended = _recommend(growth, alternating)
    return (
        growth      = growth,
        alternating = alternating,
        S           = sk.S,
        β           = sk.β,
        A           = sk.A,
        Nstar       = Nstar,
        partial_opt = partial_opt,
        εN          = εN,
        recommended = recommended,
    )
end

# Classify the tail ratios |a[k+1]/a[k]|. Geometric and factorial both have
# bounded-or-growing ratios; the distinguishing feature is whether the ratios
# themselves grow with k (factorial, r_k ∝ k) or stay flat (geometric).
function _classify_growth(a::AbstractVector{T}) where {T<:Number}
    n = length(a)
    n ≥ 5 || return :unknown
    # Use the last K ratios (at least 4 points, up to half the input).
    K = min(n - 1, max(4, n ÷ 2))
    r = Vector{Float64}(undef, K)
    @inbounds for i in 1:K
        k = n - K + i - 1            # index into a; ratio r_k = |a[k+1]|/|a[k]|
        ak = abs(a[k])
        iszero(ak) && return :unknown
        rk = float(abs(a[k+1]) / ak)
        isfinite(rk) || return :unknown
        r[i] = rk
    end
    m = sum(r) / K
    m < 0.95 && return :convergent
    # Factorial: ratios scale with k, so the ratio of ratios across the window
    # is appreciably > 1 (≈ k_end / k_start). Geometric: roughly flat (≈ 1).
    r[end] / r[1] > 1.2 ? :factorial : :geometric
end

# True iff every neighbouring pair in the tail window flips sign (using real
# part). For complex `a` this is a coarse rule of thumb, not a definition.
function _is_alternating_tail(a::AbstractVector{T}; window::Integer = 10) where {T<:Number}
    n = length(a)
    n ≥ 3 || return false
    k_start = max(1, n - window)
    @inbounds for k in k_start:(n - 1)
        x = real(a[k]); y = real(a[k+1])
        (iszero(x) || iszero(y)) && return false
        sign(x) == sign(y) && return false
    end
    return true
end

function _safe_stokes_fit(a::AbstractVector)
    length(a) < 4 && return (S = missing, β = missing, A = missing)
    try
        f = stokes_fit(a)
        return (S = f.S, β = f.β, A = f.A)
    catch
        return (S = missing, β = missing, A = missing)
    end
end

function _recommend(growth::Symbol, alternating::Bool)
    if growth === :convergent
        return [:shanks, :richardson, :pade]
    elseif growth === :factorial
        return alternating ?
            [:borel_pade, :conformal_borel_pade, :borel_leroy_pade_odm] :
            [:borel_pade_median, :borel_pade_lateral, :borel_pade_discontinuity]
    elseif growth === :geometric
        return [:shanks, :pade, :levin]
    else
        return [:optimal_truncation, :stokes_fit]
    end
end

"""
    compare(methods, a; reference = nothing) -> Vector{NamedTuple}

Run each `method` in `methods` via [`resum`](@ref) and tabulate the results.
`methods` is any iterable of [`AbstractResummation`](@ref) tags (vector,
tuple, generator, …).

Each row is a `NamedTuple`:

- `method::String`  — the tag's type name (e.g. `"BorelPade"`).
- `result`          — `resum(method, a)` or `missing` on failure.
- `residual`        — `abs(result - reference)` if `reference` is supplied
                      and the call succeeded, else `missing`.
- `error`           — `nothing` on success, an error-message `String` on
                      failure.

Per-method failures are caught and recorded so that one bad tag (e.g. a
singular Padé fit or a `length(a) < n+m+1` mismatch) doesn't abort the
survey.

# Example

```julia
rows = compare([BorelPade(10, 10), Pade(10, 10)], a; reference = 0.5963473623)
```
"""
function compare(methods, a; reference = nothing)
    rows = NamedTuple[]
    for m in methods
        label = string(nameof(typeof(m)))
        result = missing
        err = nothing
        try
            result = resum(m, a)
        catch e
            err = sprint(showerror, e)
        end
        residual = (reference === nothing || result === missing) ?
                   missing : abs(result - reference)
        push!(rows, (method = label, result = result,
                     residual = residual, error = err))
    end
    return rows
end
