"""
    optimal_truncation(a) -> (Nstar, partial_sum, smallest_term)

Smallest-term (superasymptotic) truncation of an asymptotic series with
coefficients `a`:

- `Nstar` is the index `k` minimising `|a[k]|`.
- `partial_sum` is `sum(a[1:Nstar])`.
- `smallest_term` is `|a[Nstar]|`, the standard textbook estimate of the
  superasymptotic remainder.

For a convergent or monotone series the answer is trivial (`Nstar = length(a)`
or `Nstar = 1`), but the function is still well-defined and returns the obvious
result.
"""
function optimal_truncation(a::AbstractVector{T}) where {T<:Number}
    isempty(a) && throw(ArgumentError("series must be non-empty"))
    Nstar = _argmin_abs(a)
    return Nstar, sum(@view a[1:Nstar]), abs(a[Nstar])
end

function _argmin_abs(a::AbstractVector)
    idx = firstindex(a)
    m = abs(a[idx])
    @inbounds for i in idx+1:lastindex(a)
        ai = abs(a[i])
        if ai < m
            m = ai
            idx = i
        end
    end
    return idx
end

"""
    superasymptotic_remainder(a)

Estimate the remainder of an asymptotic series at its optimal-truncation
order: returns `|a[Nstar]|`, the magnitude of the smallest-magnitude term.
This is the textbook *superasymptotic* error estimate.
"""
function superasymptotic_remainder(a::AbstractVector{T}) where {T<:Number}
    isempty(a) && throw(ArgumentError("series must be non-empty"))
    return abs(a[_argmin_abs(a)])
end

"""
    terminant(p, σ)

Berry–Howls terminant function

    T_p(σ) = exp(iπp) · Γ(p) · Γ(1 − p, σ) / (2πi)

where `Γ(s, z)` is the upper incomplete gamma. The terminant smoothly
crosses the Stokes line: for `σ ≫ p` it decays exponentially to zero, for
`σ ≪ p` it tends to one, and at the optimal-truncation regime `σ ≈ p` it
takes the value `≈ 1/2`.

Generic over real and complex inputs; promotes to a complex floating-point
type compatible with the inputs.
"""
function terminant(p::Number, σ::Number)
    R = float(real(promote_type(typeof(p), typeof(σ))))
    T = Complex{R}
    pT, σT = T(p), T(σ)
    return exp(im * R(π) * pT) * gamma(pT) * gamma(one(T) - pT, σT) / (2 * R(π) * im)
end

"""
    hyperasymptotic(a; x = 1, level = 1, kind = :subleading,
                    action = nothing, β = nothing, A = nothing,
                    c = nothing,
                    action_2 = nothing, β_2 = nothing, A_2 = nothing)

Hyperasymptotic resummation of the formal power series `Σₖ a[k+1] · xᵏ`.

`level = 0` returns the optimal-truncation partial sum evaluated at `x`
(the *superasymptotic* sum). `level = 1` adds the leading single-instanton
correction, lifting the remainder magnitude from `e^{-|S|/x}` to
`e^{-2|S|/x}` for problems with a single Stokes singularity at instanton
action `S` on the positive real `t`-axis.

The large-order asymptotics `a[n+1] ~ A · Γ(n+β) / action^{n+β}` are taken
verbatim if `(action, β, A)` are provided, or extracted via
[`stokes_fit`](@ref) otherwise. All three must be provided together or
none. The level-1 correction is the leading lateral-Borel piece,

    R_{N*}(x) ≈ -i π · A · x^{-β} · exp(-action / x),

which matches the imaginary part of the `side = +1` lateral Borel–Padé
sum (see [`borel_pade_lateral`](@ref)) — equivalent to picking the
`+iε` deformation of the Borel-Laplace contour around the leading
instanton singularity. The exported [`terminant`](@ref) function provides
the smoothed-Stokes-multiplier building block for users wanting to refine
the crossover further.

`level = 2` adds a second-order correction selected by `kind`:

- `kind = :subleading` (default) — next-to-leading 1/k term within the
  *same* Stokes singularity:

      R^{(2)}_{N*}(x) ≈ -i π · A · x^{-β} · exp(-S/x) · (1 + c₁ · x).

  The subleading coefficient `c₁` is taken from `c[1]` if `c` is provided,
  otherwise extracted via [`stokes_fit`](@ref) with `subleading = 1`. `c`
  can be a vector of any length ≥ 1 (only `c[1]` is used).

- `kind = :two_instanton` — add a *second* Stokes singularity:

      R^{(2)}_{N*}(x) ≈ -i π · A · x^{-β} · e^{-S/x}
                       - i π · A₂ · x^{-β₂} · e^{-S₂/x}.

  The second triple `(action_2, β_2, A_2)` is taken verbatim if supplied
  (all three together), otherwise extracted by [`stokes_fit`](@ref) on
  the *residual* series `a[k+1] - A · Γ(k+β) / S^{k+β}` — needs enough
  coefficients (`length(a) ≥ 12`) for the second extraction to be
  stable, and throws if the residual fit produces a near-degenerate
  `S₂ ≈ S`.

The `(action, β, A)` and `(c | action_2, β_2, A_2)` are independent
keyword groups: each enforces all-or-none internally, but you may
explicit-supply one group while auto-extracting the other.

The convention is `Re(action) > 0` (positive instanton action with the
Stokes line on the positive real `t`-axis). For Borel-summable problems
whose Stokes singularity sits on the *negative* real axis the level-1
formula is not applicable, and the function throws — pass an explicit
positive `action` only when the singularity is genuinely on the positive
real axis.

Levels `≥ 3` need late-of-late machinery or full trans-series input and
currently throw; the caller should drop down to [`resum_transseries`](@ref)
for those.
"""
function hyperasymptotic(a::AbstractVector{T}; x = 1, level::Integer = 1,
                         kind::Symbol = :subleading,
                         action::Union{Number,Nothing} = nothing,
                         β::Union{Number,Nothing} = nothing,
                         A::Union{Number,Nothing} = nothing,
                         c::Union{AbstractVector{<:Number},Nothing} = nothing,
                         action_2::Union{Number,Nothing} = nothing,
                         β_2::Union{Number,Nothing} = nothing,
                         A_2::Union{Number,Nothing} = nothing) where {T<:Number}
    isempty(a) && throw(ArgumentError("series must be non-empty"))
    level ∈ (0, 1, 2) || throw(ArgumentError("hyperasymptotic supports level ∈ {0, 1, 2}; \
        higher levels need late-of-late machinery or full trans-series input outside this \
        function's scope"))

    Nstar = _optimal_index_at(a, x)
    U = promote_type(T, typeof(x))
    pow = one(U)
    partial = U(a[1])
    @inbounds for k in 2:Nstar
        pow *= U(x)
        partial += U(a[k]) * pow
    end
    level == 0 && return partial

    nprovided = count(!isnothing, (action, β, A))
    nprovided ∈ (0, 3) ||
        throw(ArgumentError("provide either all of (action, β, A) or none; got $nprovided of 3"))
    Sv, βv, Av = if nprovided == 0
        f = stokes_fit(a)
        f.S, f.β, f.A
    else
        action, β, A
    end
    real(Sv) > 0 || throw(ArgumentError("hyperasymptotic level≥1 needs a Stokes singularity \
        with positive real `action` (got $Sv); for Borel-summable problems whose Stokes \
        singularity sits off the positive real axis the level-1 formula does not apply"))

    V = promote_type(U, typeof(Sv), typeof(βv), typeof(Av), ComplexF64)
    prefactor1 = -im * V(π) * V(Av) * V(x)^(-V(βv)) * exp(-V(Sv) / V(x))
    correction = prefactor1

    if level == 2
        if kind === :subleading
            cv = if c === nothing
                stokes_fit(a; subleading = 1).c
            else
                c
            end
            length(cv) ≥ 1 || throw(ArgumentError("hyperasymptotic level=2 :subleading needs \
                at least one subleading coefficient (c[1]); got empty `c`"))
            correction = prefactor1 * (one(V) + V(cv[1]) * V(x))
        elseif kind === :two_instanton
            S2, β2v, A2v = _resolve_two_instanton(a, Sv, βv, Av, action_2, β_2, A_2)
            prefactor2 = -im * V(π) * V(A2v) * V(x)^(-V(β2v)) * exp(-V(S2) / V(x))
            correction = prefactor1 + prefactor2
        else
            throw(ArgumentError("hyperasymptotic `kind` must be :subleading or :two_instanton \
                (got $(repr(kind)))"))
        end
    end

    return V(partial) + correction
end

# Resolve the second-singularity triple (S₂, β₂, A₂) for level=2 :two_instanton.
# Either uses the user-supplied (action_2, β_2, A_2) when all three are present,
# or extracts them by subtracting the leading large-order prediction
# A · Γ(k+β) / S^{k+β} from each a[k+1] and running stokes_fit on the residual.
function _resolve_two_instanton(a::AbstractVector{T}, S::Number, β::Number, A::Number,
                                action_2, β_2, A_2) where {T<:Number}
    nprovided = count(!isnothing, (action_2, β_2, A_2))
    if nprovided == 3
        return action_2, β_2, A_2
    elseif nprovided != 0
        throw(ArgumentError("provide either all of (action_2, β_2, A_2) or none; \
            got $nprovided of 3"))
    end
    length(a) ≥ 12 ||
        throw(ArgumentError("hyperasymptotic level=2 :two_instanton needs \
            length(a) ≥ 12 to extract the second Stokes triple from the residual \
            (got $(length(a))); supply (action_2, β_2, A_2) explicitly for shorter inputs"))
    U = promote_type(T, typeof(S), typeof(β), typeof(A))
    residual = Vector{U}(undef, length(a))
    @inbounds for k in 0:length(a)-1
        # a[k+1] ≈ A · Γ(k+β) / S^{k+β} + (subleading singularity contribution)
        residual[k+1] = U(a[k+1]) - U(A) * gamma(U(k) + U(β)) / U(S)^(U(k) + U(β))
    end
    f = stokes_fit(residual)
    abs(f.S - S) > abs(S) / 10 ||
        throw(ArgumentError("hyperasymptotic level=2 :two_instanton extracted a near-degenerate \
            second action S₂ ≈ S (got S₂=$(f.S), S=$S); the residual fit is unreliable. \
            Supply (action_2, β_2, A_2) explicitly."))
    return f.S, f.β, f.A
end

# Index k that minimises |a[k] · x^{k-1}|, i.e. optimal truncation aware
# of the coupling x. Falls back to plain |a[k]| when x = 1.
function _optimal_index_at(a::AbstractVector, x)
    idx = firstindex(a)
    U = promote_type(eltype(a), typeof(x))
    pow = one(U)
    m = abs(a[idx])
    @inbounds for i in idx+1:lastindex(a)
        pow *= U(x)
        ai = abs(U(a[i]) * pow)
        if ai < m
            m = ai
            idx = i
        end
    end
    return idx
end
