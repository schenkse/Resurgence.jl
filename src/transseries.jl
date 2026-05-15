# Trans-series Σⱼ e^{-Sⱼ/g} g^{βⱼ} (Σₖ aⱼₖ gᵏ) — the natural object of
# resurgence theory. The j = 0 sector is the (typically divergent) perturbative
# answer; j > 0 sectors are exponentially suppressed non-perturbative
# contributions (instantons, multi-instantons, …). Arithmetic encodes the
# action-additive combinatorics behind dilute-instanton-gas constructions:
# multiplying two sectors adds their actions and convolves their perturbative
# tails.

"""
    Sector{T}(S, β, a)

A single trans-series sector representing `e^{-S/g} g^β · (Σₖ a[k+1] gᵏ)`.

`S` is the action / exponent in the exponential prefactor, `β` the power of
`g`, and `a` the vector of perturbative coefficients in `g` with the constant
term first. The outer constructor `Sector(S, β, a)` promotes its three
arguments to a common element type.
"""
struct Sector{T<:Number}
    S::T
    β::T
    a::Vector{T}
end

function Sector(S, β, a::AbstractVector)
    T = promote_type(typeof(S), typeof(β), eltype(a))
    return Sector{T}(T(S), T(β), Vector{T}(a))
end

"""
    TransSeries{T}

A trans-series stored as a `Vector{Sector{T}}`. Construct via

- `TransSeries(a::AbstractVector)` — wraps `a` as the single perturbative
  sector `Sector(0, 0, a)`.
- `TransSeries([Sector(S₁, β₁, a₁), Sector(S₂, β₂, a₂), …])` — explicit
  list of sectors; mixed element types are promoted.
"""
struct TransSeries{T<:Number}
    sectors::Vector{Sector{T}}
end

function TransSeries(sectors::AbstractVector{<:Sector})
    isempty(sectors) &&
        throw(ArgumentError("TransSeries needs at least one sector; use the parametric constructor TransSeries{T}(Sector{T}[]) for an explicit empty trans-series"))
    T = mapreduce(s -> typeof(s.S), promote_type, sectors)
    return TransSeries{T}(Sector{T}[_promote_sector(T, s) for s in sectors])
end

function TransSeries(a::AbstractVector{T}) where {T<:Number}
    return TransSeries{T}([Sector{T}(zero(T), zero(T), Vector{T}(a))])
end

_promote_sector(::Type{T}, s::Sector{T}) where {T<:Number} = s
_promote_sector(::Type{T}, s::Sector) where {T<:Number} =
    Sector{T}(T(s.S), T(s.β), Vector{T}(s.a))

_ts_eltype(::TransSeries{T}) where {T} = T

function _convert_ts(::Type{T}, ts::TransSeries) where {T<:Number}
    _ts_eltype(ts) === T && return ts
    return TransSeries{T}(Sector{T}[_promote_sector(T, s) for s in ts.sectors])
end

# ---------- arithmetic ----------

function Base.:+(t1::TransSeries{T}, t2::TransSeries{T}) where {T<:Number}
    return TransSeries{T}(_merge_sectors(vcat(t1.sectors, t2.sectors)))
end

function Base.:+(t1::TransSeries, t2::TransSeries)
    T = promote_type(_ts_eltype(t1), _ts_eltype(t2))
    return _convert_ts(T, t1) + _convert_ts(T, t2)
end

function Base.:-(ts::TransSeries{T}) where {T<:Number}
    return TransSeries{T}(Sector{T}[Sector{T}(s.S, s.β, -s.a) for s in ts.sectors])
end

Base.:-(t1::TransSeries, t2::TransSeries) = t1 + (-t2)

function Base.:*(α::Number, ts::TransSeries{T}) where {T<:Number}
    R = promote_type(typeof(α), T)
    αR = R(α)
    return TransSeries{R}(Sector{R}[
        Sector{R}(R(s.S), R(s.β), αR .* s.a) for s in ts.sectors
    ])
end

Base.:*(ts::TransSeries, α::Number) = α * ts

function Base.:*(t1::TransSeries{T}, t2::TransSeries{T}) where {T<:Number}
    secs = Sector{T}[]
    for s1 in t1.sectors, s2 in t2.sectors
        push!(secs, _sector_mul(s1, s2))
    end
    return TransSeries{T}(_merge_sectors(secs))
end

function Base.:*(t1::TransSeries, t2::TransSeries)
    T = promote_type(_ts_eltype(t1), _ts_eltype(t2))
    return _convert_ts(T, t1) * _convert_ts(T, t2)
end

# ---------- private helpers ----------

function _sector_mul(s1::Sector{T}, s2::Sector{T}) where {T<:Number}
    return Sector{T}(s1.S + s2.S, s1.β + s2.β, _cauchy(s1.a, s2.a))
end

function _cauchy(a::Vector{T}, b::Vector{T}) where {T<:Number}
    (isempty(a) || isempty(b)) && return T[]
    L = length(a) + length(b) - 1
    r = zeros(T, L)
    @inbounds for i in eachindex(a)
        ai = a[i]
        iszero(ai) && continue
        for j in eachindex(b)
            r[i + j - 1] += ai * b[j]
        end
    end
    return r
end

function _pad_add(a::Vector{T}, b::Vector{T}) where {T<:Number}
    if length(a) ≥ length(b)
        r = copy(a)
        @inbounds for k in eachindex(b)
            r[k] += b[k]
        end
        return r
    else
        r = copy(b)
        @inbounds for k in eachindex(a)
            r[k] += a[k]
        end
        return r
    end
end

function _merge_sectors(secs::Vector{Sector{T}}) where {T<:Number}
    out = Sector{T}[]
    for s in secs
        idx = findfirst(t -> t.S == s.S && t.β == s.β, out)
        if idx === nothing
            push!(out, s)
        else
            old = out[idx]
            out[idx] = Sector{T}(old.S, old.β, _pad_add(old.a, s.a))
        end
    end
    return out
end

# ---------- exponentiation ----------

"""
    transseries_exp(ts::TransSeries; order = 4) -> TransSeries

Truncated Taylor expansion `1 + ts + ts²/2! + … + ts^order / order!`.

Useful for canonical "exp(small)" constructions in resurgence: e.g., taking a
single-instanton trans-series and exponentiating to obtain the
dilute-instanton-gas trans-series, which contains all multi-instanton
sectors up to `order` copies of the seed.
"""
function transseries_exp(ts::TransSeries{T}; order::Integer = 4) where {T<:Number}
    order ≥ 0 || throw(ArgumentError("order must be ≥ 0"))
    one_ts = TransSeries{T}([Sector{T}(zero(T), zero(T), T[one(T)])])
    result = one_ts
    term = one_ts
    fact = one(T)
    for k in 1:order
        term = term * ts
        fact *= T(k)
        result = result + (one(T) / fact) * term
    end
    return result
end
