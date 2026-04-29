"""
    pade(a, m, n) -> (p::Polynomial, q::Polynomial)

Return the Padé approximant `[n/m]` of the formal power series with coefficients
`a` (constant term first): a numerator polynomial `p` of degree `n` and a
denominator polynomial `q` of degree `m` such that

    p(z) / q(z) = a₀ + a₁ z + a₂ z² + … + O(z^{n+m+1})

Requires `length(a) ≥ n + m + 1`.

This is a thin wrapper around `Polynomials.fit(RationalFunction, …)` and
inherits its element-type promotion rules — `Float64` series produce `Float64`
polynomials, `BigFloat` series produce `BigFloat` polynomials.
"""
function pade(a::AbstractVector{T}, m::Integer, n::Integer) where {T<:Number}
    n ≥ 0 || throw(ArgumentError("n must be ≥ 0"))
    m ≥ 0 || throw(ArgumentError("m must be ≥ 0"))
    length(a) ≥ n + m + 1 ||
        throw(ArgumentError("Padé [$n/$m] needs length(a) ≥ $(n+m+1), got $(length(a))"))
    p, q = fit(RationalFunction, Polynomial(collect(a)), n, m)
    return p, q
end

"""
    pade_value(a, m, n, x)

Evaluate the Padé approximant `[n/m]` of `a` at `x`, i.e. return `p(x) / q(x)`
where `(p, q) = pade(a, m, n)`.
"""
function pade_value(a::AbstractVector{T}, m::Integer, n::Integer, x) where {T<:Number}
    p, q = pade(a, m, n)
    return p(x) / q(x)
end
