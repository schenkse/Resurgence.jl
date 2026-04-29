"""
    pade(a, m, n) -> (p::Polynomial, q::Polynomial)

Return the Padé approximant `[n/m]` of the formal power series with coefficients
`a` (constant term first): a numerator polynomial `p` of degree `n` and a
denominator polynomial `q` of degree `m` such that

    p(z) / q(z) = a₀ + a₁ z + a₂ z² + … + O(z^{n+m+1})

Requires `length(a) ≥ n + m + 1`.
"""
function pade(a::AbstractVector{T}, m::Integer, n::Integer) where {T<:Number}
    n ≥ 0 || throw(ArgumentError("n must be ≥ 0"))
    m ≥ 0 || throw(ArgumentError("m must be ≥ 0"))
    length(a) ≥ n + m + 1 ||
        throw(ArgumentError("Padé [$n/$m] needs length(a) ≥ $(n+m+1), got $(length(a))"))

    if m == 0
        return Polynomial(collect(a[1:n+1])), Polynomial(T[one(T)])
    end

    A = Matrix{T}(undef, m, m)
    @inbounds for i in 1:m, j in 1:m
        idx = n + i - j + 1
        A[i, j] = idx ≥ 1 ? a[idx] : zero(T)
    end
    rhs = T[-a[n + i + 1] for i in 1:m]

    qtail = try
        A \ rhs
    catch e
        e isa LinearAlgebra.SingularException || rethrow()
        pinv(A) * rhs
    end

    qcoeffs = vcat(one(T), qtail)
    pcoeffs = Vector{T}(undef, n + 1)
    @inbounds for k in 0:n
        s = zero(T)
        for j in 0:min(k, m)
            s += qcoeffs[j+1] * a[k - j + 1]
        end
        pcoeffs[k+1] = s
    end
    return Polynomial(pcoeffs), Polynomial(qcoeffs)
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
