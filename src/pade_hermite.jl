"""
    hermite_pade(a, n, m, l) -> (P::Polynomial, Q::Polynomial, R::Polynomial)

Hermite / quadratic Padé approximant of the formal power series with
coefficients `a` (constant term first): polynomials `P, Q, R` of degree at
most `n, m, l` such that

    P(z) + Q(z) f(z) + R(z) f(z)² = O(z^{n+m+l+2})

normalised so that `q_0 = 1`. Captures algebraic (square-root) branch
points of the Borel transform that linear Padé cannot represent.

The denominator-tail / numerator-then-back substitution mirrors the linear
case in [`pade`](@ref): the `(m + l + 1) × (m + l + 1)` system for
`(q_1, …, q_m, r_0, …, r_l)` is solved via LU; on
`LinearAlgebra.SingularException` it falls back to `pinv` (SVD), the same
rank-deficient strategy used by `pade`.

Requires `length(a) ≥ n + m + l + 2`.
"""
function hermite_pade(a::AbstractVector{T}, n::Integer, m::Integer, l::Integer) where {T<:Number}
    n ≥ 0 || throw(ArgumentError("n must be ≥ 0"))
    m ≥ 0 || throw(ArgumentError("m must be ≥ 0"))
    l ≥ 0 || throw(ArgumentError("l must be ≥ 0"))
    length(a) ≥ n + m + l + 2 || throw(ArgumentError(
        "Hermite-Padé [$n/$m/$l] needs length(a) ≥ $(n+m+l+2), got $(length(a))"))

    N = n + m + l + 2
    # Cauchy square: g[k] = (f²)_{k-1} = Σ_{j=1}^{k} a[j] · a[k-j+1] for k = 1..N
    g = Vector{T}(undef, N)
    @inbounds for k in 1:N
        s = zero(T)
        for j in 1:k
            s += a[j] * a[k - j + 1]
        end
        g[k] = s
    end

    rows = m + l + 1
    A = zeros(T, rows, rows)
    rhs = Vector{T}(undef, rows)
    @inbounds for r in 1:rows
        k = n + r
        # q-block columns 1..m correspond to unknowns q_1..q_m; coefficient
        # of q_j in the order-k matching equation is a_{k-j} = a[k-j+1].
        for j in 1:m
            idx = k - j + 1
            A[r, j] = idx ≥ 1 ? a[idx] : zero(T)
        end
        # r-block columns m+1..m+l+1 correspond to unknowns r_0..r_l;
        # coefficient of r_{j-1} is g_{k-(j-1)} = g[k-j+2].
        for j in 1:l+1
            idx = k - j + 2
            A[r, m + j] = idx ≥ 1 ? g[idx] : zero(T)
        end
        rhs[r] = -a[k + 1]
    end

    sol = try
        A \ rhs
    catch e
        e isa LinearAlgebra.SingularException || rethrow()
        pinv(A) * rhs
    end

    qcoeffs = vcat(one(T), sol[1:m])
    rcoeffs = collect(sol[m + 1 : m + l + 1])

    # P_k = -(Q·f)_k - (R·f²)_k for k = 0..n  (the order-k matching condition).
    pcoeffs = Vector{T}(undef, n + 1)
    @inbounds for k in 0:n
        s = zero(T)
        for j in 0:min(m, k)
            s -= qcoeffs[j + 1] * a[k - j + 1]
        end
        for j in 0:min(l, k)
            s -= rcoeffs[j + 1] * g[k - j + 1]
        end
        pcoeffs[k + 1] = s
    end

    return Polynomial(pcoeffs), Polynomial(qcoeffs), Polynomial(rcoeffs)
end

"""
    hermite_pade_value(a, n, m, l, x; branch = nothing)

Evaluate the Hermite / quadratic Padé approximant of `a` at `x` by solving
`R(x)·y² + Q(x)·y + P(x) = 0` and selecting one branch.

Branch selection:
- `branch = nothing` (default): pick the principal branch that matches
  `f(0) = a[1]`. The discriminant at `x = 0` is identically
  `(1 + 2·R(0)·a[1])²`, so the sign is determined analytically without
  evaluating the other root.
- `branch = +1` or `branch = -1`: fix the choice of the `±√(disc)` term
  manually. Useful for exploring the conjugate branch (e.g. the lateral
  Stokes copy of a square-root branch).

For real-typed `a` the discriminant `Q(x)² − 4·P(x)·R(x)` may be negative
at points beyond the radius of validity of the approximant; in that case
`sqrt` raises `DomainError`. Pass a complex `x` (e.g. `complex(x)`) to
follow the branch into the complex plane.
"""
function hermite_pade_value(a::AbstractVector{T}, n::Integer, m::Integer, l::Integer, x;
                            branch::Union{Integer,Nothing} = nothing) where {T<:Number}
    P, Q, R = hermite_pade(a, n, m, l)
    sgn = if branch === nothing
        ref = one(T) + 2 * R(zero(T)) * a[firstindex(a)]
        real(ref) ≥ 0 ? 1 : -1
    else
        (branch == 1 || branch == -1) ||
            throw(ArgumentError("branch must be +1, -1, or nothing"))
        Int(branch)
    end
    Px, Qx, Rx = P(x), Q(x), R(x)
    s = sqrt(Qx^2 - 4 * Px * Rx)
    return (-Qx + sgn * s) / (2 * Rx)
end
