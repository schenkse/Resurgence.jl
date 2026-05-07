"""
    obtain_poles_real(coeffs)

Return the real-valued roots of the polynomial whose coefficients are
`coeffs` (in `Polynomials.jl` ordering: constant term first).
"""
function obtain_poles_real(coeffs::AbstractVector{<:Number})
    rts = PolynomialRoots.roots(coeffs)
    return real.(filter(isreal, rts))
end

"""
    move_poles(roots, ε; side = +1, real_tol = 100·eps(real(T)))

Return a new vector of complex roots obtained from `roots` by shifting any root
that lies on the positive real axis off the real axis by `side · iε`. With
`side = +1` (default) the shift is `+iε` (upper half-plane); with `side = −1`
the shift is `−iε` (lower half-plane). Non-mutating and type-honest: the
result is always `Vector{Complex{R}}` where `R = real(eltype(roots))`.

Roots with `|imag(z)| ≤ real_tol` are treated as real (with their imaginary
noise zeroed) before classification. This guards against numerical noise from
`PolynomialRoots.roots` on rank-deficient Padé denominators, which can return
strictly real roots with `|imag|` of order `eps`. Clean projection preserves
the conjugate-pair symmetry of the original polynomial — required for the
lateral Borel sums `borel_pade_lateral(±1)` to be complex conjugates of each
other when the input series is real.
"""
function move_poles(roots::AbstractVector{T}, ε::Real;
                    side::Integer = +1,
                    real_tol::Real = 100 * eps(real(T))) where {T<:Number}
    side == 1 || side == -1 ||
        throw(ArgumentError("side must be +1 or -1 (got $side)"))
    R = real(T)
    εs = R(side) * R(ε)
    rtol = R(real_tol)
    out = Vector{Complex{R}}(undef, length(roots))
    @inbounds for i in eachindex(roots)
        z = roots[i]
        re, im = real(z), imag(z)
        if abs(im) ≤ rtol
            out[i] = re > zero(R) ?
                Complex{R}(re, εs) :
                Complex{R}(re, zero(R))
        else
            out[i] = Complex{R}(z)
        end
    end
    return out
end

"""
    poles_regularized(coeffs, ε; side = +1) -> Polynomial

Take a polynomial `q(t)` (given by `coeffs`, constant term first), shift any
positive-real roots by `side · iε` to push them off the integration contour,
and return the resulting polynomial as a `Polynomial`, normalized so its
constant term is `1`. `side ∈ {+1, −1}` selects the upper / lower half-plane
shift; the two sides correspond to the two lateral Borel sums.
"""
function poles_regularized(coeffs::AbstractVector{<:Number}, ε::Real;
                           side::Integer = +1)
    rts = PolynomialRoots.roots(coeffs)
    shifted = move_poles(rts, ε; side = side)
    poly = fromroots(Polynomial, shifted)
    pc = Polynomials.coeffs(poly)
    # Normalize so q(0) = 1; constant term sits at index 1.
    return Polynomial(pc ./ pc[1])
end
