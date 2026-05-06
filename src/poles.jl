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
    move_poles(roots, ε)

Return a new vector of complex roots obtained from `roots` by shifting any root
that lies on the positive real axis off into the upper half-plane by `+iε`.
Non-mutating and type-honest: the result is always `Vector{Complex{R}}` where
`R = real(eltype(roots))`.
"""
function move_poles(roots::AbstractVector{T}, ε::Real) where {T<:Number}
    R = real(T)
    out = Vector{Complex{R}}(undef, length(roots))
    @inbounds for i in eachindex(roots)
        z = roots[i]
        on_pos_real_axis = isreal(z) && real(z) > zero(R)
        out[i] = on_pos_real_axis ? Complex{R}(real(z), R(ε)) : Complex{R}(z)
    end
    return out
end

"""
    poles_regularized(coeffs, ε) -> Polynomial

Take a polynomial `q(t)` (given by `coeffs`, constant term first), shift any
positive-real roots by `+iε` to push them off the integration contour, and
return the resulting polynomial as a `Polynomial`, normalized so its constant
term is `1`.
"""
function poles_regularized(coeffs::AbstractVector{<:Number}, ε::Real)
    rts = PolynomialRoots.roots(coeffs)
    shifted = move_poles(rts, ε)
    poly = fromroots(Polynomial, shifted)
    pc = Polynomials.coeffs(poly)
    # Normalize so q(0) = 1; constant term sits at index 1.
    return Polynomial(pc ./ pc[1])
end
