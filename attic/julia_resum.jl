# some resummation techniques for julia
using Polynomials
using PyCall
using SpecialFunctions

function shanks(a::Vector{T}, n::Int) where T <: Real
    An = sum(a[1:n]); Anp = sum(a[1:(n+1)]); Anm = sum(a[1:(n-1)])
    return Anp - (Anp - An)^2 / ((Anp - An) - (An - Anm))
end

@pyimport scipy.interpolate as scinterp

# numerator of order n, denominator of order m
function pade(a::Vector{T}, m::Int, n::Int) where T <: Real
    p, q = fit(RationalFunction, Polynomial(a), n, m)
    return p, q
    #p, q = scinterp.pade(a, m, n)
    #return Polynomial(reverse(p[:c])), Polynomial(reverse(q[:c]))
end

# pade function supporting BigFloats
function pade(a::Vector{BigFloat}, m::Int, n::Int)
    p, q = fit(RationalFunction, Polynomial(a), n, m)
    return p, q
end

function resum_pade(a::Vector{T}, n::Int, m::Int, x::Real=1) where T <: Real
    # numerator of order n, denominator of order m
    p, q = pade(a, m, n)
    # WARNING: use Polynomials.jl here
    #P = fit(RationalFunction, Polynomial(a), n, m)
    #P(x)
    p(x) / q(x)
end

using QuadGK

# helper function to normalize the first non-vanishing coefficient of a series
# the series needs to have a leading term 1 + ...
function normalized_series_odd(a::Vector{T}, iind::Int=1) where T <: Real
    # remove the leading entry
    an = copy(a)
    popfirst!(an)
    an = an ./ a[2]

    an
end

function borel_transform(a::Vector{T}) where T <: Real
    # need to use BigInts if the order is too large
    if length(a) > 21
        factorial_div = [1 / convert(Float64, factorial(big(k))) for k=0:length(a)-1]
    else
        factorial_div = [1 / factorial(k) for k=0:length(a)-1]
    end
    Ba = factorial_div .* a
    
    Ba
end

# Borel-Le Roy transform
function borel_lr_transform(a::Vector{T}, b::T) where T <: Real
    # need to use BigInts if the order is too large?
    gamma_div = [1 / gamma(k + 1 + b) for k=0:length(a)-1]
    Ba = gamma_div .* a

    Ba
end

ispos(x) = x > zero(typeof(x)) ? true : false

# another helper function
function chopv!(x::Vector{T}, tol::Float64=1e-16) where T <: Number
    for (i, el) in enumerate(x)
        rel, iel = reim(el)
        if abs(rel) <= tol
            rel = 0.0
        end
        if abs(iel) <= tol
            iel = 0.0
        end
        x[i] = iel == 0.0 ? rel : rel + 1im * iel
    end
    nothing
end


# helper function to move poles around
function move_poles!(x::Vector{T}, eps::Real) where T <: Number
    # warning this is by definition not type stable
    for (i, el) in enumerate(x)
        x[i] += isreal(el) && ispos(real(el)) ? eps * 1im : zero(typeof(el))
    end
    nothing
end

function obtain_poles_real(polycoeffs::Vector{T}) where T <: Number
    # determine the roots of the polynomial defined by polycoeffs
    roots = PolynomialRoots.roots(polycoeffs)
    real_poles = real.(filter(isreal, roots))
    real_poles
end

function poles_regularized(polycoeffs::Vector{T}, eps::Real) where T <: Number
    # determine the roots of the polynomial defined by polycoeffs
    poles = PolynomialRoots.roots(polycoeffs)
    move_poles!(poles, eps)
    poly_reg = fromroots(Polynomial, poles)
    poly_reg_coeffs = coeffs(poly_reg)
    # normalize the polynomial
    poly_reg_normalized_coeffs = poly_reg_coeffs ./ real(poly_reg_coeffs[1])
    poly_reg_normalized = Polynomial(poly_reg_normalized_coeffs)
    poly_reg_normalized
end

function resum_borel_pade(a::Vector{T}, n::Int, m::Int, x::Real=1;
                          return_error::Bool=false,
                          regularize_poles::Bool=false, eps::Real=1e-16) where T <: Real
    # Borel transform
    Ba = borel_transform(a[1:n+m+1])
    # numerator of order n, denominator of order m
    # obtain Pade approximant
    Bp, Bq = pade(Ba, m, n)
    # regularize the poles on the real axis before Laplace transforming
    # Laplace transform
    if regularize_poles
        Bqr = poles_regularized(coeffs(Bq), eps)
        int, int_err = quadgk(t -> Bp(x*t) / Bqr(x*t) * exp(-t), 0, Inf, order=24)
        # WARNING: this currently works only for x = 1
        #real_poles = obtain_poles_real(coeffs(Bq))
        # filter positive roots
        #int_domain = filter(ispos, sort(real_poles))
        #prepend!(int_domain, 0.0)
        #append!(int_domain, Inf)
        #int, int_err = quadgk(t -> Bp(x*t) / Bq(x*t) * exp(-t), int_domain..., order=24)
    else
        int, int_err = quadgk(t -> Bp(x*t) / Bq(x*t) * exp(-t), 0, Inf)
    end
    if return_error
        return int, int_err
    end
    int
end

function resum_borel_lr_pade(a::Vector{T}, n::Int, m::Int, x::Real=1, b::Real=-0.5;
                             return_error::Bool=false) where T <: Real
    # Borel transform
    Ba = borel_lr_transform(a[1:n+m+1], b)
    # numerator of order n, denominator of order m
    # obtain Pade approximant
    Bp, Bq = pade(Ba, m, n)
    # regularize the poles on the real axis before Laplace transforming
    # Laplace transform
    int, int_err = quadgk(t -> t^b * Bp(x*t) / Bq(x*t) * exp(-t), 0, Inf)
    if return_error
        return int, int_err
    end
    int
end

# helper function to split the vectors to match the python mpmath format
function split_vector(x::Vector{T}, n::Vector{Int64}) where T <: Number
    result = Vector{Vector{eltype(x)}}()
    start = firstindex(x)
    for len in n
        push!(result, x[start:(start + len - 1)])
        start += len
    end
    return result
end

function split_vector(x::Vector{T}, n::Int64) where T <: Number
    result = Vector{Vector{eltype(x)}}()
    lenx = length(x)
    start = firstindex(x)
    push!(result, x[start:(start + n - 1)])
    push!(result, x[n+1:lenx])
    return result
end

# convert to pyvectors for pycall
function split_vector_py(x::Vector{T}, n::Int64) where T <: Number
    result = Vector{PyVector{eltype(x)}}()
    lenx = length(x)
    start = firstindex(x)
    push!(result, PyVector(x[start:(start + n - 1)]))
    push!(result, PyVector(x[n+1:lenx]))
    return PyVector(result)
end

sparsify!(x, eps::Float64=1e-14) = x[abs(x) .< eps] = zero(eltype(x));

@pyimport mpmath as mpmath

# Meijer-G function as a wrapper around mpmath
function meijerg(as::Vector{T}, bs::Vector{T}, m::Int, n::Int, p::Int, q::Int, z::Number; series::Int=0) where T <: Number
    # convert input to match pycall mpmath format
    as_py = split_vector_py(as, n)
    bs_py = split_vector_py(bs, m)
    # sometimes it is useful to specify series specifically
    meijerg_o = series != 0 ? mpmath.meijerg(as_py, bs_py, z, series=series) : mpmath.meijerg(as_py, bs_py, z)
    # convert back to julia types
    meijerg_value = convert(Complex, meijerg_o)
    if imag(meijerg_value) == 0.0
        meijerg_value = real(meijerg_value)
    end
    return meijerg_value
end

# uses coefficients of Borel transform as input
function borel_ratios(b::Vector{T}) where T <: Real
    r = similar(b, Float64)
    popfirst!(r)
    for (n, bn) in enumerate(b)
        if n == length(b)
            break
        end
        r[n] = b[n+1] / bn
    end
    return r
end

# helper function to compute ps and qs
function construct_A(r::Vector{T}, N::Int) where T <: Number
    l = isodd(N) ? Int((N-1)/2) : Int(N/2)
    dimA = 2*l+1
    A = Matrix{eltype(r)}(undef, dimA, dimA)
    for i=1:dimA, j=1:dimA
        if j <= l+1
            A[i,j] = (i-1)^(j-1)
        else
            A[i,j] = -r[i] * (i-1)^(j-1 - l)
        end
    end
    A
end

# Pade approximate the Borel ratios
# use Borel ratios as input
# @WARNING: length of r must be odd!!!
function determine_pq(r::Vector{T}) where T <: Real
    N = length(r); l = isodd(N) ? Int((N-1)/2) : Int(N/2)
    A = construct_A(r, N)
    # solve the linear system
    pqv = A \ r
    chopv!(pqv)
    istart = firstindex(pqv); iend = lastindex(pqv)
    p = pqv[istart:l+1]; q = l+2 == iend ? [pqv[iend]] : pqv[(l+2):iend]
    p, q
end

# determine the polynomial roots
using PolynomialRoots

function determine_xyvec(p::Vector{T}, q::Vector{T}) where T <: Real
    xs = PolynomialRoots.roots(p); ys = PolynomialRoots.roots([one(eltype(q)), q...])
    chopv!(xs); chopv!(ys)
    xvec = -xs; yvec = -ys
    # add the first entry manually
    prepend!(xvec, one(eltype(xs)))

    return xvec, yvec
end

function resum_meijerg(a::Vector{T}, n::Int, x::Number=1; series::Int=0) where T <: Real
    # Borel transform
    Ba = borel_transform(a)

    # determine Borel ratios
    if iseven(n)
        Ba_odd = normalized_series_odd(Ba)
        rBa = borel_ratios(Ba_odd)[1:n-1]
    else
        rBa = borel_ratios(Ba)[1:n]
    end
    # determine pms and qms
    pv, qv = determine_pq(rBa)
    pl = last(pv); ql = last(qv)
    # determine xvec and yvec
    xv, yv = determine_xyvec(pv, qv)
    l = length(xv) - 1
    # determine the prefactor
    if l > 1
        gamma_ratio = prod(gamma, yv[1:l]) / prod(gamma, xv[2:(l+1)])
    else
        gamma_ratio = gamma(yv[1]) / gamma(xv[2])
    end
    # extend the xvector once again
    prepend!(xv, one(eltype(xv)))
    prepend!(yv, one(eltype(yv)))

    meijerg_res = gamma_ratio * meijerg(yv, xv, l+2, 1, l+1, l+2, -ql / (pl * x); series=series)
    if iseven(n)
        meijerg_res *= Ba[2] * x
        meijerg_res += 1
    end

    meijerg_res
end