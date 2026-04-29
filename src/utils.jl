"""
    inv_factorials(::Type{T}, n) -> Vector{T}

Return `[1/0!, 1/1!, …, 1/(n-1)!]` in element type `T`, computed iteratively as
`f[k+1] = f[k] / k`. This avoids both the overflow of `factorial(::Int)` past
`k = 21` and unnecessary `BigInt` allocation, while staying generic over
`Float64`, `BigFloat`, and complex variants.
"""
function inv_factorials(::Type{T}, n::Integer) where {T<:Number}
    n ≥ 0 || throw(ArgumentError("n must be non-negative"))
    f = Vector{T}(undef, n)
    n == 0 && return f
    f[1] = one(T)
    for k in 1:n-1
        f[k+1] = f[k] / T(k)
    end
    return f
end

"""
    chop!(x, tol = eps(real(eltype(x))) * 100)

Zero out entries of `x` smaller than `tol` in absolute value. Two methods:

- For real-eltype `x`, simply zeros |x[i]| < tol.
- For complex-eltype `x`, zeros real and imaginary parts independently. The
  eltype is preserved (no silent demotion to real), keeping the function
  type-stable.
"""
function chop!(x::AbstractVector{T}, tol::Real = eps(real(T)) * 100) where {T<:Real}
    @inbounds for i in eachindex(x)
        if abs(x[i]) ≤ tol
            x[i] = zero(T)
        end
    end
    return x
end

function chop!(x::AbstractVector{T}, tol::Real = eps(real(T)) * 100) where {T<:Complex}
    R = real(T)
    @inbounds for i in eachindex(x)
        re, im = reim(x[i])
        if abs(re) ≤ tol
            re = zero(R)
        end
        if abs(im) ≤ tol
            im = zero(R)
        end
        x[i] = T(re, im)
    end
    return x
end

"""
    sparsify!(x, ε = eps(real(eltype(x))) * 100)

Element-wise variant of [`chop!`](@ref) returning the same vector. Kept as a
separate name for compatibility with the original code; behaves identically to
`chop!` for real and complex inputs.
"""
sparsify!(x::AbstractVector, ε::Real = eps(real(eltype(x))) * 100) = chop!(x, ε)

"""
    split_vector(x, n::Vector{Int}) -> Vector{Vector}
    split_vector(x, n::Int)        -> Vector{Vector}

Split `x` into contiguous chunks. With a vector of lengths, returns one chunk
per length. With a single `Int`, returns `[x[1:n], x[n+1:end]]`.

Used by Meijer-G parameter packing (see v0.2).
"""
function split_vector(x::AbstractVector{T}, n::AbstractVector{<:Integer}) where {T}
    out = Vector{Vector{T}}(undef, length(n))
    start = firstindex(x)
    for (i, len) in enumerate(n)
        out[i] = x[start:start+len-1]
        start += len
    end
    return out
end

function split_vector(x::AbstractVector{T}, n::Integer) where {T}
    return [x[firstindex(x):firstindex(x)+n-1], x[firstindex(x)+n:lastindex(x)]]
end
