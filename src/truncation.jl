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
    Nstar = argmin(abs.(a))
    return Nstar, sum(@view a[1:Nstar]), abs(a[Nstar])
end

"""
    superasymptotic_remainder(a)

Estimate the remainder of an asymptotic series at its optimal-truncation
order: returns `|a[Nstar]|`, the magnitude of the smallest-magnitude term.
This is the textbook *superasymptotic* error estimate.
"""
function superasymptotic_remainder(a::AbstractVector{T}) where {T<:Number}
    isempty(a) && throw(ArgumentError("series must be non-empty"))
    return abs(a[argmin(abs.(a))])
end
