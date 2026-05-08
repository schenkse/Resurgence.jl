# Meijer-G resummation

Pure-Julia Meijer-G resummation: the implementation collapses the G-function onto a single `HypergeometricFunctions.pFq` call via the Slater identity (the doubled `b_h = 1` from Borel–Laplace cancels two of the numerator parameters), so no Meijer-G is ever evaluated directly and there is no Python dependency.

See the [Stieltjes tutorial](../tutorials/stieltjes.md#meijer-g) for a note on degenerate drivers (the unshifted Stieltjes series is one such trap; the shifted `aₖ = (−1)ᵏ (k+1)!` is the canonical non-degenerate driver).

## Public functions

```@autodocs
Modules = [Resurgence]
Pages   = ["src/meijerg.jl"]
Private = false
Order   = [:function, :type]
```
