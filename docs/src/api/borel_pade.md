# Borel–Padé resummation

The full Borel–Padé family: the standard `borel_pade`, the Le Roy variant `borel_leroy_pade`, the conformal-map variant `conformal_borel_pade`, and the lateral / median / discontinuity sums used for non-Borel-summable inputs.
See the [lateral and median sums tutorial](../tutorials/lateral_sums.md) for a worked example of the latter.

`quadgk` keyword arguments (`rtol`, `atol`, `order`) flow through to the Laplace integration step in every method here.

## Public functions

```@autodocs
Modules = [Resurgence]
Pages   = ["src/borel_pade.jl"]
Private = false
Order   = [:function, :type]
```
