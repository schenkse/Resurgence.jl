# Padé approximants

Linear `[n/m]` Padé, continued-fraction (qd) Padé, and Hermite / quadratic
Padé for algebraic branch points. The pole-handling utilities used by all
Borel–Padé variants are also re-exposed here.

For when to reach for which variant, see the
[Methods guide](../methods_guide.md#padé-family).

## Linear Padé

```@autodocs
Modules = [Resurgence]
Pages   = ["src/pade.jl"]
Private = false
Order   = [:function, :type]
```

## Continued-fraction Padé (qd algorithm)

```@autodocs
Modules = [Resurgence]
Pages   = ["src/pade_cf.jl"]
Private = false
Order   = [:function, :type]
```

## Hermite / quadratic Padé

```@autodocs
Modules = [Resurgence]
Pages   = ["src/pade_hermite.jl"]
Private = false
Order   = [:function, :type]
```

## Pole handling

```@autodocs
Modules = [Resurgence]
Pages   = ["src/poles.jl"]
Private = false
Order   = [:function, :type]
```
