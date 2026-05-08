# References

Canonical references for the techniques implemented in Resurgence.jl, plus
the broader resurgence and large-order literature that motivates them.

## Background reading

- Bender & Orszag, *Advanced Mathematical Methods for Scientists and
  Engineers* [BenderOrszag1978](@cite). The standard textbook on
  asymptotic methods, optimal truncation, and Borel summation as applied
  to physics. Chapter on Padé approximants is the gentlest entry point.
- Costin, *Asymptotics and Borel Summability* [Costin2008](@cite).
  Modern, mathematically rigorous treatment of Borel summation, lateral
  sums, and resurgence; the standard reference for the trans-series
  framework.

## Resummation techniques

- Shanks transformation: [Shanks1955](@cite).
- Wynn ε-algorithm (planned, A1 on the [roadmap](roadmap.md)):
  [Wynn1956](@cite).
- Levin transforms: [Levin1973](@cite).
- Weniger δ-transformation and review of nonlinear sequence transformations:
  [Weniger1989](@cite).

## Borel summation and resurgence

- Sokal, "An improvement of Watson's theorem on Borel summability"
  [Sokal1980](@cite). The standard sufficient conditions for a divergent
  series to be uniquely Borel summable.
- Stokes, "On the discontinuity of arbitrary constants which appear in
  divergent developments" [Stokes1864](@cite). The original Stokes
  phenomenon paper.
- Écalle, *Les fonctions résurgentes* [Ecalle1981](@cite). The foundation
  of alien calculus and the mathematical theory of resurgence.

## Large-order behaviour and instantons

- Bender & Wu, "Anharmonic Oscillator" [BenderWu1969](@cite) and
  "Anharmonic Oscillator. II. A Study of Perturbation Theory in Large
  Order" [BenderWu1973](@cite). The starting point of large-order
  perturbation theory in physics; the canonical example of factorial
  divergence with extractable instanton action.
- Suslov, "Divergent perturbation series" [Suslov2005](@cite). Review of
  large-order behaviour in field theory, with explicit recipes for
  extracting `S`, `β`, `A` from coefficients — directly relevant to
  [`stokes_fit`](@ref).
- Zinn-Justin & Jentschura, "Multi-instantons and exact results I"
  [ZinnJustinJentschura2004](@cite). Modern resurgence-style treatment
  with multi-instanton corrections.

## Hyperasymptotics

- Berry & Howls, "Hyperasymptotics for integrals with saddles"
  [BerryHowls1990](@cite). The original hyperasymptotic terminant
  expansion; the basis for the planned `hyperasymptotic` method (B3 on
  the [roadmap](roadmap.md)).

## Bibliography

```@bibliography
```
