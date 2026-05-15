# Trans-series

Resurgence's natural object: `Σⱼ e^{−Sⱼ/g} g^{βⱼ} (Σₖ aⱼₖ gᵏ)`.
The `j = 0` sector is the perturbative answer that the Borel/Padé/Meijer-G methods produce on their own; the `j > 0` sectors carry instanton, multi-instanton, and other exponentially small contributions.
[`TransSeries`](@ref) packages them together as a list of [`Sector`](@ref) records and supports action-additive arithmetic; [`resum_transseries`](@ref) evaluates the whole object at a coupling `g` by dispatching each sector's perturbative tail through the existing [`resum`](@ref) surface and applying the prefactor `e^{−Sⱼ/g} g^{βⱼ}`.

```@autodocs
Modules = [Resurgence]
Pages   = ["src/transseries.jl"]
Private = false
Order   = [:type, :function]
```
