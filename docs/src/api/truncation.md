# Optimal truncation

Optimal-truncation / superasymptotic estimates: stop summing at the smallest term, take the partial sum, and bound the remainder by the size of the smallest term.
Free in the sense that no extra computation is needed beyond the partial sums you already have.

## Public functions

```@autodocs
Modules = [Resurgence]
Pages   = ["src/truncation.jl"]
Private = false
Order   = [:function, :type]
```
