# Unified `resum` API

A thin dispatch layer on top of the per-method functions: each resummation method gets an `AbstractResummation` subtype, and `resum(method, a; ...)` forwards to the canonical implementation.
Use this when you want to pass a method around as a value (for parameter sweeps, configurable benchmarks, etc.); use the per-method functions for one-off calls.

## Public types and functions

```@autodocs
Modules = [Resurgence]
Pages   = ["src/api.jl"]
Private = false
Order   = [:type, :function]
```
