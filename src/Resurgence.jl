module Resurgence

using LinearAlgebra
using GenericLinearAlgebra  # extends svd/pinv to BigFloat (rank-deficient pade fallback)
using Polynomials
using PolynomialRoots
using QuadGK
using SpecialFunctions

include("utils.jl")
include("poles.jl")
include("series_acceleration.jl")
include("pade.jl")
include("borel.jl")
include("conformal.jl")
include("borel_pade.jl")
include("truncation.jl")
include("api.jl")

# functional API
export shanks, richardson
export pade, pade_value
export borel_transform, borel_leroy_transform, borel_ratios
export conformal_map, inverse_conformal_map, conformal_reseries
export borel_pade, borel_leroy_pade, conformal_borel_pade
export optimal_truncation, superasymptotic_remainder
export chop!, sparsify!

# unified API
export AbstractResummation, Shanks, Richardson, Pade, BorelPade, BorelLeRoyPade, ConformalBorelPade
export resum

end # module
