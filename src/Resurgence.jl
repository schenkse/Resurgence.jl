module Resurgence

using LinearAlgebra
using GenericLinearAlgebra
using HypergeometricFunctions
using Polynomials
using PolynomialRoots
using QuadGK
using SpecialFunctions

include("utils.jl")
include("poles.jl")
include("series_acceleration.jl")
include("levin.jl")
include("pade.jl")
include("pade_cf.jl")
include("pade_hermite.jl")
include("borel.jl")
include("conformal.jl")
include("borel_pade.jl")
include("meijerg.jl")
include("truncation.jl")
include("stokes.jl")
include("api.jl")
include("transseries.jl")

# functional API
export shanks, richardson, wynn_eps, cesaro, abel
export theta_brezinski, rho_brezinski
export aitken_steffensen
export levin, weniger, sidi_s
export pade, pade_value
export pade_cf, pade_cf_value
export hermite_pade, hermite_pade_value
export borel_transform, borel_leroy_transform, borel_ratios
export conformal_map, inverse_conformal_map, conformal_reseries
export conformal_map_pair, inverse_conformal_map_pair, conformal_reseries_pair
export borel_pade, borel_leroy_pade, conformal_borel_pade, conformal_borel_pade_pair
export borel_pade_lateral, borel_pade_median, borel_pade_discontinuity
export borel_leroy_pade_lateral, borel_leroy_pade_median, borel_leroy_pade_discontinuity
export borel_leroy_pade_odm
export borel_meijerg
export optimal_truncation, superasymptotic_remainder, terminant, hyperasymptotic
export stokes_action, stokes_exponent, stokes_constant, stokes_fit
export chop!, sparsify!

# unified API
export AbstractResummation, Shanks, Richardson, WynnEps, Pade, PadeCF, HermitePade,
    BorelPade, BorelPadeLateral, BorelPadeMedian,
    BorelLeRoyPade, BorelLeRoyPadeODM, ConformalBorelPade, ConformalBorelPadePair, MeijerG,
    Cesaro, Abel, Levin, Weniger, SidiS, BrezinskiTheta, BrezinskiRho,
    Hyperasymptotic
export resum

# trans-series
export Sector, TransSeries, transseries_exp, resum_transseries

end # module
