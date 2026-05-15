using Test
using Resurgence

@testset "Resurgence.jl" begin
    include("test_utils.jl")
    include("test_poles.jl")
    include("test_shanks.jl")
    include("test_richardson.jl")
    include("test_wynn_eps.jl")
    include("test_cesaro_abel.jl")
    include("test_levin.jl")
    include("test_pade.jl")
    include("test_pade_cf.jl")
    include("test_pade_hermite.jl")
    include("test_borel.jl")
    include("test_borel_pade.jl")
    include("test_meijerg.jl")
    include("test_conformal.jl")
    include("test_truncation.jl")
    include("test_stokes.jl")
    include("test_api.jl")
    include("test_aqua.jl")
end
