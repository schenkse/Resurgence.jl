using Test
using Resurgence

@testset "Resurgence.jl" begin
    include("test_utils.jl")
    include("test_poles.jl")
    include("test_shanks.jl")
    include("test_pade.jl")
    include("test_borel.jl")
    include("test_borel_pade.jl")
    include("test_conformal.jl")
    include("test_truncation.jl")
    include("test_api.jl")
end
