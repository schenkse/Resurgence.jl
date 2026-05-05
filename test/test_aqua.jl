using Test
using Aqua
using Resurgence

@testset "Aqua" begin
    Aqua.test_all(Resurgence; ambiguities = false)
    Aqua.test_ambiguities(Resurgence)
end
