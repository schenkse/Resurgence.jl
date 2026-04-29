using Test
using Resurgence

@testset "pade" begin
    @testset "geometric series 1/(1-z) ↔ Padé [0/1] is exact" begin
        # 1 + z + z² + … should reproduce 1/(1-z) exactly with [0/1].
        a = ones(Float64, 6)
        v = pade_value(a, 1, 0, 0.5)
        @test v ≈ 1 / (1 - 0.5)
    end

    @testset "exp(z) Padé [3/3] at z = 0.5" begin
        a = [1 / factorial(k) for k in 0:7]
        v = pade_value(a, 3, 3, 0.5)
        @test v ≈ exp(0.5) atol = 1e-6
    end

    @testset "argument validation" begin
        a = ones(5)
        @test_throws ArgumentError pade(a, 3, 3)   # needs length ≥ 7
        @test_throws ArgumentError pade(a, -1, 1)
    end

    @testset "BigFloat" begin
        a = [BigFloat(1) / factorial(big(k)) for k in 0:10]
        v = pade_value(a, 4, 4, BigFloat("0.5"))
        @test v isa BigFloat
        # Padé[4/4] of exp truncates the rational approximation at finite order;
        # error scales like the next neglected coefficient, ~1e-10.
        @test abs(v - exp(BigFloat("0.5"))) < BigFloat("1e-9")
    end
end
