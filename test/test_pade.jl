using Test
using Resurgence

@testset "pade" begin
    @testset "geometric series 1/(1-z) ↔ Padé [0/1] is exact" begin
        # 1 + z + z² + … should reproduce 1/(1-z) exactly with [0/1].
        a = ones(Float64, 6)
        v = pade_value(a, 0, 1, 0.5)
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
        @test_throws ArgumentError pade(a, -1, 1)  # n < 0
        @test_throws ArgumentError pade(a, 1, -1)  # m < 0
    end

    @testset "BigFloat" begin
        a = [BigFloat(1) / factorial(big(k)) for k in 0:10]
        v = pade_value(a, 4, 4, BigFloat("0.5"))
        @test v isa BigFloat
        # Padé[4/4] of exp truncates the rational approximation at finite order;
        # error scales like the next neglected coefficient, ~1e-10.
        @test abs(v - exp(BigFloat("0.5"))) < BigFloat("1e-9")
    end

    @testset "Complex eltype" begin
        # Geometric series 1 + z + z² + … with [0/1] is exactly 1/(1-z),
        # which is well-defined for any z ≠ 1, real or complex.
        a = ones(ComplexF64, 6)
        x = 0.5 + 0.1im
        v = pade_value(a, 0, 1, x)
        @test v isa ComplexF64
        @test v ≈ 1 / (1 - x)
    end

    @testset "BigFloat rank-deficient fallback" begin
        # Borel transform of the Stieltjes series is exactly 1/(1+t), a degree-1
        # rational. Padé [n/m] with n,m ≥ 1 is therefore rank-deficient and the
        # LU fast path throws SingularException; the catch branch routes through
        # `pinv`, which only works on BigFloat once GenericLinearAlgebra brings
        # in a generic SVD. This test guards that the fallback path stays alive.
        b = BigFloat[(-1)^k for k in 0:6]    # Borel coefficients of 1/(1+t)
        # 1/(1+t) at t = 1/2 is 2/3
        v = pade_value(b, 2, 2, BigFloat("0.5"))
        @test v isa BigFloat
        @test isfinite(v)
        @test abs(v - BigFloat(2) / 3) < BigFloat("1e-30")
    end
end
