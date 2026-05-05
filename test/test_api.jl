using Test
using Resurgence

@testset "unified resum API" begin
    a_stieltjes = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:24]
    direct = borel_pade(a_stieltjes; n = 10, m = 10, x = 1)

    @testset "BorelPade dispatches to borel_pade" begin
        v = resum(BorelPade(10, 10), a_stieltjes)
        @test v ≈ direct
    end

    @testset "BorelLeRoyPade dispatches to borel_leroy_pade" begin
        # b ≠ 0 keeps the Padé linear system non-singular on this exact-rational
        # Borel transform.
        v = resum(BorelLeRoyPade(10, 10; b = -0.25), a_stieltjes)
        v_direct = borel_leroy_pade(a_stieltjes; n = 10, m = 10, b = -0.25, x = 1)
        @test v ≈ v_direct
    end

    @testset "ConformalBorelPade dispatches" begin
        v = resum(ConformalBorelPade(10, 10; sing = 1.0), a_stieltjes)
        v_direct = conformal_borel_pade(a_stieltjes; n = 10, m = 10, x = 1, sing = 1.0)
        @test v ≈ v_direct
    end

    @testset "Pade dispatches to pade_value" begin
        a = ones(Float64, 6)
        v = resum(Pade(1, 0; x = 0.5), a)
        @test v ≈ 1 / (1 - 0.5)
    end

    @testset "Shanks dispatches with depth" begin
        N = 30
        a = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        @test resum(Shanks(5; depth = 3), a) ≈ shanks(a, 5; depth = 3)
    end

    @testset "Richardson dispatches with depth" begin
        a = Float64[1 / k^2 for k in 1:30]
        @test resum(Richardson(10; depth = 3), a) ≈ richardson(a, 10; depth = 3)
    end

    @testset "kwargs forwarded to underlying call" begin
        # rtol forwarded; loose tol shouldn't make this fail catastrophically.
        v = resum(BorelPade(10, 10; rtol = 1e-3), a_stieltjes)
        @test isapprox(v, 0.5963473623; atol = 1e-3)
    end
end
