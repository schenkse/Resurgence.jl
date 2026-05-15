using Test
using Resurgence

@testset "levin" begin
    @testset "Leibniz / Gregory series for π" begin
        N = 30
        a = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        # u/t/v variants should all converge much faster than the plain sum.
        for variant in (:u, :t, :v)
            v = levin(a, 3; depth = 20, variant = variant)
            @test abs(v - π) < 1e-10
        end
        # Reference comparison: Levin u beats one Shanks step.
        s_shanks = shanks(a, 3; depth = 1)
        v_levin = levin(a, 3; depth = 20, variant = :u)
        @test abs(v_levin - π) < abs(s_shanks - π)
    end

    @testset "Stieltjes series — factorial divergence" begin
        # aₖ = (-1)ᵏ k!, Borel sum = e · E₁(1) = 0.5963473623231940743410784993...
        ref = 0.5963473623231940743410784993
        a = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:24]
        v = levin(a, 3; depth = 15, variant = :u)
        @test abs(v - ref) < 1e-8
    end

    @testset "argument validation" begin
        a = collect(1.0:6.0)
        @test_throws ArgumentError levin(a, 1; depth = 0)
        @test_throws ArgumentError levin(a, 0; depth = 3)
        @test_throws ArgumentError levin(a, 1; depth = 10, variant = :u)
        @test_throws ArgumentError levin(a, 1; depth = 6, variant = :v)  # needs a[8]
        @test_throws ArgumentError levin(a, 1; depth = 2, variant = :bogus)
    end

    @testset "BigFloat propagates" begin
        N = 25
        a = BigFloat[(-1)^k * factorial(big(k)) for k in 0:N-1]
        v = levin(a, 3; depth = 20, variant = :u)
        @test v isa BigFloat
        ref = BigFloat("0.5963473623231940743410784993")
        @test abs(v - ref) < 1e-10
    end

    @testset "Complex linearity" begin
        N = 30
        a_real = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        a_cplx = ComplexF64.(a_real) .* (1 + 0.1im)
        v_real = levin(a_real, 3; depth = 20, variant = :u)
        v_cplx = levin(a_cplx, 3; depth = 20, variant = :u)
        @test v_cplx isa ComplexF64
        @test v_cplx ≈ v_real * (1 + 0.1im)
    end

    @testset "degenerate weights return partial sum" begin
        # All-zero tail: every ω vanishes, no usable terms.
        a = zeros(Float64, 10); a[1] = 1.0
        v = levin(a, 2; depth = 5, variant = :t)
        @test isfinite(v)
        @test v ≈ 1.0
    end
end
