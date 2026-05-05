using Test
using Resurgence

@testset "richardson" begin
    @testset "Σ 1/k² → π²/6" begin
        # Partial sums Aₙ = π²/6 − 1/n + 1/(2n²) − …, the textbook 1/n
        # asymptotic expansion that Richardson is built for.
        N = 30
        a = Float64[1 / k^2 for k in 1:N]
        target = π^2 / 6
        partial = sum(a)
        r1 = richardson(a, 10)
        r3 = richardson(a, 10; depth = 3)
        r5 = richardson(a, 10; depth = 5)
        @test abs(r1 - target) < abs(partial - target)
        @test abs(r3 - target) < abs(r1 - target)
        @test abs(r5 - target) < abs(r3 - target)
        @test abs(r5 - target) < 1e-8
    end

    @testset "argument validation" begin
        a = collect(1.0:5.0)
        @test_throws ArgumentError richardson(a, 0)             # n must be ≥ 1
        @test_throws ArgumentError richardson(a, 5)             # need length ≥ n+1
        @test_throws ArgumentError richardson(a, 1; depth = 0)  # depth must be ≥ 1
        @test_throws ArgumentError richardson(a, 3; depth = 3)  # need length ≥ n+depth
    end

    @testset "depth = 1 closed form" begin
        # T⁽¹⁾ₙ = (n+1) Aₙ₊₁ − n Aₙ
        a = Float64[1 / k^2 for k in 1:10]
        n = 4
        An  = sum(a[1:n])
        Anp = sum(a[1:n+1])
        @test richardson(a, n) ≈ (n + 1) * Anp - n * An
    end

    @testset "BigFloat propagates" begin
        a = BigFloat[1 / BigFloat(k)^2 for k in 1:30]
        r = richardson(a, 10; depth = 5)
        @test r isa BigFloat
        @test abs(r - BigFloat(π)^2 / 6) < 1e-8
    end

    @testset "Complex eltype" begin
        # Richardson is linear in `a`; scaling by a complex constant scales
        # the result identically. Smoke-test type-genericity, not numerics.
        N = 30
        a_real = Float64[1 / k^2 for k in 1:N]
        a_cplx = ComplexF64.(a_real) .* (1 + 0.1im)
        r_real = richardson(a_real, 10; depth = 5)
        r_cplx = richardson(a_cplx, 10; depth = 5)
        @test r_cplx isa ComplexF64
        @test r_cplx ≈ r_real * (1 + 0.1im)
    end
end
