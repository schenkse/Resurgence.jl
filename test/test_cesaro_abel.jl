using Test
using Resurgence

@testset "cesaro" begin
    @testset "Grandi 1 − 1 + 1 − 1 + … averages to 1/2" begin
        # Canonical Cesàro example: Aₙ alternates 1, 0, 1, 0, … so C⁽¹⁾_{2k} = 1/2.
        a = Float64[(-1.0)^k for k in 0:99]
        @test cesaro(a, 100) ≈ 0.5
    end

    @testset "constant series a[k] = 1 gives C⁽¹⁾ₙ = (n+1)/2" begin
        a = ones(Float64, 10)
        # Partial sums are 1..10; their mean is 5.5.
        @test cesaro(a, 10) ≈ 5.5
    end

    @testset "Leibniz / Gregory series for π/4" begin
        # (C,1) averaging smooths out the Leibniz oscillation; depth=2 also
        # converges but is biased at finite N (iterated arithmetic means
        # weight early partial sums heavily), so we don't require it to beat
        # depth=1 here — only that both stay within a reasonable band of π.
        N = 100
        a = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        c1 = cesaro(a, N)
        c2 = cesaro(a, N; depth = 2)
        @test abs(c1 - π) < abs(sum(a) - π)
        @test abs(c2 - π) < 0.1
    end

    @testset "argument validation" begin
        a = collect(1.0:5.0)
        @test_throws ArgumentError cesaro(a, 0)
        @test_throws ArgumentError cesaro(a, 6)
        @test_throws ArgumentError cesaro(a, 3; depth = 0)
    end

    @testset "BigFloat propagates" begin
        a = BigFloat[(-1)^k * BigFloat(1) / (2k + 1) for k in 0:50]
        s = cesaro(4 * a, 51; depth = 2)
        @test s isa BigFloat
    end

    @testset "Complex eltype is linear" begin
        N = 30
        a_real = Float64[(-1.0)^k for k in 0:N-1]
        a_cplx = ComplexF64.(a_real) .* (1 + 0.1im)
        @test cesaro(a_cplx, N; depth = 2) ≈ cesaro(a_real, N; depth = 2) * (1 + 0.1im)
    end
end

@testset "abel" begin
    @testset "geometric series, closed form" begin
        # a[k+1] = 1, k = 0..N-1, Σ xᵏ = (1 − xᴺ)/(1 − x).
        N = 8
        a = ones(Float64, N)
        for x in (0.25, 0.5, 0.75, 0.99)
            @test abel(a; x = x) ≈ (1 - x^N) / (1 - x)
        end
    end

    @testset "default x = 1 reduces to the plain partial sum" begin
        a = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:9]
        @test abel(a) ≈ sum(a)
    end

    @testset "argument validation / empty input" begin
        @test abel(Float64[]) === 0.0
        @test abel(Float64[]; x = 0.5) === 0.0
    end

    @testset "BigFloat propagates" begin
        a = BigFloat[BigFloat(1) / (k + 1) for k in 0:20]
        s = abel(a; x = BigFloat("0.5"))
        @test s isa BigFloat
    end

    @testset "promotion Float64 + BigFloat → BigFloat" begin
        a = ones(Float64, 5)
        s = abel(a; x = BigFloat("0.5"))
        @test s isa BigFloat
        @test s ≈ BigFloat(1 - 0.5^5) / BigFloat("0.5")
    end

    @testset "Complex eltype" begin
        a = ComplexF64[1.0, 0.5im, -0.25, -0.125im]
        x = 0.5
        @test abel(a; x = x) ≈ 1 + 0.5im * x + (-0.25) * x^2 + (-0.125im) * x^3
    end
end
