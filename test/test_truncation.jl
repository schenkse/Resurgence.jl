using Test
using Resurgence

@testset "truncation" begin
    @testset "Stieltjes series at small z has Nstar ≈ 1/z" begin
        z = 0.1
        N = 30
        # a_k = (-1)^k k! z^k, |a_k| smallest near k ≈ 1/z.
        a = Float64[(-1.0)^k * factorial(big(k)) * z^k for k in 0:N-1]
        # truncate to keep things in Float64 range
        a = Float64.(a)
        Nstar, partial, ε = optimal_truncation(a)
        # The minimum-magnitude term sits near index k = round(1/z) = 10.
        # We use 1-based indexing, so Nstar should be around 11.
        @test 8 ≤ Nstar ≤ 12
        @test partial ≈ sum(a[1:Nstar])
        @test ε == abs(a[Nstar])
    end

    @testset "superasymptotic_remainder agrees with optimal_truncation" begin
        a = Float64[(-1.0)^k * factorial(big(k)) * 0.1^k for k in 0:20]
        a = Float64.(a)
        _, _, ε = optimal_truncation(a)
        @test superasymptotic_remainder(a) ≈ ε
    end

    @testset "empty series errors cleanly" begin
        @test_throws ArgumentError optimal_truncation(Float64[])
        @test_throws ArgumentError superasymptotic_remainder(Float64[])
    end

    @testset "monotone series ⇒ trivial answer" begin
        a = [3.0, 2.0, 1.0]
        Nstar, partial, ε = optimal_truncation(a)
        @test Nstar == 3
        @test partial == sum(a)
        @test ε == 1.0
    end
end
