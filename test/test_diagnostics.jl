using Test
using Resurgence

@testset "diagnostics" begin

    @testset "diagnose: convergent geometric" begin
        a = Float64[0.5^k for k in 0:20]
        d = diagnose(a)
        @test d.growth == :convergent
        @test d.alternating == false
        @test :shanks in d.recommended
        @test :pade in d.recommended
        # optimal_truncation for a monotone series puts N* at the end
        @test d.Nstar == length(a)
    end

    @testset "diagnose: factorial alternating (Stieltjes)" begin
        a = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:20]
        d = diagnose(a)
        @test d.growth == :factorial
        @test d.alternating == true
        @test first(d.recommended) == :borel_pade
        @test :conformal_borel_pade in d.recommended
        # stokes_action for an alternating-sign Borel-summable series is negative
        @test !ismissing(d.S)
        @test real(d.S) < 0
    end

    @testset "diagnose: factorial non-alternating" begin
        a = Float64[Float64(factorial(big(k))) for k in 0:20]
        d = diagnose(a)
        @test d.growth == :factorial
        @test d.alternating == false
        @test first(d.recommended) == :borel_pade_median
        @test :borel_pade_lateral in d.recommended
        @test !ismissing(d.S)
        @test real(d.S) > 0
    end

    @testset "diagnose: short input bypasses stokes_fit" begin
        d = diagnose([1.0, 2.0, 3.0])
        @test ismissing(d.S)
        @test ismissing(d.β)
        @test ismissing(d.A)
        @test d.Nstar ≥ 1
    end

    @testset "diagnose: empty input errors" begin
        @test_throws ArgumentError diagnose(Float64[])
    end

    @testset "diagnose: geometric divergent" begin
        # |a_k| = 2^k → ratios constant = 2, not factorial growth
        a = Float64[2.0^k for k in 0:20]
        d = diagnose(a)
        @test d.growth == :geometric
        @test :shanks in d.recommended
    end

end
