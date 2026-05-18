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

    @testset "compare: two methods with reference" begin
        a = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:24]
        ref = 0.5963473623231940
        rows = compare([BorelPade(10, 10), Pade(10, 10)], a; reference = ref)
        @test length(rows) == 2
        @test rows[1].method == "BorelPade"
        @test rows[2].method == "Pade"
        for r in rows
            @test r.result !== missing
            @test r.residual < 1.0          # both at least in the ballpark
            @test r.error === nothing
        end
        # Borel-Padé should be dramatically better than plain Padé on Stieltjes
        @test rows[1].residual < rows[2].residual
    end

    @testset "compare: failing tag is caught, others survive" begin
        a = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:24]
        # Pade(50, 50) needs length(a) ≥ 101; it will throw ArgumentError.
        rows = compare([BorelPade(10, 10), Pade(50, 50)], a; reference = 0.596)
        @test rows[1].result !== missing
        @test rows[1].error === nothing
        @test rows[2].result === missing
        @test rows[2].error isa String
        @test !isempty(rows[2].error)
        @test rows[2].residual === missing
    end

    @testset "compare: no reference → residual missing" begin
        a = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:24]
        rows = compare([BorelPade(10, 10)], a)
        @test rows[1].residual === missing
        @test rows[1].result !== missing
    end

    @testset "compare: accepts tuples" begin
        a = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:24]
        rows = compare((BorelPade(10, 10), Pade(10, 10)), a)
        @test length(rows) == 2
    end

end
