using Test
using Resurgence
using SpecialFunctions: gamma

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

@testset "terminant" begin
    @testset "decays to 0 as σ → ∞" begin
        @test abs(terminant(0.5, 50.0)) < 1e-15
        @test abs(terminant(2.5, 100.0)) < 1e-30
    end

    @testset "monotone-ish in σ at fixed p" begin
        # Magnitude decays as σ grows (the Γ(1-p, σ) factor → 0).
        @test abs(terminant(2.5, 1.0)) > abs(terminant(2.5, 5.0)) > abs(terminant(2.5, 20.0))
    end

    @testset "matches the closed-form expression" begin
        p, σ = 2.5, 1.0
        expected = exp(im * π * p) * gamma(p) * gamma(1 - p, σ) / (2 * π * im)
        @test terminant(p, σ) ≈ expected
    end

    @testset "BigFloat propagates" begin
        t = terminant(big"2.5", big"1.0")
        @test t isa Complex{BigFloat}
    end

    @testset "σ → ∞ asymptotic of Γ(1−p, σ) through two correction terms" begin
        # Independent of the function's own formula: cross-check the large-σ
        # asymptotic (NIST DLMF 8.11.2),
        #     Γ(1−p, σ) ~ σ^(−p) e^(−σ) · [1 − p/σ + p(p+1)/σ² − …],
        # which feeds into |T_p(σ)| ~ |Γ(p)| · σ^(−p) e^(−σ) · (series) / (2π).
        # Truncating after two correction terms leaves an O(p³/σ³) tail.
        p = 5.5
        for σ in (40.0, 80.0, 160.0)
            series = 1 - p / σ + p * (p + 1) / σ^2
            ref = abs(gamma(p)) * σ^(-p) * exp(-σ) * series / (2π)
            @test isapprox(abs(terminant(p, σ)) / ref, 1.0; rtol = 5 * p^2 / σ^2)
        end
    end
end

@testset "hyperasymptotic" begin
    stieltjes(::Type{T}, N) where {T} = T[(-1)^k * factorial(big(k)) for k in 0:N-1]
    factorials(::Type{T}, N) where {T} = T[factorial(big(k)) for k in 0:N-1]

    @testset "level=0 reduces to optimal-truncation partial sum at x" begin
        a = Float64.(stieltjes(BigFloat, 25))
        x = 0.1
        v = hyperasymptotic(a; x = x, level = 0)
        # Compare against the direct sum at the optimal truncation index.
        mags = [abs(a[k] * x^(k-1)) for k in eachindex(a)]
        Nstar = argmin(mags)
        direct = sum(a[k] * x^(k-1) for k in 1:Nstar)
        @test v ≈ direct
    end

    @testset "Borel-summable (action with negative real part) errors" begin
        # Stieltjes at x>0 has Stokes singularity at t = -1 (S<0); the level-1
        # formula does not apply.
        a = Float64.(stieltjes(BigFloat, 25))
        @test_throws ArgumentError hyperasymptotic(a; x = 0.1, level = 1,
                                                   action = -1.0, β = 1.0, A = 1.0)
    end

    @testset "level=1 on a_k = k! captures the lateral-sum imaginary part" begin
        # a_k = k! has Stokes singularity at t = +1: the single-instanton
        # terminant correction here adds the imaginary discontinuity.
        a = Float64.(factorials(BigFloat, 25))
        x = 0.1
        v1 = hyperasymptotic(a; x = x, level = 1, action = 1.0, β = 1.0, A = 1.0)
        # Lateral sum's imaginary part is the reference; both should be of
        # order π·exp(-1/x)/x ≈ 1.4e-3 here, with the same sign.
        Lp = borel_pade_lateral(a; n = 10, m = 10, x = x, side = +1)
        @test imag(v1) * imag(Lp) > 0           # same sign
        # within an order of magnitude
        @test 0.1 < abs(imag(v1) / imag(Lp)) < 10
    end

    @testset "level=1 imag part matches the closed-form −π·A·x^(−β)·exp(−S/x)" begin
        # The level-1 correction is exactly -i·π·A·x^(-β)·exp(-S/x). With
        # real partial sum and real (A, β, S), imag(v) equals this closed
        # form to machine precision.
        a = Float64.(factorials(BigFloat, 25))
        x = 0.1
        v = hyperasymptotic(a; x = x, level = 1, action = 1.0, β = 1.0, A = 1.0)
        expected_imag = -π * 1.0 * x^(-1.0) * exp(-1.0 / x)
        @test imag(v) ≈ expected_imag rtol = 1e-12
    end

    @testset "level=1 auto-extracts (action, β, A) via stokes_fit" begin
        # a_k = k! has stokes_fit S ≈ +1, β ≈ 1, A ≈ 1 — auto-extract should
        # succeed and the level-1 sum should track the lateral imaginary part.
        a = Float64.(factorials(BigFloat, 25))
        v = hyperasymptotic(a; x = 0.1, level = 1)
        Lp = borel_pade_lateral(a; n = 10, m = 10, x = 0.1, side = +1)
        @test imag(v) * imag(Lp) > 0
        @test 0.1 < abs(imag(v) / imag(Lp)) < 10
    end

    @testset "level ≥ 2 throws" begin
        a = Float64.(stieltjes(BigFloat, 25))
        @test_throws ArgumentError hyperasymptotic(a; x = 0.1, level = 2)
        @test_throws ArgumentError hyperasymptotic(a; x = 0.1, level = -1)
    end

    @testset "partial provision of (action, β, A) errors" begin
        a = Float64.(stieltjes(BigFloat, 25))
        @test_throws ArgumentError hyperasymptotic(a; x = 0.1, level = 1, action = 1.0)
        @test_throws ArgumentError hyperasymptotic(a; x = 0.1, level = 1,
                                                   action = 1.0, β = 1.0)
    end

    @testset "BigFloat smoke" begin
        a = factorials(BigFloat, 25)
        v = hyperasymptotic(a; x = big"0.1", level = 1,
                            action = big"1", β = big"1", A = big"1")
        @test v isa Complex{BigFloat}
    end

    @testset "empty input errors" begin
        @test_throws ArgumentError hyperasymptotic(Float64[]; x = 0.1)
    end

    @testset "Hyperasymptotic API tag dispatches" begin
        a = Float64.(stieltjes(BigFloat, 25))
        v_func = hyperasymptotic(a; x = 0.1, level = 0)
        v_api = resum(Hyperasymptotic(; x = 0.1, level = 0), a)
        @test v_func ≈ v_api
    end
end
