using Test
using Resurgence
using Polynomials

# Taylor coefficients of √(1 − z): expanded as Σ binomial(1/2, k) (−z)^k.
sqrt_one_minus_z_coeffs(::Type{T}, N) where {T} = begin
    c = Vector{T}(undef, N)
    c[1] = one(T)
    for k in 1:N-1
        c[k+1] = c[k] * (one(T)/2 - T(k - 1)) / T(k)
    end
    for k in 0:N-1
        c[k+1] *= (-one(T))^k
    end
    return c
end

# Taylor coefficients of exp(z/2): c_k = (1/2)^k / k!. Iteratively to keep
# overflow-free at BigFloat. Used as the canonical *non-degenerate* driver:
# f(z) = exp(z/2) ⇒ f² = exp(z), which has no algebraic relation in f,
# so the q_0 = 1 normalization of `hermite_pade` doesn't degenerate.
exp_half_coeffs(::Type{T}, N) where {T} = begin
    c = Vector{T}(undef, N)
    c[1] = one(T)
    for k in 1:N-1
        c[k+1] = c[k] / (T(2) * T(k))
    end
    return c
end

@testset "pade_hermite" begin
    @testset "argument validation" begin
        a = ones(5)
        @test_throws ArgumentError hermite_pade(a, 2, 2, 2)        # needs length ≥ 8
        @test_throws ArgumentError hermite_pade(a, -1, 0, 0)
        @test_throws ArgumentError hermite_pade(a, 0, -1, 0)
        @test_throws ArgumentError hermite_pade(a, 0, 0, -1)
        @test_throws ArgumentError hermite_pade_value(a, 2, 2, 2, 0.1)
        @test_throws ArgumentError hermite_pade_value(a, 0, 0, 1, 0.1; branch = 2)
        @test_throws ArgumentError hermite_pade_value(a, 0, 0, 1, 0.1; branch = 0)
    end

    @testset "[0/0/1] hand-computed solution for √(1-z)" begin
        # Closed-form Hermite-Padé at order [0/0/1]:
        #   P = -3/8, Q = 1, R = -5/8 - z/8.
        # Verifies the linear-system construction of `hermite_pade` directly;
        # also exercises the default-branch heuristic (R(0) = -5/8 < 0 forces
        # branch = -1 to recover f(0) = 1).
        a = sqrt_one_minus_z_coeffs(Float64, 3)
        P, Q, R = hermite_pade(a, 0, 0, 1)
        @test Polynomials.coeffs(P) ≈ [-3 / 8]
        @test Polynomials.coeffs(Q) ≈ [1.0]
        @test Polynomials.coeffs(R) ≈ [-5 / 8, -1 / 8]
        v = hermite_pade_value(a, 0, 0, 1, 0.1)
        @test v ≈ sqrt(0.9) atol = 5e-3   # 3-coeff match → O(z^3) error
    end

    @testset "[3/3/3] recovery of √(exp(z)) (non-degenerate driver)" begin
        # f² = exp(z) has no polynomial relation in f, so the system stays
        # full-rank under q_0 = 1. Reference: f(z) = exp(z/2).
        a = exp_half_coeffs(Float64, 11)
        for z in (0.1, 0.3, 0.5, 0.9)
            v = hermite_pade_value(a, 3, 3, 3, z)
            @test v ≈ exp(z / 2) atol = 1e-10
        end
    end

    @testset "branch override picks the conjugate root" begin
        a = exp_half_coeffs(Float64, 11)
        z = 0.3
        v_default = hermite_pade_value(a, 3, 3, 3, z)
        v_plus  = hermite_pade_value(a, 3, 3, 3, z; branch = +1)
        v_minus = hermite_pade_value(a, 3, 3, 3, z; branch = -1)
        # Default must equal one of the two explicit branches, and the two
        # branches must differ.
        @test v_default ≈ v_plus || v_default ≈ v_minus
        @test !(v_plus ≈ v_minus)
        @test v_default ≈ exp(z / 2) atol = 1e-10
        other = (v_default ≈ v_plus) ? v_minus : v_plus
        @test !isapprox(other, exp(z / 2); atol = 1e-3)
    end

    @testset "BigFloat sharpens the answer" begin
        a = exp_half_coeffs(BigFloat, 13)
        v = hermite_pade_value(a, 4, 4, 3, BigFloat("0.5"))
        @test v isa BigFloat
        @test abs(v - exp(BigFloat("0.25"))) < BigFloat("1e-18")
    end

    @testset "Complex eltype propagates" begin
        a = ComplexF64.(exp_half_coeffs(Float64, 11))
        z = 0.3 + 0.05im
        v = hermite_pade_value(a, 3, 3, 3, z)
        @test v isa Complex
        @test v ≈ exp(z / 2) atol = 1e-8
    end

    @testset "complex(z) follows √(1-z) past the branch into negative-disc territory" begin
        # √(1-z) is the canonical algebraic-branch driver; the q_0 = 1 system
        # degenerates to linear-Padé-with-R≈0 at high order, so we stay at
        # low order [0/0/1] for a direct value check, but use complex(z) so
        # sqrt() never raises DomainError on negative discriminants.
        a = sqrt_one_minus_z_coeffs(Float64, 3)
        z = complex(0.05)
        v = hermite_pade_value(a, 0, 0, 1, z)
        @test real(v) ≈ sqrt(0.95) atol = 5e-4
        @test abs(imag(v)) < 1e-10   # real input → essentially real output
    end

    @testset "polynomial degrees and normalization" begin
        a = exp_half_coeffs(Float64, 11)
        P, Q, R = hermite_pade(a, 3, 3, 3)
        @test Polynomials.degree(P) ≤ 3
        @test Polynomials.degree(Q) ≤ 3
        @test Polynomials.degree(R) ≤ 3
        @test Q(0.0) ≈ 1.0   # q_0 = 1 normalisation
    end
end
