# ModelingToolkit DAE initialization tests for DASKR.jl
# Following the pattern from Sundials.jl:
# https://github.com/SciML/Sundials.jl/blob/master/test/common_interface/initialization.jl

using ModelingToolkit, SciMLBase, DASKR, Test, DiffEqBase
using SymbolicIndexingInterface
using ModelingToolkit: t_nounits as t, D_nounits as D

@testset "MTK DAE Initialization" begin
    # Define a DAE system with:
    # - Differential equation: D(x) ~ p * y + q * t
    # - Algebraic constraint: x^3 + y^3 ~ 5
    # - Initialization equation: p^2 + q^2 ~ 3
    # - Parameters p and q are "missing" and need to be determined from initialization
    @variables x(t) [guess = 1.0] y(t) [guess = 1.0]
    @parameters p = missing [guess = 1.0] q = missing [guess = 1.0]
    @mtkcompile sys = System(
        [D(x) ~ p * y + q * t, x^3 + y^3 ~ 5],
        t;
        initialization_eqs = [p^2 + q^2 ~ 3]
    )

    # Expected values after initialization:
    # x = 1.0 (given)
    # y = cbrt(4) ≈ 1.587 (from x^3 + y^3 = 5 with x = 1)
    # p = 1.0 (given)
    # q = sqrt(2) ≈ 1.414 (from p^2 + q^2 = 3 with p = 1)
    # D(x) = p * y + q * t = 1 * cbrt(4) + sqrt(2) * 0 = cbrt(4) at t = 0
    # D(y) from differentiating constraint: 3x^2*D(x) + 3y^2*D(y) = 0
    #       => D(y) = -x^2/y^2 * D(x) = -1/cbrt(4)^2 * cbrt(4) = -1/cbrt(4)

    @testset "DAEProblem{$iip}" for iip in [true, false]
        # Create DAEProblem with initial guesses
        # Give x = 1.0 and p = 1.0, let initialization solve for y and q
        prob = DAEProblem(
            sys,
            [D(x) => cbrt(4), D(y) => -1 / cbrt(4), p => 1.0],
            (0.0, 0.4)
        )

        @testset "OverrideInit (MTK initialization)" begin
            # Test with init() to check initialization works
            integ = init(prob, daskr())
            @test integ.initializealg isa DASKR.DefaultInit
            @test integ[x] ≈ 1.0
            @test integ[y] ≈ cbrt(4)
            @test integ.ps[p] ≈ 1.0
            @test integ.ps[q] ≈ sqrt(2)

            # Test full solve
            sol = solve(prob, daskr())
            @test SciMLBase.successful_retcode(sol)
            @test sol[x, 1] ≈ 1.0
            @test sol[y, 1] ≈ cbrt(4)
            @test sol.ps[p] ≈ 1.0
            @test sol.ps[q] ≈ sqrt(2)
        end

        @testset "CheckInit" begin
            # CheckInit should fail with incomplete initial values
            @test_throws Any init(prob, daskr(); initializealg = SciMLBase.CheckInit())

            # Create a problem with correct initial values - all values specified
            # D(x) = p*y = 1*cbrt(4) = cbrt(4)
            # D(y) = -x²/y²*D(x) = -1/cbrt(4)²*cbrt(4) = -1/cbrt(4)
            prob_correct = DAEProblem(
                sys,
                [
                    D(x) => cbrt(4),
                    D(y) => -1 / cbrt(4),
                    p => 1.0,
                    x => 1.0,
                    y => cbrt(4),
                    q => sqrt(2)
                ],
                (0.0, 0.4)
            )

            # Need to convert to IIP/OOP after creation to get proper numeric arrays
            if iip
                prob_correct_typed = DAEProblem{true}(
                    prob_correct.f,
                    prob_correct.du0,
                    prob_correct.u0,
                    prob_correct.tspan,
                    prob_correct.p
                )
            else
                prob_correct_typed = DAEProblem{false}(
                    prob_correct.f,
                    prob_correct.du0,
                    prob_correct.u0,
                    prob_correct.tspan,
                    prob_correct.p
                )
            end

            # This should work since all values are correct
            @test_nowarn init(prob_correct_typed, daskr(); initializealg = SciMLBase.CheckInit())
        end
    end
end
