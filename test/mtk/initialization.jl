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

        @testset "OverrideInit (MTK initialization) via DefaultInit" begin
            # Test full solve - DefaultInit should route to OverrideInit for MTK problems
            # DASKR uses solve() directly (no integrator interface like Sundials)
            sol = solve(prob, daskr())
            @test SciMLBase.successful_retcode(sol)
            @test sol[x, 1] ≈ 1.0
            @test sol[y, 1] ≈ cbrt(4)
            @test sol.ps[p] ≈ 1.0
            @test sol.ps[q] ≈ sqrt(2)
        end

        @testset "Explicit OverrideInit" begin
            # Test with explicit OverrideInit
            sol = solve(prob, daskr(); initializealg = SciMLBase.OverrideInit())
            @test SciMLBase.successful_retcode(sol)
            @test sol[x, 1] ≈ 1.0
            @test sol[y, 1] ≈ cbrt(4)
            @test sol.ps[p] ≈ 1.0
            @test sol.ps[q] ≈ sqrt(2)
        end

        @testset "CheckInit" begin
            # CheckInit should fail with incomplete initial values
            @test_throws Any solve(prob, daskr(); initializealg = SciMLBase.CheckInit())

            # Create a problem with all correct initial values specified
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

            # This should work since all values are correct
            sol_correct = solve(prob_correct, daskr(); initializealg = SciMLBase.CheckInit())
            @test SciMLBase.successful_retcode(sol_correct)
        end
    end
end

@testset "Simple Robertson DAE with MTK" begin
    # Test the Robertson DAE from the Sundials examples
    @variables y1(t) [guess = 1.0] y2(t) [guess = 0.0] y3(t) [guess = 0.0]
    
    @mtkcompile rob_sys = System(
        [D(y1) ~ -0.04 * y1 + 1.0e4 * y2 * y3,
         D(y2) ~ 0.04 * y1 - 1.0e4 * y2 * y3 - 3.0e7 * y2^2,
         y1 + y2 + y3 ~ 1.0],  # algebraic constraint
        t
    )
    
    # Create problem with consistent initial conditions
    prob = DAEProblem(
        rob_sys,
        [D(y1) => -0.04, D(y2) => 0.04],  # derivatives
        (0.0, 100.0)
    )
    
    sol = solve(prob, daskr())
    @test SciMLBase.successful_retcode(sol)
    
    # Check mass conservation (y1 + y2 + y3 = 1)
    for i in 1:length(sol.t)
        @test sol[y1, i] + sol[y2, i] + sol[y3, i] ≈ 1.0 atol=1e-5
    end
end
