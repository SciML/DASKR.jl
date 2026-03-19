# ModelingToolkit DAE initialization tests for DASKR.jl
# These tests verify that DASKR works correctly with MTK-generated DAE problems.
#
# Note: DASKR's initialization support is simpler than Sundials. Complex MTK
# initialization features (like solving for missing parameters) are not yet
# fully supported. These tests focus on the working functionality.

using ModelingToolkit, SciMLBase, DASKR, Test, DiffEqBase
using SymbolicIndexingInterface
using ModelingToolkit: t_nounits as t, D_nounits as D

@testset "Robertson DAE with MTK" begin
    # Test the Robertson DAE - a classic DAE test problem
    # This is an index-1 DAE with one algebraic constraint
    @variables y1(t) [guess = 1.0] y2(t) [guess = 0.0] y3(t) [guess = 0.0]

    @mtkcompile rob_sys = System(
        [
            D(y1) ~ -0.04 * y1 + 1.0e4 * y2 * y3,
            D(y2) ~ 0.04 * y1 - 1.0e4 * y2 * y3 - 3.0e7 * y2^2,
            y1 + y2 + y3 ~ 1.0,
        ],  # algebraic constraint
        t
    )

    # Create problem with consistent initial conditions
    # y1=1, y2=0, y3=0 satisfies the algebraic constraint y1+y2+y3=1
    # D(y1) = -0.04*1 + 0 = -0.04
    # D(y2) = 0.04*1 - 0 - 0 = 0.04
    prob = DAEProblem(
        rob_sys,
        [D(y1) => -0.04, D(y2) => 0.04],  # derivatives
        (0.0, 100.0)
    )

    @testset "DefaultInit (consistent ICs)" begin
        sol = solve(prob, daskr())
        @test SciMLBase.successful_retcode(sol)

        # Check initial conditions
        @test sol[y1, 1] ≈ 1.0
        @test sol[y2, 1] ≈ 0.0
        @test sol[y3, 1] ≈ 0.0

        # Check mass conservation throughout (y1 + y2 + y3 = 1)
        for i in 1:length(sol.t)
            @test sol[y1, i] + sol[y2, i] + sol[y3, i] ≈ 1.0 atol = 1.0e-5
        end
    end

    @testset "NoInit" begin
        sol = solve(prob, daskr(); initializealg = DiffEqBase.NoInit())
        @test SciMLBase.successful_retcode(sol)
    end

    @testset "CheckInit (consistent ICs should pass)" begin
        sol = solve(prob, daskr(); initializealg = SciMLBase.CheckInit())
        @test SciMLBase.successful_retcode(sol)
    end
end

@testset "Simple DAE with parameter" begin
    # A simpler DAE with a constant parameter (no initialization required for parameter)
    @variables x(t) [guess = 0.0] y(t) [guess = 1.0]
    @parameters k = 0.5

    @mtkcompile simple_sys = System(
        [
            D(x) ~ -k * x + y,
            x + y ~ 1.0,
        ],  # algebraic constraint
        t
    )

    # Initial conditions: x=0, y=1 satisfies x+y=1
    # D(x) = -0.5*0 + 1 = 1
    prob = DAEProblem(
        simple_sys,
        [D(x) => 1.0],
        (0.0, 10.0)
    )

    @testset "Solve with DefaultInit" begin
        sol = solve(prob, daskr())
        @test SciMLBase.successful_retcode(sol)

        # Check constraint is satisfied
        for i in 1:length(sol.t)
            @test sol[x, i] + sol[y, i] ≈ 1.0 atol = 1.0e-5
        end
    end

    @testset "Check symbolic indexing" begin
        sol = solve(prob, daskr())
        @test SciMLBase.successful_retcode(sol)

        # Verify we can access parameters symbolically
        @test sol.ps[k] ≈ 0.5
    end
end
