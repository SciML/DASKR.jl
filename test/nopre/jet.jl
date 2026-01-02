using DASKR, DiffEqBase, Test, JET

@testset "JET static analysis" begin
    # Test that the callback functions are type-stable
    test_res = (t, y, yp, res) -> (res[1] = yp[1] - y[1]; nothing)
    @test_opt target_modules = (DASKR,) DASKR.res_c(test_res)

    # Test for actual errors (not optimizations) in the daskr constructor
    # The constructor has runtime dispatch due to linear_solver being a type parameter,
    # which is expected behavior for this design pattern
    @test_call target_modules = (DASKR,) daskr()
    @test_call target_modules = (DASKR,) daskr(linear_solver = :Banded, jac_upper = 1, jac_lower = 1)
    @test_call target_modules = (DASKR,) daskr(linear_solver = :SPIGMR)
end
