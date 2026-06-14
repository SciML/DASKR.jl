using DASKR, DiffEqBase, Test, JET

@testset "JET static analysis" begin
    # Test that res_c doesn't have type errors
    # Note: res_c has inherent runtime dispatch due to the ARM-compatible design that uses
    # unsafe_pointer_to_objref to pass user data through rpar. This loses type information,
    # causing userdata.fun to be typed as Any. This is acceptable because:
    # 1. The callback is called from C code, so Julia's type inference cannot optimize it anyway
    # 2. The alternative (closures with @cfunction) fails on ARM platforms
    test_res = (t, y, yp, res) -> (res[1] = yp[1] - y[1]; nothing)
    @test_call target_modules = (DASKR,) DASKR.res_c(test_res)

    # Test for actual errors (not optimizations) in the daskr constructor
    # The constructor has runtime dispatch due to linear_solver being a type parameter,
    # which is expected behavior for this design pattern
    @test_call target_modules = (DASKR,) daskr()
    @test_call target_modules = (DASKR,) daskr(linear_solver = :Banded, jac_upper = 1, jac_lower = 1)
    @test_call target_modules = (DASKR,) daskr(linear_solver = :SPIGMR)
end
