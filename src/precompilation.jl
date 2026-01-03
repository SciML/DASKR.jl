# Precompilation workload for DASKR.jl
using PrecompileTools

@setup_workload begin
    # Minimal setup - avoid heavy dependencies
    @compile_workload begin
        # Precompile the most common operations

        # 1. Precompile the C callback wrappers
        # These now return (callback, userdata) tuples
        test_res = (t, y, yp, res) -> (res[1] = yp[1] - y[1]; nothing)
        res_callback, res_userdata = res_c(test_res)

        test_rt = (t, y, yp, rval) -> (rval[1] = y[1] - 0.5; nothing)
        rt_callback, rt_userdata = rt_c(test_rt)

        test_jac = (t, y, yp, pd, cj) -> (pd[1, 1] = cj - 1.0; nothing)
        jac_callback, jac_userdata = jac_c(test_jac)

        # 2. Precompile common_res_c and common_jac_c with typical parameter types
        test_f! = (out, du, u, p, t) -> (out[1] = du[1] - u[1]; nothing)
        common_res_callback1, common_res_userdata1 = common_res_c(test_f!, nothing)
        common_res_callback2, common_res_userdata2 = common_res_c(test_f!, (1.0,))

        # 3. Precompile daskr algorithm constructors
        daskr()
        daskr(linear_solver = :Banded, jac_upper = 1, jac_lower = 1)
        daskr(linear_solver = :SPIGMR)

        # 4. Precompile DAE problem construction
        function precompile_residual(r, yp, y, p, t)
            r[1] = yp[1] - y[1]
            r[2] = y[1] + y[2] - 1.0
        end

        u0 = [0.5, 0.5]
        du0 = [0.5, -0.5]
        tspan = (0.0, 1.0)

        # Create DAE problem to precompile problem construction
        prob = DiffEqBase.DAEProblem(precompile_residual, du0, u0, tspan,
            differential_vars = [true, false])

        # Note: We cannot call solve() during precompilation because
        # unsafe_solve requires `lib` which is set in __init__().
        # The solve path will be compiled on first use.
    end
end
