using DASKR, BenchmarkTools

const SUITE = BenchmarkGroup()

# Robertson DAE
function dae_robertson!(resid, du, u, p, t)
    resid[1] = du[1] + 0.04 * u[1] - 1.0e4 * u[2] * u[3]
    resid[2] = du[2] - 0.04 * u[1] + 1.0e4 * u[2] * u[3] + 3.0e7 * u[2]^2
    resid[3] = u[1] + u[2] + u[3] - 1.0
    return nothing
end
u0 = [1.0, 0.0, 0.0]
du0 = [-0.04, 0.04, 0.0]
prob = DAEProblem(
    dae_robertson!, du0, u0, (0.0, 10.0);
    differential_vars = [true, true, false]
)

# =============================================================================
# DASKR solves
# =============================================================================

SUITE["solve"] = BenchmarkGroup()

SUITE["solve"]["robertson"] = @benchmarkable solve($prob, daskr())
SUITE["solve"]["robertson_dense"] = @benchmarkable solve(
    $prob, daskr(; linear_solver = :Dense)
)
SUITE["solve"]["long_tspan"] = @benchmarkable solve(
    $(
        DAEProblem(
            dae_robertson!, du0, u0, (0.0, 1000.0);
            differential_vars = [true, true, false]
        )
    ), daskr()
)
