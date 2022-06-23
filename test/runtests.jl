using DASKR
using Test
using DiffEqBase

# Test the raw interface

function vanderpol(t, y, yp, res)
    res[1] = yp[1] - y[1] * (1 - y[2]^2) + y[2]
    res[2] = yp[2] - y[1]
end

let
    y = [0.0, 1.0]
    yp = [-1.0, 0.0]
    id = Int32[1, 1]
    tstart = 0.0
    tstop = 50.0
    Nsteps = 500
    abstol = 1e-4
    reltol = 1e-4

    tstep = tstop / Nsteps
    tout = [tstep]
    idid = Int32[0]
    info = zeros(Int32, 20)

    info[11] = 0
    info[16] = 0    # == 1 to ignore algebraic variables in the error calculation
    info[17] = 0
    info[18] = 2    # more initialization info
    N = Int32[length(y)]
    t = [tstart]
    nrt = Int32[0]
    rpar = [0.0]
    rtol = [reltol]
    atol = [abstol]
    lrw = Int32[N[1]^3 + 9 * N[1] + 60 + 3 * nrt[1]]
    rwork = zeros(lrw[1])
    liw = Int32[2 * N[1] + 40]
    iwork = zeros(Int32, liw[1])
    iwork[40 .+ (1:N[1])] = id
    jroot = zeros(Int32, max(nrt[1], 1))
    ipar = Int32[length(y), nrt[1], length(y)]
    res = DASKR.res_c(vanderpol)
    rt = Int32[0]
    jac = Int32[0]
    psol = Int32[0]
    DASKR.unsafe_solve(res, N, t, y, yp, tout, info, rtol, atol, idid, rwork, lrw, iwork,
                       liw, rpar, ipar, jac, psol, rt, nrt, jroot)
    @show (t, y, yp)
    tout = [5.0]
    DASKR.unsafe_solve(res, N, t, y, yp, tout, info, rtol, atol, idid, rwork, lrw, iwork,
                       liw, rpar, ipar, jac, psol, rt, nrt, jroot)
    @show (t, y, yp)
end

# Test the JuliaDiffEq common interface
function resrob(r, yp, y, p, tres)
    r[1] = -0.04 * y[1] + 1.0e4 * y[2] * y[3]
    r[2] = -r[1] - 3.0e7 * y[2] * y[2] - yp[2]
    r[1] -= yp[1]
    r[3] = y[1] + y[2] + y[3] - 1.0
end

function testjac(res, du, u, p, t)
    res[1] = du[1] - 1.5 * u[1] + 1.0 * u[1] * u[2]
    res[2] = du[2] + 3 * u[2] - u[1] * u[2]
end

jac_called = false
function testjac_jac(J, du, u, p, gamma, t)
    global jac_called
    jac_called = true
    J[1, 1] = gamma - 1.5 + 1.0 * u[2]
    J[1, 2] = 1.0 * u[1]
    J[2, 1] = -1 * u[2]
    J[2, 2] = gamma + 3 - u[1]
    nothing
end

let
    u0 = [1.0, 0, 0]
    du0 = [-0.04, 0.04, 0.0]
    prob = DAEProblem(resrob, du0, u0, (0.0, 100000.0))
    dt = 1000
    saveat = float(collect(0:dt:100000))
    sol = solve(prob, daskr())
    @test length(sol.t) > 2
    sol = solve(prob, daskr(), save_everystep = false)
    @test length(sol.u) == length(sol.t) == 2
    prob2 = DAEProblem(resrob, du0, u0, (0.0, 100000.0),
                       differential_vars = [true, true, false])
    sol = solve(prob2, daskr(), saveat = saveat)
    @test sol.t == saveat
    sol = solve(prob2, daskr(), saveat = dt)
    @test sol.t == saveat
    sol = solve(prob2, daskr(), saveat = saveat,
                save_everystep = true)
    @test minimum([t ∈ sol.t for t in saveat])
    sol = solve(prob, daskr(), saveat = saveat, save_everystep = true)
    @test intersect(sol.t, saveat) == saveat

    # Test for callback
    @test_throws ErrorException solve(prob, daskr(), saveat = saveat,
                                      save_everystep = true,
                                      callback = (() -> true))

    # Check for warnings
    @info "Testing for Compatibility Warnings"
    sol = solve(prob, daskr(), saveat = saveat, save_everystep = true,
                verbose = true, save_idxs = true, d_discontinuities = true,
                isoutofdomain = true,
                unstable_check = true, calck = true, progress = true,
                dtmin = 1, dtmax = 2, dense = true,
                internalnorm = 0, gamma = 0.5, beta1 = 1.23, beta2 = 2.34,
                qmin = 1.0, qmax = 2.0)

    prob3 = DAEProblem(testjac, [0.5, -2.0], ones(2), (0.0, 10.0))
    sol = solve(prob3, daskr())
    @test maximum(sol[end]) < 2 #should be cyclic

    # inconsistent initial conditions
    function f!(res, du, u, p, t)
        res[1] = du[1] - 1.01
        return
    end
    u0 = [0.0]
    tspan = (0.0, 10.0)
    du0 = [0.0]
    dae_prob = DAEProblem(f!, du0, u0, tspan, differential_vars = [true])
    sol = solve(dae_prob, daskr())

    # Jacobian
    function f2!(res, du, u, p, t)
        res[1] = 1.01du[1]
        return
    end

    function f2_jac!(out, du, u, p, gamma, t)
        global jac_called
        jac_called = true
        out[1] = 1.01
    end
    u0 = [0.0]
    tspan = (0.0, 10.0)
    du0 = [0.0]
    dae_prob = DAEProblem(DAEFunction(f2!, jac = f2_jac!),
                          du0, u0, tspan, differential_vars = [true])
    sol = solve(dae_prob, daskr())
    @test jac_called
    nothing
end
