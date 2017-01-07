using DASKR
using Base.Test
using DiffEqBase

include("events.jl")

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
    liw = Int32[2*N[1] + 40]
    iwork = zeros(Int32, liw[1])
    iwork[40 + (1:N[1])] = id
    jroot = zeros(Int32, max(nrt[1], 1))
    ipar = Int32[length(y), nrt[1], length(y)]
    res = DASKR.res_c(vanderpol)
    rt = Int32[0]
    jac = Int32[0]
    psol = Int32[0]
    DASKR.unsafe_solve(res, N, t, y, yp, tout, info, rtol, atol, idid, rwork, lrw, iwork, liw, rpar, ipar, jac, psol, rt, nrt, jroot)
    @show (t, y, yp)
    tout = [5.0]
    DASKR.unsafe_solve(res, N, t, y, yp, tout, info, rtol, atol, idid, rwork, lrw, iwork, liw, rpar, ipar, jac, psol, rt, nrt, jroot)
    @show (t, y, yp)
end


# Test the JuliaDiffEq common interface
function resrob(tres, y, yp, r)
    r[1]  = -0.04*y[1] + 1.0e4*y[2]*y[3]
    r[2]  = -r[1] - 3.0e7*y[2]*y[2] - yp[2]
    r[1] -=  yp[1]
    r[3]  =  y[1] + y[2] + y[3] - 1.0
end
let
    u0 = [1.0, 0, 0]
    du0 = [-0.04, 0.04, 0.0]
    prob = DAEProblem(resrob,u0,du0,(0.0,100000.0))
    dt = 1000
    saveat = float(collect(0:dt:100000))
    sol = solve(prob, DASKR.Algorithm())
    sol = solve(prob, DASKR.Algorithm(),save_timeseries=false)
    sol = solve(prob, DASKR.Algorithm(), saveat = saveat, isdiff = [true, true, false])
    sol = solve(prob, DASKR.Algorithm(), saveat = saveat,
                      save_timeseries = false,
                      isdiff = [true, true, false])
    sol = solve(prob, DASKR.Algorithm(), saveat = saveat, save_timeseries = false)

    @test intersect(sol.t, saveat) == saveat

    @test sol.t == saveat
end
