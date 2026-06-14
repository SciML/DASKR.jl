using DASKR
using Test

# Test the raw interface

function vanderpol(t, y, yp, res)
    res[1] = yp[1] - y[1] * (1 - y[2]^2) + y[2]
    return res[2] = yp[2] - y[1]
end

let
    y = [0.0, 1.0]
    yp = [-1.0, 0.0]
    id = Int32[1, 1]
    tstart = 0.0
    tstop = 50.0
    Nsteps = 500
    abstol = 1.0e-4
    reltol = 1.0e-4

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
    rtol = [reltol]
    atol = [abstol]
    lrw = Int32[N[1]^3 + 9 * N[1] + 60 + 3 * nrt[1]]
    rwork = zeros(lrw[1])
    liw = Int32[2 * N[1] + 40]
    iwork = zeros(Int32, liw[1])
    iwork[40 .+ (1:N[1])] = id
    jroot = zeros(Int32, max(nrt[1], 1))
    ipar = Int32[length(y), nrt[1], length(y)]
    # res_c now returns (callback, userdata) - userdata is passed as rpar
    res, rpar = DASKR.res_c(vanderpol)
    rt = Int32[0]
    jac = Int32[0]
    psol = Int32[0]
    DASKR.unsafe_solve(
        res, N, t, y, yp, tout, info, rtol, atol, idid, rwork, lrw, iwork,
        liw, rpar, ipar, jac, psol, rt, nrt, jroot
    )
    @show (t, y, yp)
    tout = [5.0]
    DASKR.unsafe_solve(
        res, N, t, y, yp, tout, info, rtol, atol, idid, rwork, lrw, iwork,
        liw, rpar, ipar, jac, psol, rt, nrt, jroot
    )
    @show (t, y, yp)
end
