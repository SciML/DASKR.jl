# DASKR.jl JuliaDiffEq common algorithms

using Compat, DiffEqBase
import DiffEqBase: solve


# Abstract Types
abstract DASKRDAEAlgorithm{LinearSolver} <: AbstractDAEAlgorithm

# DAE Algorithms
immutable Algorithm{LinearSolver} <: DASKRDAEAlgorithm{LinearSolver} end
Algorithm(;linear_solver=:Dense) = Algorithm{linear_solver}()


## Solve for DAEs uses raw_solver

function solve{uType,duType,tType,isinplace,F,LinearSolver}(
    prob::AbstractDAEProblem{uType,duType,tType,isinplace,F},
    alg::DASKRDAEAlgorithm{LinearSolver},
    timeseries = [], ts = [], ks = [];
    callback = () -> nothing, abstol = 1/10^6, reltol = 1/10^3,
    saveat = Float64[], adaptive = true, maxiter = Int(1e5),
    timeseries_errors = true, save_timeseries = true,
    userdata = nothing, isdiff = fill(true, length(prob.u0)), kwargs...)

    tspan = prob.tspan
    t0 = tspan[1]
    T = tspan[end]

    save_ts = sort(unique([t0;saveat;T]))

    if T < save_ts[end]
        error("Final saving timepoint is past the solving timespan")
    end
    if t0 > save_ts[1]
        error("First saving timepoint is before the solving timespan")
    end

    if typeof(prob.u0) <: Number
        u0 = [prob.u0]
    else
        u0 = vec(deepcopy(prob.u0))
    end

    if typeof(prob.du0) <: Number
        du0 = [prob.du0]
    else
        du0 = vec(deepcopy(prob.du0))
    end

    sizeu = size(prob.u0)
    sizedu = size(prob.du0)

    ### Fix the more general function to DASKR allowed style
    if !isinplace && (typeof(prob.u0)<:Vector{Float64} || typeof(prob.u0)<:Number)
        f! = (t,u,du,out) -> (out[:] = prob.f(t,u,du); nothing)
    elseif !isinplace && typeof(prob.u0)<:AbstractArray
        f! = (t,u,du,out) -> (out[:] = vec(prob.f(t,reshape(u,sizeu),reshape(du,sizedu))); nothing)
    elseif typeof(prob.u0)<:Vector{Float64}
        f! = prob.f
    else # Then it's an in-place function on an abstract array
        f! = (t,u,du,out) -> (prob.f(t,reshape(u,sizeu),reshape(du,sizedu),out);
                          u = vec(u); du=vec(du); 0)
    end

    id = Int32[x ? 1 : -1 for x in isdiff]
    tstart = 0.0
    tstop = 50.0
    Nsteps = 500

    tstep = tstop / Nsteps
    tout = [tstep]
    idid = Int32[0]
    info = zeros(Int32, 20)

    info[3] = save_timeseries
    info[11] = 0
    info[16] = 0    # == 1 to ignore algebraic variables in the error calculation
    info[17] = 0
    info[18] = 2    # more initialization info
    N = Int32[length(u0)]
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
    ipar = Int32[length(u0), nrt[1], length(u0)]
    res = DASKR.res_c(f!)
    rt = Int32[0]
    jac = Int32[0]
    psol = Int32[0]

    ures = Vector{Vector{Float64}}()
    ts   = [t0]

    u = copy(u0)
    du = copy(du0)
    # The Inner Loops : Style depends on save_timeseries
    for k in 1:length(save_ts)
        tout = [save_ts[k]]
        while t[1] < save_ts[k]
            DASKR.unsafe_solve(res, N, t, u, du, tout, info, rtol, atol, idid, rwork, lrw, iwork, liw, rpar, ipar, jac, psol, rt, nrt, jroot)
            push!(ures,copy(u))
            push!(ts, t[1])
        end
    end
    ### Finishing Routine

    timeseries = Vector{uType}(0)
    if typeof(prob.u0)<:Number
        for i=1:length(ures)
            push!(timeseries,ures[i][1])
        end
    else
        for i=1:length(ures)
            push!(timeseries,reshape(ures[i],sizeu))
        end
    end

    build_solution(prob,alg,ts,timeseries,
                      timeseries_errors = timeseries_errors)
end
