# DASKR.jl JuliaDiffEq common algorithms

using DiffEqBase
import DiffEqBase: solve

# Abstract Types
@compat abstract type DASKRDAEAlgorithm{LinearSolver} <: AbstractDAEAlgorithm end

# DAE Algorithms
immutable daskr{LinearSolver} <: DASKRDAEAlgorithm{LinearSolver} end
Base.@pure daskr(;linear_solver=:Dense) = daskr{linear_solver}()

export daskr

## Solve for DAEs uses raw_solver

function solve{uType,duType,tType,isinplace,LinearSolver}(
    prob::AbstractDAEProblem{uType,duType,tType,isinplace},
    alg::DASKRDAEAlgorithm{LinearSolver},
    timeseries = [], ts = [], ks = [];

    verbose=true,
    callback = nothing, abstol = 1/10^6, reltol = 1/10^3,
    saveat = Float64[], adaptive = true, maxiter = Int(1e5),
    timeseries_errors = true, save_everystep = isempty(saveat), dense = save_everystep,
    save_start = true, save_timeseries = nothing,
    userdata = nothing,
    kwargs...)

    if verbose
        warned = !isempty(kwargs) && check_keywords(alg, kwargs, warnlist)
        if !(typeof(prob.f) <: AbstractParameterizedFunction)
            if has_tgrad(prob.f)
                warn("Explicit t-gradient given to this stiff solver is ignored.")
                warned = true
            end
            if has_jac(prob.f)
                warn("Explicit Jacobian given to this stiff solver is ignored.")
                warned = true
            end
        end
        warned && warn_compat()
    end

    if save_timeseries != nothing
        warn("save_timeseries is deprecated. Use save_everystep instead")
        save_everystep = save_timeseries
    end

    if callback != nothing || prob.callback != nothing
        error("DASKR is not compatible with callbacks.")
    end

    tspan = prob.tspan
    t0 = tspan[1]
    T = tspan[end]

    if typeof(saveat) <: Number
        saveat_vec = convert(Vector{tType},saveat+tspan[1]:saveat:(tspan[end]-saveat))
        # Exclude the endpoint because of floating point issues
    else
        saveat_vec =  convert(Vector{tType},collect(saveat))
    end

    if !isempty(saveat_vec) && saveat_vec[end] == tspan[2]
        pop!(saveat_vec)
    end

    if !isempty(saveat_vec) && saveat_vec[1] == tspan[1]
        save_ts = sort(unique([saveat_vec[2:end];tspan[2]]))
    else
        save_ts = sort(unique([saveat_vec;tspan[2]]))
    end

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

    if prob.differential_vars == nothing
        id = ones(Int32,length(u0))
    else
        id = Int32[x ? 1 : -1 for x in prob.differential_vars]
    end

    tstart = 0.0
    tstop = 50.0
    Nsteps = 500

    tstep = tstop / Nsteps
    tout = [tstep]
    idid = Int32[0]
    info = zeros(Int32, 20)

    info[3] = save_everystep
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
    dures = Vector{Vector{Float64}}()
    save_start ? ts = [t0] : ts = Float64[]
    save_start ? start_idx = 1 : start_idx = 2
    save_start && push!(ures, copy(u0))
    save_start && dense && push!(dures, copy(du0))

    u = copy(u0)
    du = copy(du0)
    # The Inner Loops : Style depends on save_timeseries
    for k in start_idx:length(save_ts)
        tout = [save_ts[k]]
        while t[1] < save_ts[k]
            DASKR.unsafe_solve(res, N, t, u, du, tout, info, rtol, atol, idid, rwork,
                               lrw, iwork, liw, rpar, ipar, jac, psol, rt, nrt, jroot)
            push!(ures,copy(u))
            push!(ts, t[1])
            dense && push!(dures,copy(du))
        end
    end
    ### Finishing Routine



    timeseries = Vector{uType}(0)
    if typeof(prob.u0)<:Number
        for i=start_idx:length(ures)
            push!(timeseries,ures[i][1])
        end
    else
        for i=start_idx:length(ures)
            push!(timeseries,reshape(ures[i],sizeu))
        end
    end

    build_solution(prob,alg,ts,timeseries,
                   du = dures,
                   dense = dense,
                   timeseries_errors = timeseries_errors,
                   retcode = :Success)
end
