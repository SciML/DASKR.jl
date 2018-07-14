# DASKR.jl JuliaDiffEq common algorithms

using Reexport
@reexport using DiffEqBase
import DiffEqBase: solve

# Abstract Types
@compat abstract type DASKRDAEAlgorithm{LinearSolver} <: DiffEqBase.AbstractDAEAlgorithm end

# DAE Algorithms
struct daskr{LinearSolver,NNEA,MKI,MKV} <: DASKRDAEAlgorithm{LinearSolver}
    jac_upper::Int
    jac_lower::Int
    max_order::Int
    non_negativity_enforcement::Int
    non_negativity_enforcement_array::NNEA
    max_krylov_iters::MKI
    num_krylov_vectors::MKV
    max_number_krylov_restarts::Int
    krylov_convergence_test_constant::Float64
    exclude_algebraic_errors::Bool
end
Base.@pure function daskr(;linear_solver=:Dense,
                  jac_upper=0,jac_lower=0,max_order = 5,
                  non_negativity_enforcement = 0,
                  non_negativity_enforcement_array = nothing,
                  max_krylov_iters = nothing,
                  num_krylov_vectors = nothing,
                  max_number_krylov_restarts = 5,
                  krylov_convergence_test_constant = 0.05,
                  exclude_algebraic_errors = false)
    daskr{linear_solver,typeof(non_negativity_enforcement_array),
          typeof(max_krylov_iters),typeof(num_krylov_vectors)}(
          jac_upper,jac_lower,max_order,non_negativity_enforcement,
                         non_negativity_enforcement_array,
                         max_krylov_iters,
                         num_krylov_vectors,
                         max_number_krylov_restarts,
                         krylov_convergence_test_constant,
                         exclude_algebraic_errors)
end

export daskr

## Solve for DAEs uses raw_solver

function DiffEqBase.solve(
    prob::DiffEqBase.AbstractDAEProblem{uType,duType,tupType,isinplace},
    alg::DASKRDAEAlgorithm{LinearSolver},
    timeseries = [], ts = [], ks = [];
    verbose=true,
    callback = nothing, abstol = 1/10^6, reltol = 1/10^3,
    saveat = Float64[], adaptive = true, maxiters = Int(1e5),
    timeseries_errors = true, save_everystep = isempty(saveat),
    dense = save_everystep && isempty(saveat),
    save_start = save_everystep || isempty(saveat) || typeof(saveat) <: Number ?
                 true : prob.tspan[1] in saveat,
    save_timeseries = nothing, dtmax = nothing,
    userdata = nothing, dt = nothing,
    kwargs...) where {uType,duType,tupType,isinplace,LinearSolver}

    tType = eltype(tupType)

    if verbose
        warned = !isempty(kwargs) && check_keywords(alg, kwargs, warnlist)
        if !(typeof(prob.f) <: DiffEqBase.AbstractParameterizedFunction)
            if DiffEqBase.has_tgrad(prob.f)
                @warn("Explicit t-gradient given to this stiff solver is ignored.")
                warned = true
            end
        end
        warned && warn_compat()
    end

    if callback != nothing || prob.callback != nothing
        error("DASKR is not compatible with callbacks.")
    end

    tspan = prob.tspan
    t0 = tspan[1]
    T = tspan[end]

    if typeof(saveat) <: Number
        if (tspan[1]:saveat:tspan[end])[end] == tspan[end]
          saveat_vec = convert(Vector{tType},collect(tType,tspan[1]+saveat:saveat:tspan[end]))
        else
          saveat_vec = convert(Vector{tType},collect(tType,tspan[1]+saveat:saveat:(tspan[end]-saveat)))
        end
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
        f! = (out,du,u,p,t) -> (out[:] = prob.f(du,u,p,t); nothing)
    elseif !isinplace && typeof(prob.u0)<:AbstractArray
        f! = (out,du,u,p,t) -> (out[:] = vec(prob.f(reshape(du,sizedu),reshape(u,sizeu),p,t)); nothing)
    elseif typeof(prob.u0)<:Vector{Float64}
        f! = prob.f
    else # Then it's an in-place function on an abstract array
        f! = (out,du,u,p,t) -> (prob.f(out,reshape(du,sizedu),reshape(u,sizeu),p,t); 0)
    end

    if prob.differential_vars == nothing
        id = ones(Int32,length(u0))
    else
        id = Int32[x ? 1 : -1 for x in prob.differential_vars]
    end

    tout = [0.0]
    idid = Int32[0]
    info = zeros(Int32, 20)

    info[3] = save_everystep
    info[11] = 0
    info[16] = 0    # == 1 to ignore algebraic variables in the error calculation
    info[17] = 0
    info[18] = 2    # more initialization info

    if LinearSolver == :Banded
        info[6] = 1
    end

    N = Int32[length(u0)]
    t = [prob.tspan[1]]
    nrt = Int32[0]
    rpar = [0.0]
    rtol = [reltol]
    atol = [abstol]
    lrw = Int32[N[1]^3 + 9 * N[1] + 60 + 3 * nrt[1]]
    rwork = zeros(lrw[1])

    liw = Int32[2*N[1] + 40]
    iwork = zeros(Int32, liw[1])
    iwork[1] = alg.jac_lower
    iwork[2] = alg.jac_upper
    iwork[40 .+ (1:N[1])] = id

    if dtmax != nothing
        info[7] = 1
        rwork[2] = dtmax
    end

    if dt != nothing
        info[8] = 1
        rwork[3] = dt
    end

    if alg.non_negativity_enforcement != 0
        info[10] = alg.non_negativity_enforcement
        iwork[40 + (1:N[1])] .= alg.non_negativity_enforcement_array
    end

    if alg.max_order != 5
        info[9] = 1
        iwork[3] = alg.max_order
    end

    if LinearSolver == :SPIGMR
        info[12] = 1
        if alg.max_krylov_iters != nothing
            iwork[24] = alg.max_krylov_iters
        end
        if alg.num_krylov_vectors != nothing
            iwork[25] = alg.num_krylov_vectors
        end
        iwork[26] = alg.max_number_krylov_restarts
        rwork[10] = alg.krylov_convergence_test_constant
    end

    if alg.exclude_algebraic_errors
        info[16] = 1
    end

    liw = Int32[2*N[1] + 40]

    jroot = zeros(Int32, max(nrt[1], 1))
    ipar = Int32[length(u0), nrt[1], length(u0)]
    res = DASKR.common_res_c(f!,prob.p)
    rt = Int32[0]

    if DiffEqBase.has_jac(f!)
      jac = common_jac_c(f!,prob.p)
      info[5] = 1 # Enables Jacobian
    else
      jac = Int32[0]
    end
    psol = Int32[0]

    ures = Vector{Float64}[]
    dures = Vector{Float64}[]
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
            if idid[1] < 0
                break
            end
            push!(ures,copy(u))
            push!(ts, t[1])
            dense && push!(dures,copy(du))
        end
        if idid[1] < 0
            break
        end
    end
    ### Finishing Routine

    if idid[1] == -1
        retcode = :MaxIters
    elseif idid[1] == -7 || idid[1] == -9 || idid[1] == -10 || idid[1] == -14
        retcode = :ConvergenceFailure
    elseif idid[1] == -12
        retcode = :InitialFailure
    elseif idid[1] < 0
        retcode = :Failure
    else
        retcode = :Success
    end

    timeseries = Vector{uType}(undef,0)
    if typeof(prob.u0)<:Number
        for i=start_idx:length(ures)
            push!(timeseries,ures[i][1])
        end
    else
        for i=start_idx:length(ures)
            push!(timeseries,reshape(ures[i],sizeu))
        end
    end

    DiffEqBase.build_solution(prob,alg,ts,timeseries,
                   du = dures,
                   dense = dense,
                   timeseries_errors = timeseries_errors,
                   retcode = retcode)
end
