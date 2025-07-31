# DASKR.jl JuliaDiffEq common algorithms

using Reexport
@reexport using DiffEqBase
import DiffEqBase: solve

# Abstract Types
abstract type DASKRDAEAlgorithm{LinearSolver} <: DiffEqBase.AbstractDAEAlgorithm end

# DAE Algorithms
"""
```julia
function daskr(;linear_solver=:Dense,
                  jac_upper=0,jac_lower=0,max_order = 5,
                  non_negativity_enforcement = 0,
                  non_negativity_enforcement_array = nothing,
                  max_krylov_iters = nothing,
                  num_krylov_vectors = nothing,
                  max_number_krylov_restarts = 5,
                  krylov_convergence_test_constant = 0.05,
                  exclude_algebraic_errors = false)
```

This is a wrapper for the well-known DASKR algorithm.

DASKR is a solver for systems of differential-algebraic equations (DAEs).
It includes options for both direct and iterative (Krylov) methods for the
solution of the linear systems arising at each (implicit) time step.
DASKR is a variant of the DASPK package [1].  In addition to all the
capabilities of DASPK, DASKR includes the ability to find the roots of a
given set of functions while integrating the DAE system.

In contrast to the older DASSL package, DASKR includes a procedure
for calculating consistent initial conditions for a large class of
problems (which includes semi-explicit index-1 systems) [2].  This
procedure includes options for inequality constraints on selected
components.  The package also includes an option to omit the algebraic
components from the local error control.

### Linear Solver Choice

Choices for the linear solver are:

- `:Dense`
- `:Banded`
- `:SPIGMR`, a Krylov method

### Other Keyword Arguments

- `jac_upper=0,jac_lower=0`: used for setting the upper and lower bands for a banded
  Jacobian. Defaults to 0. Ignored unless the linear solver is `:Banded`.
- `max_order = 5`: the maximum order of the BDF method.
- `non_negativity_enforcement = 0`, whether to enforce non-negativty in the solution.
  Defaults to `0` or false, can be set to `1` to enable.
- `non_negativity_enforcement_array = nothing`, an array of 0's and 1's for specifying
  non-negativity enforcement to a subset of states.
- `max_krylov_iters=nothing`: maximum number of iterations for the Krylov subspace linear
  solver before rejecting a step. Defaults to `nothing` or an automatic detection.
- `num_krylov_vectors=nothing`: maximum number of history states in the GMRES vector.
  Defaults to `nothing` or an automatic choice
- `max_number_krylov_restarts=5`: If you figure out what this is, open an issue or PR.
- `krylov_convergence_test_constant = 0.05`: Some constant in DASKR's convergence test.
- `exclude_algebraic_errors = false`: whether algebraic variables are included in the
  adaptive time stepping error check. Defaults to false.
"""
struct daskr{LinearSolver, NNEA, MKI, MKV} <: DASKRDAEAlgorithm{LinearSolver}
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
Base.@pure function daskr(; linear_solver = :Dense,
                          jac_upper = 0, jac_lower = 0, max_order = 5,
                          non_negativity_enforcement = 0,
                          non_negativity_enforcement_array = nothing,
                          max_krylov_iters = nothing,
                          num_krylov_vectors = nothing,
                          max_number_krylov_restarts = 5,
                          krylov_convergence_test_constant = 0.05,
                          exclude_algebraic_errors = false)
    daskr{linear_solver, typeof(non_negativity_enforcement_array),
          typeof(max_krylov_iters), typeof(num_krylov_vectors)}(jac_upper, jac_lower,
                                                                max_order,
                                                                non_negativity_enforcement,
                                                                non_negativity_enforcement_array,
                                                                max_krylov_iters,
                                                                num_krylov_vectors,
                                                                max_number_krylov_restarts,
                                                                krylov_convergence_test_constant,
                                                                exclude_algebraic_errors)
end

export daskr

## Solve for DAEs uses raw_solver

function DiffEqBase.__solve(prob::DiffEqBase.AbstractDAEProblem{uType, duType, tupType,
                                                                isinplace},
                            alg::DASKRDAEAlgorithm{LinearSolver},
                            timeseries = [], ts = [], ks = [];
                            verbose = true,
                            callback = nothing, abstol = 1 / 10^6, reltol = 1 / 10^3,
                            saveat = Float64[], adaptive = true, maxiters = Int(1e5),
                            timeseries_errors = true, save_everystep = isempty(saveat),
                            dense = save_everystep && isempty(saveat),
                            save_start = save_everystep || isempty(saveat) ||
                                             saveat isa Number ?
                                         true : prob.tspan[1] in saveat,
                            save_timeseries = nothing, dtmax = nothing,
                            userdata = nothing, dt = nothing, alias_u0 = false,
                            kwargs...) where {uType, duType, tupType, isinplace,
                                              LinearSolver}
    tType = eltype(tupType)

    if verbose
        warned = !isempty(kwargs) && check_keywords(alg, kwargs, warnlist)
        if !(prob.f isa DiffEqBase.AbstractParameterizedFunction)
            if DiffEqBase.has_tgrad(prob.f)
                @warn("Explicit t-gradient given to this stiff solver is ignored.")
                warned = true
            end
        end
        warned && warn_compat()
    end

    if callback !== nothing || :callback in keys(prob.kwargs)
        error("DASKR is not compatible with callbacks.")
    end

    tspan = prob.tspan
    t0 = tspan[1]
    T = tspan[end]

    if saveat isa Number
        if (tspan[1]:saveat:tspan[end])[end] == tspan[end]
            saveat_vec = convert(Vector{tType},
                                 collect(tType, (tspan[1] + saveat):saveat:tspan[end]))
        else
            saveat_vec = convert(Vector{tType},
                                 collect(tType,
                                         (tspan[1] + saveat):saveat:(tspan[end] - saveat)))
        end
    else
        saveat_vec = convert(Vector{tType}, collect(saveat))
    end

    if !isempty(saveat_vec) && saveat_vec[end] == tspan[2]
        pop!(saveat_vec)
    end

    if !isempty(saveat_vec) && saveat_vec[1] == tspan[1]
        save_ts = sort(unique([saveat_vec[2:end]; tspan[2]]))
    else
        save_ts = sort(unique([saveat_vec; tspan[2]]))
    end

    if T < save_ts[end]
        error("Final saving timepoint is past the solving timespan")
    end
    if t0 > save_ts[1]
        error("First saving timepoint is before the solving timespan")
    end

    if prob.u0 isa Number
        u0 = [prob.u0]
    else
        if alias_u0
            u0 = vec(prob.u0)
        else
            u0 = vec(deepcopy(prob.u0))
        end
    end

    if prob.du0 isa Number
        du0 = [prob.du0]
    else
        if alias_u0
            du0 = vec(prob.du0)
        else
            du0 = vec(deepcopy(prob.du0))
        end
    end

    sizeu = size(prob.u0)
    sizedu = size(prob.du0)

    ### Fix the more general function to DASKR allowed style
    if !isinplace && (prob.u0 isa Vector{Float64} || prob.u0 isa Number)
        f! = (out, du, u, p, t) -> (out[:] = prob.f(du, u, p, t); nothing)
    elseif !isinplace && prob.u0 isa AbstractArray
        f! = (out, du, u, p, t) -> (out[:] = vec(prob.f(reshape(du, sizedu),
                                                        reshape(u, sizeu), p, t));
                                    nothing)
    elseif prob.u0 isa Vector{Float64}
        f! = prob.f
    else # Then it's an in-place function on an abstract array
        f! = (out, du, u, p, t) -> (prob.f(out, reshape(du, sizedu), reshape(u, sizeu), p,
                                           t);
                                    0)
    end

    if prob.differential_vars === nothing
        id = ones(Int32, length(u0))
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

    liw = Int32[2 * N[1] + 40]
    iwork = zeros(Int32, liw[1])
    iwork[1] = alg.jac_lower
    iwork[2] = alg.jac_upper
    iwork[40 .+ (1:N[1])] = id

    if dtmax !== nothing
        info[7] = 1
        rwork[2] = dtmax
    end

    if dt !== nothing
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
        if alg.max_krylov_iters !== nothing
            iwork[24] = alg.max_krylov_iters
        end
        if alg.num_krylov_vectors !== nothing
            iwork[25] = alg.num_krylov_vectors
        end
        iwork[26] = alg.max_number_krylov_restarts
        rwork[10] = alg.krylov_convergence_test_constant
    end

    if alg.exclude_algebraic_errors
        info[16] = 1
    end

    liw = Int32[2 * N[1] + 40]

    jroot = zeros(Int32, max(nrt[1], 1))
    ipar = Int32[length(u0), nrt[1], length(u0)]
    res = DASKR.common_res_c(f!, prob.p)
    rt = Int32[0]

    if DiffEqBase.has_jac(f!)
        jac = common_jac_c(f!, prob.p)
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

    if alias_u0
        u = u0
        du = du0
    else
        u = copy(u0)
        du = copy(du0)
    end
    # The Inner Loops : Style depends on save_timeseries
    for k in start_idx:length(save_ts)
        tout = [save_ts[k]]
        while t[1] < save_ts[k]
            DASKR.unsafe_solve(res, N, t, u, du, tout, info, rtol, atol, idid, rwork,
                               lrw, iwork, liw, rpar, ipar, jac, psol, rt, nrt, jroot)
            if idid[1] < 0
                break
            end
            push!(ures, copy(u))
            push!(ts, t[1])
            dense && push!(dures, copy(du))
        end
        if idid[1] < 0
            break
        end
    end
    ### Finishing Routine

    if idid[1] == -1
        retcode = ReturnCode.MaxIters
    elseif idid[1] == -7 || idid[1] == -9 || idid[1] == -10 || idid[1] == -14
        retcode = ReturnCode.ConvergenceFailure
    elseif idid[1] == -12
        retcode = ReturnCode.InitialFailure
    elseif idid[1] < 0
        retcode = ReturnCode.Failure
    else
        retcode = ReturnCode.Success
    end

    timeseries = Vector{uType}(undef, 0)
    if prob.u0 isa Number
        for i in start_idx:length(ures)
            push!(timeseries, ures[i][1])
        end
    else
        for i in start_idx:length(ures)
            push!(timeseries, reshape(ures[i], sizeu))
        end
    end

    DiffEqBase.build_solution(prob, alg, ts, timeseries,
                              du = dures,
                              dense = dense,
                              timeseries_errors = timeseries_errors,
                              retcode = retcode)
end
