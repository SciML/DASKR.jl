# DAE Initialization support for DASKR.jl
# Following the pattern from Sundials.jl:
# https://github.com/SciML/Sundials.jl/blob/master/src/common_interface/initialize.jl

using DiffEqBase: BrownFullBasicInit, ShampineCollocationInit, DefaultInit, NoInit
using SciMLBase: CheckInit, OverrideInit, isinplace, ReturnCode, get_initial_values

"""
    perform_initialization!(prob, alg, u0, du0, p, t0, f!, abstol, reltol, initializealg, info, iwork, differential_vars)

Perform DAE initialization based on the specified initialization algorithm.
Returns `(u0, du0, p, success, init_type)` where `init_type` is the INFO(11) value for DASKR.

For DASKR, INFO(11) controls initialization:
- INFO(11) = 0: Assume initial conditions are consistent (default)
- INFO(11) = 1: Given Y_d (differential vars), calculate Y_a and Y'_d (BrownFullBasicInit)
- INFO(11) = 2: Given Y', calculate Y (ShampineCollocationInit)
"""
function perform_initialization! end

# DefaultInit - routes to OverrideInit if initialization_data exists, otherwise NoInit
# Note: We use NoInit (not CheckInit) by default because DASKR can handle some
# inconsistent initial conditions internally via its INFO(11) mechanism.
# This matches the original DASKR.jl behavior before explicit initialization support.
function perform_initialization!(
    prob, alg, u0, du0, p, t0, f!, abstol, reltol,
    initializealg::DefaultInit,
    info, iwork, differential_vars
)
    if prob.f.initialization_data !== nothing
        return perform_initialization!(
            prob, alg, u0, du0, p, t0, f!, abstol, reltol,
            OverrideInit(),
            info, iwork, differential_vars
        )
    else
        # Use NoInit by default to preserve DASKR's original behavior
        # Users can explicitly use CheckInit if they want strict validation
        return perform_initialization!(
            prob, alg, u0, du0, p, t0, f!, abstol, reltol,
            NoInit(),
            info, iwork, differential_vars
        )
    end
end

# NoInit - do nothing, assume user provided consistent initial conditions
function perform_initialization!(
    prob, alg, u0, du0, p, t0, f!, abstol, reltol,
    initializealg::NoInit,
    info, iwork, differential_vars
)
    return u0, du0, p, true, Int32(0)  # INFO(11) = 0
end

# CheckInit - verify that initial conditions satisfy the DAE constraints
function perform_initialization!(
    prob, alg, u0, du0, p, t0, f!, abstol, reltol,
    initializealg::CheckInit,
    info, iwork, differential_vars
)
    # Evaluate the DAE residual at the initial conditions
    residual = similar(u0)
    f!(residual, du0, u0, p, t0)
    
    # Check if residuals are within tolerance
    max_residual = maximum(abs.(residual))
    if max_residual >= abstol
        error(
            """
            DAE initialization failed with CheckInit: Initial conditions do not satisfy the DAE constraints.

            The residual norm is $(max_residual), which exceeds the tolerance $(abstol).

            Note that the initial conditions include both `du0` (derivatives) and `u0` (states),
            and the choice of derivatives must be compatible with the states.

            To resolve this issue, you have several options:
            1. Fix your initial conditions (both `du0` and `u0`) to satisfy the DAE constraints
            2. Use DASKR's built-in initialization: initializealg = BrownFullBasicInit()
               - This computes Y_a and Y'_d given Y_d (differential variables)
            3. Use Shampine's collocation initialization: initializealg = ShampineCollocationInit()
               - This computes Y given Y'
            4. If using ModelingToolkit, use: initializealg = OverrideInit()

            Example for automatic initialization:
            solve(prob, daskr(); initializealg = BrownFullBasicInit())
            """
        )
    end
    
    return u0, du0, p, true, Int32(0)  # INFO(11) = 0
end

# BrownFullBasicInit - use DASKR's built-in initialization (INFO(11) = 1)
# Given Y_d (differential variables), calculate Y_a (algebraic variables) and Y'_d
function perform_initialization!(
    prob, alg, u0, du0, p, t0, f!, abstol, reltol,
    initializealg::BrownFullBasicInit,
    info, iwork, differential_vars
)
    # First check if we actually need initialization
    residual = similar(u0)
    f!(residual, du0, u0, p, t0)
    
    tol = initializealg.abstol
    if any(abs.(residual) .>= tol)
        if differential_vars === nothing
            error("Must supply differential_vars argument to DAEProblem constructor to use BrownFullBasicInit.")
        end
        # Set INFO(11) = 1 for DASKR's built-in initialization
        # This tells DASKR to compute consistent Y_a and Y'_d given Y_d
        return u0, du0, p, true, Int32(1)
    end
    
    return u0, du0, p, true, Int32(0)
end

# ShampineCollocationInit - use DASKR's built-in initialization (INFO(11) = 2)
# Given Y', calculate Y
function perform_initialization!(
    prob, alg, u0, du0, p, t0, f!, abstol, reltol,
    initializealg::ShampineCollocationInit,
    info, iwork, differential_vars
)
    # First check if we actually need initialization
    residual = similar(u0)
    f!(residual, du0, u0, p, t0)
    
    if any(abs.(residual) .>= reltol)
        # Set INFO(11) = 2 for DASKR's built-in initialization
        # This tells DASKR to compute consistent Y given Y'
        return u0, du0, p, true, Int32(2)
    end
    
    return u0, du0, p, true, Int32(0)
end

# OverrideInit - use SciMLBase's get_initial_values (for MTK compatibility)
function perform_initialization!(
    prob, alg, u0, du0, p, t0, f!, abstol, reltol,
    initializealg::OverrideInit,
    info, iwork, differential_vars
)
    # Use SciMLBase's get_initial_values to compute consistent initial conditions
    # This is the pattern used by ModelingToolkit
    
    # We need to create a simple integrator-like object that SciMLBase can work with
    # Since DASKR doesn't have a proper integrator, we'll use get_initial_values directly
    
    # Check if get_initial_values is available and the problem has initialization_data
    if prob.f.initialization_data === nothing
        # No initialization data, fall back to CheckInit
        return perform_initialization!(
            prob, alg, u0, du0, p, t0, f!, abstol, reltol,
            CheckInit(),
            info, iwork, differential_vars
        )
    end
    
    # Create a minimal integrator-like object for get_initial_values
    # get_initial_values expects: prob, integrator, f, initializealg, Val(isinplace); nlsolve_alg, abstol, reltol
    integrator = DASKRInitIntegrator(prob, u0, du0, p, t0, abstol, reltol)
    
    try
        u0_new, p_new, success = get_initial_values(
            prob, integrator, prob.f, initializealg, Val(isinplace(prob));
            abstol=abstol, reltol=reltol
        )
        
        if !success
            return u0, du0, p, false, Int32(0)
        end
        
        # Update u0 in place
        if isinplace(prob)
            u0 .= vec(u0_new)
        else
            u0 = vec(u0_new)
        end
        
        return u0, du0, p_new, true, Int32(0)
    catch e
        # If get_initial_values fails, report the failure
        @warn "OverrideInit failed: $e"
        return u0, du0, p, false, Int32(0)
    end
end

"""
    DASKRInitIntegrator

A minimal integrator-like struct to interface with SciMLBase.get_initial_values.
This is needed because DASKR doesn't use a proper integrator interface.
"""
mutable struct DASKRInitIntegrator{uType, duType, pType, tType, aTol, rTol, probType}
    sol::Any  # Will hold a simple solution-like object
    u::uType
    du::duType
    p::pType
    t::tType
    abstol::aTol
    reltol::rTol
    prob::probType
end

function DASKRInitIntegrator(prob, u0, du0, p, t0, abstol, reltol)
    # Create a minimal solution-like object
    sol = (prob=prob, u=[copy(u0)], t=[t0])
    return DASKRInitIntegrator(sol, copy(u0), copy(du0), p, t0, abstol, reltol, prob)
end

# Define the opts interface that get_initial_values expects
struct DASKROpts{aTol, rTol}
    abstol::aTol
    reltol::rTol
end

Base.getproperty(integrator::DASKRInitIntegrator, s::Symbol) = begin
    if s === :opts
        return DASKROpts(getfield(integrator, :abstol), getfield(integrator, :reltol))
    else
        return getfield(integrator, s)
    end
end
