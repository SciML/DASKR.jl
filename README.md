# DASKR

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Build Status](https://travis-ci.org/JuliaDiffEq/DASKR.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/DASKR.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaDiffEq/DASKR.jl/badge.svg)](https://coveralls.io/github/JuliaDiffEq/DASKR.jl)
[![codecov.io](http://codecov.io/github/JuliaDiffEq/DASKR.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaDiffEq/DASKR.jl?branch=master)

[![DASKR](http://pkg.julialang.org/badges/DASKR_0.5.svg)](http://pkg.julialang.org/?pkg=DASKR)
[![DASKR](http://pkg.julialang.org/badges/DASKR_0.6.svg)](http://pkg.julialang.org/?pkg=DASKR)

A solver for differential algebraic equations (DAE). This wraps the original DASKR FORTRAN solver. DASKR is a derivative of the DASSL solver with root finding. 

An interface to the JuliaDiffEq common interface is also provided.

## Common Interface Example

```julia
using DASKR
u0 = [1.0, 0, 0]
du0 = [-0.04, 0.04, 0.0]
tspan = (0.0,100000.0)

function resrob(r,yp,y,p,tres)
    r[1]  = -0.04*y[1] + 1.0e4*y[2]*y[3]
    r[2]  = -r[1] - 3.0e7*y[2]*y[2] - yp[2]
    r[1] -=  yp[1]
    r[3]  =  y[1] + y[2] + y[3] - 1.0
end

prob = DAEProblem(resrob,u0,du0,tspan)    
sol = solve(prob, daskr())
```

The options for `solve` are documented [at the common solver options page](http://docs.juliadiffeq.org/latest/basics/common_solver_opts.html). For more details, see the [ODE Tutorial](http://docs.juliadiffeq.org/latest/tutorials/ode_example.html) and the [DAE Tutorial](http://docs.juliadiffeq.org/latest/tutorials/dae_example.html) pages from DifferentialEquations.jl.
