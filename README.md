# DASKR

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Build Status](https://travis-ci.com/SciML/DASKR.jl.svg?branch=master)](https://travis-ci.com/SciML/DASKR.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaDiffEq/DASKR.jl/badge.svg)](https://coveralls.io/github/JuliaDiffEq/DASKR.jl)
[![codecov.io](http://codecov.io/github/SciML/DASKR.jl/coverage.svg?branch=master)](http://codecov.io/github/SciML/DASKR.jl?branch=master)

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

The options for `solve` are documented [at the common solver options page](http://diffeq.sciml.ai/dev/basics/common_solver_opts.html). For more details, see the [ODE Tutorial](http://diffeq.sciml.ai/dev/tutorials/ode_example.html) and the [DAE Tutorial](http://diffeq.sciml.ai/dev/tutorials/dae_example.html) pages from DifferentialEquations.jl.

## Citing

Please cite DifferentialEquations.jl and the original algorithm:

```
@article{rackauckas2017differentialequations,
  title={Differentialequations. jl--a performant and feature-rich ecosystem for solving differential equations in julia},
  author={Rackauckas, Christopher and Nie, Qing},
  journal={Journal of Open Research Software},
  volume={5},
  number={1},
  year={2017},
  publisher={Ubiquity Press}
}

@article{brown1994using,
  title={Using Krylov methods in the solution of large-scale differential-algebraic systems},
  author={Brown, Peter N and Hindmarsh, Alan C and Petzold, Linda R},
  journal={SIAM Journal on Scientific Computing},
  volume={15},
  number={6},
  pages={1467--1488},
  year={1994},
  publisher={SIAM}
}

@article{brown1998consistent,
  title={Consistent initial condition calculation for differential-algebraic systems},
  author={Brown, Peter N and Hindmarsh, Alan C and Petzold, Linda R},
  journal={SIAM Journal on Scientific Computing},
  volume={19},
  number={5},
  pages={1495--1512},
  year={1998},
  publisher={SIAM}
}
```
