# API

The DASKR algorithm itself is documented on the [home page](index.md).

## Reexported SciML common interface

`using DASKR` also brings in the parts of the SciML common interface needed to build a
DAE problem, solve it, and inspect the result, so they do not have to be imported
separately. DASKR does not define these names -- they are owned and documented by
[SciMLBase](https://docs.sciml.ai/SciMLBase/stable/) and
[DiffEqBase](https://docs.sciml.ai/DiffEqDocs/stable/), and that is where their
documentation lives:

  - Problems: `DAEProblem`, `EnsembleProblem`
  - Functions: `DAEFunction`
  - Solutions: `DAESolution`, `EnsembleSolution`, `EnsembleSummary`, `DEStats`
  - Ensemble algorithms: `EnsembleSerial`, `EnsembleThreads`, `EnsembleDistributed`,
    `EnsembleSplitThreads`, and the `EnsembleAnalysis` module
  - Solving: `solve`, `remake`
  - Return status: `ReturnCode`, `successful_retcode`
  - Initialization algorithms, passed as the `initializealg` keyword to `solve`:
    `DefaultInit`, `NoInit`, `CheckInit`, `OverrideInit`, `BrownFullBasicInit`,
    `ShampineCollocationInit` (the last two are owned by DiffEqBase)
  - `NullParameters`

Anything else from SciMLBase or DiffEqBase must be imported from SciMLBase or DiffEqBase
directly. Two groups are deliberately absent:

  - **Callbacks.** DASKR errors on a `callback` keyword ("DASKR is not compatible with
    callbacks"), so `ContinuousCallback` and friends are not part of its surface.
  - **The integrator interface** (`init`, `step!`, `solve!`, `reinit!`, ...). DASKR
    implements `SciMLBase.__solve` only; it has no integrator.
