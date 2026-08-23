module DASKR

using Compat: Compat
using DASKR_jll: DASKR_jll, libdaskr
using Libdl: Libdl
using SciMLLogging: SciMLLogging, @SciMLMessage
using DiffEqBase: DEVerbosity

# The SciML common interface that DASKR reexports (see the `export` block below), so that
# `using DASKR` on its own is enough to build a DAE problem, solve it, and inspect the
# result -- the workflow the README and the `daskr` docstring document. The DAE
# initialization algorithms are the ones `perform_initialization!` in src/initialize.jl
# has methods for. Every name stays owned and documented upstream.
using SciMLBase: DAEFunction, DAEProblem, DAESolution, DEStats, EnsembleAnalysis,
    EnsembleDistributed, EnsembleProblem, EnsembleSerial, EnsembleSolution,
    EnsembleSplitThreads, EnsembleSummary, EnsembleThreads, NullParameters, ReturnCode,
    remake, solve, successful_retcode

# Reexported SciML common interface; approved via `reexports_allow` in test/qa/qa.jl.
export DAEFunction, DAEProblem, DAESolution, DEStats, EnsembleAnalysis,
    EnsembleDistributed, EnsembleProblem, EnsembleSerial, EnsembleSolution,
    EnsembleSplitThreads, EnsembleSummary, EnsembleThreads, NullParameters, ReturnCode,
    remake, solve, successful_retcode
# DAE initialization algorithms accepted by the `initializealg` keyword; `CheckInit`,
# `NoInit` and `OverrideInit` are imported in src/initialize.jl from SciMLBase, the
# other three from DiffEqBase.
export BrownFullBasicInit, CheckInit, DefaultInit, NoInit, OverrideInit,
    ShampineCollocationInit

_process_verbose_param(v::SciMLLogging.AbstractVerbosityPreset) = DEVerbosity(v)
_process_verbose_param(v::Bool) = v ? DEVerbosity() : DEVerbosity(SciMLLogging.None())
_process_verbose_param(v::DEVerbosity) = v

include("core.jl")
include("initialize.jl")
include("common.jl")

const warnkeywords = (
    :save_idxs, :d_discontinuities, :isoutofdomain, :unstable_check,
    :calck, :progress, :dtmin,
    :internalnorm, :gamma, :beta1, :beta2, :qmax, :qmin, :qoldinit,
)

function __init__()
    global lib = Libdl.dlopen(libdaskr)
    return global warnlist = Set(warnkeywords)
end

# Precompilation workload
include("precompilation.jl")

end # module
