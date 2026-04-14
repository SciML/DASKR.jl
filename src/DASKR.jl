module DASKR

using Compat
using DASKR_jll
using Libdl
using SciMLLogging: SciMLLogging, Standard, @SciMLMessage
using DiffEqBase: DEVerbosity

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
