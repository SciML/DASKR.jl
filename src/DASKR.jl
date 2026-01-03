module DASKR

using Compat
using DASKR_jll
using Libdl

include("core.jl")
include("common.jl")

const warnkeywords = (:save_idxs, :d_discontinuities, :isoutofdomain, :unstable_check,
    :calck, :progress, :dtmin,
    :internalnorm, :gamma, :beta1, :beta2, :qmax, :qmin, :qoldinit)

function __init__()
    global lib = Libdl.dlopen(libdaskr)
    global warnlist = Set(warnkeywords)
end

# Precompilation workload
include("precompilation.jl")

end # module
