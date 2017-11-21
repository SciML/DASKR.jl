__precompile__()

module DASKR

using Compat

include("core.jl")
include("common.jl")

const warnkeywords =
    (:save_idxs, :d_discontinuities, :isoutofdomain, :unstable_check,
     :calck, :progress, :timeseries_steps, :dtmin,
     :internalnorm, :gamma, :beta1, :beta2, :qmax, :qmin, :qoldinit)

const dllname = joinpath(dirname(dirname(@__FILE__)),
                         "deps",
                         is_windows() ? "daskr$(Sys.WORD_SIZE)" : "daskr")

function __init__()
    global lib = Libdl.dlopen(dllname)
    const global warnlist = Set(warnkeywords)
end
end # module
