__precompile__()

module DASKR

using Compat
using Libdl

include("core.jl")
include("common.jl")

const warnkeywords =
    (:save_idxs, :d_discontinuities, :isoutofdomain, :unstable_check,
     :calck, :progress, :timeseries_steps, :dtmin,
     :internalnorm, :gamma, :beta1, :beta2, :qmax, :qmin, :qoldinit)

const dllname = joinpath(dirname(dirname(@__FILE__)),
                         "deps",
                         Sys.iswindows() ? "daskr$(Sys.WORD_SIZE)" : "daskr")

function __init__()
    global lib = Libdl.dlopen(dllname)
    global warnlist = Set(warnkeywords)
end
end # module
