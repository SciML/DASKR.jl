__precompile__()

module DASKR

using Compat

include("core.jl")
include("common.jl")

const warnkeywords =
    (:save_idxs, :d_discontinuities, :isoutofdomain, :unstable_check,
     :calck, :progress, :timeseries_steps, :dense, :dtmin, :dtmax,
     :internalnorm, :gamma, :beta1, :beta2, :qmax, :qmin, :qoldinit)

const dllname = joinpath(dirname(dirname(@__FILE__)),
                         "deps",
                         is_windows() ? "daskr$(Sys.WORD_SIZE)" : "daskr")

function __init__()
    global lib = Libdl.dlopen(dllname)
    const global warnlist = Set(warnkeywords)
end

# support function (from DiffEqBase.jl)
export check_keywords
function check_keywords(alg, kwargs, warnlist)
    flg = false
    for (kw, val) in kwargs
        if kw in warnlist
            if val != nothing
                flg = true
                warn(string("The ", kw, " argument is ignored by ", alg, "."))
            end
        end
    end
    flg && warn("Please see http://docs.juliadiffeq.org/latest/basics/compatibility_chart.html")
    nothing
end

end # module
