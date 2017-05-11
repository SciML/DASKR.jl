__precompile__()

module DASKR

include("core.jl")
include("common.jl")

const dllname = joinpath(dirname(dirname(@__FILE__)),
                         "deps",
                         is_windows() ? "daskr$(Sys.WORD_SIZE)" : "daskr")
function __init__()
    global lib = Libdl.dlopen(dllname)
end

end # module
