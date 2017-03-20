__precompile__()

module DASKR

include("core.jl")
include("common.jl")

if is_windows()
    const dllname = joinpath(dirname(dirname(@__FILE__)), "deps", "daskr$(Sys.WORD_SIZE).dll")
else
    const dllname = joinpath(dirname(dirname(@__FILE__)), "deps", "daskr.so")
end

function __init__() 
    global lib = Libdl.dlopen(dllname)
end

end # module
