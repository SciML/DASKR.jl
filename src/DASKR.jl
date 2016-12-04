module DASKR

include("core.jl")
include("common.jl")

if is_windows()
    const dllname = Pkg.dir("DASKR", "deps", "daskr$(Sys.WORD_SIZE).dll")
else
    const dllname = Pkg.dir("DASKR", "deps", "daskr.so")
end

function __init__() 
    global lib = Libdl.dlopen(dllname)
end

end # module
