using DASKR
using Test
using SafeTestsets
using DiffEqBase
using Pkg

const GROUP = get(ENV, "GROUP", "All")

@safetestset "Raw interface" begin
    include("raw_interface.jl")
end

@safetestset "Common interface" begin
    include("common_interface.jl")
end

# NoPre group: JET static analysis tests
# Only run on non-prerelease Julia versions
if GROUP == "NoPre" && isempty(VERSION.prerelease)
    Pkg.activate(joinpath(@__DIR__, "nopre"))
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    @time include("nopre/jet.jl")
end
