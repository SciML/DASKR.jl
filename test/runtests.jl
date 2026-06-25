using SafeTestsets
using Test
using SciMLTesting

run_tests(;
    core = function ()
        @safetestset "Raw interface" begin
            include("raw_interface.jl")
        end
        return @safetestset "Common interface" begin
            include("common_interface.jl")
        end
    end,
    groups = Dict(
        # NoPre runs the JET static analysis in its own environment. The original
        # dispatcher gated this on a non-prerelease Julia (the NoPre CI lane lists
        # only "1"/"lts", so this is belt-and-suspenders); preserve the guard.
        "NoPre_JET" => (;
            env = joinpath(@__DIR__, "NoPre"),
            body = function ()
                return if isempty(VERSION.prerelease)
                    @time include(joinpath(@__DIR__, "NoPre", "jet.jl"))
                end
            end,
        ),
        "MTK_Init" => (;
            env = joinpath(@__DIR__, "MTK"),
            body = function ()
                return @time include(joinpath(@__DIR__, "MTK", "initialization.jl"))
            end,
        ),
    ),
    qa = (;
        env = joinpath(@__DIR__, "qa"),
        body = joinpath(@__DIR__, "qa", "qa.jl"),
    ),
    # The original runtests.jl ran the core interface tests unconditionally (no GROUP
    # guard) and then, for NoPre/MTK, additionally ran the group file in its sub-env.
    # Reproduce that "core (main env) + group (sub-env)" composition per GROUP value.
    umbrellas = Dict(
        "NoPre" => ["Core", "NoPre_JET"],
        "MTK" => ["Core", "MTK_Init"],
    ),
    all = ["Core"],
)
