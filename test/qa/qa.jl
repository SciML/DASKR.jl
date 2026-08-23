using SciMLTesting, DASKR, Test
using JET

# The SciML common interface DASKR deliberately reexports so that `using DASKR` is enough
# to build a DAE problem, solve it, and inspect the result. Owned and documented upstream;
# kept in sync with the reexport `export` blocks in src/DASKR.jl.
const REEXPORTS = (
    :BrownFullBasicInit, :CheckInit, :DAEFunction, :DAEProblem, :DAESolution, :DEStats,
    :DefaultInit, :EnsembleAnalysis, :EnsembleDistributed, :EnsembleProblem,
    :EnsembleSerial, :EnsembleSolution, :EnsembleSplitThreads, :EnsembleSummary,
    :EnsembleThreads, :NoInit, :NullParameters, :OverrideInit, :ReturnCode,
    :ShampineCollocationInit, :remake, :solve, :successful_retcode,
)

run_qa(DASKR; reexports_allow = REEXPORTS)

@testset "Reexport surface" begin
    # Every approved reexport must actually be reachable from `using DASKR`, so the
    # allow-list cannot drift into approving names the package no longer provides.
    # `isdefined(@__MODULE__, ...)` tests the property directly: this file's
    # `using DASKR` is what has to bring the name into scope.
    @testset "$name" for name in REEXPORTS
        @test name in names(DASKR)
        @test isdefined(@__MODULE__, name)
    end
end
