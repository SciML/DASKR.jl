using SciMLTesting, DASKR, Test
using JET

const DIFFEQBASE_FACADE = (:DAEProblem, :solve)

@testset "Intentional facade API" begin
    @test DASKR.DAEProblem === DASKR.DiffEqBase.DAEProblem
    @test DASKR.solve === DASKR.DiffEqBase.solve
end

run_qa(
    DASKR;
    reexports_allow = DIFFEQBASE_FACADE,
    api_docs_kwargs = (; rendered_ignore = DIFFEQBASE_FACADE),
)
