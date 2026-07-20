using SciMLTesting, DASKR, Test
using JET

diffeqbase_reexports = Tuple(names(DASKR.DiffEqBase))

run_qa(
    DASKR;
    api_docs_kwargs = (; rendered = true, rendered_ignore = diffeqbase_reexports),
    explicit_imports = true,
)
