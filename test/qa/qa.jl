using SciMLTesting, DASKR, Test
using JET

run_qa(
    DASKR;
    explicit_imports = true,
    # no_implicit_imports: DASKR `@reexport using DiffEqBase` and whole-module
    # `using` of Compat/DASKR_jll/Libdl/PrecompileTools/Reexport, so DiffEqBase's
    # (and through it SciMLBase's check_keywords/warn_compat) names and the JLL/Libdl
    # symbols are reached implicitly. The reexport is deliberate (DASKR re-exports the
    # DiffEqBase solve interface) and cannot be made explicit; tracked for a careful
    # follow-up rather than a risky mass refactor.
    # https://github.com/SciML/DASKR.jl/issues/102
    ei_broken = (:no_implicit_imports,),
    ei_kwargs = (;
        # `__solve` is owned by SciMLBase but DASKR extends it as `DiffEqBase.__solve`
        # (the canonical DAE-interface entry point in DASKR's declared dep DiffEqBase;
        # `SciMLBase.__solve === DiffEqBase.__solve`). It is not public in SciMLBase
        # either, so it cannot yet be migrated to the owner -- pending SciMLBase#1411.
        all_qualified_accesses_via_owners = (;
            ignore = (
                :__solve,
            ),
        ),
        # Genuinely-non-public names on the registered releases (SciMLBase 3.28.1 /
        # DiffEqBase 7.6.0): DiffEqBase.__solve (the DAE-interface entry point DASKR
        # extends -- not public in SciMLBase either) and Base.@pure (a stable Base
        # internal).
        all_qualified_accesses_are_public = (;
            ignore = (
                Symbol("@pure"),
                :__solve,
            ),
        ),
    ),
)
