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
        # SciMLBase-owned solve-interface names reached through DiffEqBase (DASKR's
        # declared dep, not SciMLBase) -- the canonical DiffEqBase DAE interface DASKR
        # subtypes/extends/calls -- plus SymbolicIndexingInterface-owned state/parameter
        # accessors reached through SciMLBase. Both go public as those base libs declare
        # their public API.
        all_qualified_accesses_via_owners = (;
            ignore = (
                :AbstractDAEAlgorithm, :AbstractDAEProblem, :AbstractParameterizedFunction,
                :__solve, :build_solution, :has_jac, :has_tgrad,
                :current_time, :parameter_values, :state_values, :symbolic_container,
            ),
        ),
        # Still non-public on the registered releases (SciMLBase 3.24.0 / DiffEqBase 7.5.7):
        # the DiffEqBase DAE-interface names DASKR subtypes/extends/calls, the
        # SymbolicIndexingInterface accessors reached through SciMLBase, and Base.@pure (a
        # stable Base internal). The SciMLBase.ReturnCode enum values DASKR uses ARE now
        # public, so they no longer need ignoring.
        all_qualified_accesses_are_public = (;
            ignore = (
                Symbol("@pure"),
                :AbstractDAEAlgorithm, :AbstractDAEProblem, :AbstractParameterizedFunction,
                :__solve, :build_solution, :has_jac, :has_tgrad,
                :current_time, :parameter_values, :state_values, :symbolic_container,
            ),
        ),
        # OverrideInit / get_initial_values are still non-public in SciMLBase 3.24.0 (the
        # DAE initialization API DASKR implements).
        all_explicit_imports_are_public = (;
            ignore = (:OverrideInit, :get_initial_values),
        ),
    ),
)
