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
        # Same names are non-public in the module they are accessed from: the
        # DiffEqBase/SciMLBase interface names above, the SciMLBase.ReturnCode enum
        # values, and Base.@pure. All become public as those modules declare their
        # public API (Base.@pure is a stable Base internal).
        all_qualified_accesses_are_public = (;
            ignore = (
                Symbol("@pure"),
                :AbstractDAEAlgorithm, :AbstractDAEProblem, :AbstractParameterizedFunction,
                :__solve, :build_solution, :has_jac, :has_tgrad,
                :current_time, :parameter_values, :state_values, :symbolic_container,
                :ConvergenceFailure, :Failure, :InitialFailure, :MaxIters, :Success,
            ),
        ),
        # OverrideInit / get_initial_values are non-public in SciMLBase (the DAE
        # initialization API DASKR implements); public once SciMLBase declares them.
        all_explicit_imports_are_public = (;
            ignore = (:OverrideInit, :get_initial_values),
        ),
    ),
)
