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
        # Owner-aware view: DiffEqBase.__solve / has_tgrad / AbstractParameterizedFunction
        # are owned by SciMLBase but accessed via DiffEqBase (DASKR's declared dep, where
        # they live as the canonical DAE interface), and the state/parameter accessors are
        # owned by SymbolicIndexingInterface but accessed via SciMLBase.
        all_qualified_accesses_via_owners = (;
            ignore = (
                :AbstractParameterizedFunction, :__solve, :has_tgrad,
                :current_time, :parameter_values, :state_values, :symbolic_container,
            ),
        ),
        # Genuinely-non-public names on the registered releases (SciMLBase 3.27.0 /
        # DiffEqBase 7.6.0): the DAE-interface hooks DASKR extends/calls via DiffEqBase
        # (__solve, has_tgrad, AbstractParameterizedFunction -- not public in DiffEqBase)
        # and Base.@pure (a stable Base internal). The SymbolicIndexingInterface accessors
        # are now extended on their owner directly, where they are public.
        all_qualified_accesses_are_public = (;
            ignore = (
                Symbol("@pure"),
                :AbstractParameterizedFunction, :__solve, :has_tgrad,
            ),
        ),
    ),
)
