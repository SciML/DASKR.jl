
"""
Return a C-style callback for the residual function `fun`. Suitable for use with `unsafe_solve`.
"""
function res_c(fun)
    newfun = function (t, y, yp, cj, delta, ires, rpar, ipar)
        n = convert(Array{Int}, unsafe_wrap(Array, ipar, (3,)))
        t = unsafe_wrap(Array, t, (1,))
        y = unsafe_wrap(Array, y, (n[1],))
        yp = unsafe_wrap(Array, yp, (n[1],))
        delta = unsafe_wrap(Array, delta, (n[1],))
        fun(first(t), y, yp, delta)
        return nothing
    end
    @cfunction($newfun, Nothing,
               # T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR
               (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                Ptr{Int32}, Ptr{Float64}, Ptr{Int32}))
end

"""
Return a C-style callback for the event-handling function `fun`. Suitable for use with `unsafe_solve`.
"""
function rt_c(fun)
    newfun = function (neq, t, y, yp, nrt, rval, rpar, ipar)
        n = convert(Array{Int}, unsafe_wrap(Array, ipar, (3,)))
        t = unsafe_wrap(Array, t, (1,))
        y = unsafe_wrap(Array, y, (n[1],))
        yp = unsafe_wrap(Array, yp, (n[1],))
        rval = unsafe_wrap(Array, rval, (n[2],))
        fun(first(t), y, yp, rval)
        return nothing
    end
    @cfunction($newfun, Nothing,
               # T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR
               (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                Ptr{Int32}, Ptr{Float64}, Ptr{Int32}))
end

"""
Return a C-style callback for the Jacobian function `fun`. Suitable for use with `unsafe_solve`.
"""
function jac_c(fun)
    newfun = function (t, y, yp, pd, cj, rpar, ipar)
        n = convert(Array{Int}, unsafe_wrap(Array, ipar, (3,)))
        _t = unsafe_wrap(Array, t, (1,))
        _y = unsafe_wrap(Array, y, (n[1],))
        _yp = unsafe_wrap(Array, yp, (n[1],))
        _pd = unsafe_wrap(Array, pd, (n[3], n[1]))
        _cj = unsafe_wrap(Array, cj, (1,))
        fun(first(_t), _y, _yp, _pd, first(_cj[1]))
        return nothing
    end
    @cfunction($newfun, Nothing,
               # T, Y, YPRIME, PD, CJ, RPAR, IPAR
               (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                Ptr{Float64}, Ptr{Int32}))
end

"""
Return a C-style callback for the residual function `fun`. Suitable for use with `unsafe_solve`.
"""
function common_res_c(fun, p)
    newfun = function (t, y, yp, cj, delta, ires, rpar, ipar)
        n = convert(Array{Int}, unsafe_wrap(Array, ipar, (3,)))
        t = unsafe_wrap(Array, t, (1,))
        y = unsafe_wrap(Array, y, (n[1],))
        yp = unsafe_wrap(Array, yp, (n[1],))
        delta = unsafe_wrap(Array, delta, (n[1],))
        fun(delta, yp, y, p, first(t))
        return nothing
    end
    @cfunction($newfun, Nothing,
               # T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR
               (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                Ptr{Int32}, Ptr{Float64}, Ptr{Int32}))
end

"""
Return a C-style callback for the Jacobian function `fun`. Suitable for use with `unsafe_solve`. For a common interface passed function.
"""
function common_jac_c(fun, p)
    newfun = function (t, y, yp, pd, cj, rpar, ipar)
        n = convert(Array{Int}, unsafe_wrap(Array, ipar, (3,)))
        _t = unsafe_wrap(Array, t, (1,))
        _y = unsafe_wrap(Array, y, (n[1],))
        _yp = unsafe_wrap(Array, yp, (n[1],))
        _pd = unsafe_wrap(Array, pd, (n[3], n[1]))
        _cj = unsafe_wrap(Array, cj, (1,))
        fun.jac(_pd, _yp, _y, p, first(_cj[1]), first(_t))
        return nothing
    end
    @cfunction($newfun, Nothing,
               # T, Y, YPRIME, PD, CJ, RPAR, IPAR
               (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                Ptr{Float64}, Ptr{Int32}))
end

"""
Direct, raw access to the DASKR solver.

```
`unsafe_solve(callback, N, t, y, yp,
              tout, info, rtol, atol,
              idid, rwork, lrw, iwork,
              liw, rpar, ipar, jac, psol,
              rt, nrt, jroot)`
```

The following shows the detailed explanations from the comments in the
FORTRAN source.

```
C      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      INTEGER NEQ, INFO(N), IDID, LRW, LIW, IWORK(LIW), IPAR(*)
C      DOUBLE PRECISION T, Y(*), YPRIME(*), TOUT, RTOL(*), ATOL(*),
C         RWORK(LRW), RPAR(*)
C      EXTERNAL RES, JAC, PSOL, RT
C
C      CALL DDASKR (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
C     *             IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC, PSOL,
C     *             RT, NRT, JROOT)
C
C  Quantities which may be altered by the code are:
C     T, Y(*), YPRIME(*), INFO(1), RTOL, ATOL, IDID, RWORK(*), IWORK(*)
C
C
C *Arguments:
C
C  RES:EXT          This is the name of a subroutine which you
C                   provide to define the residual function G(t,y,y')
C                   of the differential/algebraic system.
C
C  NEQ:IN           This is the number of equations in the system.
C
C  T:INOUT          This is the current value of the independent
C                   variable.
C
C  Y(*):INOUT       This array contains the solution components at T.
C
C  YPRIME(*):INOUT  This array contains the derivatives of the solution
C                   components at T.
C
C  TOUT:IN          This is a point at which a solution is desired.
C
C  INFO(N):IN       This is an integer array used to communicate details
C                   of how the solution is to be carried out, such as
C                   tolerance type, matrix structure, step size and
C                   order limits, and choice of nonlinear system method.
C                   N must be at least 20.
C
C  RTOL,ATOL:INOUT  These quantities represent absolute and relative
C                   error tolerances (on local error) which you provide
C                   to indicate how accurately you wish the solution to
C                   be computed.  You may choose them to be both scalars
C                   or else both arrays of length NEQ.
C
C  IDID:OUT         This integer scalar is an indicator reporting what
C                   the code did.  You must monitor this variable to
C                   decide what action to take next.
C
C  RWORK:WORK       A real work array of length LRW which provides the
C                   code with needed storage space.
C
C  LRW:IN           The length of RWORK.
C
C  IWORK:WORK       An integer work array of length LIW which provides
C                   the code with needed storage space.
C
C  LIW:IN           The length of IWORK.
C
C  RPAR,IPAR:IN     These are real and integer parameter arrays which
C                   you can use for communication between your calling
C                   program and the RES, JAC, and PSOL subroutines.
C
C  JAC:EXT          This is the name of a subroutine which you may
C                   provide (optionally) for calculating Jacobian
C                   (partial derivative) data involved in solving linear
C                   systems within DDASKR.
C
C  PSOL:EXT         This is the name of a subroutine which you must
C                   provide for solving linear systems if you selected
C                   a Krylov method.  The purpose of PSOL is to solve
C                   linear systems involving a left preconditioner P.
C
C  RT:EXT           This is the name of the subroutine for defining
C                   constraint functions Ri(T,Y,Y')) whose roots are
C                   desired during the integration.  This name must be
C                   declared external in the calling program.
C
C  NRT:IN           This is the number of constraint functions
C                   Ri(T,Y,Y').  If there are no constraints, set
C                   NRT = 0, and pass a dummy name for RT.
C
C  JROOT:OUT        This is an integer array of length NRT for output
C                   of root information.
C
C *Overview
C
C  The DDASKR solver uses the backward differentiation formulas of
C  orders one through five to solve a system of the form G(t,y,y') = 0
C  for y = Y and y' = YPRIME.  Values for Y and YPRIME at the initial
C  time must be given as input.  These values should be consistent,
C  that is, if T, Y, YPRIME are the given initial values, they should
C  satisfy G(T,Y,YPRIME) = 0.  However, if consistent values are not
C  known, in many cases you can have DDASKR solve for them -- see
C  INFO(11). (This and other options are described in detail below.)
C
C  Normally, DDASKR solves the system from T to TOUT.  It is easy to
C  continue the solution to get results at additional TOUT.  This is
C  the interval mode of operation.  Intermediate results can also be
C  obtained easily by specifying INFO(3).
C
C  On each step taken by DDASKR, a sequence of nonlinear algebraic
C  systems arises.  These are solved by one of two types of
C  methods:
C    * a Newton iteration with a direct method for the linear
C      systems involved (INFO(12) = 0), or
C    * a Newton iteration with a preconditioned Krylov iterative
C      method for the linear systems involved (INFO(12) = 1).
C
C  The direct method choices are dense and band matrix solvers,
C  with either a user-supplied or an internal difference quotient
C  Jacobian matrix, as specified by INFO(5) and INFO(6).
C  In the band case, INFO(6) = 1, you must supply half-bandwidths
C  in IWORK(1) and IWORK(2).
C
C  The Krylov method is the Generalized Minimum Residual (GMRES)
C  method, in either complete or incomplete form, and with
C  scaling and preconditioning.  The method is implemented
C  in an algorithm called SPIGMR.  Certain options in the Krylov
C  method case are specified by INFO(13) and INFO(15).
C
C  If the Krylov method is chosen, you may supply a pair of routines,
C  JAC and PSOL, to apply preconditioning to the linear system.
C  If the system is A*x = b, the matrix is A = dG/dY + CJ*dG/dYPRIME
C  (of order NEQ).  This system can then be preconditioned in the form
C  (P-inverse)*A*x = (P-inverse)*b, with left preconditioner P.
C  (DDASKR does not allow right preconditioning.)
C  Then the Krylov method is applied to this altered, but equivalent,
C  linear system, hopefully with much better performance than without
C  preconditioning.  (In addition, a diagonal scaling matrix based on
C  the tolerances is also introduced into the altered system.)
C
C  The JAC routine evaluates any data needed for solving systems
C  with coefficient matrix P, and PSOL carries out that solution.
C  In any case, in order to improve convergence, you should try to
C  make P approximate the matrix A as much as possible, while keeping
C  the system P*x = b reasonably easy and inexpensive to solve for x,
C  given a vector b.
C
C  While integrating the given DAE system, DDASKR also searches for
C  roots of the given constraint functions Ri(T,Y,Y') given by RT.
C  If DDASKR detects a sign change in any Ri(T,Y,Y'), it will return
C  the intermediate value of T and Y for which Ri(T,Y,Y') = 0.
C  Caution: If some Ri has a root at or very near the initial time,
C  DDASKR may fail to find it, or may find extraneous roots there,
C  because it does not yet have a sufficient history of the solution.
C
C *Description
C
C------INPUT - WHAT TO DO ON THE FIRST CALL TO DDASKR-------------------
C
C
C  The first call of the code is defined to be the start of each new
C  problem.  Read through the descriptions of all the following items,
C  provide sufficient storage space for designated arrays, set
C  appropriate variables for the initialization of the problem, and
C  give information about how you want the problem to be solved.
C
C
C  RES -- Provide a subroutine of the form
C
C             SUBROUTINE RES (T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR)
C
C         to define the system of differential/algebraic
C         equations which is to be solved. For the given values
C         of T, Y and YPRIME, the subroutine should return
C         the residual of the differential/algebraic system
C             DELTA = G(T,Y,YPRIME)
C         DELTA is a vector of length NEQ which is output from RES.
C
C         Subroutine RES must not alter T, Y, YPRIME, or CJ.
C         You must declare the name RES in an EXTERNAL
C         statement in your program that calls DDASKR.
C         You must dimension Y, YPRIME, and DELTA in RES.
C
C         The input argument CJ can be ignored, or used to rescale
C         constraint equations in the system (see Ref. 2, p. 145).
C         Note: In this respect, DDASKR is not downward-compatible
C         with DDASSL, which does not have the RES argument CJ.
C
C         IRES is an integer flag which is always equal to zero
C         on input.  Subroutine RES should alter IRES only if it
C         encounters an illegal value of Y or a stop condition.
C         Set IRES = -1 if an input value is illegal, and DDASKR
C         will try to solve the problem without getting IRES = -1.
C         If IRES = -2, DDASKR will return control to the calling
C         program with IDID = -11.
C
C         RPAR and IPAR are real and integer parameter arrays which
C         you can use for communication between your calling program
C         and subroutine RES. They are not altered by DDASKR. If you
C         do not need RPAR or IPAR, ignore these parameters by treat-
C         ing them as dummy arguments. If you do choose to use them,
C         dimension them in your calling program and in RES as arrays
C         of appropriate length.
C
C  NEQ -- Set it to the number of equations in the system (NEQ .GE. 1).
C
C  T -- Set it to the initial point of the integration. (T must be
C       a variable.)
C
C  Y(*) -- Set this array to the initial values of the NEQ solution
C          components at the initial point.  You must dimension Y of
C          length at least NEQ in your calling program.
C
C  YPRIME(*) -- Set this array to the initial values of the NEQ first
C               derivatives of the solution components at the initial
C               point.  You must dimension YPRIME at least NEQ in your
C               calling program.
C
C  TOUT - Set it to the first point at which a solution is desired.
C         You cannot take TOUT = T.  Integration either forward in T
C         (TOUT .GT. T) or backward in T (TOUT .LT. T) is permitted.
C
C         The code advances the solution from T to TOUT using step
C         sizes which are automatically selected so as to achieve the
C         desired accuracy.  If you wish, the code will return with the
C         solution and its derivative at intermediate steps (the
C         intermediate-output mode) so that you can monitor them,
C         but you still must provide TOUT in accord with the basic
C         aim of the code.
C
C         The first step taken by the code is a critical one because
C         it must reflect how fast the solution changes near the
C         initial point.  The code automatically selects an initial
C         step size which is practically always suitable for the
C         problem.  By using the fact that the code will not step past
C         TOUT in the first step, you could, if necessary, restrict the
C         length of the initial step.
C
C         For some problems it may not be permissible to integrate
C         past a point TSTOP, because a discontinuity occurs there
C         or the solution or its derivative is not defined beyond
C         TSTOP.  When you have declared a TSTOP point (see INFO(4)
C         and RWORK(1)), you have told the code not to integrate past
C         TSTOP.  In this case any tout beyond TSTOP is invalid input.
C
C  INFO(*) - Use the INFO array to give the code more details about
C            how you want your problem solved.  This array should be
C            dimensioned of length 20, though DDASKR uses only the
C            first 15 entries.  You must respond to all of the following
C            items, which are arranged as questions.  The simplest use
C            of DDASKR corresponds to setting all entries of INFO to 0.
C
C       INFO(1) - This parameter enables the code to initialize itself.
C              You must set it to indicate the start of every new
C              problem.
C
C          **** Is this the first call for this problem ...
C                yes - set INFO(1) = 0
C                 no - not applicable here.
C                      See below for continuation calls.  ****
C
C       INFO(2) - How much accuracy you want of your solution
C              is specified by the error tolerances RTOL and ATOL.
C              The simplest use is to take them both to be scalars.
C              To obtain more flexibility, they can both be arrays.
C              The code must be told your choice.
C
C          **** Are both error tolerances RTOL, ATOL scalars ...
C                yes - set INFO(2) = 0
C                      and input scalars for both RTOL and ATOL
C                 no - set INFO(2) = 1
C                      and input arrays for both RTOL and ATOL ****
C
C       INFO(3) - The code integrates from T in the direction of TOUT
C              by steps.  If you wish, it will return the computed
C              solution and derivative at the next intermediate step
C              (the intermediate-output mode) or TOUT, whichever comes
C              first.  This is a good way to proceed if you want to
C              see the behavior of the solution.  If you must have
C              solutions at a great many specific TOUT points, this
C              code will compute them efficiently.
C
C          **** Do you want the solution only at
C               TOUT (and not at the next intermediate step) ...
C                yes - set INFO(3) = 0 (interval-output mode)
C                 no - set INFO(3) = 1 (intermediate-output mode) ****
C
C       INFO(4) - To handle solutions at a great many specific
C              values TOUT efficiently, this code may integrate past
C              TOUT and interpolate to obtain the result at TOUT.
C              Sometimes it is not possible to integrate beyond some
C              point TSTOP because the equation changes there or it is
C              not defined past TSTOP.  Then you must tell the code
C              this stop condition.
C
C           **** Can the integration be carried out without any
C                restrictions on the independent variable T ...
C                 yes - set INFO(4) = 0
C                  no - set INFO(4) = 1
C                       and define the stopping point TSTOP by
C                       setting RWORK(1) = TSTOP ****
C
C       INFO(5) - used only when INFO(12) = 0 (direct methods).
C              To solve differential/algebraic systems you may wish
C              to use a matrix of partial derivatives of the
C              system of differential equations.  If you do not
C              provide a subroutine to evaluate it analytically (see
C              description of the item JAC in the call list), it will
C              be approximated by numerical differencing in this code.
C              Although it is less trouble for you to have the code
C              compute partial derivatives by numerical differencing,
C              the solution will be more reliable if you provide the
C              derivatives via JAC.  Usually numerical differencing is
C              more costly than evaluating derivatives in JAC, but
C              sometimes it is not - this depends on your problem.
C
C           **** Do you want the code to evaluate the partial deriv-
C                atives automatically by numerical differences ...
C                 yes - set INFO(5) = 0
C                  no - set INFO(5) = 1
C                       and provide subroutine JAC for evaluating the
C                       matrix of partial derivatives ****
C
C       INFO(6) - used only when INFO(12) = 0 (direct methods).
C              DDASKR will perform much better if the matrix of
C              partial derivatives, dG/dY + CJ*dG/dYPRIME (here CJ is
C              a scalar determined by DDASKR), is banded and the code
C              is told this.  In this case, the storage needed will be
C              greatly reduced, numerical differencing will be performed
C              much cheaper, and a number of important algorithms will
C              execute much faster.  The differential equation is said
C              to have half-bandwidths ML (lower) and MU (upper) if
C              equation i involves only unknowns Y(j) with
C                             i-ML .le. j .le. i+MU .
C              For all i=1,2,...,NEQ.  Thus, ML and MU are the widths
C              of the lower and upper parts of the band, respectively,
C              with the main diagonal being excluded.  If you do not
C              indicate that the equation has a banded matrix of partial
C              derivatives the code works with a full matrix of NEQ**2
C              elements (stored in the conventional way).  Computations
C              with banded matrices cost less time and storage than with
C              full matrices if  2*ML+MU .lt. NEQ.  If you tell the
C              code that the matrix of partial derivatives has a banded
C              structure and you want to provide subroutine JAC to
C              compute the partial derivatives, then you must be careful
C              to store the elements of the matrix in the special form
C              indicated in the description of JAC.
C
C          **** Do you want to solve the problem using a full (dense)
C               matrix (and not a special banded structure) ...
C                yes - set INFO(6) = 0
C                 no - set INFO(6) = 1
C                       and provide the lower (ML) and upper (MU)
C                       bandwidths by setting
C                       IWORK(1)=ML
C                       IWORK(2)=MU ****
C
C       INFO(7) - You can specify a maximum (absolute value of)
C              stepsize, so that the code will aNothing passing over very
C              large regions.
C
C          ****  Do you want the code to decide on its own the maximum
C                stepsize ...
C                 yes - set INFO(7) = 0
C                  no - set INFO(7) = 1
C                       and define HMAX by setting
C                       RWORK(2) = HMAX ****
C
C       INFO(8) -  Differential/algebraic problems may occasionally
C              suffer from severe scaling difficulties on the first
C              step.  If you know a great deal about the scaling of
C              your problem, you can help to alleviate this problem
C              by specifying an initial stepsize H0.
C
C          ****  Do you want the code to define its own initial
C                stepsize ...
C                 yes - set INFO(8) = 0
C                  no - set INFO(8) = 1
C                       and define H0 by setting
C                       RWORK(3) = H0 ****
C
C       INFO(9) -  If storage is a severe problem, you can save some
C              storage by restricting the maximum method order MAXORD.
C              The default value is 5.  For each order decrease below 5,
C              the code requires NEQ fewer locations, but it is likely
C              to be slower.  In any case, you must have
C              1 .le. MAXORD .le. 5.
C          ****  Do you want the maximum order to default to 5 ...
C                 yes - set INFO(9) = 0
C                  no - set INFO(9) = 1
C                       and define MAXORD by setting
C                       IWORK(3) = MAXORD ****
C
C       INFO(10) - If you know that certain components of the
C              solutions to your equations are always nonnegative
C              (or nonpositive), it may help to set this
C              parameter.  There are three options that are
C              available:
C              1.  To have constraint checking only in the initial
C                  condition calculation.
C              2.  To enforce nonnegativity in Y during the integration.
C              3.  To enforce both options 1 and 2.
C
C              When selecting option 2 or 3, it is probably best to try
C              the code without using this option first, and only use
C              this option if that does not work very well.
C
C          ****  Do you want the code to solve the problem without
C                invoking any special inequality constraints ...
C                 yes - set INFO(10) = 0
C                  no - set INFO(10) = 1 to have option 1 enforced
C                  no - set INFO(10) = 2 to have option 2 enforced
C                  no - set INFO(10) = 3 to have option 3 enforced ****
C
C                  If you have specified INFO(10) = 1 or 3, then you
C                  will also need to identify how each component of Y
C                  in the initial condition calculation is constrained.
C                  You must set:
C                  IWORK(40+I) = +1 if Y(I) must be .GE. 0,
C                  IWORK(40+I) = +2 if Y(I) must be .GT. 0,
C                  IWORK(40+I) = -1 if Y(I) must be .LE. 0, while
C                  IWORK(40+I) = -2 if Y(I) must be .LT. 0, while
C                  IWORK(40+I) =  0 if Y(I) is not constrained.
C
C       INFO(11) - DDASKR normally requires the initial T, Y, and
C              YPRIME to be consistent.  That is, you must have
C              G(T,Y,YPRIME) = 0 at the initial T.  If you do not know
C              the initial conditions precisely, in some cases
C              DDASKR may be able to compute it.
C
C              Denoting the differential variables in Y by Y_d
C              and the algebraic variables by Y_a, DDASKR can solve
C              one of two initialization problems:
C              1.  Given Y_d, calculate Y_a and Y'_d, or
C              2.  Given Y', calculate Y.
C              In either case, initial values for the given
C              components are input, and initial guesses for
C              the unknown components must also be provided as input.
C
C          ****  Are the initial T, Y, YPRIME consistent ...
C
C                 yes - set INFO(11) = 0
C                  no - set INFO(11) = 1 to calculate option 1 above,
C                    or set INFO(11) = 2 to calculate option 2 ****
C
C                  If you have specified INFO(11) = 1, then you
C                  will also need to identify  which are the
C                  differential and which are the algebraic
C                  components (algebraic components are components
C                  whose derivatives do not appear explicitly
C                  in the function G(T,Y,YPRIME)).  You must set:
C                  IWORK(LID+I) = +1 if Y(I) is a differential variable
C                  IWORK(LID+I) = -1 if Y(I) is an algebraic variable,
C                  where LID = 40 if INFO(10) = 0 or 2 and LID = 40+NEQ
C                  if INFO(10) = 1 or 3.
C
C       INFO(12) - Except for the addition of the RES argument CJ,
C              DDASKR by default is downward-compatible with DDASSL,
C              which uses only direct (dense or band) methods to solve
C              the linear systems involved.  You must set INFO(12) to
C              indicate whether you want the direct methods or the
C              Krylov iterative method.
C          ****   Do you want DDASKR to use standard direct methods
C                 (dense or band) or the Krylov (iterative) method ...
C                   direct methods - set INFO(12) = 0.
C                   Krylov method  - set INFO(12) = 1,
C                       and check the settings of INFO(13) and INFO(15).
C
C       INFO(13) - used when INFO(12) = 1 (Krylov methods).
C              DDASKR uses scalars MAXL, KMP, NRMAX, and EPLI for the
C              iterative solution of linear systems.  INFO(13) allows
C              you to override the default values of these parameters.
C              These parameters and their defaults are as follows:
C              MAXL = maximum number of iterations in the SPIGMR
C                 algorithm (MAXL .le. NEQ).  The default is
C                 MAXL = MIN(5,NEQ).
C              KMP = number of vectors on which orthogonalization is
C                 done in the SPIGMR algorithm.  The default is
C                 KMP = MAXL, which corresponds to complete GMRES
C                 iteration, as opposed to the incomplete form.
C              NRMAX = maximum number of restarts of the SPIGMR
C                 algorithm per nonlinear iteration.  The default is
C                 NRMAX = 5.
C              EPLI = convergence test constant in SPIGMR algorithm.
C                 The default is EPLI = 0.05.
C              Note that the length of RWORK depends on both MAXL
C              and KMP.  See the definition of LRW below.
C          ****   Are MAXL, KMP, and EPLI to be given their
C                 default values ...
C                  yes - set INFO(13) = 0
C                   no - set INFO(13) = 1,
C                        and set all of the following:
C                        IWORK(24) = MAXL (1 .le. MAXL .le. NEQ)
C                        IWORK(25) = KMP  (1 .le. KMP .le. MAXL)
C                        IWORK(26) = NRMAX  (NRMAX .ge. 0)
C                        RWORK(10) = EPLI (0 .lt. EPLI .lt. 1.0) ****
C
C        INFO(14) - used with INFO(11) > 0 (initial condition
C               calculation is requested).  In this case, you may
C               request control to be returned to the calling program
C               immediately after the initial condition calculation,
C               before proceeding to the integration of the system
C               (e.g. to examine the computed Y and YPRIME).
C               If this is done, and if the initialization succeeded
C               (IDID = 4), you should reset INFO(11) to 0 for the
C               next call, to prevent the solver from repeating the
C               initialization (and to aNothing an infinite loop).
C          ****   Do you want to proceed to the integration after
C                 the initial condition calculation is done ...
C                 yes - set INFO(14) = 0
C                  no - set INFO(14) = 1                        ****
C
C        INFO(15) - used when INFO(12) = 1 (Krylov methods).
C               When using preconditioning in the Krylov method,
C               you must supply a subroutine, PSOL, which solves the
C               associated linear systems using P.
C               The usage of DDASKR is simpler if PSOL can carry out
C               the solution without any prior calculation of data.
C               However, if some partial derivative data is to be
C               calculated in advance and used repeatedly in PSOL,
C               then you must supply a JAC routine to do this,
C               and set INFO(15) to indicate that JAC is to be called
C               for this purpose.  For example, P might be an
C               approximation to a part of the matrix A which can be
C               calculated and LU-factored for repeated solutions of
C               the preconditioner system.  The arrays WP and IWP
C               (described under JAC and PSOL) can be used to
C               communicate data between JAC and PSOL.
C          ****   Does PSOL operate with no prior preparation ...
C                 yes - set INFO(15) = 0 (no JAC routine)
C                  no - set INFO(15) = 1
C                       and supply a JAC routine to evaluate and
C                       preprocess any required Jacobian data.  ****
C
C         INFO(16) - option to exclude algebraic variables from
C               the error test.
C          ****   Do you wish to control errors locally on
C                 all the variables...
C                 yes - set INFO(16) = 0
C                  no - set INFO(16) = 1
C                       If you have specified INFO(16) = 1, then you
C                       will also need to identify  which are the
C                       differential and which are the algebraic
C                       components (algebraic components are components
C                       whose derivatives do not appear explicitly
C                       in the function G(T,Y,YPRIME)).  You must set:
C                       IWORK(LID+I) = +1 if Y(I) is a differential
C                                      variable, and
C                       IWORK(LID+I) = -1 if Y(I) is an algebraic
C                                      variable,
C                       where LID = 40 if INFO(10) = 0 or 2 and
C                       LID = 40 + NEQ if INFO(10) = 1 or 3.
C
C       INFO(17) - used when INFO(11) > 0 (DDASKR is to do an
C              initial condition calculation).
C              DDASKR uses several heuristic control quantities in the
C              initial condition calculation.  They have default values,
C              but can  also be set by the user using INFO(17).
C              These parameters and their defaults are as follows:
C              MXNIT  = maximum number of Newton iterations
C                 per Jacobian or preconditioner evaluation.
C                 The default is:
C                 MXNIT =  5 in the direct case (INFO(12) = 0), and
C                 MXNIT = 15 in the Krylov case (INFO(12) = 1).
C              MXNJ   = maximum number of Jacobian or preconditioner
C                 evaluations.  The default is:
C                 MXNJ = 6 in the direct case (INFO(12) = 0), and
C                 MXNJ = 2 in the Krylov case (INFO(12) = 1).
C              MXNH   = maximum number of values of the artificial
C                 stepsize parameter H to be tried if INFO(11) = 1.
C                 The default is MXNH = 5.
C                 NOTE: the maximum number of Newton iterations
C                 allowed in all is MXNIT*MXNJ*MXNH if INFO(11) = 1,
C                 and MXNIT*MXNJ if INFO(11) = 2.
C              LSOFF  = flag to turn off the linesearch algorithm
C                 (LSOFF = 0 means linesearch is on, LSOFF = 1 means
C                 it is turned off).  The default is LSOFF = 0.
C              STPTOL = minimum scaled step in linesearch algorithm.
C                 The default is STPTOL = (unit roundoff)**(2/3).
C              EPINIT = swing factor in the Newton iteration convergence
C                 test.  The test is applied to the residual vector,
C                 premultiplied by the approximate Jacobian (in the
C                 direct case) or the preconditioner (in the Krylov
C                 case).  For convergence, the weighted RMS norm of
C                 this vector (scaled by the error weights) must be
C                 less than EPINIT*EPCON, where EPCON = .33 is the
C                 analogous test constant used in the time steps.
C                 The default is EPINIT = .01.
C          ****   Are the initial condition heuristic controls to be
C                 given their default values...
C                  yes - set INFO(17) = 0
C                   no - set INFO(17) = 1,
C                        and set all of the following:
C                        IWORK(32) = MXNIT (.GT. 0)
C                        IWORK(33) = MXNJ (.GT. 0)
C                        IWORK(34) = MXNH (.GT. 0)
C                        IWORK(35) = LSOFF ( = 0 or 1)
C                        RWORK(14) = STPTOL (.GT. 0.0)
C                        RWORK(15) = EPINIT (.GT. 0.0)  ****
C
C         INFO(18) - option to get extra printing in initial condition
C                calculation.
C          ****   Do you wish to have extra printing...
C                 no  - set INFO(18) = 0
C                 yes - set INFO(18) = 1 for minimal printing, or
C                       set INFO(18) = 2 for full printing.
C                       If you have specified INFO(18) .ge. 1, data
C                       will be printed with the error handler routines.
C                       To print to a non-default unit number L, include
C                       the line  CALL XSETUN(L)  in your program.  ****
C
C   RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL)
C               error tolerances to tell the code how accurately you
C               want the solution to be computed.  They must be defined
C               as variables because the code may change them.
C               you have two choices --
C                     Both RTOL and ATOL are scalars (INFO(2) = 0), or
C                     both RTOL and ATOL are vectors (INFO(2) = 1).
C               In either case all components must be non-negative.
C
C               The tolerances are used by the code in a local error
C               test at each step which requires roughly that
C                        abs(local error in Y(i)) .le. EWT(i) ,
C               where EWT(i) = RTOL*abs(Y(i)) + ATOL is an error weight
C               quantity, for each vector component.
C               (More specifically, a root-mean-square norm is used to
C               measure the size of vectors, and the error test uses the
C               magnitude of the solution at the beginning of the step.)
C
C               The true (global) error is the difference between the
C               true solution of the initial value problem and the
C               computed approximation.  Practically all present day
C               codes, including this one, control the local error at
C               each step and do not even attempt to control the global
C               error directly.
C
C               Usually, but not always, the true accuracy of
C               the computed Y is comparable to the error tolerances.
C               This code will usually, but not always, deliver a more
C               accurate solution if you reduce the tolerances and
C               integrate again.  By comparing two such solutions you
C               can get a fairly reliable idea of the true error in the
C               solution at the larger tolerances.
C
C               Setting ATOL = 0. results in a pure relative error test
C               on that component.  Setting RTOL = 0. results in a pure
C               absolute error test on that component.  A mixed test
C               with non-zero RTOL and ATOL corresponds roughly to a
C               relative error test when the solution component is
C               much bigger than ATOL and to an absolute error test
C               when the solution component is smaller than the
C               threshold ATOL.
C
C               The code will not attempt to compute a solution at an
C               accuracy unreasonable for the machine being used.  It
C               will advise you if you ask for too much accuracy and
C               inform you as to the maximum accuracy it believes
C               possible.
C
C  RWORK(*) -- a real work array, which should be dimensioned in your
C               calling program with a length equal to the value of
C               LRW (or greater).
C
C  LRW -- Set it to the declared length of the RWORK array.  The
C               minimum length depends on the options you have selected,
C               given by a base value plus additional storage as
C               described below.
C
C               If INFO(12) = 0 (standard direct method), the base value
C               is BASE = 60 + max(MAXORD+4,7)*NEQ + 3*NRT.
C               The default value is MAXORD = 5 (see INFO(9)).  With the
C               default MAXORD, BASE = 60 + 9*NEQ + 3*NRT.
C               Additional storage must be added to the base value for
C               any or all of the following options:
C                 If INFO(6) = 0 (dense matrix), add NEQ**2.
C                 If INFO(6) = 1 (banded matrix), then:
C                    if INFO(5) = 0, add (2*ML+MU+1)*NEQ
C                                           + 2*[NEQ/(ML+MU+1) + 1], and
C                    if INFO(5) = 1, add (2*ML+MU+1)*NEQ.
C                 If INFO(16) = 1, add NEQ.
C
C               If INFO(12) = 1 (Krylov method), the base value is
C               BASE = 60 + (MAXORD+5)*NEQ + 3*NRT
C                         + [MAXL + 3 + min(1,MAXL-KMP)]*NEQ
C                         + (MAXL+3)*MAXL + 1 + LENWP.
C               See PSOL for description of LENWP.  The default values
C               are: MAXORD = 5 (see INFO(9)), MAXL = min(5,NEQ) and
C               KMP = MAXL  (see INFO(13)).  With these default values,
C               BASE = 101 + 18*NEQ + 3*NRT + LENWP.
C               Additional storage must be added to the base value for
C               the following option:
C                 If INFO(16) = 1, add NEQ.
C
C
C  IWORK(*) -- an integer work array, which should be dimensioned in
C              your calling program with a length equal to the value
C              of LIW (or greater).
C
C  LIW -- Set it to the declared length of the IWORK array.  The
C             minimum length depends on the options you have selected,
C             given by a base value plus additions as described below.
C
C             If INFO(12) = 0 (standard direct method), the base value
C             is BASE = 40 + NEQ.
C             IF INFO(10) = 1 or 3, add NEQ to the base value.
C             If INFO(11) = 1 or INFO(16) =1, add NEQ to the base value.
C
C             If INFO(12) = 1 (Krylov method), the base value is
C             BASE = 40 + LENIWP.  See PSOL for description of LENIWP.
C             If INFO(10) = 1 or 3, add NEQ to the base value.
C             If INFO(11) = 1 or INFO(16) =1, add NEQ to the base value.
C
C
C  RPAR, IPAR -- These are arrays of double precision and integer type,
C             respectively, which are available for you to use
C             for communication between your program that calls
C             DDASKR and the RES subroutine (and the JAC and PSOL
C             subroutines).  They are not altered by DDASKR.
C             If you do not need RPAR or IPAR, ignore these
C             parameters by treating them as dummy arguments.
C             If you do choose to use them, dimension them in
C             your calling program and in RES (and in JAC and PSOL)
C             as arrays of appropriate length.
C
C  JAC -- This is the name of a routine that you may supply
C         (optionally) that relates to the Jacobian matrix of the
C         nonlinear system that the code must solve at each T step.
C         The role of JAC (and its call sequence) depends on whether
C         a direct (INFO(12) = 0) or Krylov (INFO(12) = 1) method
C         is selected.
C
C         **** INFO(12) = 0 (direct methods):
C           If you are letting the code generate partial derivatives
C           numerically (INFO(5) = 0), then JAC can be absent
C           (or perhaps a dummy routine to satisfy the loader).
C           Otherwise you must supply a JAC routine to compute
C           the matrix A = dG/dY + CJ*dG/dYPRIME.  It must have
C           the form
C
C           SUBROUTINE JAC (T, Y, YPRIME, PD, CJ, RPAR, IPAR)
C
C           The JAC routine must dimension Y, YPRIME, and PD (and RPAR
C           and IPAR if used).  CJ is a scalar which is input to JAC.
C           For the given values of T, Y, and YPRIME, the JAC routine
C           must evaluate the nonzero elements of the matrix A, and
C           store these values in the array PD.  The elements of PD are
C           set to zero before each call to JAC, so that only nonzero
C           elements need to be defined.
C           The way you store the elements into the PD array depends
C           on the structure of the matrix indicated by INFO(6).
C           *** INFO(6) = 0 (full or dense matrix) ***
C               Give PD a first dimension of NEQ.  When you evaluate the
C               nonzero partial derivatives of equation i (i.e. of G(i))
C               with respect to component j (of Y and YPRIME), you must
C               store the element in PD according to
C                  PD(i,j) = dG(i)/dY(j) + CJ*dG(i)/dYPRIME(j).
C           *** INFO(6) = 1 (banded matrix with half-bandwidths ML, MU
C                            as described under INFO(6)) ***
C               Give PD a first dimension of 2*ML+MU+1.  When you
C               evaluate the nonzero partial derivatives of equation i
C               (i.e. of G(i)) with respect to component j (of Y and
C               YPRIME), you must store the element in PD according to
C                  IROW = i - j + ML + MU + 1
C                  PD(IROW,j) = dG(i)/dY(j) + CJ*dG(i)/dYPRIME(j).
C
C          **** INFO(12) = 1 (Krylov method):
C            If you are not calculating Jacobian data in advance for use
C            in PSOL (INFO(15) = 0), JAC can be absent (or perhaps a
C            dummy routine to satisfy the loader).  Otherwise, you may
C            supply a JAC routine to compute and preprocess any parts of
C            of the Jacobian matrix  A = dG/dY + CJ*dG/dYPRIME that are
C            involved in the preconditioner matrix P.
C            It is to have the form
C
C            SUBROUTINE JAC (RES, IRES, NEQ, T, Y, YPRIME, REWT, SAVR,
C                            WK, H, CJ, WP, IWP, IER, RPAR, IPAR)
C
C           The JAC routine must dimension Y, YPRIME, REWT, SAVR, WK,
C           and (if used) WP, IWP, RPAR, and IPAR.
C           The Y, YPRIME, and SAVR arrays contain the current values
C           of Y, YPRIME, and the residual G, respectively.
C           The array WK is work space of length NEQ.
C           H is the step size.  CJ is a scalar, input to JAC, that is
C           normally proportional to 1/H.  REWT is an array of
C           reciprocal error weights, 1/EWT(i), where EWT(i) is
C           RTOL*abs(Y(i)) + ATOL (unless you supplied routine DDAWTS
C           instead), for use in JAC if needed.  For example, if JAC
C           computes difference quotient approximations to partial
C           derivatives, the REWT array may be useful in setting the
C           increments used.  The JAC routine should do any
C           factorization operations called for, in preparation for
C           solving linear systems in PSOL.  The matrix P should
C           be an approximation to the Jacobian,
C           A = dG/dY + CJ*dG/dYPRIME.
C
C           WP and IWP are real and integer work arrays which you may
C           use for communication between your JAC routine and your
C           PSOL routine.  These may be used to store elements of the
C           preconditioner P, or related matrix data (such as factored
C           forms).  They are not altered by DDASKR.
C           If you do not need WP or IWP, ignore these parameters by
C           treating them as dummy arguments.  If you do use them,
C           dimension them appropriately in your JAC and PSOL routines.
C           See the PSOL description for instructions on setting
C           the lengths of WP and IWP.
C
C           On return, JAC should set the error flag IER as follows..
C             IER = 0    if JAC was successful,
C             IER .ne. 0 if JAC was unsuccessful (e.g. if Y or YPRIME
C                        was illegal, or a singular matrix is found).
C           (If IER .ne. 0, a smaller stepsize will be tried.)
C           IER = 0 on entry to JAC, so need be reset only on a failure.
C           If RES is used within JAC, then a nonzero value of IRES will
C           override any nonzero value of IER (see the RES description).
C
C         Regardless of the method type, subroutine JAC must not
C         alter T, Y(*), YPRIME(*), H, CJ, or REWT(*).
C         You must declare the name JAC in an EXTERNAL statement in
C         your program that calls DDASKR.
C
C PSOL --  This is the name of a routine you must supply if you have
C         selected a Krylov method (INFO(12) = 1) with preconditioning.
C         In the direct case (INFO(12) = 0), PSOL can be absent
C         (a dummy routine may have to be supplied to satisfy the
C         loader).  Otherwise, you must provide a PSOL routine to
C         solve linear systems arising from preconditioning.
C         When supplied with INFO(12) = 1, the PSOL routine is to
C         have the form
C
C         SUBROUTINE PSOL (NEQ, T, Y, YPRIME, SAVR, WK, CJ, WGHT,
C                          WP, IWP, B, EPLIN, IER, RPAR, IPAR)
C
C         The PSOL routine must solve linear systems of the form
C         P*x = b where P is the left preconditioner matrix.
C
C         The right-hand side vector b is in the B array on input, and
C         PSOL must return the solution vector x in B.
C         The Y, YPRIME, and SAVR arrays contain the current values
C         of Y, YPRIME, and the residual G, respectively.
C
C         Work space required by JAC and/or PSOL, and space for data to
C         be communicated from JAC to PSOL is made available in the form
C         of arrays WP and IWP, which are parts of the RWORK and IWORK
C         arrays, respectively.  The lengths of these real and integer
C         work spaces WP and IWP must be supplied in LENWP and LENIWP,
C         respectively, as follows..
C           IWORK(27) = LENWP = length of real work space WP
C           IWORK(28) = LENIWP = length of integer work space IWP.
C
C         WK is a work array of length NEQ for use by PSOL.
C         CJ is a scalar, input to PSOL, that is normally proportional
C         to 1/H (H = stepsize).  If the old value of CJ
C         (at the time of the last JAC call) is needed, it must have
C         been saved by JAC in WP.
C
C         WGHT is an array of weights, to be used if PSOL uses an
C         iterative method and performs a convergence test.  (In terms
C         of the argument REWT to JAC, WGHT is REWT/sqrt(NEQ).)
C         If PSOL uses an iterative method, it should use EPLIN
C         (a heuristic parameter) as the bound on the weighted norm of
C         the residual for the computed solution.  Specifically, the
C         residual vector R should satisfy
C              SQRT (SUM ( (R(i)*WGHT(i))**2 ) ) .le. EPLIN
C
C         PSOL must not alter NEQ, T, Y, YPRIME, SAVR, CJ, WGHT, EPLIN.
C
C         On return, PSOL should set the error flag IER as follows..
C           IER = 0 if PSOL was successful,
C           IER .lt. 0 if an unrecoverable error occurred, meaning
C                 control will be passed to the calling routine,
C           IER .gt. 0 if a recoverable error occurred, meaning that
C                 the step will be retried with the same step size
C                 but with a call to JAC to update necessary data,
C                 unless the Jacobian data is current, in which case
C                 the step will be retried with a smaller step size.
C           IER = 0 on entry to PSOL so need be reset only on a failure.
C
C         You must declare the name PSOL in an EXTERNAL statement in
C         your program that calls DDASKR.
C
C RT --   This is the name of the subroutine for defining the vector
C         R(T,Y,Y') of constraint functions Ri(T,Y,Y'), whose roots
C         are desired during the integration.  It is to have the form
C             SUBROUTINE RT(NEQ, T, Y, YP, NRT, RVAL, RPAR, IPAR)
C             DIMENSION Y(NEQ), YP(NEQ), RVAL(NRT),
C         where NEQ, T, Y and NRT are INPUT, and the array RVAL is
C         output.  NEQ, T, Y, and YP have the same meaning as in the
C         RES routine, and RVAL is an array of length NRT.
C         For i = 1,...,NRT, this routine is to load into RVAL(i) the
C         value at (T,Y,Y') of the i-th constraint function Ri(T,Y,Y').
C         DDASKR will find roots of the Ri of odd multiplicity
C         (that is, sign changes) as they occur during the integration.
C         RT must be declared EXTERNAL in the calling program.
C
C         CAUTION.. Because of numerical errors in the functions Ri
C         due to roundoff and integration error, DDASKR may return
C         false roots, or return the same root at two or more nearly
C         equal values of T.  If such false roots are suspected,
C         the user should consider smaller error tolerances and/or
C         higher precision in the evaluation of the Ri.
C
C         If a root of some Ri defines the end of the problem,
C         the input to DDASKR should nevertheless allow
C         integration to a point slightly past that root, so
C         that DDASKR can locate the root by interpolation.
C
C NRT --  The number of constraint functions Ri(T,Y,Y').  If there are
C         no constraints, set NRT = 0 and pass a dummy name for RT.
C
C JROOT -- This is an integer array of length NRT, used only for output.
C         On a return where one or more roots were found (IDID = 5),
C         JROOT(i) = 1 or -1 if Ri(T,Y,Y') has a root at T, and
C         JROOT(i) = 0 if not.  If nonzero, JROOT(i) shows the direction
C         of the sign change in Ri in the direction of integration:
C         JROOT(i) = 1  means Ri changed from negative to positive.
C         JROOT(i) = -1 means Ri changed from positive to negative.
C
C
C  OPTIONALLY REPLACEABLE SUBROUTINE:
C
C  DDASKR uses a weighted root-mean-square norm to measure the
C  size of various error vectors.  The weights used in this norm
C  are set in the following subroutine:
C
C    SUBROUTINE DDAWTS (NEQ, IWT, RTOL, ATOL, Y, EWT, RPAR, IPAR)
C    DIMENSION RTOL(*), ATOL(*), Y(*), EWT(*), RPAR(*), IPAR(*)
C
C  A DDAWTS routine has been included with DDASKR which sets the
C  weights according to
C    EWT(I) = RTOL*ABS(Y(I)) + ATOL
C  in the case of scalar tolerances (IWT = 0) or
C    EWT(I) = RTOL(I)*ABS(Y(I)) + ATOL(I)
C  in the case of array tolerances (IWT = 1).  (IWT is INFO(2).)
C  In some special cases, it may be appropriate for you to define
C  your own error weights by writing a subroutine DDAWTS to be
C  called instead of the version supplied.  However, this should
C  be attempted only after careful thought and consideration.
C  If you supply this routine, you may use the tolerances and Y
C  as appropriate, but do not overwrite these variables.  You
C  may also use RPAR and IPAR to communicate data as appropriate.
C  ***Note: Aside from the values of the weights, the choice of
C  norm used in DDASKR (weighted root-mean-square) is not subject
C  to replacement by the user.  In this respect, DDASKR is not
C  downward-compatible with the original DDASSL solver (in which
C  the norm routine was optionally user-replaceable).
C
C
C------OUTPUT - AFTER ANY RETURN FROM DDASKR----------------------------
C
C  The principal aim of the code is to return a computed solution at
C  T = TOUT, although it is also possible to obtain intermediate
C  results along the way.  To find out whether the code achieved its
C  goal or if the integration process was interrupted before the task
C  was completed, you must check the IDID parameter.
C
C
C   T -- The output value of T is the point to which the solution
C        was successfully advanced.
C
C   Y(*) -- contains the computed solution approximation at T.
C
C   YPRIME(*) -- contains the computed derivative approximation at T.
C
C   IDID -- reports what the code did, described as follows:
C
C                     *** TASK COMPLETED ***
C                Reported by positive values of IDID
C
C           IDID = 1 -- A step was successfully taken in the
C                   interval-output mode.  The code has not
C                   yet reached TOUT.
C
C           IDID = 2 -- The integration to TSTOP was successfully
C                   completed (T = TSTOP) by stepping exactly to TSTOP.
C
C           IDID = 3 -- The integration to TOUT was successfully
C                   completed (T = TOUT) by stepping past TOUT.
C                   Y(*) and YPRIME(*) are obtained by interpolation.
C
C           IDID = 4 -- The initial condition calculation, with
C                   INFO(11) > 0, was successful, and INFO(14) = 1.
C                   No integration steps were taken, and the solution
C                   is not considered to have been started.
C
C           IDID = 5 -- The integration was successfully completed
C                   by finding one or more roots of R(T,Y,Y') at T.
C
C                    *** TASK INTERRUPTED ***
C                Reported by negative values of IDID
C
C           IDID = -1 -- A large amount of work has been expended
C                     (about 500 steps).
C
C           IDID = -2 -- The error tolerances are too stringent.
C
C           IDID = -3 -- The local error test cannot be satisfied
C                     because you specified a zero component in ATOL
C                     and the corresponding computed solution component
C                     is zero.  Thus, a pure relative error test is
C                     impossible for this component.
C
C           IDID = -5 -- There were repeated failures in the evaluation
C                     or processing of the preconditioner (in JAC).
C
C           IDID = -6 -- DDASKR had repeated error test failures on the
C                     last attempted step.
C
C           IDID = -7 -- The nonlinear system solver in the time
C                     integration could not converge.
C
C           IDID = -8 -- The matrix of partial derivatives appears
C                     to be singular (direct method).
C
C           IDID = -9 -- The nonlinear system solver in the integration
C                     failed to achieve convergence, and there were
C                     repeated  error test failures in this step.
C
C           IDID =-10 -- The nonlinear system solver in the integration
C                     failed to achieve convergence because IRES was
C                     equal  to -1.
C
C           IDID =-11 -- IRES = -2 was encountered and control is
C                     being returned to the calling program.
C
C           IDID =-12 -- DDASKR failed to compute the initial Y, YPRIME.
C
C           IDID =-13 -- An unrecoverable error was encountered inside
C                     the user's PSOL routine, and control is being
C                     returned to the calling program.
C
C           IDID =-14 -- The Krylov linear system solver could not
C                     achieve convergence.
C
C           IDID =-15,..,-32 -- Not applicable for this code.
C
C                    *** TASK TERMINATED ***
C                reported by the value of IDID=-33
C
C           IDID = -33 -- The code has encountered trouble from which
C                   it cannot recover.  A message is printed
C                   explaining the trouble and control is returned
C                   to the calling program.  For example, this occurs
C                   when invalid input is detected.
C
C   RTOL, ATOL -- these quantities remain unchanged except when
C               IDID = -2.  In this case, the error tolerances have been
C               increased by the code to values which are estimated to
C               be appropriate for continuing the integration.  However,
C               the reported solution at T was obtained using the input
C               values of RTOL and ATOL.
C
C   RWORK, IWORK -- contain information which is usually of no interest
C               to the user but necessary for subsequent calls.
C               However, you may be interested in the performance data
C               listed below.  These quantities are accessed in RWORK
C               and IWORK but have internal mnemonic names, as follows..
C
C               RWORK(3)--contains H, the step size h to be attempted
C                        on the next step.
C
C               RWORK(4)--contains TN, the current value of the
C                        independent variable, i.e. the farthest point
C                        integration has reached.  This will differ
C                        from T if interpolation has been performed
C                        (IDID = 3).
C
C               RWORK(7)--contains HOLD, the stepsize used on the last
C                        successful step.  If INFO(11) = INFO(14) = 1,
C                        this contains the value of H used in the
C                        initial condition calculation.
C
C               IWORK(7)--contains K, the order of the method to be
C                        attempted on the next step.
C
C               IWORK(8)--contains KOLD, the order of the method used
C                        on the last step.
C
C               IWORK(11)--contains NST, the number of steps (in T)
C                        taken so far.
C
C               IWORK(12)--contains NRE, the number of calls to RES
C                        so far.
C
C               IWORK(13)--contains NJE, the number of calls to JAC so
C                        far (Jacobian or preconditioner evaluations).
C
C               IWORK(14)--contains NETF, the total number of error test
C                        failures so far.
C
C               IWORK(15)--contains NCFN, the total number of nonlinear
C                        convergence failures so far (includes counts
C                        of singular iteration matrix or singular
C                        preconditioners).
C
C               IWORK(16)--contains NCFL, the number of convergence
C                        failures of the linear iteration so far.
C
C               IWORK(17)--contains LENIW, the length of IWORK actually
C                        required.  This is defined on normal returns
C                        and on an illegal input return for
C                        insufficient storage.
C
C               IWORK(18)--contains LENRW, the length of RWORK actually
C                        required.  This is defined on normal returns
C                        and on an illegal input return for
C                        insufficient storage.
C
C               IWORK(19)--contains NNI, the total number of nonlinear
C                        iterations so far (each of which calls a
C                        linear solver).
C
C               IWORK(20)--contains NLI, the total number of linear
C                        (Krylov) iterations so far.
C
C               IWORK(21)--contains NPS, the number of PSOL calls so
C                        far, for preconditioning solve operations or
C                        for solutions with the user-supplied method.
C
C               IWORK(36)--contains the total number of calls to the
C                        constraint function routine RT so far.
C
C               Note: The various counters in IWORK do not include
C               counts during a prior call made with INFO(11) > 0 and
C               INFO(14) = 1.
C
C
C------INPUT - WHAT TO DO TO CONTINUE THE INTEGRATION  -----------------
C              (CALLS AFTER THE FIRST)
C
C     This code is organized so that subsequent calls to continue the
C     integration involve little (if any) additional effort on your
C     part.  You must monitor the IDID parameter in order to determine
C     what to do next.
C
C     Recalling that the principal task of the code is to integrate
C     from T to TOUT (the interval mode), usually all you will need
C     to do is specify a new TOUT upon reaching the current TOUT.
C
C     Do not alter any quantity not specifically permitted below.  In
C     particular do not alter NEQ, T, Y(*), YPRIME(*), RWORK(*),
C     IWORK(*), or the differential equation in subroutine RES.  Any
C     such alteration constitutes a new problem and must be treated
C     as such, i.e. you must start afresh.
C
C     You cannot change from array to scalar error control or vice
C     versa (INFO(2)), but you can change the size of the entries of
C     RTOL or ATOL.  Increasing a tolerance makes the equation easier
C     to integrate.  Decreasing a tolerance will make the equation
C     harder to integrate and should generally be aNothinged.
C
C     You can switch from the intermediate-output mode to the
C     interval mode (INFO(3)) or vice versa at any time.
C
C     If it has been necessary to prevent the integration from going
C     past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
C     code will not integrate to any TOUT beyond the currently
C     specified TSTOP.  Once TSTOP has been reached, you must change
C     the value of TSTOP or set INFO(4) = 0.  You may change INFO(4)
C     or TSTOP at any time but you must supply the value of TSTOP in
C     RWORK(1) whenever you set INFO(4) = 1.
C
C     Do not change INFO(5), INFO(6), INFO(12-17) or their associated
C     IWORK/RWORK locations unless you are going to restart the code.
C
C                    *** FOLLOWING A COMPLETED TASK ***
C
C     If..
C     IDID = 1, call the code again to continue the integration
C                  another step in the direction of TOUT.
C
C     IDID = 2 or 3, define a new TOUT and call the code again.
C                  TOUT must be different from T.  You cannot change
C                  the direction of integration without restarting.
C
C     IDID = 4, reset INFO(11) = 0 and call the code again to begin
C                  the integration.  (If you leave INFO(11) > 0 and
C                  INFO(14) = 1, you may generate an infinite loop.)
C                  In this situation, the next call to DDASKR is
C                  considered to be the first call for the problem,
C                  in that all initializations are done.
C
C     IDID = 5, call the code again to continue the integration in the
C                  direction of TOUT.  You may change the functions
C                  Ri defined by RT after a return with IDID = 5, but
C                  the number of constraint functions NRT must remain
C                  the same.  If you wish to change the functions in
C                  RES or in RT, then you must restart the code.
C
C                    *** FOLLOWING AN INTERRUPTED TASK ***
C
C     To show the code that you realize the task was interrupted and
C     that you want to continue, you must take appropriate action and
C     set INFO(1) = 1.
C
C     If..
C     IDID = -1, the code has taken about 500 steps.  If you want to
C                  continue, set INFO(1) = 1 and call the code again.
C                  An additional 500 steps will be allowed.
C
C
C     IDID = -2, the error tolerances RTOL, ATOL have been increased
C                  to values the code estimates appropriate for
C                  continuing.  You may want to change them yourself.
C                  If you are sure you want to continue with relaxed
C                  error tolerances, set INFO(1) = 1 and call the code
C                  again.
C
C     IDID = -3, a solution component is zero and you set the
C                  corresponding component of ATOL to zero.  If you
C                  are sure you want to continue, you must first alter
C                  the error criterion to use positive values of ATOL
C                  for those components corresponding to zero solution
C                  components, then set INFO(1) = 1 and call the code
C                  again.
C
C     IDID = -4  --- cannot occur with this code.
C
C     IDID = -5, your JAC routine failed with the Krylov method.  Check
C                  for errors in JAC and restart the integration.
C
C     IDID = -6, repeated error test failures occurred on the last
C                  attempted step in DDASKR.  A singularity in the
C                  solution may be present.  If you are absolutely
C                  certain you want to continue, you should restart
C                  the integration.  (Provide initial values of Y and
C                  YPRIME which are consistent.)
C
C     IDID = -7, repeated convergence test failures occurred on the last
C                  attempted step in DDASKR.  An inaccurate or ill-
C                  conditioned Jacobian or preconditioner may be the
C                  problem.  If you are absolutely certain you want
C                  to continue, you should restart the integration.
C
C
C     IDID = -8, the matrix of partial derivatives is singular, with
C                  the use of direct methods.  Some of your equations
C                  may be redundant.  DDASKR cannot solve the problem
C                  as stated.  It is possible that the redundant
C                  equations could be removed, and then DDASKR could
C                  solve the problem.  It is also possible that a
C                  solution to your problem either does not exist
C                  or is not unique.
C
C     IDID = -9, DDASKR had multiple convergence test failures, preceded
C                  by multiple error test failures, on the last
C                  attempted step.  It is possible that your problem is
C                  ill-posed and cannot be solved using this code.  Or,
C                  there may be a discontinuity or a singularity in the
C                  solution.  If you are absolutely certain you want to
C                  continue, you should restart the integration.
C
C     IDID = -10, DDASKR had multiple convergence test failures
C                  because IRES was equal to -1.  If you are
C                  absolutely certain you want to continue, you
C                  should restart the integration.
C
C     IDID = -11, there was an unrecoverable error (IRES = -2) from RES
C                  inside the nonlinear system solver.  Determine the
C                  cause before trying again.
C
C     IDID = -12, DDASKR failed to compute the initial Y and YPRIME
C                  vectors.  This could happen because the initial
C                  approximation to Y or YPRIME was not very good, or
C                  because no consistent values of these vectors exist.
C                  The problem could also be caused by an inaccurate or
C                  singular iteration matrix, or a poor preconditioner.
C
C     IDID = -13, there was an unrecoverable error encountered inside
C                  your PSOL routine.  Determine the cause before
C                  trying again.
C
C     IDID = -14, the Krylov linear system solver failed to achieve
C                  convergence.  This may be due to ill-conditioning
C                  in the iteration matrix, or a singularity in the
C                  preconditioner (if one is being used).
C                  Another possibility is that there is a better
C                  choice of Krylov parameters (see INFO(13)).
C                  Possibly the failure is caused by redundant equations
C                  in the system, or by inconsistent equations.
C                  In that case, reformulate the system to make it
C                  consistent and non-redundant.
C
C     IDID = -15,..,-32 --- Cannot occur with this code.
C
C                       *** FOLLOWING A TERMINATED TASK ***
C
C     If IDID = -33, you cannot continue the solution of this problem.
C                  An attempt to do so will result in your run being
C                  terminated.
C
C  ---------------------------------------------------------------------
C
C***REFERENCES
C  1.  L. R. Petzold, A Description of DASSL: A Differential/Algebraic
C      System Solver, in Scientific Computing, R. S. Stepleman et al.
C      (Eds.), North-Holland, Amsterdam, 1983, pp. 65-68.
C  2.  K. E. Brenan, S. L. Campbell, and L. R. Petzold, Numerical
C      Solution of Initial-Value Problems in Differential-Algebraic
C      Equations, Elsevier, New York, 1989.
C  3.  P. N. Brown and A. C. Hindmarsh, Reduced Storage Matrix Methods
C      in Stiff ODE Systems, J. Applied Mathematics and Computation,
C      31 (1989), pp. 40-91.
C  4.  P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Using Krylov
C      Methods in the Solution of Large-Scale Differential-Algebraic
C      Systems, SIAM J. Sci. Comp., 15 (1994), pp. 1467-1488.
C  5.  P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Consistent
C      Initial Condition Calculation for Differential-Algebraic
C      Systems, SIAM J. Sci. Comp. 19 (1998), pp. 1495-1512.
```
"""
function unsafe_solve(callback, N, t, y, yp,
                      tout, info, rtol, atol,
                      idid, rwork, lrw, iwork,
                      liw, rpar, ipar, jac, psol,
                      rt, nrt, jroot)
    ccall(Libdl.dlsym(lib, :ddaskr_), Nothing,
          (Ptr{Nothing}, Ref{Int32}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, # RES, NEQ, T, Y, YPRIME
           Ptr{Float64}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},            # TOUT, INFO, RTOL, ATOL
           Ref{Int32}, Ptr{Float64}, Ref{Int32}, Ptr{Int32},                # IDID, RWORK, LRW, IWORK
           Ref{Int32}, Any, Ptr{Int32}, Ptr{Nothing}, Ptr{Nothing},               # LIW, RPAR, IPAR, JAC, PSOL
           Ptr{Nothing}, Ptr{Int32}, Ptr{Int32}),                              # RT, NRT, JROOT
          callback, N, t, y, yp, tout, info, rtol, atol,
          idid, rwork, lrw, iwork, liw, rpar, ipar, jac, psol,
          rt, nrt, jroot)
end
