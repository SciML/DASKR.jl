module EventTest

using DASKR, Base.Test
using DiffEqBase

function f(t,u,du,r)
    r[1] = du[1] - u[2]
    r[2] = du[2] + 9.81
end

function condition(t,u,du,r) # Event when event_f(t,u,k) == 0
    r[1] = u[1]
end

function affect!(integrator)
    println("Event!")
    integrator.u[2] = -integrator.u[2]
end

callback = DASKR.Callback(condition, [affect!], [x -> nothing])

u0 = [50.0,0.0]
du0 = [0.0,0.0]
tspan = (0.0,15.0)
prob = DAEProblem(f, u0, du0, tspan)

integrator = init(prob, DASKR.Algorithm(), callback = callback)

while integrator.daskr_args.t[] < 15
    println(integrator.daskr_args.t[], " ", integrator.daskr_args.u[1])
    step!(integrator)
end

end   # module