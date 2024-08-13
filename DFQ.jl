using OrdinaryDiffEq
using BenchmarkTools

##

#du/dt = f(u,t,p), definition of a differential equation

function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8/3) * u[3]
    return du
end

u0 = [1.0; 0; 0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz!, u0, tspan)
sol = solve(prob, Tsit5())

using Plots;
plot(sol, vars = (1,2,3))