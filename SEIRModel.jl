# import packages
import Pkg; Pkg.add("StaticArrays")
Pkg.add("CSV"); Pkg.add("DataFrames")

using DifferentialEquations
using StaticArrays
using Plots; gr()
using CSV
using DataFrames

# SEIR function of the model
function SEIR(dx,x0,params,t)

    S, E, I, R = x0
    beta, U_lambda, mu, gamma, delta, alpha = params

    dx[1] = U_lambda - (beta * S * I) - (mu * S)       # dS/dt
    dx[2] = (beta * S * I) - (gamma + mu + delta) * E  # dE/dt
    dx[3] = (gamma * E) - (alpha + mu) * I             # dI/dt
    dx[4] = (delta * E) + (alpha * I) - (mu * R)       # dR/dt

end

# Initial values for S, E, I, R respectively
x0 = [34218200.0,1000.0,157.0,100.0]

# Parameters of the SEIR model
beta     = 0.00000000118   # transmission rate from susceptible to infected
U_lambda = 0.00007         # new births per population
mu       = 0.00003         # rate of natural death
gamma    = 0.2             # transmission rate of infected from exposed
delta    = 0.1             # rate of recovery from exposed
alpha    = 0.03            # rate of recovery from infected

# The timespan for the model
tspan = (0.0, 260.0)

# ODE problem
params   = [beta, U_lambda, mu, gamma, delta, alpha]
problem  = ODEProblem(SEIR,x0,tspan,params)
solution = solve(problem)

# plot the Susceptibles result from the model
plot(solution,label="S", vars=(0,1))

# plot the Exposed result from the model
plot!(solution,label="E", vars=(0,2))

# plot the Infected result from the model
plot!(solution,label="I", vars=(0,3))

# plot the Recovered result from the model
plot!(solution,label="R", vars=(0,4), legend=:topright)

# plot real world infected data
df = DataFrame(CSV.File("ksa data filtered.csv"))
plot(df.date,df.total_cases, xlabel="days", ylabel ="Total cases", legend = false, )

# Save the graph of real world data
savefig("real world data.png")
