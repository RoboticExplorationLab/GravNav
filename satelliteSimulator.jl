using Statistics
using DifferentialEquations
using Plots 
using LinearAlgebra
using Random, Distributions
using JLD2
using StaticArrays
include("functions.jl")
include("parameters.jl")


## Satellite 1 ##
r1 = _r0 + rand(Normal(1/10000, 1/100000)) # 10 m = 1/10,000 * (100 km); sigma = 1 m
v1 = _v0 + rand(Normal(9/2500, 9/25000)) # 10 cm/s = 1/1,000,000 (100 km) / s = 9/2500 (100km/h); sigma = 1 cm/s
x0 = [r1, 0, 0, 0, v1*cos(_theta), v1*sin(_theta)]
w = [_mu, _Re, r1, _J2] # Parameters for getOrbit function

prob = ODEProblem(dynamics, x0, tspan, w) 
sol1 = solve(prob, Vern7(), reltol = 1e-8); #, saveat = saveRate); 
##################

## Satellite 2 ##
r2 = _r0 + rand(Normal(1/10000, 1/100000)) # 10 m = 1/10,000 * (100 km); sigma = 1 m
v2 = _v0 + rand(Normal(9/2500, 9/25000)) # 10 cm/s = 1/1,000,000 (100 km) / s = 9/2500 (100km/h); sigma = 1 cm/s
x0 = [r2, 0, 0, 0, v2*cos(_theta), v2*sin(_theta)]
w = [_mu, _Re, r2, _J2] # Parameters for getOrbit function

prob = ODEProblem(dynamics, x0, tspan, w) 
sol2 = solve(prob, Vern7(), reltol = 1e-8);# , saveat = saveRate);
#################

## Satellite 3 ## 
r3 = _r0 + rand(Normal(1/10000, 1/100000)) # 10 m = 1/10,000 * (100 km); sigma = 1 m
v3 = _v0 + rand(Normal(9/2500, 9/25000)) # 10 cm/s = 1/1,000,000 (100 km) / s = 9/2500 (100km/h); sigma = 1 cm/s
x0 = [r3, 0, 0, 0, v3*cos(_theta), v3*sin(_theta)]
w = [_mu, _Re, r3, _J2] # Parameters for getOrbit function

prob = ODEProblem(dynamics, x0, tspan, w) 
sol3 = solve(prob, Vern7(), reltol = 1e-8);# , saveat = saveRate);
#################


### PLOTS ###

# Positions 
pos1 = getPositionPlot(sol1)
pos2 = getPositionPlot(sol2)
pos3 = getPositionPlot(sol3) 
display(plot(pos1, pos2, pos3, layout = (3,1), reuse = false))


# Get differences of radii
y = getDifferences(sol1, sol2, sol3)
# simData = Dict("y"=>y, "sol1"=>sol1, "sol2"=>sol2, "sol3"=>sol3)
@save "Testing.jld2" y=y sol1=sol1 sol2=sol2 sol3=sol3
println("Saved!")

readline() # Pauses before closing terminal
