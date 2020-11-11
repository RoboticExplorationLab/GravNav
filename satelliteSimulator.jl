using Statistics
using DifferentialEquations
using Plots 
using LinearAlgebra
using Random, Distributions
using JLD
include("functions.jl")
include("parameters.jl")


## Satellite 1 ##
r1 = r0 + rand(Normal(1/10000, 1/100000)) # 10 m = 1/10,000 * (100 km); sigma = 1 m
v1 = v0 + rand(Normal(9/2500, 9/25000)) # 10 cm/s = 1/1,000,000 (100 km) / s = 9/2500 (100km/h); sigma = 1 cm/s
x0 = [r1, 0, 0, 0, v1*cos(theta), v1*sin(theta)]
w = [mu, Re, r1, J2] # Parameters for getOrbit function

prob = ODEProblem(getOrbit, x0, tspan, w) 
sol1 = solve(prob, Vern9(), saveat = saveRate); # Vern9() loses 3% energy, Feagin12() gains 7% per year
##################

## Satellite 2 ##
r2 = r0 + rand(Normal(1/10000, 1/100000)) # 10 m = 1/10,000 * (100 km); sigma = 1 m
v2 = v0 + rand(Normal(9/2500, 9/25000)) # 10 cm/s = 1/1,000,000 (100 km) / s = 9/2500 (100km/h); sigma = 1 cm/s
x0 = [r2, 0, 0, 0, v2*cos(theta), v2*sin(theta)]
w = [mu, Re, r2, J2] # Parameters for getOrbit function

prob = ODEProblem(getOrbit, x0, tspan, w) 
sol2 = solve(prob, Vern9(), saveat = saveRate);
#################

## Satellite 3 ## 
r3 = r0 + rand(Normal(1/10000, 1/100000)) # 10 m = 1/10,000 * (100 km); sigma = 1 m
v3 = v0 + rand(Normal(9/2500, 9/25000)) # 10 cm/s = 1/1,000,000 (100 km) / s = 9/2500 (100km/h); sigma = 1 cm/s
x0 = [r3, 0, 0, 0, v3*cos(theta), v3*sin(theta)]
w = [mu, Re, r3, J2] # Parameters for getOrbit function

prob = ODEProblem(getOrbit, x0, tspan, w) 
sol3 = solve(prob, Vern9(), saveat = saveRate);
#################


### PLOTS ###

# Positions 
pos1 = getPositionPlot(sol1)
pos2 = getPositionPlot(sol2)
pos3 = getPositionPlot(sol3) 
display(plot(pos1, pos2, pos3, layout = (3,1)))

# Differences
difPos_21 = getPositionPlot(sol2 - sol1)
difPos_32 = getPositionPlot(sol3 - sol2)
difPos_13 = getPositionPlot(sol1 - sol3)
display(plot(difPos_21, difPos_32, difPos_13, layout = (3,1), title = ["||S2 - S1||" "||S3 - S2||" "||S1 - S3||"])) #, xticks = ([0:(orbitTime/saveRate):length(sol1);]) ))
png("Differences")

# || Difference in radii ||
difPos_21 = zeros(length(sol1))
difPos_31 = zeros(length(sol1))
difPos_32 = zeros(length(sol1))

for i = 1:length(sol1)
    difPos_21[i] = norm(sol2[1:3,i] - sol1[1:3,i]) 
    difPos_31[i] = norm(sol3[1:3,i] - sol1[1:3,i])
    difPos_32[i] = norm(sol3[1:3,i] - sol2[1:3,i])
end 

y = hcat(difPos_21, difPos_31, difPos_32);
simData = Dict("y"=>y, "sol1"=>sol1, "sol2"=>sol2, "sol3"=>sol3)
println(size(y))
println(size(sol1))


save("test.jld", "y", y)
save("truth.jld", "sol1", sol1, "sol2", sol2, "sol3", sol3)

# readline() # Pauses before closing terminal
