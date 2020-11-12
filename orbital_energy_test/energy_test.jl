# this is how we activate our virtual environment
using Pkg
cd(dirname(@__FILE__))
Pkg.activate(".")

# now that we are in our virtual environment, we can load all our packages
# NOTE: OrdinaryDiffEq is just an ODE specific version of DifferentialEquations
using LinearAlgebra, OrdinaryDiffEq, StaticArrays, MATLAB


# units of km for position, and km/s for velocity
const GM_EARTH = 3.986004415e14/(1e9)
const J2_EARTH = 0.0010826358191967
const R_EARTH = 6.3781363e6/1000

function specific_energy(rv_eci)
    """Specific orbital energy from position and velocity"""
    r_eci = @views rv_eci[1:3]
    v_eci = @views rv_eci[4:6]

    return norm(v_eci)^2/2 - GM_EARTH/norm(r_eci)
end

function dynamics(dx,x,p,t)
    """ODE for orbital motion"""
    r1 = @views x[1:3]
    v1 = @views x[4:6]

    dx[1:3] = v1
    dx[4:6] =FODE_J2(r1)
end

function FODE_J2(r_eci)
    """Acceleration caused by spherical earth and J2"""
    # eci position stuff
    r = norm(r_eci)
    x,y,z = r_eci

    # precompute repeated stuff
    Re_r_sqr = 1.5*J2_EARTH*(R_EARTH/r)^2
    five_z_sqr = 5*z^2/r^2

    return  (-GM_EARTH/r^3)*SVector(x*(1 - Re_r_sqr*(five_z_sqr - 1)),
                             y*(1 - Re_r_sqr*(five_z_sqr - 1)),
                             z*(1 - Re_r_sqr*(five_z_sqr - 3)))
end


# initial conditions [r_eci (km); v_eci (km/s)]
x0 = [-1861.7015559490976, -6645.09702340011, 0.032941584155793194,
      -0.9555092604431692, 0.29391099680436356, 7.541418280028347]


# run it for a year
tspan = (0.0,365*86400.0)
prob = ODEProblem(dynamics,x0,tspan)

# Vern7 is probably the best for orbits, reltol
# <https://diffeq.sciml.ai/stable/basics/common_solver_opts/>
sol = solve(prob,Vern7(),reltol = 1e-8)

# get specific energy for each point in the trajectory
energy_vec = [specific_energy(sol.u[i]) for i = 1:length(sol.u)]

# plot using the ELITE matlab plotting api
t = sol.t
mat"
figure
hold on
title('Specific Mechanical Energy for J2 Perturbed Orbit')
plot($t/(86400),$energy_vec)
ylabel('Specific Mechanical Energy')
xlabel('Days')
xlim([0 365])
hold off
"
