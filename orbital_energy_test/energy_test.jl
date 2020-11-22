# this is how we activate our virtual environment
# using Pkg
# cd(dirname(@__FILE__))
# Pkg.activate(".")
# Pkg.instantiate()

# now that we are in our virtual environment, we can load all our packages
# NOTE: OrdinaryDiffEq is just an ODE specific version of DifferentialEquations
using LinearAlgebra, OrdinaryDiffEq, StaticArrays, MATLAB


# units of km for position, and km/s for velocity
GM_EARTH = 3.986004415e14/(1e9)
J2_EARTH = 0.0010826358191967
R_EARTH = 6.3781363e6/1000

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
tspan = (0.0,60000.0)
prob = ODEProblem(dynamics,x0,tspan)

# Vern7 is probably the best for orbits, reltol
# <https://diffeq.sciml.ai/stable/basics/common_solver_opts/>
sol = solve(prob,Vern7(),reltol = 1e-12)

# get specific energy for each point in the trajectory
energy_vec = [specific_energy(sol.u[i]) for i = 1:length(sol.u)]

X = sol.u

Xmat = mat_from_vec(X)

r_eci1 = Xmat[1:3,:]
v_eci1 = Xmat[4:6,:]

# plot using the ELITE matlab plotting api
t1 = sol.t

## part two
dscale = 1e6
tscale = 3600/5
GM_EARTH = 3.986004415e14*tscale^2/dscale^3
J2_EARTH = 0.0010826358191967
R_EARTH = 6.3781363e6/dscale

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
x0 = [-1.8617015559490975e6, -6.645097023400109e6, 32.941584155793194,
             -955.5092604431692, 293.91099680436355, 7541.418280028347]
x0[1:3] /= dscale
x0[4:6] /= (dscale/tscale)


# run it for a year
tspan = (0.0,60000.0/tscale)
prob = ODEProblem(dynamics,x0,tspan)

# Vern7 is probably the best for orbits, reltol
# <https://diffeq.sciml.ai/stable/basics/common_solver_opts/>
sol = solve(prob,Vern7(),reltol = 1e-12)

# get specific energy for each point in the trajectory
energy_vec = [specific_energy(sol.u[i]) for i = 1:length(sol.u)]

X = sol.u

Xmat = mat_from_vec(X)

r_eci2 = Xmat[1:3,:]
v_eci2 = Xmat[4:6,:]

# plot using the ELITE matlab plotting api
t2 = sol.t


t2 = tscale*t2
mat"
figure(1)
hold on
plot($t2,$r_eci2(1,:))
plot($t1,$r_eci1(1,:)/1000,'*')
legend('test1','test2')
hold off
"
# mat"
# figure(2)
# hold on
# plot(1:10,randn(10,1))
# plot(1:10,randn(10,1))
# plot($t1,$r_eci1(1,:),'*')
# hold off
# "

# mat"
# [X,Y,Z] = sphere(30)
# R = $R_EARTH
# figure
# hold on
# surf(R*X,R*Y,R*Z)
# plot3($r_eci(1,:),$r_eci(2,:),$r_eci(3,:))
# hold off
# "


mat"
figure
hold on
title('Specific Mechanical Energy for J2 Perturbed Orbit')
plot($energy_vec)
ylabel('Specific Mechanical Energy')
xlabel('Days')
xlim([0 365])
hold off
"
