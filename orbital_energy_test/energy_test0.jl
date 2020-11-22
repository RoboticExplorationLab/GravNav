# this is how we activate our virtual environment
using Pkg
cd(dirname(@__FILE__))
Pkg.activate(".")
# Pkg.instantiate()

# now that we are in our virtual environment, we can load all our packages
# NOTE: OrdinaryDiffEq is just an ODE specific version of DifferentialEquations
using LinearAlgebra, OrdinaryDiffEq, StaticArrays, MATLAB
using Attitude
using SatelliteDynamics

epc = Epoch("2018-12-20")

function accel_perturbations(epc::Epoch, x::Array{<:Real} ;
             mass::Real=1.0, area_drag::Real=0.01, coef_drag::Real=2.3,
             area_srp::Real=1.0, coef_srp::Real=1.8,
             n_grav::Integer=10, m_grav::Integer=10)

    # Extract position and velocity
    r = x[1:3]
    v = x[4:6]

    # Compute ECI to ECEF Transformation -> IAU2010 Theory
    PN = bias_precession_nutation(epc)
    E  = earth_rotation(epc)
    W  = polar_motion(epc)
    R  = W * E * PN

    # Compute sun and moon position
    r_sun  = sun_position(epc)
    r_moon = moon_position(epc)

    # Compute acceleration
    a = zeros(eltype(x), 3)

    # spherical harmonic gravity
    a += accel_gravity(x, R, n_grav, m_grav)

    # atmospheric drag
    ρ = density_harris_priester(epc,r)
    a += accel_drag([r;v],ρ,mass, area_drag, coef_drag, Array{Real, 2}(PN))

    # SRP
    nu = eclipse_conical(x, r_sun)
    a += nu*accel_srp(x, r_sun, mass, area_srp, coef_srp)

    # third body sun
    a += accel_thirdbody_sun(x, r_sun)

    # third body moon
    a += accel_thirdbody_moon(x, r_moon)

    return a
end
function dynamics(dx,x,p,t)
    """ODE for orbital motion"""
    r1 = x[1:3]
    v1 = x[4:6]

    # ECI_Q_ECEF = rECEFtoECI(epc + t)
    # a_grav = accel_gravity(r1,ECI_Q_ECEF,10,10)
    # a_sun = accel_thirdbody_sun(epc + t,r1)
    # a_moon = accel_thirdbody_moon(epc+t,r1)
    #
    # ρ = density_harris_priester(epc+t,r1)
    # PN = bias_precession_nutation(epc+t)
    # a_drag = accel_drag([r1;v1],ρ,1.0, 0.01, 2.3, Array{Real, 2}(PN))
    dx[1:3] = v1
    dx[4:6] = accel_perturbations(epc+t,x)
end


# initial conditions [r_eci (km); v_eci (km/s)]
x0 = [-1861.7015559490976, -6645.09702340011, 0.032941584155793194,
      -0.9555092604431692, 0.29391099680436356, 7.541418280028347]
x0*=1000

# run it for a year
tspan = (0.0,600000.0)
prob = ODEProblem(dynamics,x0,tspan)

# Vern7 is probably the best for orbits, reltol
# <https://diffeq.sciml.ai/stable/basics/common_solver_opts/>
sol = solve(prob,Vern7(),reltol = 1e-12)

# get specific energy for each point in the trajectory
# energy_vec = [specific_energy(sol.u[i]) for i = 1:length(sol.u)]

X = sol.u

Xmat = mat_from_vec(X)

r_eci1 = Xmat[1:3,:]
v_eci1 = Xmat[4:6,:]

# plot using the ELITE matlab plotting api
t1 = sol.t

## part two
dscale = 1e6
tscale = 3600/5
# GM_EARTH = 3.986004415e14*tscale^2/dscale^3
# J2_EARTH = 0.0010826358191967
# R_EARTH = 6.3781363e6/dscale

function specific_energy(rv_eci)
    """Specific orbital energy from position and velocity"""
    r_eci = @views rv_eci[1:3]
    v_eci = @views rv_eci[4:6]

    return norm(v_eci)^2/2 - GM_EARTH/norm(r_eci)
end

function dynamics(dx,x,p,t)
    """ODE for orbital motion"""
    r1 = dscale * x[1:3]
    v1 = (dscale/tscale) * x[4:6]

    dx[1:3] = v1 / (dscale/tscale)
    dx[4:6] = accel_perturbations(epc+t*tscale,[r1;v1]) / (dscale/tscale^2)
end

# initial conditions [r_eci (km); v_eci (km/s)]
x0 = [-1.8617015559490975e6, -6.645097023400109e6, 32.941584155793194,
             -955.5092604431692, 293.91099680436355, 7541.418280028347]
x0[1:3] /= dscale
x0[4:6] /= (dscale/tscale)


# run it for a year
tspan = (0.0,600000.0/tscale)
prob = ODEProblem(dynamics,x0,tspan)

# Vern7 is probably the best for orbits, reltol
# <https://diffeq.sciml.ai/stable/basics/common_solver_opts/>
sol = solve(prob,Vern7(),reltol = 1e-14, abstol = 1e-10)

# get specific energy for each point in the trajectory
# energy_vec = [specific_energy(sol.u[i]) for i = 1:length(sol.u)]

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
plot($t1,$r_eci1(1,:)/$dscale,'*')
plot($t2,$r_eci2(1,:))
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


# mat"
# figure
# hold on
# title('Specific Mechanical Energy for J2 Perturbed Orbit')
# plot($energy_vec)
# ylabel('Specific Mechanical Energy')
# xlabel('Days')
# xlim([0 365])
# hold off
# "
