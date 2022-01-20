using LinearAlgebra
using SatelliteDynamics

# distance and time scaling, these clean up the scaling issues
dscale = 1e6
tscale = 3600/5

"""
These functions heavily rely on SatelliteDynamics.jl for all the forces. Feel
free to read through the documentation:
https://github.com/sisl/SatelliteDynamics.jl           # repo
https://sisl.github.io/SatelliteDynamics.jl/latest/    # docs
The general gist is that it relies on an Epoch type for the time, and you can
add time to this epoch by simply adding numbers in seconds to it. To avoid any
scaling issues, I have included a scaled version of the dynamics below. All
this does is store scaled versions of position and velocity in the integrator,
then when it comes time to calculate the dynamics, the units are converted back
to meters and meters/second. Below I have included an accel_perturbations
function that calculates the full acceleration on the spacecraft when the input
units are [m; m/s], and the outpout units are [m/s;m/s²].
"""


function accel_perturbations(epc::Epoch, x::Array{<:Real} ;
             mass::Real=1.0, area_drag::Real=0.01, coef_drag::Real=2.3,
             area_srp::Real=1.0, coef_srp::Real=1.8,
             n_grav::Integer=10, m_grav::Integer=10)
    """Accelerations for spacecraft in LEO, ForwardDiff friendly. Units are
    meters for position and meters/second for velocity. x is [r;v]."""

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

    # Compute acceleration (eltype(x) makes this forward diff friendly)
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


function scaled_dynamics(x,u,t)
    """ODE for orbital motion"""

    # take the scaled units and convert them to m and m/s
    r = x[1:3]*dscale           # m
    v = x[4:6]*(dscale/tscale)  # m/s

    # the current time is (epc + t*dscale)
    return [v / (dscale/tscale);
            accel_perturbations(epc+t*tscale,[r;v]) / (dscale/tscale^2)]
end

function rk4(f, u, x_n, h,t_n)
    """Vanilla RK4"""
    k1 = h*f(x_n,u,t_n)
    k2 = h*f(x_n+k1/2, u,t_n + h/2)
    k3 = h*f(x_n+k2/2, u, t_n + h/2)
    k4 = h*f(x_n+k3, u, t_n + h)
    return (x_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4))
end
