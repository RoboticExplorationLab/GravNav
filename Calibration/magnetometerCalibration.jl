# get the virtual environment all set up
using Pkg
using Plots
# using PyPlot
cd(@__DIR__)
Pkg.activate(".")
Pkg.instantiate()


using LinearAlgebra, ForwardDiff, Attitude, StaticArrays
using SparseArrays, IterativeSolvers, Infiltrator
const FD = ForwardDiff
using SatelliteDynamics
include("mag_field.jl")


using Random
# Random.seed!(6345)

function accel_perturbations(epc::Epoch, x::Array{<:Real} ;
             mass::Real=1.0, area_drag::Real=0.01, coef_drag::Real=2.3,
             area_srp::Real=1.0, coef_srp::Real=1.8,
             n_grav::Integer=10, m_grav::Integer=10)
    """Accelerations for spacecraft in LEO, ForwardDiff friendly"""

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

function measurement(x,t)
    """Measurement function for our system. Comprises magnetic field estimate along each axis"""

    r = SVector(x[1],x[2],x[3])
    v = SVector(x[4],x[5],x[6])

    # the current time is (epc + t*dscale)
    t_current = epc+t*tscale

    # sun position
    r_sun  = sun_position(t_current)

    # Orthogonal components of the geomagnetic field
    # - why are they scaled by 1e6....? <- IDK
    bx, by, bz = (IGRF13(r*dscale,t_current)) 

    b_field = [bx; by; bz];
    b_hat = (T_matrix*b_field) + bias 

    # Teslas are kg per amp s^2, so we need to scale by tScale^2 ?
    bx_hat = b_hat[1];
    by_hat = b_hat[2];
    bz_hat = b_hat[3];

    return (tscale^2)*SVector(bx_hat, by_hat, bz_hat, bx, by, bz)
    # return SVector(bx_hat, by_hat, bz_hat)
    # return SVector(bx, by, bz)
end

function dynamics(x,u,t)
    """ODE for orbital motion"""
    r1 = x[1:3]*dscale
    v1 = x[4:6]*(dscale/tscale)

    # this is my debugging stuff for if one of them hits the earth
    if norm(r1)<R_EARTH
        error("Impact 1")
    end

    # the current time is (epc + t*dscale)
    t_current = epc+t*tscale

    return [v1 / (dscale/tscale);
            accel_perturbations(t_current,[r1;v1]) / (dscale/tscale^2)]
end

function rk4(f, u, x_n, h,t_n)
    """Vanilla RK4"""
    x_n = SVector{nx}(x_n)

    k1 = h*f(x_n,u,t_n)
    k2 = h*f(x_n+k1/2, u,t_n + h/2)
    k3 = h*f(x_n+k2/2, u, t_n + h/2)
    k4 = h*f(x_n+k3, u, t_n + h)


    return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))

end

function generate_data(x0,T,dt,R)
    """Forward rollout of the sim, no process noise"""

    X = fill(zeros(nx),T)
    Y = fill(zeros(m),T)

    X[1] = x0

    u = SVector(0,0,0)
    for i = 1:(T-1)
        t = (i-1)*dt
        X[i+1] = rk4(dynamics,u,X[i],dt,t) #+ sqrt(Q)*randn(nx)
        Y[i] = measurement(X[i],t) + sqrt(R)*randn(m)
    end
    Y[T] = measurement(X[T],(T-1)*dt) +  sqrt(R)*randn(m)

    return X,Y
end

function plotEarth()
    Re = (6371000 / dscale)

    phi = 0:pi/50:2*pi;
    theta = 0:pi/100:pi
    x = [cos(t)*sin(p) for t in theta, p in phi];
    y = [sin(t)*sin(p) for t in theta, p in phi];
    z = [cos(p) for t in theta, p in phi];

    return Re*x, Re*y, Re*z

end



"""   SIMULATION PARAMETERS   """
T = 193 # number of knot points (~193 = 1 orbit)

dscale, tscale = (1e6), (3600/5)

nx = 6      # Size of state matrix [pos, vel]
m  = 6      # Size of measurement matrix [pred, meas] 

dt = 30.0/tscale # sample time is 30 seconds

""" Noise """
measurement_noise = 0.0001*(tscale^2); # Teslas (kg/ (amp s^2))  

R = Diagonal((measurement_noise*ones(m)).^2)
cholR = sqrt(R)
invcholR = inv(cholR)

# Q = (1e-3)*Diagonal(@SVector ones(nx))
# cholQ = sqrt(Q)
# invcholQ = inv(cholQ)

""" Model Parameters """
s_a, s_b, s_c = 1 .+ 0.05*randn(3)                 # Scale Factors 
bias_x, bias_y, bias_z = 1.0 * randn(3);           # Bias (Teslas)
ρ, ϕ, λ = 5.0*randn(3);                            # Non-orthogonality angles  (in DEGREES)
ρ, ϕ, λ = deg2rad(ρ), deg2rad(ϕ), deg2rad(λ)

# T does NOT include the time-varying current (yet)
T_matrix = [s_a 0 0; s_b*sin(ρ) s_b*cos(ρ) 0; s_c*sin(λ) s_c*sin(ϕ)*cos(λ) s_c*cos(ϕ)*cos(λ)]
bias = [bias_x; bias_y; bias_z]
 
""" Initial Conditions """
epc = Epoch(2019, 1, 1, 12, 0, 0, 0.0) # initial time for sim

# Declare initial state in terms of osculating orbital elements
oe0  = [R_EARTH + 550e3, 0.01, 45.0, 12, 25, 0]
eci0 = sOSCtoCART(oe0, use_degrees=true)

eci0[1:3] /= dscale
eci0[4:6] /= (dscale/tscale)
x0 = eci0 


# GENERATE DATA
X,Y = generate_data(x0,T,dt,R);
xVals = mat_from_vec(X)
yVals = mat_from_vec(Y) / (tscale^2)


# SOLVE FOR T, BIAS
d = 3
I3 = Matrix(1I, d, d)
A = [yVals[4,1]*I3 yVals[5,1]*I3 yVals[6,1]*I3 I3]
for i = 2:size(yVals)[2]
    row = [yVals[4,i]*I3 yVals[5,i]*I3 yVals[6,i]*I3 I3]
    global A = [A; row]
end

y_meas = yVals[1:3, :];
y_meas = reshape(y_meas, length(y_meas));

params = A \ y_meas;

T_hat = reshape(params[1:9], 3, 3)
b_hat = params[10:12]

println("T Error: ", sum(abs.(T_matrix - T_hat)))
println("Bias Error: ", sum(abs.(bias - b_hat)))



# Ex, Ey, Ez = plotEarth()

# pygui(true)

# plt = plot3D(xVals[1,:], xVals[2,:], xVals[3,:])
# plt = surf(Ex, Ey, Ez) 
# plt = surface(Ex[:], Ey[:], Ez[:], color = :black, legend = false)

# plt = scatter3d!(xVals[1,:], xVals[2,:], xVals[3,:])
# display(plt)


