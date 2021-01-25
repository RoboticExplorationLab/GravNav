# get the virtual environment all set up
using Pkg
using Plots
cd(@__DIR__)
Pkg.activate(".")
Pkg.instantiate()

using LinearAlgebra, ForwardDiff, Attitude, StaticArrays
using SparseArrays, IterativeSolvers, Infiltrator
const FD = ForwardDiff
using Geodesy, Optim, DifferentialEquations, Statistics
using SatelliteDynamics


using Random
# Random.seed!(645)

# distance and time scaling
dscale = 1e6
tscale = 3600/5
numMinutes = 5 # Number of minutes to propagate data (this many before AND after)

##### DATA GENERATION #####
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

function measurement(x, epc) 
    """ Measurement function for our system. 
         Consists of differences between the distances from the satellite to each of three groundstations."""

    gs1_eci = sECEFtoECI(epc, dscale * gs1_ecef / tscale) 
    gs2_eci = sECEFtoECI(epc, dscale * gs2_ecef / tscale)
    gs3_eci = sECEFtoECI(epc, dscale * gs3_ecef / tscale)

    # Get distances from each groundstation to satellite
    r1 = norm(x[1:3]-(gs1_eci / (dscale/tscale)))
    r2 = norm(x[1:3]-(gs2_eci / (dscale/tscale)))
    r3 = norm(x[1:3]-(gs3_eci / (dscale/tscale)))


    # Return difference of distances
    return SVector( (r2-r1), 
                    (r3-r1),
                    (r3-r2) )
end

function dynamics(x,u,t)
    """ODE for orbital motion"""
    r = x[1:3]*dscale
    v = x[4:6]*(dscale/tscale)

    # Ensure satellite doesn't collide with Earth
    if norm(r)<R_EARTH
        error("Impact 1")
    end

    # Return [v; a]
    return [v / (dscale/tscale);
            accel_perturbations(epc+t*tscale,[r;v]) / (dscale/tscale^2)]
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

function generate_data(x0)
    """ Generates data for a specified amount of minutes before and after the satellite is directly overhead.
        No process noise is added. No measurement noise is added.
        Takes in an initial state and timestep.
        Returns a vector of times, along with associated states and measurements.
    """

    tHist = Array(range(-(numMinutes * 60.0 / tscale), (numMinutes * 60.0 / tscale), step = dt)) # numMinutes before -> numMinutes after
    T = length(tHist)
    midPoint = convert(Int, ceil(T/2))
    
    X = zeros(nx, T)
    Y = zeros(m, T)
    u = SVector(0,0,0)
    X[:,midPoint] = x0
    
    # Go backward from midpoint
    for i = midPoint:-1:2
        t = tHist[i] 
        X[:,i-1] = rk4(dynamics, u, X[:,i], -dt, t)
        Y[:,i] = measurement(X[:,i], (epc+(t*tscale))) + (randn(3) * 240.0  / dscale)
    end
    t = tHist[1]
    Y[:,1] = measurement(X[:,1], (epc+(t*tscale))) + (randn(3) * 240.0  / dscale)
    
    # Go forward from midpoint
    for i = midPoint:(T-1) 
        t = tHist[i]
        X[:, i+1] = rk4(dynamics, u, X[:,i], dt, t)
        Y[:,i] = measurement(X[:,i], (epc+(t*tscale))) + (randn(3) * 240.0  / dscale)
    end
    t = tHist[T]
    Y[:,T] = measurement(X[:,T], (epc+(t*tscale))) + (randn(3) * 240.0  / dscale)
    
    return X, Y, tHist
end


## Set locations of the three groundstations
# Groundstation 1: Wean Hall
weanHall_lla = LLA(40.4427217, -79.9479335, 17.0) # Wean Hall, lat-lon-alt
weanHall_ecef = ECEF(weanHall_lla, wgs84)         # Wean Hall, Earth-Centered Earth-Fixed 
gs1_ecef = Array(weanHall_ecef) ./ dscale 

# Groundstation 2: National Robotics Engineering Center 
# nrec_lla = LLA(40.4712452, -79.9661841, 16.6)
# nrec_ecef = ECEF(nrec_lla, wgs84)
# gs2_ecef = Array(nrec_ecef) ./ dscale

hartwood_lla = LLA(40.558423, -79.883541, 378.0) # Hartwood Acres Park
hartwood_ecef = ECEF(hartwood_lla, wgs84)
gs2_ecef = Array(hartwood_ecef) ./ dscale


# Groundstation 3: Gascola
gascola_lla = LLA(40.4652298, -79.7889967, 268.0)
gascola_ecef = ECEF(gascola_lla, wgs84)
gs3_ecef = Array(gascola_ecef) ./ dscale

## Initial Conditions (polar orbit passing directly over GS at t = 0)
μ  = (3.9860044188)*(10^(14))     # m^3 / s^2
μ  = (μ * (tscale^2))/(dscale^3)  # Convert to appropriate units
Re = 6371000 / dscale             # Radius of Earth, in appropriate units
r0 = Re + (550000 / dscale) #+ (100000*randn()/dscale)
r0 = r0*(gs1_ecef/norm(gs1_ecef))   # Initial radius of orbit (directly above GS1)
v0 = cross(r0, cross(r0,[0;0;1]))
v0 = (sqrt(μ/norm(r0))*v0/norm(v0)) # Initial velocity (magnitude of velocity * unit vector in appropriate direction)
x0 = [r0; v0];  # Initial state of satellite, adjusted units 

epc = Epoch(2020, 3, 15, 12, 0, 0, 0.0) # Initial time for sim. Chosen arbitrarily 
dt = 0.1 / tscale

nx = 6  # Size of x vector
m =  3  # Size of measurement vector


# run sim
x, yhat, tHist = generate_data(x0);
T = length(tHist)
println("----------Finished Generating Data----------")


##### GAUSS NEWTON SOLVER #####
# Set initial guess for Gauss-Newton
xhat = zeros(size(x))
xhat[1:3, :] = x[1:3,:] + 100000 * randn(3, size(x)[2]) / dscale          # Position 
xhat[4:6, :] = x[4:6,:] + 100 * randn(3, size(x)[2]) / (dscale/tscale)    # Velocity 

x = reshape(x, length(x), 1)
xhat = reshape(xhat, length(xhat), 1) 
yhat = reshape(yhat, length(yhat), 1)

# Noise parameters (Set fairly arbitrarily)
posNoise_σ = (1e-3) / dscale # 1e-3 (adjusted) m   (Comes from drag force)
velNoise_σ = (1e-5) / (dscale/tscale) # 1e-5 (adjusted) m/s
measurementNoise_σ = (240) / dscale # This number was selected after some monte carlo testing was performed by Dr. Manchester

Q = Diagonal([posNoise_σ*ones(3); velNoise_σ*ones(3)] .^ 2)
cholQ = sqrt(Q)
invcholQ = inv(cholQ)

R = Matrix( (measurementNoise_σ^2)I, m, m)
cholR = sqrt(R)
invcholR = inv(cholR)

# indexing for x and y within the residual vector
idx_x = [(t-1)*nx .+ (1:nx) for t = 1:T]
idx_y = [(T-1)*nx + (t-1)*m .+ (1:m) for t = 1:T]


function residual(x)
    """residual vector for Gauss-Newton. rᵀr = MAP cost function"""

    u = 0.0 #SVector(0,0,0)
    r = zeros(eltype(x),nx*(T-1) + m*(T))
    midPoint = convert(Int, ceil(T/2))
    
    # NOTE - At this point, I could probably just use the forward propatation from the start point...
    ## States
    for i = 1:(T-1)
        xt = @view x[idx_x[i]]
        xtp1 = @view x[idx_x[i+1]]
        t = tHist[i]
        r[idx_x[i]] = invcholQ*(xtp1 - rk4(dynamics, u, xt, dt,t))
    end
    
    
    ## Measurements
    for i = 1:T
        xt = @view x[idx_x[i]]
        t = tHist[i] 
        r[idx_y[i]] = invcholR*(yhat[(3*(i-1)+1):(3*i)] - measurement(xt, epc + (t * tscale) ))
    end
    
    return r
end

function sparse_jacobian!(J,x)
    """Modify sparse jacobian in place"""
    #NOTE: this is just a fancy version of FD.jacobian(residual,x)
    u = 0.0
    for i = 1:T
        if i < T
            k = idx_x[i]
            xt = @view x[k]
            t = tHist[i]
            A_fx(x) = rk4(dynamics, u, x, dt,t)
            J[k,k] = -invcholQ*FD.jacobian(A_fx,xt)
            J[k,k .+ nx] = invcholQ
            _measurementClosure(x) = measurement(x, epc + (t*tscale) )
            J[idx_y[i],k] = -invcholR*FD.jacobian(_measurementClosure, xt)
        else
            k = idx_x[i]
            xt = @view x[k]
            t = tHist[i]
            _measurement_closure(x) = measurement(x, epc + (t*tscale) )
            J[idx_y[i],k] = -invcholR*FD.jacobian(_measurement_closure, xt)
        end
    end
    return nothing
end

function gauss_newton(x0)
    """Gauss-Newton for batch estimation"""

    # copy initial guess
    x = copy(x0)

    # create sparse jacobian
    J = spzeros(nx*(T-1) + m*(T),nx*T)

    Ds = 0.0
    v = zeros(length(x))

    # run Gauss-Newton for 100 iterations max
    for i = 1:25

        # ∂r/∂x
        sparse_jacobian!(J,x)
        # this is the same as:
        # J = FD.jacobian(residual,x)

        # calculate residual at x
        r = residual(x)

        # solve for Gauss-Newton step (direct, or indirect)
        v = -J\r
        # lsqr!(v,-J,r)

        # calculate current cost
        S_k = dot(r,r)
        # @show S_k

        # step size (learning rate)
        α = 1.0

        # run a simple line search
        for ii = 1:20
            x_new = x + α*v
            S_new = norm(residual(x_new))^2

            # this could be updated for strong frank-wolfe conditions
            if S_new < S_k
                x = copy(x_new)
                Ds = S_k - S_new
                # @show ii
                break
            else
                α /= 2
            end
            if ii == 20
                @warn "line search failed"
                Ds = 0

            end
        end

        # depending on problems caling, termination criteria should be updated
        if Ds < 1e-8
            break
        end

        # ----------------------------output stuff-----------------------------
        if rem((i-1),4)==0
            println("iter      α           S          dS")
        end
        S_display = round(S_k,sigdigits = 3)
        dS_display = round(Ds,sigdigits = 3)
        alpha_display = round(α,sigdigits = 3)
        println("$i         $alpha_display      $S_display    $dS_display")

    end
    return x
end

x_gn = gauss_newton(xhat);
println("----------Finished Gauss Newton----------")

##### PLOTTING #####
x_gn = reshape(x_gn, nx, :)
xhat = reshape(xhat, nx, :)
yhat = reshape(yhat, m, :)
x = reshape(x, nx, :)
yhat = reshape(yhat,m, :)

initError = zeros(size(x))
finalError = zeros(size(x))

for i = 1:T 
    initError[:,i]  = abs.(x[:,i] - xhat[:,i]) * dscale
    finalError[:,i] = abs.(x[:,i] - x_gn[:,i]) * dscale
end

tHist *= tscale
#Compare initial estimate and final estimate errors 
# err = plot( tHist, initError[1,:], label = false,linestyle = :dash)
# err = plot!(tHist, initError[2,:], label = false,linestyle = :dash)
# err = plot!(tHist, initError[3,:], label = false,linestyle = :dash)

# Final Error
err = plot( tHist, finalError[1,:], label = "x", title = "Error")
err = plot!(tHist, finalError[2,:], label = "y", xguidefontsize = 8, yguidefontsize = 8)
err = plot!(tHist, finalError[3,:], label = "z", xlabel = "time\n(s)", ylabel = "Error\n(m)")

# # Satellite Position (actual and estimate)
# pos = plot( tHist, x_gn[1,:], label = "̂x")
# pos = plot!(tHist, x_gn[2,:], label = "̂y")
# pos = plot!(tHist, x_gn[3,:], label = "̂z")
# pos = plot!(tHist, x[1,:], label = false, linestyle = :dash, title = "Satellite Position")
# pos = plot!(tHist, x[2,:], label = false, linestyle = :dash, xlabel = "time\n($tscale s)", ylabel = "Position\n($dscale m)")
# pos = plot!(tHist, x[3,:], label = false, linestyle = :dash, xguidefontsize = 8, yguidefontsize = 8)

# Measurements
meas= plot( tHist, yhat[1,:], label = "r2-r1", title = "Measurements")
meas= plot!(tHist, yhat[2,:], label = "r3-r1", xlabel = "time\n(s)", ylabel = "Relative Distance\n(m)")
meas= plot!(tHist, yhat[3,:], label = "r3-r2", xguidefontsize = 8, yguidefontsize = 8)


display(plot(meas, err, layout = (2,1)))




