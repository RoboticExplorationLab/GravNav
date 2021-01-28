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
function rk4_orbital(f::Function, t_n, x_n, u, h::Real)
    """Runge-Kutta 4th order integration. Epoch for time.
    Args:
        ODE:           f(t,x,u)        | function
        initial time:  t_n             | epoch or float64
        initial cond:  x_n             | vec
        control input: u               | vec
        time step:     h               | scalar
    Returns:
        x_{n+1}
    """

    k1 = h*f(t_n,x_n,u)
    k2 = h*f(t_n+h/2,x_n+k1/2,u)
    k3 = h*f(t_n+h/2,x_n+k2/2,u)
    k4 = h*f(t_n+h,x_n+k3,u)

    return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))

end
function J2_accel(r_eci)

    # relavent parameters
    J2 = J2_EARTH
    # J2 = 0.0
    μ = GM_EARTH

    # eci position stuff
    r = norm(r_eci)
    x,y,z = r_eci

    # precompute repeated stuff
    Re_r_sqr = 1.5*J2*(R_EARTH/r)^2
    five_z_sqr = 5*z^2/r^2

    return  SA[(-μ/r^3)*x*(1 - Re_r_sqr*(five_z_sqr - 1)),
                      (-μ/r^3)*y*(1 - Re_r_sqr*(five_z_sqr - 1)),
                      (-μ/r^3)*z*(1 - Re_r_sqr*(five_z_sqr - 3))]

end
function dynamics(t,x,u)
    r = SA[x[1],x[2],x[3]]
    # a = -GM_EARTH*r/(norm(r)^3)
    a =J2_accel(r)
    return SA[x[4],x[5],x[6],a[1],a[2],a[3]]
end

function measurement(x, t)
    """ Measurement function for our system.
         Consists of differences between the distances from the satellite to each of three groundstations."""
    epct = epc + t
    gs1_eci = sECEFtoECI(epct, gs_ecef[1])
    gs2_eci = sECEFtoECI(epct, gs_ecef[2])
    gs3_eci = sECEFtoECI(epct, gs_ecef[3])

    # Get distances from each groundstation to satellite
    r1 = norm(x[1:3]-(gs1_eci))
    r2 = norm(x[1:3]-(gs2_eci))
    r3 = norm(x[1:3]-(gs3_eci))


    # Return difference of distances
    return SA[(r2-r1),(r3-r1),(r3-r2)]
end

function generate_measurements(dt,N,eci0,M,sensor_σ)
    # dt = 10
    # N = 40
    X = [@SVector zeros(6) for i = 1:N]
    X[1] = copy(eci0)

    for i = 1:N-1
        X[i+1] = rk4_orbital(dynamics,0,X[i],0,dt)
    end

    Xm = mat_from_vec(X)
    eci_hist = Xm[1:3,:]
    ecef_hist = [rECItoECEF(epc + (i-1)*dt)*eci_hist[:,i] for i = 1:N]
    ecef_hist = mat_from_vec(ecef_hist)
    t_hist = 0:dt:(N-1)*dt
    X_r = [SA[X[i][1],X[i][2],X[i][3]] for i = 1:length(X)]
    eci_hist_cubic = CubicSplineInterpolation(t_hist,X_r)


    Y_ts = [rand_in_range(0,(N-1)*dt) for i = 1:M]
    # M = 25
    Y = [@SVector zeros(3) for i = 1:M]
    for i = 1:M
        ty = Y_ts[i]
        Y[i]=measurement(eci_hist_cubic(ty),ty) + sensor_σ*randn(3)
    end

    idx_r = [((i-1)*3 .+ (1:3)) for i = 1:M]

    return X,eci_hist,ecef_hist,Y,Y_ts, idx_r, t_hist
end

function rollout_guess(x_gn)
    X = [@SVector zeros(6) for i = 1:N]
    X[1] = copy(x_gn)
    for i = 1:N-1
        X[i+1] = rk4_orbital(dynamics,0,X[i],0,dt)
    end
    X_r = [SA[X[i][1],X[i][2],X[i][3]] for i = 1:length(X)]
    ecef_hist_gn = [rECItoECEF(epc + (i-1)*dt)*X_r[i] for i = 1:N]
    ecef_hist_gn = mat_from_vec(ecef_hist_gn)

    return X, X_r, ecef_hist_gn
end
function rollout(eci_initial)
    X = [@SVector zeros(6) for i = 1:N]
    X[1] = copy(eci_initial)
    for i = 1:N-1
        X[i+1] = rk4_orbital(dynamics,0,X[i],0,dt)
    end
    X_r = [SA[X[i][1],X[i][2],X[i][3]] for i = 1:length(X)]
    eci_hist_cubic = CubicSplineInterpolation(t_hist,X_r)
    return eci_hist_cubic
end
