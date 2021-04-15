using SatelliteDynamics, LinearAlgebra, Attitude, MATLAB
using StaticArrays, Interpolations
using Interpolations
using DelimitedFiles, ProgressMeter
# include(joinpath(dirname(@__FILE__),"mag_field.jl"))


function sim_driver()

    # Declare simulation initial Epoch (NOTE: this must be after 2020)
    epc0 = Epoch(2020, 1, 1, 12, 0, 0, 0.0)

    # Declare initial state in terms of osculating orbital elements
    sma = R_EARTH + 550e3  # semi-major axis     (m)
    ecc = 0.01             # eccentricity        ()
    inc = 96.0             # inclination         (deg)
    Ω = 75.0               # RAAN                (deg)
    ω = 30.0               # argument of perigee (deg)
    M = 0.0                # mean anomaly        (deg)
    oe0  = [sma, ecc, inc, Ω, ω, M]

    # Convert osculating elements to Cartesean state
    eci0 = sOSCtoCART(oe0, use_degrees=true)

    # Set the propagation end time to 5 orbit periods after the start
    dt = 10.0
    # T = 3*orbit_period(oe0[1])
    T = 86400+100
    T = round(T - rem(T,dt))
    epcf = epc0 + T

    # Initialize State Vector
    orb  = EarthInertialState(epc0, eci0, dt=dt,
            mass=100.0, n_grav=10, m_grav=10,
            drag=true, srp=true,
            moon=true, sun=true,
            relativity=false
    )

    # simulate and process
    sim_data = postprocess_sim(sim!(orb, epcf)...)

    return sim_data

end

function postprocess_sim(t_hist, epc_hist, eci_hist)
    N = size(eci_hist,2)

    r_hist = [eci_hist[1:3,i] for i = 1:N]
    v_hist = [eci_hist[4:6,i] for i = 1:N]
    eci_hist = [[r_hist[i];v_hist[i]] for i = 1:length(r_hist)]

    # magnetic field
    B_eci_T = [IGRF13(r_hist[i],epc_hist[i]) for i = 1:N]

    # sun location
    r_sun_eci = [sun_position(epc_hist[i]) for i = 1:N]

    # ecef stuff
    ecef_hist = [rECItoECEF(epc_hist[i])*r_hist[i] for i = 1:N]

    ## LLA hist
    lla_hist = [sECEFtoGEOD(ecef_hist[i]) for i = 1:N]

    # create our sim data struct
    sim_data = sim_results(t_hist,epc_hist,eci_hist,r_hist,v_hist,B_eci_T,
                           r_sun_eci,ecef_hist,lla_hist)

    return sim_data
end

struct sim_results
    t_hist     ::Array{Float64,1}
    epc_hist   ::Array{Epoch,1}
    eci_hist   ::Array{Array{Float64,1},1}
    r_hist     ::Array{Array{Float64,1},1}
    v_hist     ::Array{Array{Float64,1},1}
    B_eci_T    ::Array{Array{Float64,1},1}
    r_sun_eci  ::Array{Array{Float64,1},1}
    ecef_hist  ::Array{Array{Float64,1},1}
    lla_hist   ::Array{Array{Float64,1},1}
end

function create_interps(sim_results)
    # create named tuple with interpolation objects
    t = sim_results.t_hist
    interp_t = 0:(t[2]-t[1]):t[end]
    eci  = LinearInterpolation(interp_t,sim_results.eci_hist)
    B    = LinearInterpolation(interp_t,sim_results.B_eci_T)
    sun  = LinearInterpolation(interp_t,sim_results.r_sun_eci)
    ecef = LinearInterpolation(interp_t,sim_results.ecef_hist)
    return interps = (eci = eci, B=B,sun = sun, ecef = ecef)

end

struct SAT
    epc0::Epoch
    interps::NamedTuple
    J::Diagonal{Float64,Array{Float64,1}}
end
struct SENSORS
    gyro::Array{Array{Float64,1},1}
    B_b::Array{Array{Float64,1},1}
    ss::Array{Array{Float64,1},1}
end
function dynamics(model::SAT,x,u,t)
    q = x[1:4]
    ω = x[5:7]
    q̇ = .5* (q ⊙ [0;ω])
    α = model.J\(u - cross(ω,model.J*ω))
    return [q̇;α]
end
function rk4(model,x_n,u,t_n,h)

    k1 = h*dynamics(model,x_n,u,t_n)
    k2 = h*dynamics(model,x_n+k1/2,u,t_n+h/2)
    k3 = h*dynamics(model,x_n+k2/2,u,t_n+h/2)
    k4 = h*dynamics(model,x_n+k3,u,t_n+h)

    return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))
end
function attitude_sim(init,interps,epc0,dt,tf)
    N = Int(round(tf/dt))

    X = [zeros(length(init)) for i = 1:N]
    sensors = SENSORS([zeros(3) for i = 1:N],[zeros(3) for i = 1:N],[zeros(6) for i = 1:N])
    X[1] = copy(init)
    model = SAT(epc0,interps,Diagonal([1,2,3.0]))

    @showprogress "simulating attitude" for i = 1:N-1
        t = dt*(i-1)
        X[i+1] = rk4(model,X[i],zeros(3),t,dt)
        normalize!(X[i+1][1:4])
        sensors.gyro[i],sensors.B_b[i],sensors.ss[i] = generate_measurements(X[i],interps,t,epc0)
    end
    sensors.gyro[N],sensors.B_b[N],sensors.ss[N] = generate_measurements(X[N],interps,dt*(N-1),epc0)
    return X, sensors
end
function sun_body_normalized(r_sun_eci,r_eci,ᴺQᴮ)
    # positon vector from spacecraft to sun
    sc_r_sun = r_sun_eci - r_eci[1:3]
    # normalize and express in the body frame
    return transpose(ᴺQᴮ)*normalize(sc_r_sun)
end
function sun_flux(r_sun_eci,r_eci,ᴺQᴮ,faces,ν)
    # normalize and express in the body frame
    s = sun_body_normalized(r_sun_eci,r_eci,ᴺQᴮ)

    # Get dot products
    I_vec = faces*s

    # only keep the positive ones (negative ones go to 0)
    I_vec = I_vec .* (I_vec .> 0.0)

    return ν*I_vec
end
function generate_measurements(x,interps,t,epc0)
    q = x[1:4]
    ω = x[5:7]
    Q = dcm_from_q(q)
    faces = [I(3);-I(3)]

    r_sun = interps.sun(t)
    r_eci = interps.eci(t)
    B = interps.B(t)
    B_body = Q'*B

    ν = eclipse_conical(r_eci,r_sun)

    I_flux = sun_flux(r_sun,r_eci,Q,faces,ν)

    # add noise
    ω += (deg2rad(1)*randn(3))
    B_body += (.05*norm(B_body)*randn(3))
    I_flux += (0.05*randn(6))

    return ω, B_body, I_flux
end
function runsim()

    sim_results = sim_driver()
    @info "done simulating orbit"
    interps = create_interps(sim_results)

    q0 = [1;0;0;0]
    ω0 = deg2rad.([1.5;1.3;-1.2])
    # ω0 = [.2;7;.2]
    X, sensors = attitude_sim([q0;ω0],interps,sim_results.epc_hist[1],.1,86400.0)
    # xm = mat_from_vec(X)
    # mat"
    # figure
    # hold on
    # plot($xm(5:7,:)')
    # hold off
    # "

    gyro = rad2deg.(mat_from_vec(sensors.gyro))
    mat"
    figure
    hold on
    plot($gyro')
    hold off
    "

    ss = mat_from_vec(sensors.ss)
    mat"
    figure
    hold on
    plot($ss')
    hold off
    "

    B_b = 1e6*mat_from_vec(sensors.B_b)
    mat"
    figure
    hold on
    plot($B_b')
    hold off
    "

    csvmat = [gyro' ss' B_b']
    writedlm("sensor_data.csv", csvmat ,',')


    return nothing
end
runsim()
# sd2 =  sim_driver()

# # plot eci position
# eci_hist = copy(mat_from_vec(sd2.eci_hist))
# mat"
# [x,y,z] = sphere(40);
# imgRGB = imread('earth.jpg');
# figure
# hold on
# title('Orbit in ECI')
# warp($R_EARTH*x,$R_EARTH*y,$R_EARTH*z,circshift(rot90(imgRGB,2),569,2))
# plot3($eci_hist(1,:),$eci_hist(2,:),$eci_hist(3,:))
# view(150,34)
# xlabel('ECI X (m)')
# ylabel('ECI Y (m)')
# zlabel('ECI Z (m)')
# axis equal
# hold off
# "
#
# # plot the ground tracks
# llam = mat_from_vec(sd2.lla_hist)
# mat"
# figure
# hold on
# title('Ground Tracks for Orbit')
# load('topo.mat', 'topo');
# topoplot = [topo(:, 181:360), topo(:, 1:180)];
# contour(-180:179, -90:89, topoplot, [0, 0], 'black');
# plot(rad2deg($llam(1,:)),rad2deg($llam(2,:)),'b.')
# xlabel('Longitude (deg)')
# ylabel('Latitude (deg)')
# axis equal
# grid on
# hold off
# "
#
# # plot magnetic field
# B_plot = mat_from_vec(sd2.B_eci_T)
# mat"
# figure
# title('Magnetic Field in ECI')
# hold on
# plot($sd2.t_hist/3600,$B_plot')
# xlabel('Time (hours)')
# ylabel('Magnetic Moment (Tesla)')
# legend('B_x','B_y','B_z')
# hold off
# "
