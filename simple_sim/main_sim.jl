using SatelliteDynamics, LinearAlgebra, Attitude, MATLAB

include(joinpath(dirname(@__FILE__),"mag_field.jl"))


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
    T    = 3*orbit_period(oe0[1])
    epcf = epc0 + T

    # Initialize State Vector
    orb  = EarthInertialState(epc0, eci0, dt=10.0,
            mass=100.0, n_grav=10, m_grav=10,
            drag=true, srp=true,
            moon=true, sun=true,
            relativity=false
    )

    # Propagate the orbit
    # t_hist, epc_hist, eci_hist = sim!(orb, epcf)

    # post process the results of the sim
    sim_data = postprocess_sim(sim!(orb, epcf)...)

    return sim_data

end

function postprocess_sim(t_hist, epc_hist, eci_hist)
    N = size(eci_hist,2)

    r_hist = [eci_hist[1:3,i] for i = 1:N]
    v_hist = [eci_hist[4:6,i] for i = 1:N]

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
    t_hist::Array{Float64,1}
    epc_hist::Array{Epoch,1}
    eci_hist::Array{Float64,2}
    r_hist::Array{Array{Float64,1},1}
    v_hist::Array{Array{Float64,1},1}
    B_eci_T::Array{Array{Float64,1},1}
    r_sun_eci::Array{Array{Float64,1},1}
    ecef_hist::Array{Array{Float64,1},1}
    lla_hist::Array{Array{Float64,1},1}
end

sd2 =  sim_driver()

# plot eci position
eci_hist = copy(sd2.eci_hist)
mat"
[x,y,z] = sphere(40);
imgRGB = imread('earth.jpg');
figure
hold on
title('Orbit in ECI')
warp($R_EARTH*x,$R_EARTH*y,$R_EARTH*z,circshift(rot90(imgRGB,2),569,2))
plot3($eci_hist(1,:),$eci_hist(2,:),$eci_hist(3,:))
view(150,34)
xlabel('ECI X (m)')
ylabel('ECI Y (m)')
zlabel('ECI Z (m)')
axis equal
hold off
"

# plot the ground tracks
llam = mat_from_vec(sd2.lla_hist)
mat"
figure
hold on
title('Ground Tracks for Orbit')
load('topo.mat', 'topo');
topoplot = [topo(:, 181:360), topo(:, 1:180)];
contour(-180:179, -90:89, topoplot, [0, 0], 'black');
plot(rad2deg($llam(1,:)),rad2deg($llam(2,:)),'b.')
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
axis equal
grid on
hold off
"

# plot magnetic field
B_plot = mat_from_vec(sd2.B_eci_T)
mat"
figure
title('Magnetic Field in ECI')
hold on
plot($sd2.t_hist/3600,$B_plot')
xlabel('Time (hours)')
ylabel('Magnetic Moment (Tesla)')
legend('B_x','B_y','B_z')
hold off
"
