using Attitude, SatelliteDynamics, MATLAB
using StaticArrays
using LinearAlgebra
using FiniteDiff
using Printf
using Interpolations


# dynamics, measurement, rollout functions
include(joinpath(dirname(@__FILE__),"dynamics.jl"))

# residual, gauss-newton functions
include(joinpath(dirname(@__FILE__),"gauss_newton_stuff.jl"))

# function for logging solver status
include(joinpath(dirname(@__FILE__),"logging.jl"))


# canonical units for scaling
tscale = 806.80415
dscale = copy(R_EARTH)


# ground stations
# this Longitude (deg), Latitude (deg), Altitude (m)
gs_lla = [[-79.9479335;40.4427217;17],   # Wean Hall
          [-79.883541;40.558423;378],    # Hartwood Acres Park
          [-79.7889967;40.4652298;268]]  # Gascola
gs_ecef = [SVector{3}(sGEODtoECEF(gs_lla[i],use_degrees = true)) for i = 1:3]

# initial orbital conditions
oe0  = [R_EARTH + 532e3, 0.001, 97.5, 78.7, 153.037, 32+210]
eci0 = sOSCtoCART(oe0, use_degrees=true)
epc = Epoch(2020, 3, 15, 12, 0, 0, 0.0) + 38000

# simulation
dt = 10 # seconds
N = 40  # number of steps for simulation
M = 4000  # number of measurements
sensor_σ = 240.0 # meters
X,eci_hist,ecef_hist,Y,Y_ts, idx_r, t_hist = generate_measurements(dt,N,eci0,M,sensor_σ)

# initial guess for Gauss-Newton
x_guess = eci0+ [100e3*randn(3);100*randn(3)]
x_guess[1:3] /= dscale
x_guess[4:6] /= (dscale/tscale)

# run Gauss-Newton
x_gn = GN(x_guess)
x_gn[1:3] *= dscale
x_gn[4:6] *= (dscale/tscale)


println(" ")
printstyled("Position error (km)\n", bold = true, color =:green)
@printf "%.4E\n" norm(x_gn[1:3] - eci0[1:3])/1000
printstyled("Velocity error (m/s)\n", bold = true, color =:green)
@printf "%.4E\n" norm(x_gn[4:6] - eci0[4:6])

X, X_r, ecef_hist_gn = rollout_guess(x_gn)




# LLA hist
lla = [sECEFtoGEOD(ecef_hist[:,i]) for i = 1:N]
llam = rad2deg.(mat_from_vec(lla))
lla_gn = [sECEFtoGEOD(ecef_hist_gn[:,i]) for i = 1:N]
llam_gn = rad2deg.(mat_from_vec(lla_gn))
gs1 = gs_lla[1]
gs2 = gs_lla[2]
gs3 = gs_lla[3]
mat"
figure
hold on
load('topo.mat', 'topo');
topoplot = [topo(:, 181:360), topo(:, 1:180)];
contour(-180:179, -90:89, topoplot, [0, 0], 'black');
h1(1) = plot($llam(1,:),$llam(2,:),'b','linewidth',2);
h2(1) = plot($llam_gn(1,:),$llam_gn(2,:),'r');
plot($gs1(1),$gs1(2),'r*')
plot($gs2(1),$gs2(2),'g*')
plot($gs3(1),$gs3(2),'b*')
legend([h1(1),h2(1)],'True Trajectory','Estimated Trajectory')
xlim([-140 -50])
ylim([0 90])
%axis equal
grid on
hold off
"






# gs1 = Array(gs_ecef[1])
# gsr = r_ecef
# nv = 1.5*gs1
# mat"
# [x,y,z] = sphere(40);
# imgRGB = imread('earth.jpg');
# figure
# hold on
# warp($R_EARTH*x,$R_EARTH*y,$R_EARTH*z,circshift(rot90(imgRGB,2),569,2))
# plot3($gs1(1),$gs1(2),$gs1(3),'g.','MarkerSize',20)
# plot3($gsr(1),$gsr(2),$gsr(3),'b.','MarkerSize',20)
# %plot3($eci_hist(1,:),$eci_hist(2,:),$eci_hist(3,:))
# plot3($ecef_hist(1,:),$ecef_hist(2,:),$ecef_hist(3,:),'r','linewidth',4)
# view(150,34)
# axis equal
# hold off
# "
