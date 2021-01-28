using Attitude, SatelliteDynamics, MATLAB
using StaticArrays
using Geodesy
using LinearAlgebra
using FiniteDiff

# canonical units (if needed)
tscale = 806.80415
dscale = copy(R_EARTH)


# ground stations
gs_ecef = [@SVector zeros(3) for i = 1:3]
gs_lla = [[-79.9479335;40.4427217;17],
          [-79.883541;40.558423;378],
          [-79.7889967;40.4652298;268]]
# Wean Hall
gs_ecef[1] = sGEODtoECEF( [-79.9479335;40.4427217;17],use_degrees = true  )

# Hartwood Acres Park
gs_ecef[2] = sGEODtoECEF( [-79.883541;40.558423;378],use_degrees = true  )

# Gascola
gs_ecef[3] = sGEODtoECEF( [-79.7889967;40.4652298;268],use_degrees = true  )


# initial orbital conditions
oe0  = [R_EARTH + 532e3, 0.001, 97.5, 78.7, 153.037, 32+210]
eci0 = sOSCtoCART(oe0, use_degrees=true)


epc = Epoch(2020, 3, 15, 12, 0, 0, 0.0) + 38000
r_eci = eci0[1:3]
r_ecef = rECItoECEF(epc)*r_eci

# @show eci0[3] - gs_ecef[1][3]
@show norm(r_ecef - gs_ecef[1])



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

function dynamics(t,x,u)
    r = SA[x[1],x[2],x[3]]
    a = -GM_EARTH*r/(norm(r)^3)
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

dt = 10
N = 25
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

# generate measurements
measurement_ts = copy(t_hist)
M = 25
Y = [@SVector zeros(3) for i = 1:M]
for i = 1:M
    Y[i]=measurement(X[i],t_hist[i])
end
idx_r = [((i-1)*3 .+ (1:3)) for i = 1:M]
# now lets optimize something
function rollout(eci_initial)
    X = [@SVector zeros(6) for i = 1:N]
    X[1] = copy(eci_initial)

    for i = 1:N-1
        X[i+1] = rk4_orbital(dynamics,0,X[i],0,dt)
    end
    X_r = [SA[X[i][1],X[i][2],X[i][3]] for i = 1:length(X)]
    eci_hist_cubic = CubicSplineInterpolation(t_hist,X_r)
    # eci_hist_cubic = LinearInterpolation(t_hist,X_r)
    return eci_hist_cubic
end
function residual(eci_initial)
    eci_hist_cubic = rollout(eci_initial)
    outp = zeros(M*3)
    for i = 1:M
        t = measurement_ts[i]
        outp[idx_r[i]] = measurement(eci_hist_cubic(t),t) - Y[i]
    end
    return outp
end


rt = residual(eci0+ 100*randn(6))

# cub = rollout(eci0)



function GN(v0)

# v0 = v1 + 100*randn(3)

new_S = 1e15
for i = 1:10

    res = residual(v0)
    S = dot(res,res)
    # @show S

    J = FiniteDiff.finite_difference_jacobian(residual,v0)
    newton_step = -J\res
    @infiltrate
    error()
    α = 1.0
    for i = 1:10
        new_v = v0 + α*newton_step
        new_res = residual(new_v)
        new_S = dot(new_res,new_res)
        if new_S < S
            v0 = copy(new_v)
            break
        else
            α /=2
        end
    end
    res = residual(v0)
    S = dot(res,res)
    @show S
    if norm(newton_step)<0.1
        break
    end

end

end


GN(eci0+ [100000*randn(3);100*randn(3)])

#TODO tomorrow or tonight:
# add random timing for the measurements and make sure it still converges on
# the correct orbit




# ## LLA hist
# lla = [sECEFtoGEOD(ecef_hist[:,i]) for i = 1:N]
# llam = rad2deg.(mat_from_vec(lla))
# gs1 = gs_lla[1]
# gs2 = gs_lla[2]
# gs3 = gs_lla[3]
# mat"
# figure
# hold on
# load('topo.mat', 'topo');
# topoplot = [topo(:, 181:360), topo(:, 1:180)];
# contour(-180:179, -90:89, topoplot, [0, 0], 'black');
# plot($llam(1,:),$llam(2,:),'b')
# plot($gs1(1),$gs1(2),'r*')
# plot($gs2(1),$gs2(2),'r*')
# plot($gs3(1),$gs3(2),'r*')
# %xlim([-82 -78])
# %ylim([38 42])
# %axis equal
# grid on
# hold off
# "
#





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
