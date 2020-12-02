using LinearAlgebra, MATLAB, ForwardDiff, Attitude, StaticArrays
using SparseArrays, IterativeSolvers, Infiltrator
const FD = ForwardDiff
using SatelliteDynamics

# distance and time scaling
dscale = 1e6
tscale = 3600/5

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
function relv(rv1,rv2)

        x1,y1,z1,dx1,dy1,dz1 = rv1
        x2,y2,z2,dx2,dy2,dz2 = rv2

        return ((dx1*x1 - dx1*x2 - dx2*x1 + dx2*x2 + dy1*y1 - dy1*y2 -
                dy2*y1 + dy2*y2 + dz1*z1 - dz1*z2 - dz2*z1 + dz2*z2)/((x1
                - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)^(1/2))

end
function measurement(x,t)
    """Measurement function for our system. Fill ranging for the first sat
    and relative positions between all 3"""
    r1 = x[1:3]*dscale
    v1 = x[4:6]*(dscale/tscale)
    r2 = x[7:9]*dscale
    v2 = x[10:12]*(dscale/tscale)
    r3 = x[13:15]*dscale
    v3 = x[16:18]*(dscale/tscale)
    r4 = x[19:21]*dscale
    v4 = x[22:24]*(dscale/tscale)

    x1 = [r1;v1]
    x2 = [r2;v2]
    x3 = [r3;v3]
    x4 = [r4;v4]

    return SVector(norm(r1-r2),norm(r1-r3),
                   norm(r1-r4),norm(r2-r4),
                   norm(r2-r3),norm(r3-r4),
                   relv(x1,x2),relv(x1,x3),
                   relv(x1,x4),relv(x2,x4),
                   relv(x2,x3),relv(x3,x4))
end
function dynamics(x,u,t)
    """ODE for orbital motion"""
    r1 = x[1:3]*dscale
    v1 = x[4:6]*(dscale/tscale)
    r2 = x[7:9]*dscale
    v2 = x[10:12]*(dscale/tscale)
    r3 = x[13:15]*dscale
    v3 = x[16:18]*(dscale/tscale)
    r4 = x[19:21]*dscale
    v4 = x[22:24]*(dscale/tscale)

    # the current time is (epc + t*dscale)
    return [v1 / (dscale/tscale);
            accel_perturbations(epc+t*tscale,[r1;v1]) / (dscale/tscale^2);
            v2 / (dscale/tscale);
            accel_perturbations(epc+t*tscale,[r2;v2]) / (dscale/tscale^2);
            v3 / (dscale/tscale);
            accel_perturbations(epc+t*tscale,[r3;v3]) / (dscale/tscale^2);
            v4 / (dscale/tscale);
            accel_perturbations(epc+t*tscale,[r4;v4]) / (dscale/tscale^2)]
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

function generate_data(x0,T,dt,Q,R)
    """Forward rollout of the sim, no process noise"""

    X = fill(zeros(nx),T)
    Y = fill(zeros(m),T)


    X[1] = copy(x0)

    u = SVector(0,0,0)
    for i = 1:(T-1)
        t = (i-1)*dt
        # @infiltrate
        # error()
        # @infiltrate
        # error()
        X[i+1] = rk4(dynamics,u,X[i],dt,t) #+ sqrt(Q)*randn(nx)
        Y[i] = measurement(X[i],t) #+ sqrt(R)*randn(m)
    end
    Y[T] = measurement(X[T],(T-1)*dt) #+  sqrt(R)*randn(m)

    return X,Y
end


# problem size
nx = 24
nu = 1
m = 12

# sample time is 30 seconds
dt = 30.0/tscale

# number of knot points
T = 500

# initial time for sim
epc = Epoch(2019, 1, 1, 12, 0, 0, 0.0)

# Declare initial state in terms of osculating orbital elements
oe0  = [R_EARTH + 550e3, 0.01, 45.0, 12, 25, 0]
eci0 = sOSCtoCART(oe0, use_degrees=true)
eci1 = eci0 + [10*normalize(randn(3));.1*normalize(randn(3))]
eci2 = eci0 + [10*normalize(randn(3));.1*normalize(randn(3))]
eci3 = eci0 + [10*normalize(randn(3));.1*normalize(randn(3))]
#scale everything
eci0[1:3] /= dscale
eci0[4:6] /= (dscale/tscale)
eci1[1:3] /= dscale
eci1[4:6] /= (dscale/tscale)
eci2[1:3] /= dscale
eci2[4:6] /= (dscale/tscale)
eci3[1:3] /= dscale
eci3[4:6] /= (dscale/tscale)
x0 = copy([eci0;eci1;eci2;eci3])

# run sim
Q = zeros(3,3)
R = zeros(3,3)
X,Y = generate_data(x0,T,dt,Q,R)

# new Q and R for gauss-newton stuff
Q = (1e-2)*1*Diagonal(@SVector ones(nx))
cholQ = sqrt(Q)
invcholQ = inv(cholQ)
R= (1e-2)*.1*Diagonal(@SVector ones(m))
cholR = sqrt(R)
invcholR = inv(cholR)

# indexing for x and y within the residual vector
idx_x = [(t-1)*nx .+ (1:nx) for t = 1:T]
idx_y = [(T-1)*nx + (t-1)*m .+ (1:m) for t = 1:T]

function residual(x)
    """residual vector for Gauss-Newton. rᵀr = MAP cost function"""
    u = SVector(0,0,0)
    r = zeros(eltype(x),nx*(T-1) + m*(T))
    for i = 1:(T-1)
        xt = @view x[idx_x[i]]
        xtp1 = @view x[idx_x[i+1]]
        t = (i-1)*dt
        # dynamics residual
        # @infiltrate
        # error()
        r[idx_x[i]] = invcholQ*(xtp1 - rk4(dynamics, u, xt, dt,t))
    end
    for i = 1:T
        xt = @view x[idx_x[i]]
        # sensor residual
        t = (i-1)*dt
        r[idx_y[i]] = invcholR*(Y[i] - measurement(xt,t))
    end

    return r
end

function sparse_jacobian!(J,x)
    """Modify sparse jacobian in place"""
    u = 0.0
    for i = 1:T
        if i < T
            k = idx_x[i]
            xt = @view x[k]
            t = (i-1)*dt
            A_fx(x) = rk4(dynamics, u, x, dt,t)
            J[k,k] = -invcholQ*FD.jacobian(A_fx,xt)
            J[k,k .+ nx] = invcholQ
            C_fx(x) = measurement(x,t)
            J[idx_y[i],k] = -invcholR*FD.jacobian(C_fx,xt)
        else
            k = idx_x[i]
            xt = @view x[k]
            t = (i-1)*dt
            _C_fx(x) = measurement(x,t)
            J[idx_y[i],k] = -invcholR*FD.jacobian(_C_fx,xt)
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

    # run Gauss-Newton for 30 iterations max
    for i = 1:60

        # ∂r/∂x
        sparse_jacobian!(J,x)

        # calculate residual at x
        r = residual(x)

        # solve for Gauss-Newton step (direct, or indirect)
        v = -J\r
        # lsqr!(v,-J,r)

        # calculate current cost
        S_k = dot(r,r)
        @show S_k

        # step size (learning rate)
        α = 1.0

        # run a simple line search
        for ii = 1:20
            x_new = x + α*v
            S_new= norm(residual(x_new))^2

            # this could be updated for strong frank-wolfe conditions
            if S_new < S_k
                x = copy(x_new)
                Ds = S_k - S_new
                @show ii
                break
            else
                α /= 2
            end
            if ii == 20
                error("line search failed")
                Ds = 0.0
                break
            end
        end

        # depending on problems caling, termination criteria should be updated
        if Ds < 1e-8
            break
        end

    end
    return x
end

# true x
x_real = vec(mat_from_vec(X))

# x1 = X[1]
# x2 = X[2]
# x3 = X[3]
# u = 0.0
# t = 0.0
# x2_t = rk4(dynamics,u,x1,dt,t)
#
# x3_t = rk4(dynamics,u,x2,dt,)
# # add noise with 50km meter σ to true solution
x_guess = x_real + (1/dscale)*randn(length(x_real))
x_gn = gauss_newton(x_guess)
#
#
X_gn = reshape(x_gn,24,:)

# this is our estimated state in the vector of vectors format
μs = vec_from_mat(X_gn)

position_error = zeros(4,length(X))
for i = 1:length(X)
    position_error[1,i] = norm(X[i][1:3] - μs[i][1:3])
    position_error[2,i] = norm(X[i][7:9] - μs[i][7:9])
    position_error[3,i] = norm(X[i][13:15] - μs[i][13:15])
    position_error[4,i] = norm(X[i][19:21] - μs[i][19:21])
end

dt *= tscale
t_vec = 0:dt:dt*(length(X)-1)
t_vec /= 3600
mat"
figure
hold on
title('Relative Position Errors')
plot($t_vec,$position_error'*$dscale)
legend('Target 1','Target 2','Target 3','Target 4')
xlabel('Time (hours)')
ylabel('Error (meters)')
hold off
%saveas(gcf,'relerror_gn.png')
"

# now get the relative distances in the chief RTN frame
r1 = [X[i][1:3]*dscale for i = 1:length(X)]
r2 = [X[i][7:9]*dscale for i = 1:length(X)]
r3 = [X[i][13:15]*dscale for i = 1:length(X)]
r4 = [X[i][19:21]*dscale for i = 1:length(X)]
v1 = [X[i][4:6]*(dscale/tscale) for i = 1:length(X)]
# v2 = [X[i][10:12]*(dscale/tscale) for i = 1:length(X)]
# v3 = [X[i][15:18]*(dscale/tscale) for i = 1:length(X)]

rp1 = fill(zeros(3),length(X))
rp2 = fill(zeros(3),length(X))
rp3 = fill(zeros(3),length(X))

for i = 1:length(X)
    ECI_Q_RTNchief = rRTNtoECI([r1[i];v1[i]])

    rp1[i] = ECI_Q_RTNchief'*(r2[i]-r1[i])
    rp2[i] = ECI_Q_RTNchief'*(r3[i]-r1[i])
    rp3[i] = ECI_Q_RTNchief'*(r4[i]-r1[i])
end

rp1 = mat_from_vec(rp1)/1000
rp2 = mat_from_vec(rp2)/1000
rp3 = mat_from_vec(rp3)/1000

mat"
figure
hold on
title('Relative Position from Chief')
plot3($rp1(1,:),$rp1(2,:),$rp1(3,:))
plot3($rp2(1,:),$rp2(2,:),$rp2(3,:))
plot3($rp3(1,:),$rp3(2,:),$rp3(3,:))
legend('Target 1','Target 2','Target 3')
xlabel('Chief R (km)')
ylabel('Chief T (km)')
zlabel('Chief N (km)')
view([23,39])
hold off
%saveas(gcf,'relp_gn.png')
"
