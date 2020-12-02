# cd(@__DIR__)
# Pkg.activate(".")

using SatelliteDynamics, LinearAlgebra, MATLAB

# Declare simulation initial Epoch
epc0 = Epoch(2020, 5, 1, 1, 0, 0, 0.0)

# Declare initial state in terms of osculating orbital elements
oe0  = [R_EARTH + 550*1e3, 0.01,
        45.6, 0, 0, 110]

# Convert osculating elements to Cartesean state
eci0 = sOSCtoCART(oe0, use_degrees=true)

# Set the propagation end time to one orbit period after the start
T    = 1*orbit_period(oe0[1])
epcf = epc0 + T

dt = 1.0
# Create an EarthInertialState orbit propagagator
orb  = EarthInertialState(epc0, eci0, dt=dt,
            area_drag = 0.01, coef_drag = 2.0, area_srp = 0.01, coef_srp = 2.0,
            mass=1.0, n_grav=2, m_grav=2,
            drag=false, srp=false,
            moon=false, sun=false,
            relativity=false
)

# Propagate the orbit
t, epc, eci_1 = sim!(orb, epcf)

# Declare initial state in terms of osculating orbital elements
oe0  = [R_EARTH + 550*1e3, 0.03,
        45.6, 0, 0, 112]

# Convert osculating elements to Cartesean state
eci0 = sOSCtoCART(oe0, use_degrees=true)

# Set the propagation end time to one orbit period after the start
T    = 1*orbit_period(oe0[1])
epcf = epc0 + T

dt = 1.0
# Create an EarthInertialState orbit propagagator
orb  = EarthInertialState(epc0, eci0, dt=dt,
            area_drag = 0.01, coef_drag = 2.0, area_srp = 0.01, coef_srp = 2.0,
            mass=1.0, n_grav=2, m_grav=2,
            drag=false, srp=false,
            moon=false, sun=false,
            relativity=false
)

# Propagate the orbit
t, epc, eci_2 = sim!(orb, epcf)


# relative ranging stuff

rel_r = [norm([eci_1[1:3,i] - eci_2[1:3,i] ]) for i = 1:size(eci_1,2)]

function relative_velocity(rv1,rv2)

        x1,y1,z1,dx1,dy1,dz1 = rv1
        x2,y2,z2,dx2,dy2,dz2 = rv2

        relv = (dx1*x1 - dx1*x2 - dx2*x1 + dx2*x2 + dy1*y1 - dy1*y2 -
                dy2*y1 + dy2*y2 + dz1*z1 - dz1*z2 - dz2*z1 + dz2*z2)/((x1
                - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)^(1/2)

        return relv
end

rel_v = [relative_velocity(eci_1[:,i],eci_2[:,i]) for i = 1:size(eci_1,2)]

diffed_rel_v = [(rel_r[i+1] - rel_r[i])/dt for i = 1:(size(eci_1,2)-1)]

mat"
figure
hold on
plot($rel_r)
hold off
"

mat"
figure
hold on
plot($rel_v)
plot($diffed_rel_v)
legend('sym','diff')
hold off
"
