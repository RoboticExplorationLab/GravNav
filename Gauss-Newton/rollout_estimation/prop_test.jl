using SatelliteDynamics
using Test
using Printf

"""In this script, I see how big a position difference we can see when
comparing propagation of a full dynamics model with just full gravity.

spoiler:  it's extremely small"""

# Declare simulation initial Epoch
epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0)

# Declare initial state in terms of osculating orbital elements
oe0  = [R_EARTH + 530e3, 0.01, 97.0, 45.0, 30.0, 0.0]

# Convert osculating elements to Cartesean state
eci0 = sOSCtoCART(oe0, use_degrees=true)

# Set the propagation end time to one orbit period after the start
# T    = orbit_period(oe0[1])
T = 15*60 # 15 minutes
epcf = epc0 + T

# Initialize State Vector
orb  = EarthInertialState(epc0, eci0, dt=1.0,
        mass=100.0, n_grav=20, m_grav=20,
        drag=true, srp=true,
        moon=true, sun=true,
        relativity=true
)

# Propagate the orbit
t1, epc, eci_hist_full = sim!(orb, epcf)

# Initialize State Vector
orb  = EarthInertialState(epc0, eci0, dt=1.0,
        mass=100.0, n_grav=10, m_grav=10,
        drag=false, srp=false,
        moon=false, sun=false,
        relativity=false
)

# Propagate the orbit
t2, epc, eci_hist_kep = sim!(orb, epcf)

println(" ")
println("Position error (m)")
@printf "%.4E\n" norm(eci_hist_full[1:3,end] - eci_hist_kep[1:3,end])
println("Velocity error (m/s)")
@printf "%.4E\n" norm(eci_hist_full[4:6,end] - eci_hist_kep[4:6,end])
