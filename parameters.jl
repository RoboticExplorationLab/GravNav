### SYSTEM PARAMETERS ###
## Constants ##
_Re = 6378.1363*(10^(-2)) # Radius of the earth (100*km)
_J2 = 0.0010826269 # Unitless (normalized), taken care of with the root(mu)
_mu = (3600^2)*(3.9860044188)*(10^(-1)); # (100*km^3/hr^2) Standard Gravitational Parameter

### SYSTEM VARIABLES ###
_r0 = (550+6371)*(10^(-2)) # Distance from two center of masses (100*km) 
_a = _r0 # Semi-major axis of (elliptical) orbit (circular orbit if a = r0)
_v0 = sqrt(_mu*( (2/_r0) - (1/_a))) 
_time = 365. * 24. * 60. * 60. # Seconds per year
_theta = getOrbitAngle(_time)

# Sim Params 
orbitTime = ((2*pi*_r0)/(_v0))
day = 24.0 #* 60.0 * 60. # in seconds
tspan = (0, 5. * orbitTime)  # Should take ~5730 for one cycle
saveRate = orbitTime/30;
