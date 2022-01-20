using ForwardDiff
using LinearAlgebra
using SatelliteDynamics

include(joinpath(@__DIR__,"scaled_dynamics.jl"))


# initial time for sim (global variable for convienence)
epc = Epoch(2019, 1, 1, 12, 0, 0, 0.0)

let
    # Declare initial state in terms of osculating orbital elements
    oe0  = [R_EARTH + 550e3, 0.01, 45.0, 12, 25, 0]
    eci0 = sOSCtoCART(oe0, use_degrees=true)

    # now we scale these units (scaling is in scaled_dynamics.jl)
    eci0[1:3] /= dscale
    eci0[4:6] /= (dscale/tscale)

    # let's get a continuous time jacobian
    x = eci0
    u = zeros(3)
    t = 4.6 # this is random
    A = ForwardDiff.jacobian(_x -> scaled_dynamics(_x,u,t),x)
    @show A

    # let's get a discrete time jacobian
    dt = 10/tscale # scaled version of 10 seconds
    Ad = ForwardDiff.jacobian(_x -> rk4(scaled_dynamics,u,_x,dt,t),x)
    @show Ad

end
