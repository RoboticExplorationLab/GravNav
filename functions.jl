function getOrbitAngle(rotationPeriod, mu = ((3.9860044188)*(10^5)),  Re = (6371), J2 = 0.0010826269, r = (550+6371))
    # Determine the angle required to get a desired J2 orbit 
    #   Note that default units are km and s, but this can be adjusted by passing in appropriate variables
    #   Returns theta (radians)

    rotRate = 2*pi/time
    theta = acos( (rotRate * -2.0 * r^(7/2) ) / (3.0*sqrt(mu)*J2*(Re^2)) ) # Radians

    return theta
end

function getPositionPlot(sol)
    # Positions
    pos = plot(sol[1,:], title = "Position", label = "x")
    pos = plot!(sol[2,:], label = "y")
    pos = plot!(sol[3,:], label = "z")
    pos = ylabel!("Position (100km)")
    return pos
end

function getVelocityPlot(sol)
    # Velocities
    vel = plot(sol[4,:], title = "Velocity", label = "Vx")
    vel = plot!(sol[5,:], label = "Vy")
    vel = plot!(sol[6,:], label = "Vz")
    vel = ylabel!("Velocity (100km/hour)")
    return vel
end

function getEnergyPlot(sol)
    # Get relative total energy of the system
    radius = (sol[1,:].^2 + sol[2,:].^2 + sol[3,:].^2).^(1/2)
    radPlot = plot(radius, title = "Radius", label = false)

    velocities = (sol[4,:].^2 + sol[5,:].^2 + sol[6,:].^2).^(1/2)

    E_k = 0.5* (sol[4,:].^2 + sol[5,:].^2 + sol[6,:].^2)
    E_p = (mu / radius)[1,:]
    E_t = E_k + E_p
    E_t = E_t ./ E_t[1] # Normalize 

    energy = plot(E_t, title = "Total Energy", ylabel = "Relative energy (%)")
    return energy
end

function getOrbit(xdot, x, w, t)
    #= Simulates an orbit. Includes J2
    x = [px, py, pz, vx, vy, vz]
    xdot = [vx, vy, vz, ax, ay, az]
    =#
    mu, Re, r, J2 = w
    xdot[1:3] = x[4:6]

    normr = norm(x[1:3])

    # Acceleration due to gravity
    a = ((-mu)/(normr^3)) * x[1:3]

    # Acceleration due to J2
    a_j2 = -((3.0 * J2 * mu * (Re^2)) / (2 * (normr^5))) * [(1 - 5 * ((x[3]/normr)^2))*x[1]; (1 - 5 * ((x[3]/normr)^2))*x[2]; (3 - 5 * ((x[3]/normr)^2))*x[3]]

    a_tot = a +  a_j2
    xdot[4:6] = a_tot 
end
