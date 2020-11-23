# Practice file for getting the Gauss-Newton method working

using ForwardDiff
using LinearAlgebra 
using Plots
using Random, Distributions

function measurement(x) # g(x)
    b1 = [10;0]
    b2 = [0;10]
    b3 = [0;0]
    b4 = [10;10]
    p = [x[1];x[2]]
    return [norm(b1-p);norm(b2-p);norm(b3-p)]
end

function dynamics(x,u,t) # f(x)
    px = x[1]
    py = x[2]
    θ  = x[3]
    return [cos(θ);sin(θ);sin(t)]
end

function rk4(f, u, x_n, h, t_n)
    k1 = h*f(x_n,u,t_n)
    k2 = h*f(x_n+k1/2, u, t_n + h/2)
    k3 = h*f(x_n+k2/2, u, t_n + h/2)
    k4 = h*f(x_n+k3, u, t_n + h)
    return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))
end

function getResidual(x)
    r = zeros(eltype(x), (2*numSteps-1),3)
    t = 0
    # Process noise section 
    for i = 1:(numSteps-1)
        t = (i-1) * 0.05
        r[i,:] = invQ * ((x[i+1,:] - rk4(dynamics, 0.0, x[i,:], 0.05, t)))
    end
    # Measurement noise section 
    for i = (numSteps):(2*numSteps-1)
        j = i - numSteps + 1
        r[i,:] = invR * (y[j,:] - measurement(x[j,:]))
    end
    return reshape(r, 3*(2*numSteps-1))
end


numSteps = 300
dt = 0.2
u = 0.0
Q = 1.0e-1 * [1 0 0; 0 1 0; 0 0 1]
R = 1.0e-2 * [1 0 0; 0 1 0; 0 0 1]
invR = inv(R)
invQ = inv(Q)


## GENERATE DATA
x_truth = zeros(Float64, numSteps,3)
x_hat = zeros(Float64, numSteps,3)
y = zeros(Float64, numSteps,3)

for i = 1:(numSteps-1)
    t = dt * (i - 1)
    x_truth[i+1,:] = rk4(dynamics, u, x_truth[i,:], dt, t) 
    x_hat[i+1,:] = x_truth[i+1,:] + sqrt(Q) * randn(3)
    y[i,:] = measurement(x_truth[i,:]) + sqrt(R)*randn(3)
end
y[numSteps,:] = measurement(x_truth[numSteps,:]) + sqrt(R)*randn(3)

initialPlot = plot(x_hat[:,1], x_hat[:,2], seriestype = :scatter, title = "Initial Estimate", label = "Estimated")
initialPlot = plot!(x_truth[:,1], x_truth[:,2], label = "Truth", xlim = [0,6], ylim = [0, 6])


## NON-LINEAR LEAST SQUARES (GAUSS NEWTON METHOD)
for i = 1:1500
    J = ForwardDiff.jacobian(getResidual, x_hat) 
    step = inv(transpose(J)*J)*transpose(J) * getResidual(x_hat)
    step = reshape(step, Int(size(step)[1]/size(x_hat)[2]), size(x_hat)[2]) # convert back to the right shape
    global x_hat = x_hat .- step 
end

finalPlot = plot(x_hat[:,1], x_hat[:,2], seriestype = :scatter, title = "Optimized Plot", label = "Estimated")
finalPlot = plot!(x_truth[:,1], x_truth[:,2], label = "Truth", xlim = [0,6], ylim = [0, 6])
display(plot(initialPlot, finalPlot, layout = (2,1)))

