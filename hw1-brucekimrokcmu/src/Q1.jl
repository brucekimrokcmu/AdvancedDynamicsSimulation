import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()
# NOTE: if GR_jll fails to precompile, don't worry

using LinearAlgebra
using ForwardDiff
using OrdinaryDiffEq
using Test
using Plots
plotly()

#Implement this function in spherical coordinates
function gravitational_potential(s)
    # input: position in spherical coordinates 
    # s = [r, θ, ϕ]
    # output: gravitational potential
    
    #Constants
    μ = 398600.44 #km^3/sec^2
    J2 = 1.75553e10 #km^5/sec^2
    
    # unpack input
    r = s[1]
    θ = s[2]
    ϕ = s[3]
    
    # TODO: input the gravitational potential
    U = 0.0
    m = 1.0
    U = m*(-μ/r + J2*(3sin(θ)^2 -1)/2r^3)
#     print(U)
    return U
end


@testset "potential" begin 
    s1 = [7e6;0.6;0.2]
    @test isapprox(gravitational_potential(s1),-0.05694292000111414,atol=1e-10)
    s2 = [7e6;0.6;0.8]
    @test isapprox(gravitational_potential(s2),-0.05694292000111414,atol=1e-10)
end

# conversion from cartesian coordinates to spherical coordinates
function cartesian_to_spherical(x)
    r = sqrt(x[1:3]'*x[1:3])
    θ = atan(x[3],sqrt(x[1:2]'*x[1:2]))
    ϕ = atan(x[2],x[1])
    
    return [r; θ; ϕ]
end

#Implement this function in cartesian coordinates using the potential, the conversion between cartesian and spherical
#coordinates above, and the chain rule. You can write out the derivatives analytically or use ForwardDiff. You can
#also call ForwardDiff on the composite of multiple functions in the following way:
# if y = f(g(x))
# ∇fₓ = ForwardDiff.gradient(_x -> f(g(_x)), xk) # evaluated at xk

function gravitational_acceleration(x)
    # input: position in cartesian coordiantes 
    # output: acceleration in cartesian coordiantes 
    
    # TODO: output the gravitational acceleration 
    a = zeros(3)
    
    s = cartesian_to_spherical(x)
    U = gravitational_potential(s)

    ∇fₓ = -ForwardDiff.gradient(x -> gravitational_potential(cartesian_to_spherical(x)), x) 
    a = ∇fₓ
#     print(a)
    return a 
end

@testset "acceleration" begin 
    x1 = [6e3,5e3,-1.2e3]
    a1 = [-0.004851778939177607,-0.004043149115981339,0.0009724072028047602]
    if isapprox(gravitational_acceleration(x1),-a1,atol=1e-10)
        error("you are just missing a negative sign")
    end
    @test isapprox(gravitational_acceleration(x1),a1,atol=1e-10)
    
    x2 = [-4e3,8e3,-1.2e3]
    a2 = [0.0021710116269053233,-0.0043420232538106465,0.0006523593723992336]
    @test isapprox(gravitational_acceleration(x2),a2,atol=1e-10)
end

function orbit_dynamics(x)
    q = x[1:3]
    v = x[4:6]
    
    # TODO: put the cartesian equations of motion here 
    ẋ = zeros(length(x))
    a = gravitational_acceleration(q)
    ẋ = [v; a]
   
    return ẋ
end


@testset "ode" begin 
    x1 = [1.0e3,-4e3,2e3,10,20,30.0]
    ẋ1 = [10.0,20.0,30.0,-0.004142608441419876, 0.016570433765679505,-0.008337337706356213]
    @test isapprox(orbit_dynamics(x1),ẋ1,atol = 1e-10)
end

#Initial conditions for ISS orbit
r0 = [6791.0; 0; 0] #km
v0 = [0; cosd(51.5)*7.66; sind(51.5)*7.66] #km/s
x0 = [r0; v0]

#Ground Truth Simulation

Tf = 24*3600 #one day
tspan = (0.0, Tf)

#DifferentialEquations.jl needs the ODE to be defined with the following signature
function orbit_dynamics!(ẋ,x,p,t)
    # this is editing the input (ẋ) in place 
    ẋ .= orbit_dynamics(x) 
end

# this is how we setup ODEs in Julia
prob = ODEProblem(orbit_dynamics!,x0,tspan)

# TODO: adjust abstol and reltol to get the desired energy performance
soln = solve(prob,Tsit5();abstol=1e-17,reltol=1e-17)

plot(soln,vars=(1,2,3),label="")

#Implement a function to calculate the total energy of the satellite
function energy(x)
    # this should take in the cartesian state and return the total energy
    # kinetic energy + gravitational potential
    
    # TODO: output the total energy 
    E = 0.0 
    m = 1.0 

    s = cartesian_to_spherical(x)
    U = gravitational_potential(s)
    v = x[4:6]
    T = 0.5m*v'*v
   
    E = U + T


    return E
end

@testset "energy" begin 
    x1 = [1.0e3,-4e3,2e3,10,20,30.0]
    E1 = 612.9791623184244
    @test isapprox(energy(x1),E1,atol=1e-10)
end


#Plot the energy along the trajectory
Etraj = zeros(length(soln.t))
for k = 1:length(soln.t)
    Etraj[k] = energy(soln.u[k])
end

plot(soln.t,Etraj,xlabel="Time (sec)",ylabel="Energy (Jules)",lw=2)

maximum(Etraj)-minimum(Etraj) #This should be < 1e-12

@testset "energy performance" begin 
    @test maximum(Etraj)-minimum(Etraj) < 1e-12
end

#Implement the classic RK4 integrator
function rk4_step(f,xk,h)
    
    #TODO: implement rk4
    xn = zeros(length(xk))
    
    k1 = f(xk)
    k2 = f(xk + h/2*k1)
    k3 = f(xk + h/2*k2)
    k4 = f(xk + h*k3)

    xn = xk + (h/6).*(k1 + 2k2 + 2k3 + k4)
    
    return xn
end

#Implement the explicit Euler integrator
function explicit_euler_step(f,xk,h)    
    #TODO: implement explicit Euler 
    xn = zeros(length(xk))
    xn = xk + h*f(xk)
        
    return xn 
end

#Implement the explicit midpoint integrator
function explicit_midpoint_step(f,xk,h)
    
    # TODO: implement explicit midpoint 
    xn = zeros(length(xk))
    xn = xk + h*f(xk + 0.5*h*f(xk))
    
    return xn 
end

#Implement the implicit midpoint integrator
#You can use ForwardDiff to compute Jacobians if you want
function implicit_midpoint_step(f,xk,h)
    # TODO: implement implicit midpoint 
    # use maximum(abs.(residual)) < 1e-12 as termination criteria
    xn = zeros(length(xk)) 
    xn = xk
    residual = xn - xk - h*f(0.5xk + 0.5xn)       
    while maximum(abs.(residual)) >= 1e-12

        Δxn = -inv(ForwardDiff.jacobian(x -> x - xk - h*f(0.5xk + 0.5xn), xn))*residual
        xn += Δxn
        residual = xn - xk - h*f(0.5xk + 0.5xn) 
    end    

    return xn
end

h_test = [1.0; 3.0; 10.0; 30.0; 100.0] #time steps to test

# allocate vectors to hold errors for each h
error_rk4 = zeros(length(h_test))
error_mid1 = zeros(length(h_test))
error_mid2 = zeros(length(h_test))
error_eul = zeros(length(h_test))

# iterate through the step sizes 
for j = 1:length(h_test)
    
    # get the current step size (h), and number of simulation steps (N)
    h = h_test[j]
    N = Int(floor(Tf/h + 1))
    thist = h.*Array(0:(N-1))
    
    # create empty arrays 
    xtraj_rk4 = zeros(6,N)
    xtraj_rk4[:,1] = x0;
    xtraj_mid1 = zeros(6,N)
    xtraj_mid1[:,1] = x0;
    xtraj_mid2 = zeros(6,N)
    xtraj_mid2[:,1] = x0;
    xtraj_eul = zeros(6,N)
    xtraj_eul[:,1] = x0;
    
    # TODO: simulate the orbit using the given time step, for all 4 integrators
    for k = 1:(N-1)
        xtraj_rk4[:,k+1] .= rk4_step(orbit_dynamics,xtraj_rk4[:,k],h)
        xtraj_mid1[:,k+1] .= explicit_midpoint_step(orbit_dynamics,xtraj_mid1[:,k],h)
        xtraj_mid2[:,k+1] .= implicit_midpoint_step(orbit_dynamics,xtraj_mid2[:,k],h)
        xtraj_eul[:,k+1] .= explicit_euler_step(orbit_dynamics,xtraj_eul[:,k],h)    
    end
    
    # TODO: calculate the error (norm) between the final state (xtraj_""[:,N]), and the
    # final true solution (soln(Tf))

    error_rk4[j] = norm(xtraj_rk4[:,N] - soln(Tf))
    error_mid1[j] = norm(xtraj_mid1[:,N] - soln(Tf))
    error_mid2[j] = norm(xtraj_mid2[:,N] - soln(Tf))
    error_eul[j] = norm(xtraj_eul[:,N] - soln(Tf))
    
end


#Plot error on a log-log scale
plot(h_test,error_eul, xlabel="Step Size (sec)", ylabel="Error", xaxis=:log, yaxis=:log, label="Euler", lw=2)
plot!(h_test,error_mid1, xaxis=:log, yaxis=:log, label="Explicit Mid", lw=2)
plot!(h_test,error_mid2, xaxis=:log, yaxis=:log, label="Implicit Mid", lw=2)
plot!(h_test,error_rk4, xaxis=:log, yaxis=:log, label="RK4", lw=2)

# answer here
# 1) About 1 or 2 steps are needed when using midpoint integrators to achieve the accuracy(=error of 1) RK4 achieves.
# 2) while RK calls dynamics 4 times, and its step size is around 40.
# Explicit Mid calls dynamics 2 times, and its step size is around 1.5
# In terms of step size, RK is 26.66 times larger, but it calculates twice more.
# Therefore RK is about 13.33 times faster.


@testset "integrators" begin 
    @test 2000 < error_eul[1] < 2100
    @test 30000 < error_eul[end] < 40000
    
    @test 0 < error_mid1[1] < 0.6
    @test 10e3< error_mid1[end] < 12e3
    
    @test 0 < error_mid2[1] < 100
    @test 1 < error_mid2[end] < 1e4
    
    @test 0 < error_rk4[1] < 1e-5
    @test 1< error_rk4[end] < 50
end
