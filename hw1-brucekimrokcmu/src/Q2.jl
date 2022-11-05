import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()

using LinearAlgebra
using ForwardDiff
using OrdinaryDiffEq
using Test
using Random
using Plots
plotly()

#Some standard functions for dealing with rotation matrices and quaternions from the class notes
function hat(ω)
    # TODO: implement the hat function 
    ω̂ = I(3)
    ω̂ = [0 -ω[3] ω[2]; ω[3] 0 -ω[1]; -ω[2] ω[1] 0]    
    return ω̂
end

function L(q)
    # TODO: implement L 
    _L = zeros(4,4)
    _L[1,:] = [q[1] -q[2] -q[3] -q[4]]
    _L[2:4,1] = q[2:4]'
    _L[2:4,2:4]= q[1]*I(3)+hat(q[2:4])    
    
    return _L
end

function R(q)
    # TODO: implement R 
    _R = zeros(4,4)
    _R[1,:] = [q[1] -q[2] -q[3] -q[4]]
    _R[2:4,1] = q[2:4]'
    _R[2:4,2:4]= q[1]*I(3)-hat(q[2:4])    
    
    return _R
end

# TODO: implement H function 
H = zeros(4,3)
H = [0 0 0; I(3)]

# add tests here 
@testset "rotation functions" begin 
    x = randn(3)
    y = randn(3)
    # test hat 
    @test isapprox(cross(x,y),hat(x)*y,atol=1e-12)
    @test isapprox(cross(y,x),hat(y)*x,atol=1e-12)
    @test isapprox(tr(hat(x)),0,atol=1e-12)
    @test isapprox(tr(hat(y)),0,atol=1e-12)
    @test isapprox(norm(hat(x) + hat(x)'),0,atol=1e-12)
    
    # test orthogonality of L and R 
    L1 = L(normalize(randn(4)))
    R1 = R(normalize(randn(4)))
    @test isapprox(norm(L1*L1'-I),0,atol=1e-10)
    @test isapprox(norm(R1*R1'-I),0,atol=1e-10)

    # make sure the rotation is correct, and H does what it's supposed to
    Random.seed!(1234)
    r = normalize(randn(3))
    θ = 0.3 
    ϕ = r*θ
    q = [cos(θ/2);r*sin(θ/2)]
    Q = exp(hat(ϕ)) 
    x = randn(3)
    x2 = Q*x
    # this tests H 
    x3 = H'*L(q)*R(q)'*H*x
    @test isapprox(x2,x3,atol = 1e-10)
end

#Inertia matrix we'll be using
J = Diagonal([1.0; 2.0; 3.0])

#Implement Euler's equation for a rigid body
function dynamics_ω(ω)
    
    # TODO: add Euler's equation 
    ω̇ = zeros(length(ω))
    τ = zeros(3)
    ω̇ = inv(J)*τ - inv(J)*hat(ω)*J*ω
#     print(ω̇ )
    return ω̇
end

@testset "Euler's" begin
    ω1 = [.2;7;.2]
    α1 = [-1.4;0.04;-0.46666666666666]
    @test isapprox(α1, dynamics_ω(ω1),atol=1e-10)
end

#Implement the full attitude dynamics using a rotation matrix
function dynamics_Q(x)
    # input: x = [vec(Q);ω]
    
    # TODO: implement 
    ẋ = zeros(length(x))
    ω = x[10:12]
    ω̂ = hat(ω)
    Q = reshape(x[1:9],(3,3))
    
    Q̇ = Q*ω̂ 
    ω̇ = dynamics_ω(ω)
    ẋ = [vec(Q̇);ω̇ ]
    
    return ẋ
end

function dynamics_Q!(ẋ,x,p,t)
    # make the ODE in-place so it works with OrdinaryDiffEq
    ẋ .= dynamics_Q(x)
    
end

@testset "attitude dynamics (DCM)" begin 
    Q1 = exp(hat([1;2;3]))
    ω1 = [0.2;0.3;0.7]
    x1 = [vec(Q1);ω1]
    ẋ1 = [0.47267683571087754,-0.4926072371847232,0.337512546219523,
          0.5043029621213005,0.3210433517191289, -0.41546322185318607,
         -0.35118036539795094,0.003154917030294274,0.08162351044578747,
         -0.2099999999999999,0.14,-0.02]
    
    @test isapprox(ẋ1,dynamics_Q(x1),atol=1e-10)
end

#Implement the full attitude dynamics using a quaternion

function dynamics_q(x)
    
    # TODO: implement attitude dynamics with quaternion 
    ẋ = zeros(length(x))
    q = x[1:4]
    ω = x[5:7]
    
    q̇ = 0.5L(q)*H*ω
    ω̇ = dynamics_ω(ω)
    
    ẋ = [q̇; ω̇ ]
#     print(ẋ)
    return ẋ
end

function dynamics_q!(ẋ,x,p,t)
    # make the ODE in-place so it works with OrdinaryDiffEq
    ẋ .= dynamics_q(x)
end

@testset "attitude dynamics (quaternion)" begin
    r1 = normalize([1;2;3])
    θ1 = deg2rad(45)
    q1 = [cos(θ1/2);r1*sin(θ1/2)]
    x1 = [q1;0.2;7;0.3]
    ẋ1 = [-0.7721871929186838;
         -0.9508318305595436;
          3.2489198311984837;
          0.47609421287426346;
         -2.0999999999999996;
          0.06000000000000001;
         -0.46666666666666673]
    @test isapprox(dynamics_q(x1),ẋ1,atol=1e-10)
    
end

#Initial Conditions
Q0 = Array(I(3))
q0 = [1.0; 0; 0; 0]
Random.seed!(1234)
ω0 = [0; 2*pi; 0] + 1e-2*randn(3) #small perturbation from intermediate-axis spin

x0Q = [Q0[:]; ω0]
x0q = [q0; ω0];

#Set up and solve with OrdinaryDiffEq.jl with default settings
Tf = 60.0
tspan = (0.0, Tf)

prob_q = ODEProblem(dynamics_q!, x0q, tspan)
prob_Q = ODEProblem(dynamics_Q!, x0Q, tspan)

sol_q = solve(prob_q,Tsit5());
sol_Q = solve(prob_Q,Tsit5());


@testset "integration of attitude dynamics" begin
    q1 = sol_q(0)[1:4]
    qf = sol_q(Tf)[1:4]
    norm(q1)
    @test isapprox(1.0,norm(q1),atol=1e-10)
    @test 1 <= norm(qf) <= 1.1
    
    h1 = norm(J*sol_q(0)[5:7])
    hf = norm(J*sol_q(Tf)[5:7])
    @test 12.53 <= h1 <= 12.6
    @test 12.53 <= hf <= 13.0
    
    Q1 = reshape(sol_Q(0)[1:9],3,3)
    Qf = reshape(sol_Q(Tf)[1:9],3,3)
    
    @test norm(Q1'*Q1 - I) < 1e-10
    @test 0 < norm(Qf'*Qf - I) < 0.2
    
    h1 = norm(J*sol_Q(0)[10:12])
    hf = norm(J*sol_Q(Tf)[10:12])
    @test 12.53 <= h1 <= 12.6
    @test 12.53 <= hf <= 13.0
    
end

plot(sol_q, vars=(5,6,7))

using MeshCat
using CoordinateTransformations
using Rotations
using GeometryBasics: HyperRectangle, Vec, Point, Mesh
using Colors: RGBA, RGB

# Create a new visualizer instance
vis = Visualizer()
render(vis)


# meshcat stuff
# delete!(vis)
t_vec = 1:0.01:60
q_hist = sol_q(t_vec).u
box = HyperRectangle(Vec(-0.1, -1, -0.1), Vec(0.2, 2, 0.2))
box2 = HyperRectangle(Vec(-0.1 - 1, -1, -0.1), Vec(2, 0.2, 0.2))
green_material = MeshPhongMaterial(color=RGBA(1, 1, 0.3, 0.5))
setobject!(vis["group1"]["greenbox"], box, green_material)
blue_material = MeshPhongMaterial(color=RGBA(0, 0, 1, 0.5))
setobject!(vis["group1"]["bluebox"], box2, blue_material)
settransform!(vis["group1"], Translation(0, 0, 0))
# animation loop
for i = 1:length(t_vec)
    settransform!(vis, LinearMap(UnitQuaternion(q_hist[i][1:4])  ))
    sleep(0.005)
end

#Plot the error in the quaternion norm vs. time
qnorm = zeros(length(sol_q.t))

for k = 1:length(sol_q.t)
    
    # TODO: calculate the norm of each quaternion  
#     norm = 0
#     for j = 1:4
#         norm += sol_q[k][j]*sol_q[k][j]
#         norm = sqrt(norm)
#     end
    
    qnorm[k] = norm(sol_q[k][1:4])

# q-1 ≜ q̄/|q|^2 = (q0−q1⋅i−q2⋅j−q3⋅k)/(q0^2+q1^2+q2^2+q3^2)
#q*q−1=1
    
#     q = sol_q[k][1:4]
#     q̄ = conj(q)
#     q1 = q̄/norm(q)^2
#     qnorm[k] = norm(q*q1') - 1


end

plot(sol_q.t,qnorm,xlabel="Time (sec)",ylabel="||q||",label="",lw=2)

#Plot the orthogonality error vs. time for the rotation matrix
Qerr = zeros(length(sol_Q.t))
for k = 1:length(sol_Q.t)
    Qk = reshape(sol_Q.u[k][1:9],3,3)
    
    # TODO: calculate the orthogonality error of the rotation matrix 
    Qerr[k] = norm(transpose(Qk)*Qk - I(3))

end

plot(sol_Q.t,Qerr,xlabel="Time (sec)",ylabel="Orthogonality Error",label="",lw=2)

function dynamics_q_stabilized(x)
    # TODO: implement stabilized quaternion attitude dynamics 
    # x = [q; ω]
    ẋ = zeros(length(x))
    
    q = x[1:4]
    ω = x[5:7]

    error = 1-norm(q)
    
    Kp = 12
    
    q̇ = 0.5L(q)*H*ω + Kp*error*(q) #/norm(q)) #why not divided by norm(q)
    ω̇ = dynamics_ω(ω)
    
    ẋ = [q̇; ω̇ ]
    
    return ẋ
end

function dynamics_q_stabilized!(ẋ,x,p,t)
    # make the ODE in-place so it works with OrdinaryDiffEq
    ẋ .= dynamics_q_stabilized(x)
end

prob_q2 = ODEProblem(dynamics_q_stabilized!, x0q, tspan)
sol_q2 = solve(prob_q2,Tsit5());

#Plot norm(q) with and without stabilization
qnorm2 = zeros(length(sol_q2.t))

for k = 1:length(sol_q2.t)
    
    # TODO: calculate quaternion norm 
    qnorm2[k] = norm(sol_q2[k][1:4]) 
    
end
@show abs(1-mean(qnorm2))

plot(sol_q.t,qnorm,xlabel="Time (sec)",ylabel="||q||",label="No Stabilization",lw=2)
plot!(sol_q2.t,qnorm2,label="Stabilization",lw=2)

@testset "stabilized attitude dynamics (quaternion)" begin
    r1 = normalize([1;2;3])
    θ1 = deg2rad(45)
    q1 = [cos(θ1/2);r1*sin(θ1/2)]
    ω1 = [0.2;7;0.3]
    ẋ1 = dynamics_q([q1;ω1])
    ẋ2 = dynamics_q([q1*1.1;ω1])    
    @test isapprox(ẋ1[5:7],ẋ2[5:7],atol=1e-12)
    @test norm(ẋ1 - dynamics_q_stabilized([q1*1.000000001;ω1]))<norm(ẋ1 - dynamics_q_stabilized([q1*1.01;ω1]))
    @test qnorm2[end]<qnorm[end]
end

function lie_midpoint_step_Q(xk,h)
    
    Qk = reshape(xk[1:9],3,3)
    ωk = xk[10:12]
    
    # TODO: implement the lie midpoint step 
    ẋ = zeros(length(xk))
    
    ωk̂ = hat(ωk)
    ωk̇ = dynamics_ω(ωk)
    ωm = ωk + 0.5*h*ωk̇
    ωm̂ = hat(ωm)
    ωṁ = dynamics_ω(ωm)
    ωn = ωk + h*ωṁ    
    ωṅ = ωn
#     ωn̂ = hat(ωn)
    
#     explicit midpoint  
#     Qn = Qk + h*f(Qk + 0.5*h*f(Qk))
#     Qm = Qk + 0.5*h*f(Qk)
#     Qm = Qk + 0.5*h*Qk̇
#     Qn = Qk + h*Qṁ    
    
    Qk̇ = Qk*ωk̂
    # exponential mapping was the answer.. to the lie algebra and then to map to exp for nonlinear space
    Qm = Qk*exp(0.5*h*ωk̂)
    Qṁ = Qm*ωm̂
    Qn = Qk*exp(h*ωm̂)
    Qṅ = Qn

    ẋ = [vec(Qṅ);ωṅ]

    return ẋ
end

@testset "lie midpoint (DCM)" begin 
    Q = exp(hat([1,2,3.0]))
    ω = [.2; 7; -.3]
    h = 0.1 
    ẋ = [-0.6027631402683421, -0.7420787949334725, 0.2932501642692329,
          0.6865916167045675, -0.2951187547149963, 0.6644523101671616,
         -0.40653234632050517, 0.6018504654194051, 0.6873917868796442,
          0.42623633333333333, 6.990138333333333, -0.37113616666666666]
    @test isapprox(lie_midpoint_step_Q([vec(Q);ω],h),ẋ,atol=1e-10)
end

#Now let's simulate the same example and plot the orthogonality error again

h = 0.1
N = N = Int(floor(Tf./h + 1))
thist = h.*Array(0:(N-1))

xtrajQ = zeros(12,N)
xtrajQ[:,1] .= x0Q

for k = 1:(N-1)
    xtrajQ[:,k+1] .= lie_midpoint_step_Q(xtrajQ[:,k], h)
end

Qerr2 = zeros(N)
for k = 1:N
    Qk = reshape(xtrajQ[1:9,k],3,3)
    
    # TODO: calculate the orthogonality error of Q
    Qerr2[k] = norm(transpose(Qk)*Qk-I(3))
    
end

plot(sol_Q.t,Qerr,xlabel="Time (sec)",ylabel="Orthogonality Error",label="RK",lw=2)
plot!(thist,Qerr2,label="Lie Group",lw=2)

@testset "Lie group DCM" begin 
    @test isapprox(maximum(Qerr2),0.0,atol = 1e-8)
    @test isapprox(minimum(Qerr2),0.0,atol = 1e-8)
end

function Expq(ϕ)
    #map from R^3 into SO(3)
    
    # TODO: implement the quaternion exponential map ϕ → q 
    q = zeros(4)
    θ = norm(ϕ)
    
    if θ == 0
        r = [1;1;1]
    else    
        r = ϕ/θ
    end    
    q = [cos(0.5*θ); r[1]*sin(0.5*θ);r[2]*sin(0.5*θ);r[3]*sin(0.5*θ)]
    
    return q
end

function lie_midpoint_step_q(xk,h)
    
    qk = xk[1:4]
    ωk = xk[5:7]
    
    # TODO: implement lie midpoint step 
    xn = zeros(length(xk))
    
#     ωk̂ = hat(ωk)
    ωk̇ = dynamics_ω(ωk)
    ωm = ωk + 0.5*h*ωk̇
#     ωm̂ = hat(ωm)
    ωṁ = dynamics_ω(ωm)
    ωn = ωk + h*ωṁ    

    
    qn = L(qk)*Expq(h*ωm)
    xn = [qn;ωn] 

#     print(xn)
    
    return xn

end


@testset "Expq" begin 
    q = Expq(zeros(3))
    @test isapprox(norm(q),1,atol = 1e-10)
    @test isapprox(q,[1;0;0;0],atol=1e-10)
    for i = 1:10 
        q = Expq(randn(3))
        @test isapprox(norm(q),1,atol = 1e-10)
    end
    
    # test if the rotation is correct
    Random.seed!(1234)
    ϕ = randn(3)
    q = Expq(ϕ)
    Q = exp(hat(ϕ))
    x = randn(3)
    x2 = Q*x
    x3 = H'*L(q)*R(q)'*H*x
    @test isapprox(x2,x3,atol = 1e-10)
end

@testset "lie midpoint (quaternion)" begin 
    q = Expq([1,2,3.0])
    ω = [.2; 7; -.3]
    h = 0.1 
    ẋ = [-0.4442718458039811, -0.035227218052985225, 0.39378058569262325,
          0.8039393139197423, 0.42623633333333333, 6.990138333333333,
         -0.37113616666666666]
    @test isapprox(lie_midpoint_step_q([q;ω],h),ẋ,atol=1e-10)
end

#Now let's simulate and plot the quaternion norm again

xtrajq = zeros(7,N)
xtrajq[:,1] .= x0q

for k = 1:(N-1)
    xtrajq[:,k+1] .= lie_midpoint_step_q(xtrajq[:,k], h)
end

qnorm3 = zeros(N)
for k = 1:N
    # TODO: calculate quaternion norm 
    qnorm3[k] = norm(xtrajq[:,k][1:4])
    
end

plot(sol_q.t,qnorm,label="No Stabilization",xlabel="Time (sec)",ylabel="||q||",lw=2)
plot!(sol_q2.t,qnorm2,label="Stabilization",lw=2)
plot!(thist,qnorm3,label="Lie Group",lw=2)

using Statistics
# test that lie group is perfect 
@testset "lie group stabilization" begin 
    @test mean(qnorm3) < mean(qnorm)
    @test mean(qnorm2) < mean(qnorm)
    @test abs(1-mean(qnorm3)) < 1e-6
    @test abs(1-mean(qnorm2)) < 1e-4
end


