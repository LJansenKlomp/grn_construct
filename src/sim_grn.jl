using DifferentialEquations
using Plots
using MAT
using Statistics

pth = "C:/Data/phdCode/data/"
tr = 0.1
net = "chondrogenic"

netfile = matread(pth*net*"_GENIEnet.mat")
nmfile = matread(pth*net*"_GENIEgenes.mat")
datfile = matread(pth*net*"_GENIEdata.mat")
A = float(netfile["x"] .> tr)
mat = float(datfile["m"])
nm = nmfile["g"]

c = cor(transpose(mat))

# ODE function for model using Hill functions
function grn_ode(du,u,p,t)
    a,b,n,Aa,Ai = p
    actH(x) = a.*x.^n./(b.+x.^n)
    inhH(x) = a.*b./(b.+x.^n)

    du[:] = -u + transpose(Aa)*actH.(u)+transpose(Ai)*inhH(u)
end

# number of nodes
N = length(A[1,:])

# Hill coefficient and other parameters
n = 3
a = 0.37
S = 0.5

# adjacency matrices
Aa = zeros(N,N)
Ai = zeros(N,N)

# if corr(A,B) > 0 then the connection is assumed to be activating,
# otherwise it is assumed to be inhibotory
for i in 1:N
    for j in 1:N
        if A[i,j] > 0
            if c[i,j] > 0
                Aa[i,j] = 1.0
            else
                Ai[i,j] = 1.0
            end
        end
    end
end


# creating ODE problem
p = (a,S^n,n,Aa,Ai)
u0 = 3*rand(N)
tspan = (0.0, 10.0)
prob = ODEProblem(grn_ode,u0,tspan,p)

# solving ODE problem
sol = solve(prob)

plot(sol)