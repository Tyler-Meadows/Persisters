using OrdinaryDiffEq, Parameters
using Plots, LaTeXStrings
using Revise
theme(:wong,
     size = (600,300),
     fontfamily = "computer modern",
     lw = 2,
     #thickness_scaling = 2,
     tickfontsize = 10,
     guidefontsize = 14,
     legendfontsize = 12,
     )

## Parameters and model setup
K = 100
pars = (θ = 1.0,
        η = 0.3,
        d = 0.03, 
        b = 0.6,
        v = 1.0,
        μ = 0.4,
        α = 80,
        λ = 60,
        m = 1e-2
        )

# Movement matrices
M = zeros(K,K) # Diffusion
for i in 1:K
    for j in 1:K
        if abs(i-j) == 1
            M[i,j] = K^2 |> Float64
        end
    end
    M[i,i] = -sum(M[i,:])
end
## Impose Boundary Conditions
M[1,1] = M[end,end] *= 2
M[1,2] = M[end,end-1] *= 2 


## First choice of p
# p(x,y) = 1
P1 = zeros(K,K)
for i in 1:pars.α
    P1[:,i] .= 1.0
end
P1 = P1./sum(P1,dims=2)/K # Normalize sum over Y

## Second choice of p
# p(x,y) = exp(-(x-y)^2/σ^2)
P2 = zeros(K,K)
for i in 1:K
    for j in 1:pars.α
        P2[i,j] = exp(-(i-j)^2/0.3^2)
    end
end
P2 = P2./ sum(P2, dims = 2)/K


v(x) = (x-1)*(pars.λ-x)*(K-x)/K^3 # Advection
V = zeros(K,K)
for i in 1:K
    for j in 1:K
        if i-j == -1
            V[j,i] = v(i)*K
        end  
    end
    V[i,i] = -sum(V[:,i])
end
pars = merge(pars,(p = P1,))

function persisters!(du,u,pars,t)
    @unpack θ,η,d,b,α,μ,v,p,m = pars
    R = u[end]
    N = u[1:K]
    χ = zeros(K,K)
    for i in 1:α
        χ[i,i] = 1.0
    end
    du[end] = θ - η*R - sum(b*R*N[1:α])
    du[1:K] .= (m*M+v*V)*N +(b*R*(1.0-μ))*χ*N - d*χ*N + (μ*R*b).*P*N
end

tspan = (0.0,100.0)


