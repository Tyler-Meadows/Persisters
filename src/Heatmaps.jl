include("Model.jl")

## Persisters only
n2 = 0.5*ones(K-pars.α); # Initial normal cells
n1 = zeros(pars.α);# Initial persisters
U₀=[n1...,n2...,0.1];


prob = ODEProblem(persisters!,U₀,tspan,pars);
sol = solve(prob, Rosenbrock23(autodiff=false),saveat=0.1);

N_sol = sol[1:K,:];
R_sol = sol[end,:];
x =range(0,1,K);

R_sol = sol[end,:]
x =range(0,1,K)
heat = heatmap(1:length(R_sol),range(0,1,K),N_sol,xlabel=L"t",ylabel=L"x")
savefig("Figures/Uniform_Kernel_heat_persisters.pdf")

## Normal Cells only
n2 = 1*zeros(K-pars.α) # Initial normal cells
n1 = 0.5*ones(pars.α)# Initial persisters
U₀=[n1...,n2...,0.1]


prob = ODEProblem(persisters!,U₀,tspan,pars)
sol = solve(prob, Rosenbrock23(autodiff=false),saveat=0.1)

N_sol = sol[1:K,:]
R_sol = sol[end,:]
x =range(0,1,K)
heat = heatmap(1:length(R_sol),range(0,1,K),N_sol,xlabel=L"t",ylabel=L"x")
savefig("Figures/Uniform_redistribution_heat_normal.pdf")
init = plot(x,N_sol[:,1],color=QueensBlue,label=L"n(x,0)")
plot!(x,N_sol[:,end], color = QueensGold,label=L"n(x,T)",xlabel = L"x")
savefig("Figures/Uniform_redistribution_normal.pdf")




## Local p
σ = 0.05
p(i,j) = exp(-(i-j)^2/K^2/σ^2)
P = zeros(K,K)
for i in 1:K
    for j in 1:pars.α
        P[i,j] = p(i,j)
    end
end
P = P./ sum(P, dims = 2)
pars = merge(pars, (p = P,v = 1.0))

prob = ODEProblem(persisters!,U₀,tspan,pars)
sol = solve(prob, Rosenbrock23(autodiff=false),saveat=1.0)

N_sol = sol[1:K,:]
R_sol = sol[end,:]
x =range(0,1,K)
heat = heatmap(1:length(R_sol),range(0,1,K),N_sol,xlabel=L"t",ylabel=L"x",labelfontsize=6)
init = plot(x,N_sol[:,1],color=QueensBlue,label=L"n(x,0)",legendfontsize=5)
plot!(x,N_sol[:,end], color = QueensGold,label=L"n(x,T)",xlabel = L"x", labelfontsize=6)
savefig("Figures/local_redistribution_persisters.pdf")
heat
savefig("Figures/local_redistribution_heat_persisters.pdf")



## How does the equilibrium distribution change as σ changes


σ_range = 10.0.^(-1.5:0.1:0)
eq = zeros(K,length(σ_range))

for (l,σ) in enumerate(σ_range)
    p(i,j) = exp(-(i-j)^2/K^2/σ^2)
    P = zeros(K,K)
    for i in 1:K
        for j in 1:pars.α
            P[i,j] = p(i,j)
        end
    end
    P = P./ sum(P, dims = 2)
    pars = merge(pars, (p = P,v=0.05))

    prob = ODEProblem(persisters!,U₀,tspan,pars)
    sol = solve(prob, Rosenbrock23(autodiff=false))
    eq[:,l] = sol[1:K,end]
    U₀ = [sol[1:K,end]...,sol[end,end]]
end
heatmap(σ_range, range(0,1,K),eq)
x = range(0,1,K)
plot(x,eq[:,1], label = "σ = $(σ_range[1])")
plot!(x,eq[:, 100], label = "σ = $(σ_range[100])")
plot!(x,eq[:,end], label = "σ = $(σ_range[end])")