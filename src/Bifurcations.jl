include("Model.jl")


U₀=[ones(K)...,1.0];
pars = (pars...,α = 20)
tspan = (0.0,10000.0)
prob = ODEProblem(persisters!,U₀,tspan,pars);
sol = solve(prob, Rosenbrock23(autodiff=false),saveat=0.1);
len = length(0.05:0.01:2.0)
eqs = zeros(K,len)
tspan = (0.0,1000.0)

# Equilibrium changes with α
for (k,th) in enumerate(0.05:0.01:2.0)
    pars = merge(pars,(θ = th,))
    prob = ODEProblem(persisters!,U₀,tspan,pars);
    sol = solve(prob, Rosenbrock23(autodiff=false),saveat=0.1);
    eqs[:,k] .= sol[1:K,end]
end

plot(0.05:0.01:2.0,eqs[:,end],label = "R",ylims = (0.0,0.1))
plot!(xlabel = L"\theta")