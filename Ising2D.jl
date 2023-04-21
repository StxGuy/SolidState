using PyPlot
using Statistics

# Find neighbors imposing periodic
# boundary conditions.
function neighborhood(M,idy,idx)
    Ny,Nx = size(M)
    
    if (idx > 1)
        Nl = M[idy,idx-1]
    else
        Nl = M[idy,Ny]
    end
    if (idx < Nx)
        Nr = M[idy,idx+1]
    else
        Nr = M[idy,1]
    end
    if (idy > 1)
        Nu = M[idy-1,idx]
    else
        Nu = M[idy,Nx]
    end
    if (idy < Ny)
        Nd = M[idy+1,idx]
    else
        Nd = M[idy,1]
    end
    
    return Nl + Nr + Nu + Nd
end

# Solve Ising model using Monte Carlo
function Ising(N,β,B)
    s = rand([-1,1],N,N)
    Energy = []
    Magnetization = []
    J = 1.0
        
    # Initial magnetization
    M = sum(s)
        
    # Initial energy
    ε = B*M
    for j in 1:N
        for i in 1:N
            z = neighborhood(s,i,j)
            ε += -0.5*J*z
        end
    end
    
    NTerm = 2E7
    for step in 1:NTerm
        # Pick a site
        idx = rand(1:N)
        idy = rand(1:N)

        # Find change in energy due to spin flipping
        # and perform Metropolis update.
        z = neighborhood(s,idy,idx)
        Δε = 2*s[idy,idx]*(J*z + B)

        if (rand() < exp(-β*Δε))
            s[idy,idx] = -s[idy,idx]
            ε += Δε
            M += 2*s[idy,idx]
        end
        
        # Save statistics
        if (mod(step,10) == 0)
            push!(Energy,ε)
            push!(Magnetization,M)
        end
    end

    return Energy,Magnetization
end

# Collect statistics
function stats(N,β,B)
    εₙ,Mₙ = Ising(N,β,B)
    
    M = abs(mean(Mₙ))/N^2
    χ = β*var(Mₙ)/N^2
    
    return M,χ
end


# Plot
Tsp = LinRange(1,4,100)
χsp = []
Msp = []
for T in Tsp
    β = 1.0/T
    M,χ = stats(20,β,0.0)
    push!(χsp,χ)
    push!(Msp,M)
end

plot(Tsp,Msp,"*")
ax1 = gca()
ax2 = ax1.twinx()
ax2.semilogy(Tsp,χsp,"s",color="orange")
show()
