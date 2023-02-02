using PyPlot
using LinearAlgebra

# Elastic constants
K = 1
G = 1.5

# Lattice constant
a = 1

# Ion masses
m1 = 1.2
m2 = 2.3

# Crystal vectors
k_space = LinRange(-1,1,100)

# Main loop
Bands = []
for k in k_space
    ki = k*π/a
    
    # Dynamical Matrix
    D = [(K+G)/m1                       -(K+G*exp(-im*ki*a));
         -(K+G*exp(im*ki*a))             (K+G)/m2]
     
    # Eigenvalues
    E = eigvals(D)
    push!(Bands,E)
end

# Plot results
plot(k_space,Bands)
xlabel("Crystal momentum [π/a]")
ylabel("Eigenvalue")
show()
    

