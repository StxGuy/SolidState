using PyPlot
using LinearAlgebra

#------------------------------
# Fundamental Constants
#------------------------------
ħ = 1.054571817E-34             # Planck's constant [J.s]
q = 1.602176634E-19             # Electron charge [C]
mo = 9.1093837015E-31           # Rest electron mass [kg]
mr = 1                          # Effective mass []

w = 1E-6                        # Physical length [m]
N = 100                         # Number of discretization points []
Δ = w/N                         # Lattice constant [m]
t = ((ħ/Δ^2)*(ħ/(mr*mo)))/q     # Hopping energy [eV]

println("Hopping energy: ",t*1E3," meV")

#------------------------------
# Create Hamiltonian
#------------------------------
Ho = SymTridiagonal(fill(2t,N),fill(-t,N-1))
V = Diagonal([t*((j-N÷2-1)/(N÷2))^2 for j in 1:N])
H = Ho + V

#------------------------------
# Find Eigenstates
#------------------------------
ea = eigvals(H)    # Eigenvalues
ev = eigvecs(H)    # Eigenvectors

φ = abs.(ev[:,2]).^2

println("E[n+1] - E[n]")
for i in 1:5
    println(ea[i+1]-ea[i])
end

#------------------------------
# Plot
#------------------------------
plot(φ)
show()
