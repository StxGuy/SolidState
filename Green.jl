using StatsBase
using Statistics
using LinearAlgebra
using PyPlot

#
# Quantum Transmission
#
function QTran()
    # Physical Constants
    ħ = 1.054571817E-34                         # J.s
    q = 1.602176634E-19                         # C
    mo = 9.1093837015E-31                       # kg
    iη = 1E-36*im                                # Differential

    # Parameters
    N = 200                                     # Number of points
    w = 10E-9                                   # Physical width of system
    mr = 1                                      # Effective mass

    Σ1 = zeros(ComplexF64,N,N)                  # Self-energies
    Σ2 = zeros(ComplexF64,N,N)
    Γ1 = zeros(ComplexF64,N,N)                  # Gamma matrices
    Γ2 = zeros(ComplexF64,N,N)

    # Derived parameters
    ao = w/N                                    # Lattice constant
    t  = 0.5*((ħ/ao^2)*(ħ/(mr*mo)))/q           # Hopping energy [eV]

    println("Hopping energy = ",t," eV")
    println("Lattice constant = ",ao/1E-10," A")
    
    # Potential
    xi = 30
    xa = 79
    U = zeros(N)
    U[xi:xa] .= 1.0
    Δ = (xa-xi+1)*ao
    
    println("Barrier = ",Δ/1E-9," nm")

    # Hamiltonian
    H = SymTridiagonal(2*t*ones(N) + U, -t*ones(N-1))
    EnergySpace = LinRange(0,5,200)
        
    # Main loop
    transmission = []
    for Energy in EnergySpace
        x = Energy/(2t)
        ka = acos(1.0 - x)
        
        # Self-energies
        Σ1[1,1] = -t*exp(im*ka)
        Σ2[N,N] = -t*exp(im*ka)

        # Gamma functions
        Γ1 = im*(Σ1 - Σ1')
        Γ2 = im*(Σ2 - Σ2')

        # Green's function
        G = inv(Energy*I - H - Σ1 - Σ2)

        # Transmission
        T = Γ1*G*Γ2*G'
        x = tr(real(T))
        
        push!(transmission,x)
    end
    
    return EnergySpace,transmission,Δ
end 

function theory(ε,a)
    ħ = 1.054571817E-34     # [J.s]
    mo = 9.1093837015E-31   # [kg]
    q = 1.602176634E-19     # [C]
    U = 1.0                 # [eV]
    
    T = []
    for εi in ε
        k1 = sqrt(Complex(2*mo*εi*q))/ħ
        k2 = sqrt(Complex(2*mo*(εi-U)*q))/ħ
        k3 = sqrt(Complex(2*mo*εi*q))/ħ
    
        Ti = 4*k1*k2*k2*k3/((k2*(k1+k3)*cos(k2*a))^2 + ((k1*k3+k2*k2)*sin(k2*a))^2)
        push!(T,real(Ti))
    end
       
    return T
end

#==============================================================#
#                             MAIN                             #
#==============================================================#
x,y,a = QTran()
z = theory(x,a)

# Plotting
plot(x,y)
plot(x,z)
xlabel("Energy [eV]")
ylabel("Transmission")
legend(["Green","Exact"])
show()





