using LinearAlgebra
using PyPlot


mutable struct Quantum
    ħ   :: Float64          # Reduced Planck constant [J.s]
    q   :: Float64          # Fundamental electron charge [C]
    mo  :: Float64          # Electron rest mass [kg]
    iη  :: ComplexF64
    t   :: Float64          # Hopping energy
    N   :: Int              # Number of lattice points
    W   :: Float64          # Width of device [m]
    U   :: Vector{Float64}  # Potential [eV]
    H   :: Matrix{Float64}  # Hamiltonian matrix
    ao  :: Float64          # Lattice constant [m]
end

function Quantum(W,U)
    ħ = 1.054571817E-34                         # J.s
    q = 1.602176634E-19                         # C
    mo = 9.1093837015E-31                       # kg
    iη = 1E-36*im                               
    
    N = length(U)
    ao = W/N
    t  = 0.5*((ħ/ao^2)*(ħ/mo))/q

    H = SymTridiagonal(2*t*ones(N) + U, -t*ones(N-1))
    
    println("Hopping energy = ",t," eV")
    println("Lattice constant = ",ao/1E-10," A")
        
    return Quantum(ħ,q,mo,iη,t,N,W,U,H,ao)
end

function getPotential(N :: Integer, a :: Integer)
    xi = (N-a)÷2 
    xa = (N+a)÷2
    
    U = zeros(N)
    U[xi:xa] .= 1.0

    return U
end    

function Σ(Q :: Quantum, n :: Int, ε :: Float64)
    Σf = zeros(ComplexF64,Q.N,Q.N)
       
    x = ε/(2Q.t)
    ka = acos(1.0 - x)
    
    Σf[n,n] = -Q.t*exp(im*ka)
    
    return Σf
end

function Γ(Σ :: Matrix{ComplexF64})
    return im*(Σ - Σ')
end

function Green(Q :: Quantum, ε :: Float64, Σ1 :: Matrix{ComplexF64}, Σ2 :: Matrix{ComplexF64})
    return inv(ε*I - Q.H - Σ1 - Σ2)
end

function transmission(Q :: Quantum, Γ1 :: Matrix{ComplexF64}, Γ2 :: Matrix{ComplexF64}, G :: Matrix{ComplexF64})
    T = Γ1*G*Γ2*G'
    x = tr(real(T))

    return x
end

function SpectralFunction(G)
    return im*(G-G')
end

# Quantum Transmission
function QTran(Q :: Quantum, EnergySpace :: LinRange{Float64,Int})
    # Main loop
    tran = []
    for ε in EnergySpace
        # Self-energies
        Σ1 = Σ(Q,1,ε)
        Σ2 = Σ(Q,Q.N,ε)

        # Gamma functions
        Γ1 = Γ(Σ1)
        Γ2 = Γ(Σ2)

        # Green's function
        G = Green(Q,ε,Σ1,Σ2)

        # Transmission
        push!(tran,transmission(Q,Γ1,Γ2,G))
    end
    
    return tran
end 

function LDOS(Q :: Quantum, ε :: Float64)
    Σ1 = Σ(Q,1,ε)
    Σ2 = Σ(Q,Q.N,ε)
        
    # Green's function
    G = Green(Q,ε,Σ1,Σ2)
    A = SpectralFunction(G)
    return real(diag(Diagonal(A)))
end

function DOS(Q :: Quantum, EnergySpace :: LinRange{Float64,Int})
    D = []
    for ε in EnergySpace
        Σ1 = Σ(Q,1,ε)
        Σ2 = Σ(Q,Q.N,ε)
    
        # Green's function
        G = Green(Q,ε,Σ1,Σ2)
        A = SpectralFunction(G)
        d = tr(real(A))
        
        push!(D,d)
    end
    
    return D
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
U = getPotential(200,40)
qSys = Quantum(10E-9,U)

if (true)
    EnergySpace = LinRange(0.0,5.0,200)
    t = QTran(qSys,EnergySpace)

    # Plotting
    plot(EnergySpace,t)
    xlabel("Energy [eV]")
    ylabel("Transmission")
    show()
end

if (false)
    y = LDOS(qSys,2.5)
    plot(y)
    xlabel("Position")
    ylabel("LDOS")
    show()
end

if (false)
    EnergySpace = LinRange(0.0,5.0,200)
    D = DOS(qSys,EnergySpace)
    
    # Plotting
    plot(EnergySpace,D)
    xlabel("Energy [eV]")
    ylabel("DOS")
    show()
end
