using PyPlot
using LinearAlgebra

# Points of high symmetry
Γ = [0,   0,   0]    # 0,   0,   0
Χ = [0,   1,   0]    # 0,   1/2, 1/2
L = [1/2, 1/2, 1/2]  # 1/2, 1/2, 1/2
U = [1/4, 1,   1/4]  # 1/4, 5/8, 5/8
K = [3/4, 3/4, 0]    # 3/8, 3/4, 3/8
W = [1/2, 1,   0]    # 1/4, 3/4, 1/2

# Paths along Brillouin zone
# - first find lenghts
nΛ = norm(L - Γ)
nΔ = norm(Γ - Χ)
nS = norm(Χ - U)
nZ = norm(Χ - W)
nQ = norm(W - L)
nΣ = norm(U - Γ)
nΞ = norm(U - W)
nΕ = norm(W - K)

ni = maximum([nΛ,nΔ,nS,nZ,nQ,nΣ,nΞ,nΕ])

N = 100
nΛ = Int(floor(N*nΛ/ni))
nΔ = Int(floor(N*nΔ/ni))
nS = Int(floor(N*nS/ni))
nZ = Int(floor(N*nZ/ni))
nQ = Int(floor(N*nQ/ni))
nΣ = Int(floor(N*nΣ/ni))
nΞ = Int(floor(N*nΞ/ni))
nΕ = Int(floor(N*nΕ/ni))

Δ = [(1-α)*Γ + α*Χ for α in LinRange(0,1,nΔ)]
S = [(1-α)*Χ + α*U for α in LinRange(0,1,nS)]
Ξ = [(1-α)*U + α*W for α in LinRange(0,1,nΞ)]
Ε = [(1-α)*W + α*K for α in LinRange(0,1,nΕ)]
Σ = [(1-α)*K + α*Γ for α in LinRange(0,1,nΣ)]
Λ = [(1-α)*Γ + α*L for α in LinRange(0,1,nΛ)]
Z = [(1-α)*Χ + α*W for α in LinRange(0,1,nZ)]
Q = [(1-α)*W + α*L for α in LinRange(0,1,nQ)]

path = vcat(Δ,S,Σ,Λ)
    
ao = 0.543E-9         # [m]
kb = 1
ks = 0.25
m1 = 1
m2 = 1


Rs = [0.0 0.0 0.0;
      0.5 0.5 0.0;
      0.5 0.0 0.5;
      0.0 0.5 0.5]
a4 = [0.25, 0.25, 0.25]

Bands = []
for p in path
    k = 2π*p

    Dii = zeros(3,3)
    Dij = zeros(3,3)
    Dji = zeros(3,3)

    for R in eachrow(Rs)
        r = reshape(R-a4,3,1)
        r /= norm(r)
                            
        Dij += (kb*I + (kb-ks)*r*r')*exp(-im*(k⋅R))
        Dji += (kb*I + (kb-ks)*r*r')*exp(im*(k⋅R))
        Dii += -kb*I - (kb-ks)*r*r'
    end

    D = [Dii/m1 Dij;Dji Dii/m2]
    E = sqrt.(-eigvals(D))
    push!(Bands,E)
end    


# Plotting
plot(Bands)
xlabel("Crystal wavevector")
ylabel("Frequency ω [rad/s]")
show()
