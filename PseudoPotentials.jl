using LinearAlgebra
using PyPlot


function Pot(g²,G)
    τ = [1 1 1]/8
    Ry = 13.605703976   # eV
    
    # Silicon
    if (g² == 3)
        Vs = -0.211*Ry
        Va = 0
    elseif (g² == 8)
        Vs = 0.04*Ry
        Va = 0
    elseif (g² == 11)
        Vs = 0.08*Ry
        Va = 0
    else
        Vs = 0
        Va = 0
    end

    x = 2π*(G⋅τ)
    return Vs*cos(x) - im*Va*sin(x)
end

function L2(v)
    return sum(v.^2)
end

function createGs()
    G = []
    for h in -3:3
        for k in -3:3
            for l in -3:3
                push!(G,[h,k,l])
            end
        end
    end
    
    return G
end

function Hamiltonian(k)
    ħ = 1.054571817E-34     # J.s
    h = 6.62607015E-34      # J.s
    m = 9.1093837015E-31    # kg
    q = 1.602176634E-19     # C
    ao = 5.428E-10          # m

    t = (h^2/(2*m*ao^2))/q  # eV

    G = createGs()
    N = size(G)[1]

    basis = [-1 1 1;
             1 -1 1;
             1 1 -1]
    
    H = zeros(N,N)

    for i in 1:N
        for j in 1:N
            if (i == j)
                g = G[i,:][1]'
                realg = g*basis
                H[i,j] = t*L2(k' + realg)
            else
                gi = G[i,:][1]'
                gj = G[j,:][1]'
                Δg = gi - gj
                realΔg = Δg*basis
                g² = L2(realΔg)
                H[i,j] = Pot(g²,realΔg)
            end
        end
    end

    return H
end

function walk()
    # Points of high symmetry
    Γ = [0,     0,   0]
    Χ = [0,     0,   1]
    L = [1/2, 1/2, 1/2]
    U = [1/4, 1/4,   1]
    K = [3/4, 3/4,   0]
    W = [1,   1/2,   0]

    # Paths along Brillouin zone
    # - first find lenghts
    nΛ = norm(L - Γ)
    nΔ = norm(Γ - Χ)
    nΞ = norm(Χ - U)
    nΣ = norm(K - Γ)

    ni = maximum([nΛ,nΔ,nΞ,nΣ])

    nΛ = Int(floor(100*nΛ/ni))
    nΔ = Int(floor(100*nΔ/ni))
    nΞ = Int(floor(100*nΞ/ni))
    nΣ = Int(floor(100*nΣ/ni))

    Λ = [(1-α)*L + α*Γ for α in LinRange(0,1,nΛ)]
    Δ = [(1-α)*Γ + α*Χ for α in LinRange(0,1,nΔ)]
    Ξ = [(1-α)*Χ + α*U for α in LinRange(0,1,nΞ)]
    Σ = [(1-α)*K + α*Γ for α in LinRange(0,1,nΣ)]
        
    # Compute LCAO along path
    bands = []
    path = vcat(Λ,Δ,Ξ,Σ)
    for k in path
        H = Hamiltonian(k)
        ε = sort(eigvals(H))[1:8]
        push!(bands,ε)
    end

    return bands
end

# Plotting
b = walk()
plot(b)
ylabel("E(k) [eV]")
xlabel("k-path")
show()
