# Si: [Ne] 3s² 3p²
# Ge: [Ar] 3d¹⁰ 4s² 4p²

# J. C. Slater and G. F. Koster,
# "Simplified LCAO Method for the Periodic Potential Problem
# Phys. Rev. 94-6 (1954) 1498-1524.
#
# Parameters for Si,C,Ge:
# D. J. Chadi and M. L. Cohen
# Phys. Stat. Sol. (b) 68 (1975) 405.

using LinearAlgebra
using PyPlot


# Points of high symmetry
Γ = [0,    0,   0]
Χ = [1,    0,   0]
L = [1/2,  1/2, 1/2]
U = [1,   1/4, 1/4]
K = [3/4,  3/4, 0]
W = [1,    1/2, 0]

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

# Slater-Koster parameters
# Vss = 4*Vssσ
# Vsp = -4*Vspσ/sqrt(3)
# Vxx = 4*(Vppσ/3 + 2*Vppπ/3)
# Vxy = 4*(Vppσ/3 - Vppπ/3)

# Si
Es = 0        # eV
Ep = 7.20     # eV
Vss = -8.13   # eV
Vsp = 5.88    # eV
Vxx = 1.71    # eV
Vxy = 7.51    # eV

# # C
# Es = 0        # eV
# Ep = 7.40     # eV
# Vss = -15.2   # eV
# Vsp = 10.25   # eV
# Vxx = 3.0     # eV
# Vxy = 8.3     # eV

# # Ge
# Es = 0        # eV
# Ep = 8.41     # eV
# Vss = -6.78   # eV
# Vsp = 5.31    # eV
# Vxx = 1.62    # eV
# Vxy = 6.82    # eV


a = 0.543E-9         # [m]

d1 = (a/4)*[ 1, 1, 1]
d2 = (a/4)*[ 1,-1,-1]
d3 = (a/4)*[-1, 1,-1]
d4 = (a/4)*[-1,-1, 1]

# Energy space 
εsp = LinRange(-9,16,100)
η = (16+9)/100
DOS = zeros(100)

# Compute LCAO along path
path = vcat(Λ,Δ,Ξ,Σ)
Bands = []
for p in path
       k = (2π/a)*p

       g1 = 0.25*(exp(im*k⋅d1) + exp(im*k⋅d2) + exp(im*k⋅d3) + exp(im*k⋅d4))
       g2 = 0.25*(exp(im*k⋅d1) + exp(im*k⋅d2) - exp(im*k⋅d3) - exp(im*k⋅d4))
       g3 = 0.25*(exp(im*k⋅d1) - exp(im*k⋅d2) + exp(im*k⋅d3) - exp(im*k⋅d4))
       g4 = 0.25*(exp(im*k⋅d1) - exp(im*k⋅d2) - exp(im*k⋅d3) + exp(im*k⋅d4))

       # ZincBlende
       H = [Es       0        0        0       Vss*g1 Vsp*g2 Vsp*g3 Vsp*g4;
            0        Ep       0        0      -Vsp*g2 Vxx*g1 Vxy*g4 Vxy*g3;
            0        0        Ep       0      -Vsp*g3 Vxy*g4 Vxx*g1 Vxy*g2;
            0        0        0        Ep     -Vsp*g4 Vxy*g3 Vxy*g2 Vxx*g1;
            Vss*g1' -Vsp*g2' -Vsp*g3' -Vsp*g4' Es     0      0      0;
            Vsp*g2'  Vxx*g1'  Vxy*g4'  Vxy*g3' 0      Ep     0      0;
            Vsp*g3'  Vxy*g4'  Vxx*g1'  Vxy*g2' 0      0      Ep     0;
            Vsp*g4'  Vxy*g3'  Vxy*g2'  Vxx*g1' 0      0      0      Ep]

       E = eigvals(H)
       push!(Bands,E)
       
       # Density of states
       for (idx,ε) in enumerate(εsp)
          G = inv((ε + im*η)*I - H)
          D = -imag(tr(G)/π)
          DOS[idx] += D
       end
end

# Plotting
subplot(121)
plot(Bands)
xlabel("k")
ylabel("Energy [eV]")
axis([0,300,-9,16])

subplot(122)
plot(DOS,εsp)
xlabel("DOS [states/eV]")
show()

