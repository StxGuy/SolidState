using DFTK
using Unitful, UnitfulAtomic
using LinearAlgebra
using Plots

# Parameters
L = 20
kgrid = [6,6,1]
Ecut = 15
temperature = 1E-3

# Geometry and pseudopotential
aₒ = 4.66
a₁ = aₒ*[1/2,-sqrt(3)/2,0]
a₂ = aₒ*[1/2, sqrt(3)/2,0]
a₃ =  L*[  0,         0,1]
lattice = [a₁ a₂ a₃]

C₁ = [1/3,-1/3,0.0]
C₂ = -C₁
positions = [C₁, C₂]
C = ElementPsp(:C,load_psp("hgh/pbe/c-q4"))
atoms = [C, C]

# Run SCF
model = model_DFT(lattice, atoms, positions; functionals=PBE(), temperature)
basis = PlaneWaveBasis(model; Ecut, kgrid)
scfres = self_consistent_field(basis)

# Construct 2D path through Brillouin zone
kpath = irrfbz_path(model; dim=2, space_group_number=13)    # Graphene
bands = compute_bands(scfres, kpath; kline_density=20
plot_bandstructure(bands)
