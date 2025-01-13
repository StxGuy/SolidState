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
aₒ = 10.26  # Bohrs
a₁ = 0.5*aₒ*[0,1,1]
a₂ = 0.5*aₒ*[1,0,1]
a₃ = 0.5*aₒ*[1,1,0]

lattice = [a₁ a₂ a₃]

Si₁ = [0.00,0.00,0.00]
Si₂ = [0.25,0.25,0.25]
positions = [Si₁, Si₂]
Si = ElementPsp(:Si,load_psp("hgh/pbe/Si-q4"))
atoms = [Si, Si]

# Run SCF
model = model_DFT(lattice, atoms, positions; functionals=PBE(), temperature)
basis = PlaneWaveBasis(model; Ecut, kgrid)
scf = self_consistent_field(basis)

# Construct 2D path through Brillouin zone
kpath = irrfbz_path(model; dim=3, space_group_number=227)
bands = compute_bands(scf, kpath; kline_density=20)
plot_bandstructure(bands)
#plot_dos(bands;temperature=1E-3,smearing=Smearing.FermiDirac())
