# Quantum-simulator
Functions for building and simulating quantum systems in Julia

For example:

# Number of lattice sites
N=10 

# First create a operator

Z_5 = Op_create(N,"Z",5,"pauli") # Pauli Z operator on site 5
Sz_5 = Op_create(N,"Z",5,"spin") # Spin Z operator on site 5 (Sz= 0.5Z)

#Other options are "X","Y","Z","+","-" (final two are \sigma^+ = X+iY and \sigma^- = X-iY)

# We can create a product state of up and down spins:
psi = Product_state(N,"Neel")

## Next we can create a variety of built-in Hamiltonians ##

# Isotropic Heisenberg model
H_XXX = H_XYZ(N,[1,1,1],zeros(N),zeros(N),zeros(N),"pauli")

# Chaotic Ising model (chaotic for nonzero hx)
hx=1.0
H_Chaotic_Ising = H_Quantum_Ising(N,hx)

# Anderson Model (3d)
d=3;W=20;t=1
H_Anderson = H_Anderson(d,N,W,t)

# Random GOE matrix
H_GOE = H_GOE(N)

# We can check if our operator commutes with one of our H's:
iszero(comm(Z_5,H_XXX))

# We can now time evolve our state by diagonalizing our H
t_array=collect(0:10)
psi_t = Diag_TE(psi,H)

# Or do the same plus measuring our operator Z_5

t_array=collect(0:10)
Ops=[Z_5]
Ops_t=Diag_TE_measure(psi,H,t_array,Ops)


