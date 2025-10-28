# Quantum-simulator

Functions for building and simulating quantum systems in Julia with exact diagonalization and sparse matrices. For more advanced simulation methods check out [PauliStrings.jl](https://paulistrings.org/dev/), [ITensors.jl](https://docs.itensor.org/Overview/), and [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl)

## Example Usage

### Create a Pauli Z operator on site 5
Operator options are "X","Y","Z","+","-" (final two are \sigma^+ = X+iY and \sigma^- = X-iY)
```julia
N = 10 #Number of sites
Op = "Z"
site = 5 # site on lattice
spin_flag = "pauli" # options are "pauli" or "spin" for \sigma^a and S^a = 0.5 \sigma^a respectively

Z_5 = Op_create(N, Op, site, spin_flag)
```


### Next create a product state, in this case the Neel state:
Options are "Neel", "Up", "Dn", "Rand", for alternating up and down, all up, all down, and random up/down spins respectively
```julia
psi = Product_state(N,"Neel") 
```
## Built-in Hamiltonians

### Isotropic Heisenberg model
```julia
Jx,Jy,Jz = 1,1,1 # Define tunneling probabilities in each direction
X_fields, Y_fields, Z_fields = zeros(N),zeros(N),zeros(N) # Set potentials in each direction
H_XXX = H_XYZ(N,[Jx,Jy,Jz],X_fields, Y_fields, Z_fields,"pauli");
```
### Chaotic Ising model (chaotic for nonzero hx)
```julia
X_field=1.0; # X-field
H_Chaotic_Ising = H_Quantum_Ising(N,hx);
```
### Anderson Model (3d)
```julia
d=3; # Spatial dimension
W=20; # Disorder
t=1; # Tunneling probability
H_Anderson = H_Anderson(d,N,W,t) #Anderson Localization Model
```
### Random GOE matrix
```julia
H_GOE = H_GOE(N);
```
## Time evolution

### We can now time evolve our state by diagonalizing our H
```julia
N = 10
t_array=collect(0:10) # Define array of times to calculate
psi =  Product_state(N,"Neel")
H = H_XXX
psi_t = Diag_TE(psi,H,t_array) # Diagonalizes H to evolve psi
```

Or do the same plus measuring Z on each site
```julia
N=10
hx=1.0
t_array=collect(0:10);
Ops=[Op_create(N,"Z",site,"pauli") for site in 1:N];
H = H_Quantum_Ising(N,hx)
Ops_t=Diag_TE_measure(psi,H,t_array,Ops);
```

