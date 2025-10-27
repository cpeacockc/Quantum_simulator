#Example code for using Quantum Simulator, a Julia code for for simulating quantum mechanics
cd(@__DIR__);include("Pauli_generator.jl");include("Models.jl")
N=10 #Number of lattice sites

Z_5 = Op_create(N,"Z",5,"pauli") #Pauli Z operator on site 5
Sz_5 = Op_create(N,"Z",5,"spin") #spin Z operator on site 5 (Sz= 0.5Z)

#Other options are "X","Y","Z","+","-" (final two are \sigma^+ = X+iY and \sigma^- = X-iY)

#Now let's create a simple product state of up and down spins
#Options: "Neel" (alternating up & down), "Down" (all down), "Up (all up), "Rand" (random up and down)
state = Product_state(N,"Neel")

