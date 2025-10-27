# Quantum-simulator
Functions for building and simulating quantum systems in Julia

For example:

N=10 #Number of lattice sites

Z_5 = Op_create(N,"Z",5,"pauli") #Pauli Z operator on site 5
Sz_5 = Op_create(N,"Z",5,"spin") #spin Z operator on site 5 (Sz= 0.5Z)

#Other options are "X","Y","Z","+","-" (final two are \sigma^+ = X+iY and \sigma^- = X-iY)
