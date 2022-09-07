# Adam Rilatt
# 15 June 2022
# drycycle.gap

n_qubits := 2;

Read("rhobasis.gap");

dryer := GenDryer(n_qubits, true);

# the Dryer matrix seems to have (2n choose n) basis vectors. 
