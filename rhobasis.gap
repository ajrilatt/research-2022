# Adam Rilatt
# 09 June 2022

#n_qubits := 2;

# ==================================== #

Read("orbdim.gap");
Read("werner.gap");

# Prints the given matrix without commas or brackets,
# and replacing zeros with empty space for easier visual
# pattern recognition. This output tends to be physically smaller
# than Display(mat) and can be used to view higher-dimension matrices.
PrettifyMat := function(mat)
	local row, item;
	Print("[ ");
	for row in mat do
		Print("\n");
		for item in row do
			if item = 0 then
				Print("   ");
			else
				# Negative entries need one less space to align them.
				if item > 0 then
					Print("  ", item);
				else
					Print(" ", item);
				fi;
			fi;
		od;
	od;
	Print("]\n");
end;

# Similarly, condense the matrix so that it shows only non-zero entries with a placeholder.
CondenseMat := function(mat)
	local row, item;
	for row in mat do
		Print("\n");
		for item in row do
			if item = 0 then
				Print(" ");
			else
				Print("#");
			fi;
		od;
	od;
end;

# Given a single integer bit (0/1), returns
# the opposite bit.
BitComp := function(bit)
	return 1 - bit;
end;

#matrix_dims := 2 ^ (2 * n_qubits);

# Generates the Dryer matrix for n qubits, which obeys the first of two
# Werner state restrictions. Call as GenDryer(n_qubits, 1) to include
# a Dryer matrix printout.
GenDryer := function(n_qs, verbose...)

	local m_dims, rho, row, col, item, i, j, jflip, k;
	
	m_dims := 2 ^ (2 * n_qs);

	rho := NullMat(m_dims, m_dims);

	for col in [0..m_dims - 1] do

		# Actual row i and column j indices
		i := dec2binlist(col, 2 * n_qs);
		j := i;	
	
		for k in [1..(2 * n_qs)] do
	
			jflip := ShallowCopy(j);
			jflip[k] := BitComp(jflip[k]);
	
			if k <= n_qs then
				rho[Bin(i) + 1][Bin(jflip) + 1] := n_qs;
			else
				rho[Bin(i) + 1][Bin(jflip) + 1] := -n_qs;
			fi;
	
		od;	
	od;
	
	if Length(verbose) <> 0 then

		for row in rho do
			Print("\n");
			for item in row do
				if item = 0 then
					Print("  ");
				else
					Print(item, " ");
				fi;
			od;
		od;

	fi;

	return rho;

end;

#GenDryerSparse := function(n_qs, verbose...)
#
#	m_dims := 2 ^ n_qs;
#	rho := SparseMatrix()
#	
#end;

#rho := GenDryer(n_qubits, true);

# Experimenting with the eigenvalues and eigenvectors of the Dryer matrix.
#Print(Eigenvalues(Rationals, rho), "\n");
#eigenspaces := Eigenspaces(Rationals, rho);

#evectors := Eigenvectors(Rationals, rho);
#Print("\n");

#for vector in evectors do
#
#	Print("\n");
#
#	for item in vector do
#		
#		if item > 0 then
#			Print("  ", item);
#		else
#			if item < 0 then
#				Print(" ", item);
#			else 
#				Print("   ");
#			fi;
#		fi;
#
#	od;	
#
#od;

#for space in eigenspaces do
#	Print(Length(Basis(space)), "\n");
#	Print("\n"); 
#od;

 
#Print("\n");

# An experiment with the basis vectors of the Dryer matrix.
#a := NullspaceMat(rho);
#Print("Before weedout...\n");
#Display(a);
#Print("\nHow many of these basis vectors have weight 2 ^ n?");

AbsoluteVectorWeight := function(v)
	local count, vi;
	count := 0;
	for vi in v do
		if vi <> 0 then
			count := count + 1;
		fi;
	od;

	return count;

end;

NonzeroBinaryIndices := function(v, n_qs)
	local index_list, i;
	index_list := [];
	for i in [1..4 ^ n_qs] do
		if v[i] <> 0 then
			Add(index_list, [dec2binlist(i - 1, Length(v)), v[i]]);
		fi;
	od;
	
	return index_list;

end;

# Heck yeah brother
WeedoutReadout := function(a)

	local element_counter, indices, vector, i, bit;

	element_counter := 1;
	indices := [];
	for vector in a do
		Print("\nBasis element ", element_counter, " (", AbsoluteVectorWeight(vector), " weight)\n=====\n");
		indices := NonzeroBinaryIndices(vector, 2 ^ (2 * Length(a)));
		for i in indices do
			Print("    Index ");
			for bit in i[1] do
				Print(bit);
			od;
			Print(" with value ", i[2], "\n");  	
		od;
		Print("\n");
		element_counter := element_counter + 1;
	od;

end;

#WeedoutReadout(a);

#for row in a do
#	Print("\n");
#	for col in row do
#		if col <> 0 then
#			if col < 0 then
#				#Print(col, " ");
#			else
#				Print(col, "  ");
#			fi;
#		else
#			Print("   ");
#		fi;
#	od;
#od;
#Print("\n");

BitWeight := function(bitstring)
	
	local sum, bit;
	sum := 0; 
	for bit in bitstring do
		if bit = 1 then
			sum := sum + 1;
		fi;
	od;
	return sum;

end;

# Generates the Dryer matrix and then places the Washer matrix
# underneath it (thus the washer-dryer combo name). The Washer matrix
# nulls out any columns where the indexing bitstrings IJ and KL are unequal,
# fulfilling the second requirement on a Werner basis.
GenWasherDryer := function(n_qs)

	local dryer, washer, I, Ibin, J, Jbin, IJbin, IJ;
	
	dryer := GenDryer(n_qs);
	washer := NullMat(2 ^ (2 * n_qs), 2 ^ (2 * n_qs));	

	for I in [0..2 ^ (2 * n_qs) - 1] do
		
		Ibin := dec2binlist(I, n_qs);

		for J in [0..2 ^ (2 * n_qs) - 1] do
			
			Jbin := dec2binlist(J, n_qs);

			IJbin := ShallowCopy(Ibin);
			Append(IJbin, Jbin);
			IJ := Bin(IJbin);

			if BitWeight(Ibin) <> BitWeight(Jbin) then
				washer[IJ + 1][IJ + 1] := 1;
			fi;
		od;
	od;

	Append(dryer, washer);

	return dryer;

end;

# Using GenWasherDryer, generate the null space and obtain vectors
# forming (we hope!) a Werner basis for n qubits.
#basis := GenWasherDryer(n_qubits);
#Print("Found ", Length(basis), " basis vectors.\n");
#Display(basis);
#Print("\n");

# Given a vector, attempt to create a matrix with r rows and c columns.
Unflatten := function(vector, r, c)

	local ri, ci, rowvector, m;

	m := NullMat(r, c);

	for ri in [1..r] do
		for ci in [1..c] do		
			m[ri][ci] := vector[r * (ri - 1) + ci];
		od;
	od;

	return m;		
	
end;

#for b in basis do
#	Print("\n");
#	Display(Unflatten(b, 2 ^ n_qubits, 2 ^ n_qubits));
#od;

# Manually creating these would be horrifyingly inefficient. The
# NCPartitionsSet() function found in orbdim.gap generates these states.
#threegonstates := [
#	[[1], [2], [3]],
#	[[1, 2], [3]],
#	[[1, 3], [2]],
#	[[2, 3], [1]],
#	[[1, 2, 3]]
#];

#Print("Good to this point.\n");

#linear_combo_passed := true;

# Given a test Werner basis and a number of qubits, generates all Werner diagram states
# with n qubits and verifies that they can be written as a linear combination of the test Werner basis.
# Writing VerifyWernerBasis(test_basis, n_qubits, true) will include an output of the coefficients to put
# on the supplied basis to generate each Werner state.
VerifyWernerBasis := function(test_basis, n_qs, verbose...)

	local state, qubit_list, solution, nullmat;
	
	qubit_list := [1..n_qs];

	# GAP irritatingly does not supply an isNullMat() function, so we perfom the check manually.
	nullmat := NullMat(2 ^ n_qs, 2 ^  n_qs);

	for state in NCPartitionsSet(qubit_list) do
		
		solution := Unflatten(SolutionMat(test_basis, flatten(WernerDiagram(state))) * test_basis, 2 ^ n_qs, 2 ^ n_qs);
		
		if solution - WernerDiagram(state) <> nullmat then
			
			if Length(verbose) <> 0 then
				Print("Failed to represent Werner Diagram ", state, " in this basis.\n");
			fi;
			return false;
		fi;
		
		if Length(verbose) <> 0 then
			Print("Represented Werner Diagram ", state, " with basis coefficients ", solution, ".\n");
		fi;

	od;
	
	return true;

end;

#VerifyWernerBasis(basis, n_qubits, true);
