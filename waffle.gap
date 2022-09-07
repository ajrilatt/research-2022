# Adam Rilatt
# 20 June 2022
# Werner Waffle Procedure

# Our goal is to construct the n-qubit mixed Werner states using
# the 2n-qubit pure Werner states. We posit that there exists a matrix
# that, by multiplication, will map between the two. Our guess is
# named the 'waffle state' because it looks a bit like a waffle:

#  --------
# /\ \ \ \ \
#|\ \ \ \ \ |
#| \ \ \ \ \|
# \ \ \ \ \/
#  --------

# This pure state can be made with singlet pairs (1, 2n), (2, 2n - 1),
# ... (n, n + 1). We additionally hypothesize that the resulting
# unitary permutation matrix will have the form

# sum(bitstrings I of length n) ( (-1) ^ (weight of I) * |I> <(tI)^c|
# where (tI) is bitstring I in reverse order, and ^c represents a bitwise
# complement.

# 21 June 2022 -- this method appears to work. It would appear that writing
# |v>|v> as |v><v| and then multiplying by waffle (or a similar matrix, 'pizza',
# which is the sum of all |I><I^c| and looks like symmetric crossing singlets)
# generates a Werner invariant state.

# 30 June 2022 -- It's worth noting that the 'waffle' and 'pizza' are not unique;
# any chord diagram that has 180 degree symmetry (as in it maps {1..n} to {n+1..2n})
# will yield a linear combination of mixed states. The waffle and pizza are simply
# convenvient.

# ===== TODO =====
# - (DONE) Check that Waffle * Mixed is a linear combo of pure states
# - (DONE) Track down and eliminate extra negative in Waffle function, if it exists
# - (DONE) Implement Pizza
# - (DONE) Fix GenPureState
# - Generate N = 4 change-of-basis matrix, generalize for large N

Read("rhobasis.gap");

compose := f -> g -> x -> f(g(x));

# Confusingly, this seems to work... but it's not the iY matrix we expect.
WaffleMaker := function(n_qubits)

	local mat, bitstrings, reversed_bitcomp, i;

	mat := NullMat(2 ^ n_qubits, 2 ^ n_qubits);

	bitstrings := Tuples([0, 1], n_qubits);
	reversed_bitcomp := List(bitstrings, compose(BitComplement)(Reversed));

	for i in [1..2 ^ n_qubits] do
		mat[Bin(bitstrings[i]) + 1][Bin(reversed_bitcomp[i]) + 1] := (-1) ^ AbsoluteVectorWeight(bitstrings[i]);
	od;

	return mat;

end;

BitCycle := blist -> Permuted(blist, CycleFromList([Length(blist), Length(blist)-1..1]));

# What about a rotation on the Waffle?
WaffleCousinMaker := function(n_qubits)

	local mat, bitstrings, cycled_reversed_bitcomp, i;
	
	mat := NullMat(2 ^ n_qubits, 2 ^ n_qubits);
	
	bitstrings := Tuples([0, 1], n_qubits);
	cycled_reversed_bitcomp := List(bitstrings, BitComplement);
	cycled_reversed_bitcomp := List(cycled_reversed_bitcomp, BitCycle);
	cycled_reversed_bitcomp := List(cycled_reversed_bitcomp, Reversed);

	for i in [1..2 ^ n_qubits] do
		mat[Bin(bitstrings[i]) + 1][Bin(cycled_reversed_bitcomp[i]) + 1] := (-1) ^ AbsoluteVectorWeight(bitstrings[i]);
	od;

	return WaffleMaker(3) - mat;

end;

Waffle_MixedtoPure := function(nc_partition)

	return WernerDiagram(nc_partition) * WaffleMaker(Length(flatten(nc_partition)));

end;

Waffle_PuretoMixed := function(nc_partition)
		
	return WernerDiagram(nc_partition) * TransposedMat(WaffleMaker(Length(flatten(nc_partition))));

end;

PizzaMaker := function(n_qubits)
	
	local mat, bitstrings, bitcomp, i;

	mat := NullMat(2 ^ n_qubits, 2 ^ n_qubits);

	bitstrings := Tuples([0, 1], n_qubits);
	bitcomp := List(bitstrings, BitComplement);
	
	for i in [1..2 ^ n_qubits] do
		mat[Bin(bitstrings[i]) + 1][Bin(bitcomp[i]) + 1] := (-1) ^ AbsoluteVectorWeight(bitstrings[i]);
	od;

	return mat;
	
end;

Pure_12_34_56 := [
		[0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 1, -1, 0],
		[0, -1, 1, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, -1, 1, 0],
		[0, 1, -1, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0]
	];;
	Pure_16_23_45 := [
		[0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 1, 0, -1, 0, 0],
		[0, 0, 0, -1, 0, 1, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, -1, 0, 1, 0, 0, 0],
		[0, 0, 1, 0, -1, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0]
	];;
	Pure_16_25_34 := [
		[0, 0, 0, 0, 0, 0, 0, 1],
		[0, 0, 0, -1, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, -1, 0, 0],
		[0, 1, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, -1, 0],
		[0, 0, 1, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 1, 0, 0, 0],
		[-1, 0, 0, 0, 0, 0, 0, 0]
	];;
	Pure_12_36_45 := [
		[0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 1, 0, -1, 0, 0],
		[0, 0, -1, 0, 1, 0, 0, 0],
		[0, 0, 0, -1, 0, 1, 0, 0],
		[0, 0, 1, 0, -1, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0]
	];;

	Pure_14_23_56 := [
		[0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 1, -1, 0],
		[0, 0, 0, 0, 0, -1, 1, 0],
		[0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0],
		[0, -1, 1, 0, 0, 0, 0, 0],
		[0, 1, -1, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0]
	];

pure_states := [Pure_16_25_34, Pure_12_36_45, Pure_14_23_56, Pure_16_23_45, Pure_12_34_56];

Verify3QubitWaffle := function()

	local Pure_12_34_56, Pure_16_23_45, Pure_16_25_34, Pure_12_36_45, Pure_14_23_56,
	      mixed_states, c, waffle, i;

	mixed_states := [
		[[1],[2],[3]],
		[[1], [2, 3]],
		[[1, 2], [3]],
		[[1, 3], [2]],
		[[1,  2,  3]]
	];

	c := [
		[1/8, 0, 0, 0, 0],
		[0, 0, 0, 0, -1/4],
		[0, 0, 0, -1/4, 0],
		[0, -1/4, -1/4, -1/4, -1/4],
		[0, E(3)/6, (E(3) ^ 2)/6, -1/6, -1/6]
	];

	waffle := WaffleMaker(3);
	
	for i in [1..5] do

		Print("\nConfirming linear combo of ", mixed_states[i], "...\n");

		Display(waffle * WernerDiagram(mixed_states[i]));
		Display(c[i][1] * pure_states[1] +
			c[i][2] * pure_states[2] +
			c[i][3] * pure_states[3] +
			c[i][4] * pure_states[4] +
			c[i][5] * pure_states[5]);
		Print("\n");
	od;

end;

VerifySymmetry := function(mat)

	local i, j;

	for i in [1..Length(mat)] do
		for j in [1..Length(mat[1])] do

			if mat[i][j] <> mat[j][i] then
				return false;
			fi;

		od;
	od;

	return true;

end;

GenPureState := function(singlets)

	# The pure states for n qubits can be generated algorithmically.
	# Given a list of integer pairs [[a, b], ... ] describing each
	# singlet in a non-crossing chord diagram, the matrix version of the
	# state is given as
	# sum(I)( (-1) ^ (weight I) * |I> PT <I ^ c| )
	# where PT is a positional tensor operation and ^c denotes bit complement.

	local n_qs, mat, position_lookup, i, j, bitstrings, bitcomp, bitstring, reconstructed_bitstring;

	n_qs := Length(singlets);
	mat := NullMat(2 ^ n_qs, 2 ^ n_qs);

	# Map the desired point number on the unitary circle to
	# its actual location in the bitstring.
	# In lieu of a proper hashmap, use indices as keys.
	position_lookup := 0 * [1..2 * n_qs];
	for i in [1..Length(singlets)] do
		position_lookup[singlets[i][1]] := i;
		position_lookup[singlets[i][2]] := i + n_qs;
	od;

	bitstrings := Tuples([0, 1], n_qs);
	bitcomp := List(bitstrings, BitComplement);

	for i in [1..2 ^ n_qs] do

		# Due to the positional tensor product, the bits are listed out-of-order
		# and must be re-ordered before they are interpreted as row and column
		# indexes.
		bitstring := [];
		Append(bitstring, bitstrings[i]);
		Append(bitstring, bitcomp[i]);
		reconstructed_bitstring := 0 * [1..2 * n_qs];

		for j in [1..2 * n_qs] do
			reconstructed_bitstring[j] := bitstring[position_lookup[j]];
		od;

		mat[Bin(reconstructed_bitstring{[1..n_qs]}) + 1][Bin(reconstructed_bitstring{[n_qs + 1..2 * n_qs]}) + 1] := (-1) ^ AbsoluteVectorWeight(bitstrings[i]);

	od; 

	return mat;

end;

# n = 3 change-of-basis matrix from pure to mixed states
#mat := NullMat(5, 5);
#pizza := PizzaMaker(3);
#bigboi := List([GenPureState([[1, 2], [3, 6], [4, 5]]),
#	   	GenPureState([[1, 4], [2, 3], [5, 6]]),
#	   	GenPureState([[1, 6], [2, 5], [3, 4]]),
#	       (GenPureState([[1, 2], [3, 4], [5, 6]]) + GenPureState([[1, 6], [2, 3], [4, 5]])) / 2, # pair component 1
#	       (GenPureState([[1, 2], [3, 4], [5, 6]]) - GenPureState([[1, 6], [2, 3], [4, 5]])) / (2 * E(4))  # pair component 2
#], flatten);
#mat[1] := SolutionMat(bigboi, flatten(pizza * WernerDiagram([[1],[2],[3]])));
#mat[2] := SolutionMat(bigboi, flatten(pizza * WernerDiagram([[1, 2], [3]])));
#mat[3] := SolutionMat(bigboi, flatten(pizza * WernerDiagram([[1, 3], [2]])));
#mat[4] := SolutionMat(bigboi, flatten(pizza * WernerDiagram([[1], [2, 3]])));
#mat[5] := SolutionMat(bigboi, flatten(pizza * WernerDiagram([[1,  2,  3]])));

#Display(mat);

#Print("\n");

#recreated_tri := mat[5][1] * GenPureState([[1, 2], [3, 6], [4, 5]]) + mat[5][2] * GenPureState([[1, 4], [2, 3], [5, 6]])
#	       + mat[5][3] * GenPureState([[1, 6], [2, 5], [3, 4]]) + mat[5][4] * (GenPureState([[1, 2], [3, 4], [5, 6]]) + GenPureState([[1, 6], [2, 3], [4, 5]])) / 2
#	       + mat[5][5] * (GenPureState([[1, 2], [3, 4], [5, 6]]) - GenPureState([[1, 6], [2, 3], [4, 5]])) / (2 * E(4));

# Generalized change-of-basis matrix.
# Step 1: Generate all pure states. How? Good question. We want only the non-crossing diagram states, of which
# 	  there are Catalan(n), out of (2n - 1)!! diagram states. We have an efficient weeding process; the
#	  question is one of generating all diagram states efficiently.
# Step 2: Identify the symmetric pairs and sigmafy them. Alternately, apply the sigmify process to all matrices
# 	  and their transposes, normalize them, then weed out any zero matrices or duplicates.
# Step 3: ???
# Step 4: Profit


# Generates all crossing and non-crossing chord diagram as a list of lists of singlets.
# NOTE: High memory usage and a fair amount of time is required. N = 5 used half a gigabyte of memory;
# N = 6 did not complete in under a minute and consumed worrying amounts of memory.
DiagramsPairs := function(n_qs)
	
	return Filtered(Combinations(Combinations([1..2 * n_qs], 2), n_qs), compose(IsDuplicateFree)(flatten));

end;

# Generates all non-crossing chord diagrams.
# NOTE: High memory usage and a fair amount of time is required. N = 5 used half a gigabyte of memory;
# N = 6 did not complete in under a minute and consumed worrying amounts of memory. A recursive version that
# generates only valid states will be implemented.
NCDiagramsPairsSlow := function(n_qs)

	local diagrams;
	diagrams := DiagramsPairs(n_qs);
	
	return Filtered(diagrams, x -> ForAll(Combinations(x, 2), y -> 
			(y[2][1] - y[1][1]) * (y[1][2] - y[2][1]) * (y[2][2] - y[1][1]) * (y[1][2] - y[2][2]) > 0));

end;


# Input: [1..2n] indices of the unitary circle, sorted.
# Output: Catalan(n) non-crossing chord diagram descriptions as lists of tuples.
NCDiagramsPairs := function(circle_indices)

	local diagrams, index, current_move, circle_m, circle_k,
	      combo, combos_m, combos_k, c1, c2, temp, circle_len;
	
	circle_len := Length(circle_indices);

	# With two indices remaining, there is only one NC state-- a single singlet.
	if circle_len = 2 then
		return [[ [circle_indices[1], circle_indices[2]] ]];
	fi;

	diagrams := [];

	# We construct each singlet pair beginning with the lowest-value index in
	# the unity circle (i.e. the one that appears at the top of the circle).
	# To enforce the non-crossing condition, we select other points on the condition
	# that they are an odd distance away from the first point.
	for index in [2, 4..circle_len] do
	
		current_move := [circle_indices[1], circle_indices[index]];

		# Placing a valid edge partitions the unitary circle into two sets
		# which can be considered smaller unitary circles.
		circle_m := Difference(circle_indices{[1..index]}, current_move);
		circle_k := Difference(circle_indices{[index + 1..circle_len]}, current_move);

		# Edges of the form (k, k + 1) result in only one subgraph. In these cases
		# we'll shift the nonempty subgraph to be subgraph m for convenience.
		if Length(circle_m) = 0 or Length(circle_k) = 0 then

			if Length(circle_m) = 0 then
				circle_m := circle_k;
			fi;

			for combo in NCDiagramsPairs(circle_m) do
				temp := [current_move];
				Append(temp, combo);
				Add(diagrams, temp);
			od;

		# Most valid edges will create two non-empty subgraphs, though.
		else

			combos_m := NCDiagramsPairs(circle_m);
			combos_k := NCDiagramsPairs(circle_k);

			for c1 in combos_m do
				for c2 in combos_k do
					temp := [current_move];
					Append(temp, c1);
					Append(temp, c2);
					Add(diagrams, temp);
				od;
			od;

		fi;

	od;
	
	return diagrams;	

end;
