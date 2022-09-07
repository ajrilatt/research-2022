# Adam Rilatt
# 7 July 2022
# Representing Top Polygons as Linear Combinations of Permutation Matrices

Read("lincomb.gap");

# Unfortunately incorrect. Uses integer indexing, not binary indexing.
#GenPermutationMats := n -> List(List(PermutationsList([1..n]), PermList), x -> PermutationMat(x, 2 ^ n, Rationals));

GenPermutationMat := function(n, perm)

	local m, bindex, i, permuted;

	m := NullMat(2 ^ n, 2 ^ n);
	
	for i in [0..2 ^ n - 1] do
		
		bindex := dec2binlist(i, n);
		permuted := Permuted(bindex, perm);

		m[Bin(permuted) + 1][Bin(bindex) + 1] := 1;

	od;

	return m;	

end;

PermutationsListN := n -> List(PermutationsList([1..n]), PermList);

GenPermutationMats := function(n)

	local ms, perm;

	ms := [];

	for perm in PermutationsListN(n) do
		Add(ms, GenPermutationMat(n, perm));	
	od;

	return ms;

end;

# Life is pain and I hate everything
NCPermutationsListN := function(n)

	local partitions, cycles, partition, cycle, term;

	partitions := NCPartitionsSet([1..n]);
	cycles := [];
	
	for partition in partitions do
		cycle := ();
		for term in partition do
			cycle := cycle * CycleFromList(term);		
		od;
		Add(cycles, cycle ^ (-1));
	od;
	
	cycles := Unique(cycles);
	
	return cycles;

end;

GenNCPermutationMats := function(n)

	local ms, perm;
	
	ms := [];
	
	for perm in NCPermutationsListN(n) do
		Add(ms, GenPermutationMat(n, perm));
	od;
	
	return ms;

end;

#carpenter := WernerDiagram([[1, 2, 5], [3, 4]]);
#recreated_carpenter := (1 / 12) * (GenPermutationMat(5, ()) - GenPermutationMat(5, (3, 4)) + E(3) * GenPermutationMat(5, (1, 2, 5)) + (E(3) ^ 2) * GenPermutationMat(5, (1, 5, 2)));
#PrettifyMat(carpenter);
#Print("\n\n");
#PrettifyMat(recreated_carpenter);
#Print("\n\n");

BasisAsPermMats := function(n, verbose...) 

	local basis, perms, flat_perms, permlist, pizza, solution, bvector, i, j;

	# Generate the sigmafied basis for some n.
	basis := GenSigmafiedBasis(n);

	# Generate the permutation matrices.
	perms := GenNCPermutationMats(n);
	flat_perms := List(perms, flatten);

	# Generate the permutations in the corresponding order to the matrices.
	permlist := NCPermutationsListN(n);

	# TODO: does this need to be transposed or not?
	pizza := TransposedMat(PizzaMaker(n));

	solution := [];
	for bvector in basis do
		Add(solution, SolutionMat(flat_perms, flatten(pizza * bvector)));
	od;

	if Length(verbose) <> 0 then
		
		Print("N = ", n , " as a Linear Combination of Permutation Matrices\n");
		Print("======================================================\n");
		
		for i in [1..Length(solution)] do
			
			Print("\nBasis vector ", i, ": ");
			
			for j in [1..WernerBasisCheck(n)] do
				
				# === TODO ===
				# - The coefficient solution[i][j] is not recognized as <> 0 if
				#   it is imaginary / cyclotomic. Write an alternate version of
				#   this nonzero check to catch imaginary / cyclotomic numbers.
				
				if RealPart(solution[i][j]) <> 0 or ImaginaryPart(solution[i][j]) <> 0 then	
					if j > 1 and (RealPart(solution[i][j]) > 0 or ImaginaryPart(solution[i][j]) > 0) then
						Print("+");
					fi;
				
					Print(solution[i][j], " * V_", permlist[j], " "); 

				fi;			

			od;

		od;
		
	fi;

	return solution;

end;

MatAsNCPermMats := function(n, mat) 

	local flat_perms, solution;

	# Generate the permutation matrices.
	flat_perms := List(GenNCPermutationMats(n), flatten);

	solution := SolutionMatDestructive(flat_perms, flatten(mat));
	return solution;

end;



CrossingPermutationFilter := p -> (p = CycleFromList(Reversed(MovedPoints(p))));

PermAsPermMats := function(n, perm)

	local perm_mat, perms, perm_mats, perms_flat, coeffs, solution_vec;
	
	solution_vec := [];

	# Given a permutation matrix...
	perm_mat := GenPermutationMat(n, perm);

	# ... what linear combination of equally-sized permutation matrices
	# forms it?
	perms := Difference(NCPermutationsListN(n), [perm]);
	perm_mats := List(perms, x -> GenPermutationMat(n, x));
	perms_flat := List(perm_mats, flatten);

	coeffs := SolutionMat(perms_flat, flatten(perm_mat));	

	for i in [1..Length(coeffs)] do

		Add(solution_vec, [coeffs[i], perms[i]]);

	od;

	return solution_vec;

end;

# Because GenPermutationMat(n, ...) takes too long to type out when
# doing things like this by hand. 
V := x -> GenPermutationMat(4, x);

# Testing 3-gon construction from permutation matrices
n := WernerDiagram([[1, 2, 3]]);
m := (1 - (E(3) ^ 2)) * V(()) + (E(3) ^ 2) * (V((1, 2)) + V((2, 3)) + V((1, 3))) + (E(3) - (E(3) ^ 2)) * V((1, 3, 2));

# Testing 4-gon construction from permutation matrices
#n := WernerDiagram([[1, 2, 3, 4]]);
# === TODO ===
#m := (3 + 2 * E(4)) * V(()) - E(4) * (V((1, 2)) + V((1, 4)) + V((2, 3)) + V((3, 4))) + c * () + (E(4) - 1) * (V((1, 3, 2)) + V((1, 4, 2)) + V((2, 4, 3)) + V((1, 4, 3))) + V((1, 4, 3, 2));

#n := V((1, 2, 3, 4));
#m := -2 * V(()) + V((1, 2)) + V((1, 3)) + V((1, 4)) + V((2, 3)) + V((2, 4)) + V((3, 4)) - V((1, 3, 2)) - V((1, 4, 2)) - V((1, 4, 3)) - V((2, 4, 3)) + V((1, 4, 3, 2));

# Fun new conjecture:
# The rho matrix for the n-gon is a linear combination of the permutation matrices for all possible diagram states.

VerifyNonZeroTopPolygon := function(ns)

	local n, perms, perm_coeffs, index;

	for n in ns do
		
		Print("Starting n = ", n, "... ");
		perm_coeffs := MatAsNCPermMats(n, WernerDiagram([[1..n]]));
		perms := NCPermutationsListN(n);		
		index := Position(perms, CycleFromList([1..n]) ^ -1, 1);
		Print(perm_coeffs[index], "\n");
		
	od;	


end;

VerifyTopPolyExclusive := function(n)

	local diagrams, perm_coeffs, perms, i, index;

	diagrams := NCPartitionsSet([1..n]);
	perm_coeffs := List(List(diagrams, WernerDiagram), x -> MatAsNCPermMats(n, x));
	perms := NCPermutationsListN(n);
	
	for i in [1..Length(perm_coeffs)] do

		index := Position(perms, CycleFromList([1..n]) ^ -1, 1);
		if perm_coeffs[i][index] <> 0 then
			Print("Nonzero coefficient on ", diagrams[i], "\n");
		fi;	
	od;

end;

# Task: for each Werner diagram, ascertain whether the permutation matrix for the 'waffle state'
# -- representable as (1, 2n)(2, 2n - 1)...(n, n + 1) -- appears in the linear combination of
# permutation matrices representing that diagram.

WaffleContribution := function(diagram)

	local n, diagram_mat, coeffs, perms, waffle_perm;
	
	n := Length(flatten(diagram));
	diagram_mat := WernerDiagram(diagram);
	coeffs := MatAsNCPermMats(n, diagram_mat);
	perms := NCPermutationsListN(n);

	# What does the waffle state look like for any given n? We generate it so we can search for it.
	waffle_perm := ();
	for i in [1..Int(n - (n mod 2)) / 2] do
		waffle_perm := waffle_perm * (i, n + 1 - i);
	od;

	return coeffs[Position(perms, waffle_perm)];

end;

WernerStatesWithWaffles := function(n)
	
	local diagrams, waffle_coeffs, i;

	diagrams := NCPartitionsSet([1..n]);
	waffle_coeffs := List(diagrams, WaffleContribution);

	for i in [1..Length(diagrams)] do
		if waffle_coeffs[i] <> 0 then
			Print(waffle_coeffs[i], " contribution of Waffle state in Werner diagram ", diagrams[i], "\n");
		fi;
	od;

end;

# Generate the change-of-basis matrix that describes the transformation from
# the mixed Werner diagrams to the linear combination of our mixed NCC basis (sigmafied basis).
# TODO: what happens if you pizzafy both sides?
ChangeOfBasisMatrix := function(n, verbose...)

	local werner_diagrams, sig_basis, solution, diagram, pizza, partitions, i;
	partitions := NCPartitionsSet([1..n]);
	werner_diagrams := List(partitions, WernerDiagram);
	sig_basis := List(GenSigmafiedBasis(n), flatten);
	pizza := TransposedMat(PizzaMaker(n));
	solution := [];

	for diagram in werner_diagrams do		
		Add(solution, SolutionMat(sig_basis, flatten(pizza * diagram)));
	od;

	if Length(verbose) <> 0 then

		Print("The following basis vectors have these IDS:\n");

		for i in [1..Length(sig_basis)] do
			Print(i, ":\n");
			Display(Unflatten(sig_basis[i], 2 ^ n, 2 ^ n));
			Print("\n");
		od;

		for i in [1..Length(solution)] do
			Print(partitions[i], " uses coefficients ");
			Display(solution[i]);
		od;
	fi;

	return solution;

end;

TopPolygonPermBasis := function(n)

	local basis, permutes, top_poly, solution, index;

	basis := GenNCPermutationMats(n);
	permutes := NCPermutationsListN(n);
	top_poly := WernerDiagram([[1..n]]);	

	solution := SolutionMatDestructive(List(basis, flatten), flatten(top_poly));
	index := Position(permutes, CycleFromList([n, n-1..1]));
	Print(solution, "\n");
	Print(permutes, "\n");
	return solution[index];

end;

# What happens when you take the V_n contribution out of rho_n?
# How does it behave differently, and does that lend itself to a proof?
# (This is just the linear combination of the v_lowers.)
GenToplessRho := function(n)
	local rho, perm_coeffs, perms, top_v_index, mat, i, mats;

	rho := WernerDiagram([[1..n]]);
	perm_coeffs := MatAsNCPermMats(n, rho);
	
	perms := NCPermutationsListN(n);
	top_v_index := Position(perms, CycleFromList([1..n]) ^ -1, 1);
	Remove(perms, top_v_index);
	Remove(perm_coeffs, top_v_index);

	mat := NullMat(2 ^ n, 2 ^ n);
	for i in [1..Length(perm_coeffs)] do
		if perm_coeffs[i] <> 0 then
			mat := mat + perm_coeffs[i] * GenPermutationMat(n, perms[i]);
		fi;
	od;

	# As an extra little experiment...
	# Can you build rho_ngon as a linear combination of these lower permutation matrices?
	mats := GenNCPermutationMats(n);
	Remove(mats, top_v_index);
	Print(SolutionMat(List(mats, flatten), flatten(rho)), "\n");
	
	return mat;

end;

SharedOrUnique := function(n)

    local permutations_set, shared_or_unique, unique_entries, i, j, k,
	  entries, nonzero_index, h, shared_diagrams, shared_entries, g, diagrams;

    permutations_set := NCPermutationsListN(n);
    shared_or_unique := [];

    for i in [1..Length(permutations_set)] do

        unique_entries := [];
        diagrams := List(permutations_set, x -> GenPermutationMat(n, x));
        
        for j in [1..2 ^ n] do
            for k in [1..2 ^ n] do
                    
                # Lord forgive me...
                entries := flatten(flatten(diagrams{[1..Length(diagrams)]}{[j]}{[k]}));
                if Sum(entries) = 1 then

		    nonzero_index := Position(entries, 1);

		#    if not(permutations_set[nonzero_index] in unique_entries) then
		   	Print("UNZ ", permutations_set[nonzero_index], "\t in position K = ", dec2binlist(k - 1, n), ", J = ", dec2binlist(j - 1, n));
			Print("(", k, ", ", j, ") dependent on ");

			# We're curious what nonzero entries this shares with the previous tiers.
			for h in [1..Length(shared_or_unique)] do
				
				shared_diagrams := List(shared_or_unique[h], x -> GenPermutationMat(n, x));
				shared_entries := flatten(flatten(shared_diagrams{[1..Length(shared_diagrams)]}{[j]}{[k]}));

				for g in [1..Length(shared_entries)] do
					if shared_entries[g] = 1 then
						Print("\n\t", shared_or_unique[h][g], " from tier ", h, ",");
					fi;	 
				od;

			od;
			
			Print("\n");
		
                   # fi;

		    Add(unique_entries, permutations_set[nonzero_index]);

                fi;

            od;

        od;

        Add(shared_or_unique, Unique(unique_entries));
        permutations_set := Difference(permutations_set, unique_entries);
    
        if Length(permutations_set) = 0 then
            break;
        fi;

    od;

    return shared_or_unique;

end;

GenMask := function(tierlist)

	local mat, perm, i, j;
	
	mat := NullMat(2 ^ 5, 2 ^ 5);
	for perm in tierlist do

		mat := mat + GenPermutationMat(5, perm);

	od;

	for i in [1..Length(mat)] do
		for j in [1..Length(mat[1])] do
			if mat[i][j] <> 0 then
				mat[i][j] := 1;
			else
				mat[i][j] := 0;
			fi;
		od;
	od;

	return mat;

end;

CombineMasks := function(m1, m2)

	local i, j, m1c;

	m1c := ShallowCopy(m1);
	
	for i in [1..Length(m1c)] do
		for j in [1..Length(m1c[1])] do
			if m1[i][j] <> 1 or m2[i][j] <> 1 then
				m1c[i][j] := 0;
			fi;
				
		od;
	od;
		
	return m1c;

end;
