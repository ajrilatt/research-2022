# Adam Rilatt
# 6 July 2022
# Top NCPDs in the Pizza Basis

# === TODO ===
# - There are too many basis vectors; determine which one is extraneous
#   and figure out what part of the code will remove it. (Current guess
#   is a failure of the filtering function.)

Read("waffle.gap");

# GAP's absolute value function throws a fit when faced with imaginary numbers,
# so we have to define our own.
FlexibleAbsoluteValue := function(x)
	
	if RealPart(x) = 0 then

		if ImaginaryPart(x) >= 0 then
			return x; 
		else
			return -x;
		fi;
	
	else
		
		if RealPart(x) > 0 then
			return x;
		else
			return -x;
		fi;
	
	fi;

end;
 
MatrixAbsoluteValue := m -> List(m, x -> List(x, FlexibleAbsoluteValue));

GenSigmafiedBasis := function(n_qs)

	local pure_states, pure_state, sigmafied_states, abs_val_states,
	      zero_mat, sigmafied_state_1, sigmafied_state_2, abs_sigma_1,
	      abs_sigma_2;

	pure_states := List(NCDiagramsPairs([1..2 * n_qs]), GenPureState);
	
	# The pure states are either symmetric or come in symmetric pairs.
	# Instead of manually finding each symmetric pair and creating a
	# symmetric matrix by adding, subtracting, and scaling them, we can
	# apply that same transform across all pure states and then remove any
	# redundant states later.
	zero_mat := NullMat(2 ^ n_qs, 2 ^ n_qs);
	sigmafied_states := [];
	abs_val_states := [];
	for pure_state in pure_states do

		sigmafied_state_1 := (pure_state + (-1 ^ n_qs) * TransposedMat(pure_state)) / 2;
		sigmafied_state_2 := (pure_state - (-1 ^ n_qs) * TransposedMat(pure_state)) / (2 * E(4));
		abs_sigma_1 := MatrixAbsoluteValue(sigmafied_state_1);
		abs_sigma_2 := MatrixAbsoluteValue(sigmafied_state_2);
		
		if sigmafied_state_1 <> zero_mat and not (abs_sigma_1 in abs_val_states) then
			Add(abs_val_states, abs_sigma_1);
			Add(sigmafied_states, sigmafied_state_1);
		fi;
	
		if sigmafied_state_2 <> zero_mat and not (abs_sigma_2 in abs_val_states) then
			Add(abs_val_states, abs_sigma_2);
			Add(sigmafied_states, sigmafied_state_2);
		fi;

	od;

	return sigmafied_states;

end;

# N-gons to represent in the permutation basis.
#ns := [2, 3, 4, 5, 6, 7];
#for n in ns do
#
#	Print("N = ", n, "\n======\n\n");
#	top_polygon := WernerDiagram([[1..n]]);
#	basis := GenNCPermutationMats(n);
#	Print(SolutionMat(List(basis, flatten), flatten(top_polygon)));
#	Print("\n\n");
#
#od;
