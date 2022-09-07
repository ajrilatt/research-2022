# Adam Rilatt
# 2 July 2022
# Multithreaded Helper Functions

Read("rhobasis.gap");

# === TODO ===
# - I wrote multithreaded code when I really should have run multiprocessed code.
#   Edit to perform multiprocessing.
# - Rename function and file to reflect this change.
# ============

# MultiPROCESSED version of NCDiagramsPairs.
NCDiagramsPairsMT := function(circle_indices)

	local diagrams, index, current_move, circle_m, circle_k,
	      combo, combos_m, combos_k, c1, c2, temp, circle_length,
	      t1;
	
	# GAP probably caches this, but I prefer to do so explicitly to be sure.
	circle_length := Length(circle_indices);

	# With two indices remaining, there is only one NC state-- a single singlet.
	if circle_length = 2 then
		return [[ [circle_indices[1], circle_indices[2]] ]];
	fi;

	diagrams := [];

	for index in [2, 4..circle_length] do
	
		current_move := [ circle_indices[1], circle_indices[index] ];

		circle_m := Difference(circle_indices{[1..index]}, current_move);
		circle_k := Difference(circle_indices{[index + 1..circle_length]}, current_move);

		if Length(circle_m) = 0 then
			circle_m := circle_k;
		fi;

		if Length(circle_m) = circle_length - 2 then
			
			# === TODO ===
			# - MultiPROCESS the NCDiagramsPairs call
			# ============
			
			t1 := RunTask(NCDiagramsPairsMT, circle_m); 

			for combo in TaskResult(t1) do #NCDiagramsPairsMT(circle_m) do

				# Appending-Adding is faster than list concatenation and flattening.
				temp := [current_move];
				Append(temp, combo);
				Add(diagrams, temp);

			od;

		else

			# === TODO ===
			# - MultiPROCESS the NCDiagramsPairs calls
			# ============

			#combos_m := NCDiagramsPairsMT(circle_m);
			#combos_k := NCDiagramsPairsMT(circle_k);
			combos_m := RunTask(NCDiagramsPairsMT, circle_m);
			combos_k := RunTask(NCDiagramsPairsMT, circle_k);			


			for c1 in TaskResult(combos_m) do #combos_m do
				for c2 in TaskResult(combos_k) do #combos_k do

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
