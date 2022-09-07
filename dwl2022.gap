# 6/2022 Werner invariance constraint for dens mats

Read("orbdim.gap");

BitComp := function(bit);
  return 1-bit;
end;

Wt := function(bitstring)
  local i, weight;
  weight := 0;
  for i in [1..Length(bitstring)] do
    weight := weight + bitstring[i];
  od;
  return weight;
end;

## turn an array vec with n^2 entries into an n x n matrix
UnFlatten := function(vec)
  local r,i,j,mat,currentrow;
  r := Sqrt(Length(vec));
  mat := [];
  for i in [1..r] do
    currentrow := [];
    for j in [1..r] do
      Append(currentrow,[vec[(i-1)*r+j]]);
    od;
    Print(currentrow,"\n");
    Append(mat,[currentrow]);
  od;
  return mat;
end;

MakeMat := function()
  local i,k,l,binlist,I,J,Itemp,Jtemp,OutMat;
  OutMat := NullMat(16,16);
  for i in [0..15] do
    binlist := dec2binlist(i,4);
    I := [binlist[1],binlist[2]];
    Itemp := [binlist[1],binlist[2]];
    J := [binlist[3],binlist[4]];
    Jtemp := [binlist[3],binlist[4]];
    for k in [1..2] do
      Itemp[k] := BitComp(I[k]);
      OutMat[Bin(binlist)+1][Bin([Itemp[1],Itemp[2],J[1],J[2]])+1] := 2;
      Itemp[k] := BitComp(Itemp[k]);
    od;
    for l in [1..2] do
      Jtemp[l] := BitComp(J[l]);
      OutMat[Bin(binlist)+1][Bin([I[1],I[2],Jtemp[1],Jtemp[2]])+1] := -2;
      Jtemp[l] := BitComp(Jtemp[l]);
    od;
    
  od;
  
  Print(OutMat);
  return NullspaceMat(OutMat);
end;

## make matrix for equations [sum C^{(k)},rho] = 0
## return basis for nullspace
# n = number of qubits
MakeMatn := function(n)
  local i,j,k,l,binlist,I,J,K,Itemp,Jtemp,OutMat,Mat2;
  OutMat := NullMat(2^(2*n),2^(2*n));
  for i in [0..2^(2*n)-1] do
    I := [];
    J := [];
    Itemp := [];
    Jtemp := [];        
    binlist := dec2binlist(i,2*n);
    for j in [1..n] do
      Append(I,[binlist[j]]);
      Append(Itemp,[binlist[j]]);
      Append(J,[binlist[n+j]]);
      Append(Jtemp,[binlist[n+j]]);            
    od;
#    Print(I,J,"\n");
    for k in [1..n] do
      Itemp[k] := BitComp(I[k]);
      K:= [];
      Append(K,Itemp);
      Append(K,J);      
#      Print(K,"\n");
      OutMat[Bin(binlist)+1][Bin(K)+1] := n;
#      if Wt(I) <> Wt(J) then Print(I,J,"mismatch \n"); fi;
      Itemp[k] := BitComp(Itemp[k]);
    od;
    for l in [1..n] do
      Jtemp[l] := BitComp(J[l]);
      K:= [];
      Append(K,I);
      Append(K,Jtemp);

      OutMat[Bin(binlist)+1][Bin(K)+1] := -n;
#      Print(I,J,Jtemp," ",K," printed minus n \n");
      Jtemp[l] := BitComp(Jtemp[l]);
    od;
  od;

## add rows to constraint matrix to zero out rho_{I,J} with wtI != wtJ
  Mat2 := NullMat(2^(2*n),2^(2*n));
  for i in [0..2^(2*n)-1] do
    I := [];
    J := [];
    binlist := dec2binlist(i,2*n);
    for j in [1..n] do
      Append(I,[binlist[j]]);
      Append(J,[binlist[n+j]]);
    od;
    if Wt(I) <> Wt(J) then Mat2[i+1][i+1] := 1; fi;
  od;
  Append(OutMat,Mat2);

#   Print(Display(OutMat));
   return NullspaceMat(TransposedMat(OutMat));
end;

### make the Werner basis for n=4
4outmat := [];
4out := MakeMatn(4);
for i in [1..14] do
  Append(4outmat,[UnFlatten(4out[i])]);
  Print(i,"\n");
  UnFlatten(4out[i]);
od;

for i in [1..14] do
  for j in [1..14] do
  if 4outmat[i] = TransposedMat(4outmat[j]) then Print(i," ",j,"\n");
    fi;
  od;
od;


