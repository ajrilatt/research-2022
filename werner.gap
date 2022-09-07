### code mostly by SNW, summer 2013
#Read("orbdim20130726.gap");

m := normalizeHS(I);

mix := n -> normalizeHS(Kron(replicate(n)(I)));
ExpandIn := ExpandInOrthogonalBasis(HilbertSchmidtInnerProduct);

# Two-qubit Werner states

werner2I := [WernerDiagram([[1],[2]])
            ,WernerDiagram([[1,2]])
            ];

werner2S := [WernerDiagram([[1,2]])
            ,WernerDiagram([[1],[2]])
            ];

werner2Iov := Overlap(werner2I,HilbertSchmidtInnerProduct);
werner2Sov := Overlap(werner2S,HilbertSchmidtInnerProduct);

werner2Ion := GramSchmidt(werner2I,HilbertSchmidtInnerProduct);
werner2Son := GramSchmidt(werner2S,HilbertSchmidtInnerProduct);

werner2orthog1 := [dmsinglet,UniformDickeMixture(2)];
werner2orthog2 := [dm2mixed,dm2mixed-dmsinglet];

# gap> VectorSpace(Rationals,werner2orthog1)=VectorSpace(Rationals,werner2orthog2);
# true

ir2 := IrrepBasis([1,0])[1];
irrepBasis2 := fmap(normalizeHS)(Concatenation([Kron([I,I])],IrrepBasis([1,0])));

ExpandIn2 := ExpandInOrthogonalBasis(HilbertSchmidtInnerProduct)(irrepBasis2);

# Three-qubit Werner states

werner3orig := [WernerDiagram([[1],[2],[3]])
               ,WernerDiagram([[1,2],[3]])
               ,WernerDiagram([[1,3],[2]])
               ,WernerDiagram([[2,3],[1]])
               ,WernerDiagram([[1,2,3]])
               ];

overlaporig:=Overlap(werner3orig,HilbertSchmidtInnerProduct);

werner3evil := [WernerGon(3)
               ,WernerDiagram([[1,2],[3]])
               ,WernerDiagram([[1,3],[2]])
               ,WernerDiagram([[2,3],[1]])
               ,WernerDiagram([[1,2,3]])
               ];

overlapevil:=Overlap(werner3evil,HilbertSchmidtInnerProduct);

werner3 := [WernerDiagram([[1,2,3]])
           ,WernerGon(3)
           ,WernerDiagram([[1],[2],[3]])
	   ,WernerDiagram([[1,2],[3]])
	   ,WernerDiagram([[1,3],[2]])
	   ];

state := GramSchmidt(werner3,HilbertSchmidtInnerProduct);

overlapnew:=Overlap(state,HilbertSchmidtInnerProduct);

# This is orthogonal
werner3ortho := [WernerDiagram([[1,2,3]])
                ,WernerGon(3)
#                ,SymmetrizeRho(SymRhoFromPoly(3,1+x^2+y^2+z^2))  This is uniform Dicke mixture
                ,UniformDickeMixture(3)
	        ,WernerDiagram([[1,2],[3]])-WernerDiagram([[1,3],[2]])
                ,WernerDiagram([[1,2],[3]])+WernerDiagram([[1,3],[2]])-2*WernerDiagram([[2,3],[1]])
                ];

overlaportho:=Overlap(werner3ortho,HilbertSchmidtInnerProduct);

werner3orthog1 := [dm3mixed
                  ,Kron([rr,I])
                  ,OrderedKron([0,1,0],rr,I)
                  ,Kron([I,rr])
                  ,Sqrt(3)*AntiSymmetrizeRho(Kron([sx,sy,sz]))
                  ];
# Without the square root of 3, werner3orig and werner3orthog1 are
# not equal as vector spaces over the rationals.

# gap> VectorSpace(Rationals,werner3orig)=VectorSpace(Rationals,werner3orthog1);                
# true

werner3SubtractMixed := [WernerDiagram([[1],[2],[3]])
               ,WernerDiagram([[1,2],[3]])-WernerDiagram([[1],[2],[3]])
               ,WernerDiagram([[1,3],[2]])-WernerDiagram([[1],[2],[3]])
               ,WernerDiagram([[2,3],[1]])-WernerDiagram([[1],[2],[3]])
               ,WernerDiagram([[1,2,3]])-WernerDiagram([[1],[2],[3]])
               ];

ir3 := IrrepBasis([1,1,0])[1];
irrepBasis3 := fmap(normalizeHS)([Kron([I,I,I])
                                 ,OrderedKron([0,1,1],I,ir2)
                                 ,OrderedKron([1,0,1],I,ir2)
                                 ,OrderedKron([1,1,0],I,ir2)
                                 ,-i*IrrepBasis([1,1,0])[1]]);

ExpandIn3 := ExpandIn(irrepBasis3);

ghz13 := (e000+e111)/Sqrt(2);
ghz23 := (e000-e111)/Sqrt(2);
w13   := (e100+e010+e001)/Sqrt(3);
w23   := (e110+e011+e101)/Sqrt(3);
T13  := (e100+E(3)*e010+E(3)^2*e001)/Sqrt(3);
T13s := (e100+E(3)^2*e010+E(3)*e001)/Sqrt(3);
T23  := (e110+E(3)*e101+E(3)^2*e011)/Sqrt(3);
T23s := (e110+E(3)^2*e101+E(3)*e011)/Sqrt(3);

# This is an orthonormal basis
cui3pure := [ghz13
            ,ghz23
            ,w13
            ,w23
            ,T13
            ,T13s
            ,T23
            ,T23s
            ];

T13T13   := OuterProduct(T13,T13);
T13T13s  := OuterProduct(T13,T13s);
T13sT13  := OuterProduct(T13s,T13);
T13sT13s := OuterProduct(T13s,T13s);
T23T23   := OuterProduct(T23,T23);
T23T23s  := OuterProduct(T23,T23s);
T23sT23  := OuterProduct(T23s,T23);
T23sT23s := OuterProduct(T23s,T23s);

# possible 3-qubit mixed state Werner basis

# w3a is the uniform Dicke mixture for 3 qubits
# w3b = ComplexConjugate(w3c);
# w3c is WernerDiagram([[1,2,3]]);
w3a := 1/4 * DM(ghz13) + 1/4 * DM(ghz23) + 1/4 * DM(w13) + 1/4 * DM(w23);
w3b := 1/2 * DM(T13) + 1/2 * DM(T23s);
w3c := 1/2 * DM(T13s) + 1/2 * DM(T23);
# w3d := 

# This is an orthonormal basis for mixed states in the Hilbert-Schmidt inner product
cui3OuterProduct := function()
  local result,vec1,vec2;
  result := [];
  for vec1 in cui3pure do
    for vec2 in cui3pure do
      result := Concatenation(result,[OuterProduct(vec1,vec2)]);
    od;
  od;
  return result;
end;

# Four-qubit Werner states

# werner4 := [WernerDiagram([[1,2,3,4]])
#            ,WernerGon(4)
#            ,UniformDickeMixture(4)
#            ,DMn(m4)
#            ,DMn(m4bar)
#            ];

# werner4overlap := Overlap(werner4,HilbertSchmidtInnerProduct);

# werner4on := GramSchmidt(werner4,HilbertSchmidtInnerProduct);

# werner4overlapon := Overlap(werner4on,HilbertSchmidtInnerProduct);

werner4 := [WernerDiagram([[1,2,3,4]])
           ,WernerDiagram([[1,2],[3,4]])
           ,WernerDiagram([[1,4],[2,3]])
           ,WernerDiagram([[1,2],[3],[4]])
           ,WernerDiagram([[2,3],[4],[1]])
           ,WernerDiagram([[3,4],[1],[2]])
           ,WernerDiagram([[4,1],[2],[3]])
           ,WernerDiagram([[1,3],[2],[4]])
           ,WernerDiagram([[2,4],[1],[3]])
           ,WernerDiagram([[1,2,3],[4]])
           ,WernerDiagram([[2,3,4],[1]])
           ,WernerDiagram([[3,4,1],[2]])
           ,WernerDiagram([[4,1,2],[3]])
           ,WernerDiagram([[1],[2],[3],[4]])
           ];

werner4Overlap := Overlap(werner4,HilbertSchmidtInnerProduct);

xyz := Sqrt(3)*AntiSymmetrizeRho(Kron([sx,sy,sz]));

r4 := 1/16*(Kron([sx,sx,sx,sx]) + Kron([sy,sy,sy,sy]) + Kron([sz,sz,sz,sz]));
r3 := FieldExtension(Rationals,x^2-3);

werner4orthog1 := [dm4mixed
                  ,OrderedKron([0,0,1,1],rr,dm2mixed)
                  ,OrderedKron([0,1,0,1],rr,dm2mixed)
                  ,OrderedKron([0,1,1,0],rr,dm2mixed)
                  ,OrderedKron([1,0,0,1],rr,dm2mixed)
                  ,OrderedKron([1,0,1,0],rr,dm2mixed)
                  ,OrderedKron([1,1,0,0],rr,dm2mixed)
                  ,OrderedKron([0,0,0,1],xyz,I)
                  ,OrderedKron([0,0,1,0],xyz,I)
                  ,OrderedKron([0,1,0,0],xyz,I)
                  ,OrderedKron([1,0,0,0],xyz,I)
                  ,PauliWeight(4)(WernerDiagram([[1,2,3,4]]))
                  ,PauliWeight(4)(WernerDiagram([[1,2],[3,4]])-r4)
                  ,PauliWeight(4)(WernerDiagram([[1,4],[2,3]])-r4)
                  ];

werner4low := [WernerDiagram([[1],[2],[3],[4]])
           ,WernerDiagram([[1,2],[3],[4]])
           ,WernerDiagram([[2,3],[4],[1]])
           ,WernerDiagram([[3,4],[1],[2]])
           ,WernerDiagram([[4,1],[2],[3]])
           ,WernerDiagram([[1,3],[2],[4]])
           ,WernerDiagram([[2,4],[1],[3]])
           ,WernerDiagram([[1,2,3],[4]])
           ,WernerDiagram([[2,3,4],[1]])
           ,WernerDiagram([[3,4,1],[2]])
           ,WernerDiagram([[4,1,2],[3]])
           ,WernerDiagram([[1,2,3,4]])
           ,WernerDiagram([[1,2],[3,4]])
           ,WernerDiagram([[1,4],[2,3]])
           ];

werner4w := [PauliWeight(0)(WernerDiagram([[1],[2],[3],[4]]))
            ,PauliWeight(2)(WernerDiagram([[1,2],[3],[4]]))
            ,PauliWeight(2)(WernerDiagram([[2,3],[4],[1]]))
            ,PauliWeight(2)(WernerDiagram([[3,4],[1],[2]]))
            ,PauliWeight(2)(WernerDiagram([[4,1],[2],[3]]))
            ,PauliWeight(2)(WernerDiagram([[1,3],[2],[4]]))
            ,PauliWeight(2)(WernerDiagram([[2,4],[1],[3]]))
            ,PauliWeight(3)(WernerDiagram([[1,2,3],[4]]))
            ,PauliWeight(3)(WernerDiagram([[2,3,4],[1]]))
            ,PauliWeight(3)(WernerDiagram([[3,4,1],[2]]))
            ,PauliWeight(3)(WernerDiagram([[4,1,2],[3]]))
            ,PauliWeight(4)(WernerDiagram([[1,2,3,4]]))
            ,PauliWeight(4)(WernerDiagram([[1,2],[3,4]]))
            ,PauliWeight(4)(WernerDiagram([[1,4],[2,3]]))
            ];

ir4a := IrrepBasis([1,0,1,0])[1];
ir4b := IrrepBasis([1,1,1,0])[1];
ir4c := IrrepBasis([1,2,1,0])[1];
irrepBasis4 := fmap(normalizeHS)([Kron([I,I,I,I])
                                 ,OrderedKron([0,0,1,1],ir2,mix(2))
                                 ,OrderedKron([0,1,0,1],ir2,mix(2))
                                 ,OrderedKron([0,1,1,0],ir2,mix(2))
                                 ,OrderedKron([1,0,0,1],ir2,mix(2))
                                 ,OrderedKron([1,0,1,0],ir2,mix(2))
                                 ,OrderedKron([1,1,0,0],ir2,mix(2))
                                 ,-i*OrderedKron([0,0,0,1],ir3,I)
                                 ,-i*OrderedKron([0,0,1,0],ir3,I)
                                 ,-i*OrderedKron([0,1,0,0],ir3,I)
                                 ,-i*OrderedKron([1,0,0,0],ir3,I)
                                 ,ir4a
                                 ,ir4b
                                 ,ir4c
                                 ]);

ExpandIn4 := ExpandIn(irrepBasis4);

# Five-qubit Werner states

# sum5 := OrderedKron([1,0,1,1,0],PauliWeight(2)(dmsinglet),PauliWeight(3)(WernerDiagram([[1,2,3]])))
#       + OrderedKron([0,1,0,1,1],PauliWeight(2)(dmsinglet),PauliWeight(3)(WernerDiagram([[1,2,3]])))
#       + OrderedKron([1,1,0,1,0],PauliWeight(2)(dmsinglet),PauliWeight(3)(WernerDiagram([[1,2,3]])))
#       + OrderedKron([1,0,1,0,1],PauliWeight(2)(dmsinglet),PauliWeight(3)(WernerDiagram([[1,2,3]])))
#       + OrderedKron([0,1,1,0,1],PauliWeight(2)(dmsinglet),PauliWeight(3)(WernerDiagram([[1,2,3]])));

# sum5a := PAF(PauliWeight(5)(normalizeDM(
#            WernerDiagram([[2,5],[1,3,4]])
#          + WernerDiagram([[1,3],[2,4,5]])
#          + WernerDiagram([[3,5],[1,2,4]])
#          + WernerDiagram([[2,4],[1,3,5]])
#          + WernerDiagram([[1,4],[2,3,5]]))));

triangleLineNoCross := [WernerDiagram([[1,2],[3,4,5]])
                       ,WernerDiagram([[2,3],[4,5,1]])
                       ,WernerDiagram([[3,4],[5,1,2]])
                       ,WernerDiagram([[4,5],[1,2,3]])
                       ,WernerDiagram([[5,1],[2,3,4]])
                       ];

triangleLineCross := [WernerDiagram([[1,3],[2,4,5]])
                     ,WernerDiagram([[2,4],[3,5,1]])
                     ,WernerDiagram([[3,5],[4,1,2]])
                     ,WernerDiagram([[4,1],[5,2,3]])
                     ,WernerDiagram([[5,2],[1,3,4]])
                     ];

triangleLine10 := Concatenation(triangleLineNoCross,triangleLineCross);

fivegon := WernerDiagram([[1,2,3,4,5]]);

############################################
# Constructing non-crossing diagram states #
############################################

# cross checking
# See NonCrossing.hs from summer2013

adjacentPairs := xs -> zip(xs,Concatenation(cdr(xs),[car(xs)]));

filterByPair := pr -> function(xs)
  if pr[1] < pr[2] then
    return filter(xs,x -> pr[1] <= x and x < pr[2]);
  else
    return filter(xs,x -> pr[1] <= x or  x < pr[2]);
  fi;
end;

partitionByList := ms -> xs -> List(adjacentPairs(SortedList(ms)),pr -> filterByPair(pr)(xs));

# GAP's not is not a function
Not := b -> not b;

nonCrossing := ms -> xs -> Length(filter(partitionByList(ms)(xs),compose(Not)(IsEmpty))) = 1;

nonCrossingPair := pr -> nonCrossing(pr[1])(pr[2]);

NCPartition := part -> And(fmap(nonCrossingPair)(Combinations(part,2)));

# Non-crossing partitions
NCPartitionsSet := xs -> filter(PartitionsSet(xs),NCPartition);

########################################
# Checking the Werner-basis conjecture #
########################################

WernerBasisCheck := n -> RankMat(List(NCPartitionsSet([1..n]),compose(flatten)(WernerDiagram)));

# gap> WernerBasisCheck(3);
# 5
# gap> WernerBasisCheck(4);
# 14
# gap> WernerBasisCheck(5);
# 42
# gap> WernerBasisCheck(6);
# 132
# gap> WernerBasisCheck(7);
# 429

####################
# Dotless diagrams #
####################

dotless := nss -> Count(1,List(nss,Length)) = 0;

# Dotless Non-crossing partitions
DotlessNCPartitionsSet := xs -> filter(PartitionsSet(xs),nss -> NCPartition(nss) and dotless(nss));

DotlessNCPartitionsSetImperative := function(xs)
  local ps,result;
  ps := PartitionsSet(xs);
  result := [];
  for part in ps do
    if dotless(part) and NCPartition(part) then
      Add(result,part);
    fi;
  od;
  return result;
end;

############################
# Dotless Diagram Checking #
############################

DotlessWernerBasisCheck := n -> RankMat(List(DotlessNCPartitionsSet([1..n]),compose(flatten)(WernerDiagram)));

# gap> DotlessWernerBasisCheck(2);
# 1
# gap> DotlessWernerBasisCheck(3);
# 1
# gap> DotlessWernerBasisCheck(4);
# 3
# gap> DotlessWernerBasisCheck(5);
# 6
# gap> DotlessWernerBasisCheck(6);
# 15
# gap> DotlessWernerBasisCheck(7);
# 36
# gap> DotlessWernerBasisCheck(8);
# 91

DotlessWernerBasisCheckImperative := n -> RankMat(List(DotlessNCPartitionsSetImperative([1..n]),compose(flatten)(WernerDiagram)));

