# Note about the representation of shifted tableaux:
#       
# SHIFTED TABLEAUX here are represented as lists of lists,
# so we need to be slightly careful about indexing.
# t[i][1] for i in [1..l] is the main diagonal (not t[i][i] as you
# might expect).

# Returns true if this is a strict partition in the sense of
# Nazarov's paper; that is, that it has distinct parts. These
# are the only types of partition that can create valid
# shifted tableaux.
IsStrictPartition := function(p)
  return Size(Set(p)) = Size(p);
end;

# Prints a shifted tableau.
PrintShiftedTableau := function(t)
  local i, j;
  for i in [1..Size(t)] do
    for j in [1..i-1] do
      Print("   ");
    od;
    for j in t[i] do
      if j < 10 then Print(" "); fi;
      Print(j);
      Print(" ");
    od;
    Print("\n");
  od;
end;

# Returns true if this is a standard shifted tableau, or if
# it could be a standard shifted tableau with some of the items
# replaced by 0 (the idea being that we can use this to find
# all the real shifted tableau by progressively filling in
# an empty diagram).
IsStandardShiftedTableau := function(tableau)
  local i, j;
    
  for i in [1..Size(tableau)] do
    for j in [1..Size(tableau[i])] do
      if tableau[i][j] <> 0 then
        if i > 1 and tableau[i-1][j+1] >= tableau[i][j] then
          return false;
        fi;
        if j > 1 and tableau[i][j-1] >= tableau[i][j] then
          return false;
        fi;
      fi;
    od;
  od;
  return true;
end;        

# Calculates the standard shifted tableau for a partition p.
StandardShiftedTableaux := function(p)
  local emptytableau, tableaulist, newtableaulist, tableau, newtableau,
        rownum, colnum, validtableau, i, n;
          
  if not IsStrictPartition(p) then
    return fail;
  fi;
  n := Sum(p);
  emptytableau := List(p, x->List([1..x], y->0));
    
  tableaulist := [ emptytableau ];

  for i in [1..n] do
    newtableaulist := [ ];
    for tableau in tableaulist do            
      for rownum in [1..Size(tableau)] do
        for colnum in [1..Size(tableau[rownum])] do
          validtableau := false;
          if tableau[rownum][colnum] = 0 then
            tableau[rownum][colnum] := i;
            if IsStandardShiftedTableau(tableau) then
              newtableau := StructuralCopy(tableau);
              tableau[rownum][colnum] := 0;
              Add(newtableaulist, newtableau);
            fi;
          fi;
        od;
      od;
    od;
    tableaulist := newtableaulist;
  od;
    
  return tableaulist;
end;        

# Computes the dimension of the representation of 2.S_n indexed
# by the partition p.
Dimension2SnRep := function(p)
  local n, l, d, m;
  if not IsStrictPartition(p) then return fail; fi;
  n := Sum(p);
  l := Size(p);
  d := (n - l) mod 2;
  m := (n - l - d) / 2;
  return Size(StandardShiftedTableaux(p)) * 2^m;
end;

NazarovMatrix := function(lambda, epsilon, k)
  local n, l, d, m,
        phi, rho, h,
        I, J, K,
        M, MM,
        i, q,
        A, B, C,
        perm, tableaux, permtableau, j;

  # p438
  n := Sum(lambda);
  l := Size(lambda);
  d := (n - l) mod 2;
  m := (n - l - d) / 2;

  # p439 (top): computation of h(\Lambda, k)
  h := function(t, k)
    local i, j;
    for i in [1..Size(t)] do;
      for j in [1..Size(t[i])] do;
        # N.B. not j - i, because main diag is t[i][1].
        if t[i][j] = k then return j - 1; fi;
      od;
    od;
  end;

  # p439 (bottom): construction of the matrices M_i, which
  # we put in the array MM[]
  I := [[0, -E(4)], [E(4), 0]];
  J := [[1,     0], [0,   -1]];
  K := [[0,     1], [1,    0]];

  MM := [];

  # M_1
  M := [[epsilon]];
  for i in [1..m] do;
    M := Tensor(M, J);
  od;
  Add(MM, M);  

  for q in [1..m] do;
    # M_2q, q in 1..m
    M := epsilon * Tensor(IdentityMat(2*(q-1)), K);
    if q = 1 then M := epsilon * K; fi;
    for i in [1..m-q] do;
      M := Tensor(M, J);
    od;
    Add(MM, M);

    # M_{2q+1}, q in 1..m
    M := epsilon * Tensor(IdentityMat(2*(q-1)), I);
    if q = 1 then M := epsilon * I; fi;
    for i in [1..m-q] do;
      M := Tensor(M, J);
    od;
    Add(MM, M);
  od;

  # p440
  phi := function(p, q)
    return Sqrt(2*q*q+1) / ((p-q) * (p+q+1));
  end;

  rho := function(p, q)
    return Sqrt(1/2 * (1 - 1/(p-q)^2) * (1-1/(p+q+1)^2));
  end;

  # Basis for U_\Lambda (stanard shifted tableaux). We need to be able to
  # decompose into subspaces depending on the action of the
  # permutation (k, k+1).
  tableaux := StandardShiftedTableaux(lambda);
  perm := (k, k+1);
  for i in [1..Size(tableaux)] do;
    permtableau := ActionOnTableau(tableaux[i], perm);
    if IsStandardShiftedTableau(permtableau) then
      j := Position(tableaux, permtableau);
      Print(i, " -> ", j, "\n");
    fi;
  od;

end;
