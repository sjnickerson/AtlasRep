###########################################################################
# Young Natural Representation of the Symmetric Group
###########################################################################

# Note:
# Partitions are represented by lists of numbers.
# Tableaux are represented by lists of rows, where each row is a list.
# Few checks are made to test whether tableaux and partitions are
# actually valid data.

###########################################################################

# -------------------------------------------------------------------------
# Test whether a tableau is standard, i.e. the entries always increase
# as row and column number increase.
IsStandardTableau := function(tableau)
    
    local i, j;
    
    for i in [1..Size(tableau)] do
        for j in [1..Size(tableau[i])] do
            if i > 1 and tableau[i-1][j] >= tableau[i][j] then
                return false;
            fi;
            if j > 1 and tableau[i][j-1] >= tableau[i][j] then
                return false;
            fi;
        od;
    od;
    return true;
    
end;

# -------------------------------------------------------------------------
# List up to [max] standard tableaux corresponding to the [partition]
StandardTableaux := function(partition, max)

    local emptytableau, tableaulist, newtableaulist, tableau, newtableau,
          rownum, colnum, validtableau, i, n;
    
    n := Sum(partition);
    
    emptytableau := List(partition, x->List([1..x], y->0));
    
    tableaulist := [ emptytableau ];
    
    for i in [1..n] do
        newtableaulist := [ ];
        for tableau in tableaulist do
            
            if max > 0 and Size(newtableaulist) >= max then
                break;
            fi;
            
            for rownum in [1..Size(tableau)] do
                for colnum in [1..Size(tableau[rownum])] do
                    validtableau := false;
                    if tableau[rownum][colnum] = 0 then
                        
                        if rownum = 1 then
                            if colnum = 1 then
                                validtableau := true;
                            elif tableau[1][colnum-1] in [1..i-1] then
                                validtableau := true;
                            fi;
                        else
                            if colnum = 1 then
                                if tableau[rownum-1][1] in [1..i-1] then
                                    validtableau := true;
                                fi;
                            else
                                if tableau[rownum-1][colnum] in [1..i-1] and
                                   tableau[rownum][colnum-1] in [1..i-1] then
                                    
                                    validtableau := true;
                                fi;
                            fi;
                        fi;
                    fi;
                    
                    if validtableau then
                        newtableau := StructuralCopy(tableau);
                        newtableau[rownum][colnum] := i;
                        Add(newtableaulist, newtableau);
                    fi;
                od;
            od;
        od;
        tableaulist := newtableaulist;
    od;
    
    return tableaulist;
    
end;

# -------------------------------------------------------------------------
# Calculate the size of the Specht module corresponding to a given
# partition by means of the hook length formula.
HookLengthFormula := function(partition)
    
    local rownum, colnum, hooklength, product;
    
    product := 1;
    
    for rownum in [1..Size(partition)] do
        for colnum in [1..partition[rownum]] do
            hooklength := partition[rownum] - colnum +
                          Number([rownum..Size(partition)],
                                 x->partition[x] >= colnum);
            product := product * hooklength;
        od;
    od;
    
    return Factorial(Sum(partition)) / product;
    
end;

# -------------------------------------------------------------------------
# Calculate the dual of a partition; i.e. interchange rows and
# columns in the Ferrers diagram
DualPartition := function(partition)
    
    return List([1..Maximum(partition)], x->Number(partition, y->y >= x));
    
end;

# -------------------------------------------------------------------------
# Calculate the dual of a tableau; i.e. interchange rows and columns
DualTableau := function(tableau)
    
    local dual, i, j;
    
    dual := List([1..Size(tableau[1])], x->[]);
    
    for i in [1..Size(tableau)] do
        for j in [1..Size(tableau[i])] do
            dual[j][i] := tableau[i][j];
        od;
    od;
    
    return dual;
    
end;

# -------------------------------------------------------------------------
# Given two partitions, decide whether partition1 >= partition2 in
# the partial order of row domination
PartitionRowDominates := function(partition1, partition2)
    
    local i, sum1, sum2;
    
    sum1 := 0;
    sum2 := 0;
    
    for i in [1..Size(partition1)] do
        sum1 := sum1 + partition1[i];
        sum2 := sum2 + partition2[i];
        if sum2 > sum1 then
            return false;
        fi;
    od;
    
    return true;
    
end;

# -------------------------------------------------------------------------
# Given two tableaux, let the corresponding tabloids be T1 and T2.
# Decide whether T1 >= T2 in the partial order of row domination
# (for tabloids)
TabloidRowDominates := function(tableau1, tableau2)
    
    local i, n, partition1, partition2;
    
    n := Maximum(Flat(tableau1));
    
    for i in [1..n] do
        partition1 := List(tableau1, x->Number(x, y->y <= i));
        partition2 := List(tableau2, x->Number(x, y->y <= i));
        
        if not PartitionRowDominates(partition1, partition2) then
            return false;
        fi;
        
    od;
    
    return true;
    
end;

# -------------------------------------------------------------------------
# Given two partitions, decide whether partition1 >= partition2 in
# the partial order of column domination
PartitionColumnDominates := function(partition1, partition2)
    
    return PartitionRowDominates(
                   DualPartition(partition1),
                   DualPartition(partition2));
    
end;

# -------------------------------------------------------------------------
# Given two tableaux, let the corresponding tabloids be T1 and T2.
# Decide whether T1 >= T2 in the partial order of column domination
# (for tabloids)
TabloidColumnDominates := function(tableau1, tableau2)
    
    return TabloidRowDominates(
                   DualTableau(tableau1),
                   DualTableau(tableau2));
    
end;

# -------------------------------------------------------------------------
# Print a tableau to the screen, lining up the elements neatly
PrintTableau := function(tableau)
    local row, entry;
    for row in tableau do
        for entry in row do
            if entry < 10 then Print(" "); fi;
            Print(entry);
            Print(" ");
        od;
        Print("\n");
    od;
end;

# -------------------------------------------------------------------------
# Given a tableau, sort the columns of the tableau, and return
# the sorted tableau along with a permutation which performs
# the sorting.
SortTableauColumns := function(tableau)
    
    local perm, dual, sort, i, n;
    
    dual := DualTableau(tableau);
    sort := List(dual, SortedList);
    
    perm := MappingPermListList(Flat(dual), Flat(sort));
    
    return rec(tableau := DualTableau(sort), perm := perm);
    
end;

# -------------------------------------------------------------------------
# Find a row descent for a tableau T, i.e. co-ordinates i and j
# such that T(i,j) > T(i+1, j)
FindRowDescent := function(tableau)
    
    local i, j;
    
    for i in [1..Size(tableau)] do
        for j in [1..Size(tableau[i]) - 1] do
            if tableau[i][j] > tableau[i][j+1] then
                return rec(i := i, j := j);
            fi;
        od;
    od;
    return fail;
    
end;

# -------------------------------------------------------------------------
# Compute the action of a permutaton on a tableau
ActionOnTableau := function(tableau, perm)
    return List(tableau, x->List(x, y->y^perm));
end;

# -------------------------------------------------------------------------
# Calculate a (standard) Garnir element g_AB = \sum_p sgn(p) p where
# p is a transversal of Sym(A) x Sym(B) in Sym(A union B).
# This element is used to remove row descents, because
# g_AB(e_t) = 0 if A, B are tableau elements below/above the row descent
# on the left/right respectively. See Sagan, The Symmetric Group (1st ed)
# chapter 2, Proposition 2.6.3
GarnirPermutations := function(A, B)
    
    local C, perm, permlist, AA, BB, CC;
    
    C := Concatenation(A, B);
    permlist := [];
    
    for AA in Combinations(C, Size(A)) do
        BB := Difference(C, AA);
        CC := Concatenation(AA, BB);
        perm := MappingPermListList(C, CC);
        if not perm = () then
            Add(permlist, perm);
        fi;
    od;
    
    return permlist;
end;

# -------------------------------------------------------------------------
# Given a tableau with a row descent, calculate a sum of tableaux which
# give the same sum of polytabloids and which do not have the
# row descent. This is achieved by acting on the tableau by a Garnir
# permutation.
RemoveRowDescentFromTableau := function(tableau, rowdescent)
    
    local A, B, height, garnir, tableaulist, coefficientlist;
    
    height := Number(tableau, x->Size(x) >= rowdescent.j);
    
    A := List([rowdescent.i..height], x->tableau[x][rowdescent.j]);
    B := List([1..rowdescent.i], x->tableau[x][rowdescent.j + 1]);
    garnir := GarnirPermutations(A, B);
    tableaulist := List(garnir, perm->ActionOnTableau(tableau, perm));
    coefficientlist := List(garnir, x-> -SignPerm(x));
    
    return rec(tableaux := tableaulist, coeffs := coefficientlist);
    
end;

# -------------------------------------------------------------------------
# Find a sum of standard tableaux which gives the same sum of
# polytabloids as a given sum of (not necessarily standard) tableaux.
# Sums of tableaux are handled as lists of tableaux and lists
# of coefficients.
FindStandardTableauSum := function(tableaulist, coefficientlist)
    local i, j, sortedtableau, newtableau, newcoefficient,
          newtableaulist, newcoefficientlist, tableau, coefficient,
          rowdescent, answer, rowdescenttableaulist,
          rowdescentcoefficientlist,
          recursiveanswer, pos, recursivetableaulistsize, newtableaulistsize;
    
    newtableaulist := [ ];
    newcoefficientlist := [ ];
    
    for i in [1..Size(tableaulist)] do
        tableau := tableaulist[i];
        coefficient := coefficientlist[i];
        sortedtableau := SortTableauColumns(tableau);
        newtableau := sortedtableau.tableau;
        newcoefficient := coefficient * SignPerm(sortedtableau.perm);
        if IsStandardTableau(newtableau) then
            pos := Position(newtableaulist, newtableau);
            if pos = fail then
                Add(newtableaulist, newtableau);
                Add(newcoefficientlist, newcoefficient);
            else
                newcoefficientlist[pos] := newcoefficientlist[pos] +
                                              newcoefficient;
            fi; 
        else
            rowdescent := FindRowDescent(newtableau);
            answer := RemoveRowDescentFromTableau(newtableau, rowdescent);
            rowdescenttableaulist := answer.tableaux;
            rowdescentcoefficientlist := newcoefficient * answer.coeffs;
            recursiveanswer := FindStandardTableauSum(
                                       rowdescenttableaulist,
                                       rowdescentcoefficientlist);
            recursivetableaulistsize := Size(recursiveanswer.tableaulist);
            for j in [1..recursivetableaulistsize] do
                pos := Position(newtableaulist,
                                recursiveanswer.tableaulist[j]);
                
                if pos = fail then
                    Add(newtableaulist, recursiveanswer.tableaulist[j]);
                    Add(newcoefficientlist, recursiveanswer.coefficientlist[j]);
                else
                    newcoefficientlist[pos] := newcoefficientlist[pos] +
                                                  recursiveanswer.
                                                  coefficientlist[j];
                                                
                fi;
            od;
        fi;
        
    od;
    
    return rec(tableaulist := newtableaulist,
               coefficientlist := newcoefficientlist);

end;

# -------------------------------------------------------------------------
# Calculate the images of some permutations in the Young representation of
# indexed by a given partition.
YoungRepresentation := function(partition, permlist)
    
    local size, k, row, gen, gens, n, tableaux, i, j, m, tableausum, pos;
    
    n := Sum(partition);
    size := HookLengthFormula(partition);
    tableaux := StandardTableaux(partition, 0);
    m := Size(tableaux);
    gens := [];
    
    for k in [1..Size(permlist)] do
        gen := [];
        for i in [1..m] do
            tableausum := FindStandardTableauSum(
                                  [ ActionOnTableau(tableaux[i], permlist[k]) ],
                                  [ 1 ] );
            row := List([1..m], x->0);            
            for j in [1..Size(tableausum.tableaulist)] do
                pos := Position(tableaux, tableausum.tableaulist[j]);
                row[pos] := tableausum.coefficientlist[j];
            od;
            Add(gen, row);
        od;
        Add(gens, gen);
    od;
    
    return gens;
    
end;

# --------------------------------------------------------------------------
# Calculate a Young representation of the symmetric group S_n, presented
# by the Dynkin diagram for the lattice A_(n-1).
YoungRepresentationByTranspositions := function(partition)
    
    local n, permlist;
    
    n := Sum(partition);
    permlist := List([1..n-1], x->(x,x+1));
    
    return YoungRepresentation(partition, permlist);
    
end;

# --------------------------------------------------------------------------
# Calculate a Young representation of the symmetric group S_n, generated
# by the elements (1 2) and (2 3 ... n) (which are ATLAS standard generators).
YoungRepresentationByStandardGeneratorsSn := function(partition)
    
    local n, perm;
    
    n := Sum(partition);
    perm := PermList(Concatenation([1], [3..n], [2]));
    
    return YoungRepresentation(partition, [ (1,2), perm ] );
end;

# --------------------------------------------------------------------------
# Test whether a set of generators satisfies the A_(n-1) Dynkin diagram
# presentation for a symmetric group.
SatisfiesSnDynkinDiagramPresentation := function(gens)
    local i, j, answer, id;
    answer := true;
    id := IdentityMat(Size(gens[1]));
    for i in [1..Size(gens)-1] do
        if not gens[i]^2 = id then
            answer := false;
            Print("Generators ", i, " does not square to 1\n");
        fi;
    od;
    for i in [1..Size(gens)-1] do
        if not (gens[i]*gens[i+1])^3 = id then
            answer := false;
            Print("The product of generator ", i, " with generator ", i+1,
                  "does not cube to 1\n");
        fi;
    od;
    for i in [1..Size(gens)] do
        for j in [i+2..Size(gens)] do
            if not gens[i] * gens[j] = gens[j] * gens[i] then
                answer := false;
                Print("Generator ", i, " does not commute with generator ", j);
            fi;
        od;
    od;
    return answer;
end;



PolytabloidForTableau := function(tableau)
    
    local elt, col, cols, gens, gps,
          eltlist, tableaulist, signlist;
    
    cols := DualTableau(tableau);
    gps := [];
    for col in cols do
        Add(gps, SymmetricGroup(col));
    od;
    eltlist := Cartesian(gps);
    tableaulist := List(eltlist, elt->ActionOnTableau(tableau, Product(elt)));
    signlist := List(eltlist, elt->SignPerm(Product(elt)));
    
    return rec(tableaux := tableaulist, coeffs := signlist);
    
end;
