CliffordEven := function(n)

  local epsilon, x, y, one, k, allxi, xi, i, j, m, mat;

  if n mod 2 <> 0 then
   return fail;
  fi;
  k := n / 2;

  epsilon := [[1,0], [0,-1]];
  x := [[0,1], [1,0]];
  y := [[0, E(4)], [-E(4), 0]];
  one := IdentityMat(2);

  allxi := [];
  for m in [1..n] do;
    xi := [[1]];
    j := Int((m + 1)/2);
    for i in [1..k] do;
      if i < j then
        mat := epsilon;
      elif i = j and 2*j - 1 = m then
        mat := x;
      elif i = j and 2*j = m then
        mat := y;
      elif i > j then
        mat := one;
      else
        return fail;
      fi;
      xi := Tensor(xi, mat);
    od;
    Add(allxi, xi);
  od;
  return allxi;
end;

SpinRepresentationEven := function(n)
  local xi;

  xi := CliffordEven(n);
  return List([1..n-1], j -> E(4)/Sqrt(2)*(xi[j] - xi[j+1]));
end;       

SpinRepresentationEvenStandardGenerators := function(n)
  local t, a, b, r, ans;
  t := SpinRepresentationEven(n);
  ans := rec();
  ans.a := t[2]*t[1];
  ans.b := Product(List(Reversed([2..n-1]), i -> t[i]));
  ans.r := t[1];
  return ans;
end;
