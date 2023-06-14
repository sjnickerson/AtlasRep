Tensor := function(a, b)
  local c, v, arow, aelt, brow, belt;
  c := [];
  for arow in a do;
    for brow in b do;
      v := [];
      for aelt in arow do;
        for belt in brow do;
          Add(v, aelt * belt);
        od;
      od;
      Add(c, v);
    od;
  od;     
  return c;   
end;
