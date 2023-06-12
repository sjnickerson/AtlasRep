SymReps := function(n, directory)
  local partitions, p, comment, splits, partitioninfo, makestdgens, y, gens,
    file, dim, traces, bclass, getnext, nextletter;

  makestdgens := function()
    local a, b, r;
    a := (1,2,3);
    if n mod 2 = 0 then;
      b := PermList(Concatenation([1], [3..n], [2]));
    else;
      b := PermList(Concatenation([1], [2], [4..n], [3]));
    fi;
    r := (1,2);
    return [a, b, r];
  end;

  partitioninfo := function(p)
    local res;
    res := rec();
    res.sndim := HookLengthFormula(p);
    res.parts := Set([p, DualPartition(p)]);
    res.splits := Size(res.parts) = 2;
    if res.splits then;
      res.andim := res.sndim;
    else;
      res.andim := res.sndim / 2;
    fi;
    return res;
  end;

  partitions := Set(Filtered(List(Partitions(n), partitioninfo), x -> x.andim > 1));
  SortBy(partitions, x -> x.andim);

  traces := function(dim, classname)
    local ct, chars, charindex;
    ct := CharacterTable(StringFormatted("A{1}", n));
    chars := Filtered(Irr(ct), x -> x[1] = dim);
    charindex := Position(ClassNames(ct), classname);
    return List(chars, char -> char[charindex]);
  end;

  gens := makestdgens();

 
  nextletter := [];
  getnext := function(dim, num)
    local i, str;
    str := [];
    if not IsBound(nextletter[dim]) then;
      nextletter[dim] := 'a';
    fi;
    for i in [1..num] do;
      Add(str, nextletter[dim]);
      nextletter[dim] := CharInt(IntChar(nextletter[dim]) + 1);
    od;
    return str;
  end;

  CreateDir(directory);
  for p in partitions do;
    if not p.splits then;
      comment := StringFormatted("fusing to size {1} in S{2}", p.sndim, n);
      file := StringFormatted("Z{2}{3}", n, p.andim, getnext(p.andim, 2));
    else;
      comment := "";
      file := StringFormatted("Z{2}{3}", n, p.andim, getnext(p.andim, 1));
    fi;

    y := YoungRepresentation(p.parts[1], gens);
    PrintFormatted("{1} - partition {2} {3}, Tr(A)={4}, Tr(B)={5}\n", file, p.parts[1], comment, Trace(y[1]), Trace(y[2])); 
    if p.splits then;
      if Trace(y[3]) < 0 then;
        PrintFormatted("  Tr(R)={1} < 0; replacing with -R\n", Trace(y[3]));
        y[3] := -y[3];
      elif Trace(y[3]) = 0 then;
        Print("  WARNING: Tr(R)=0; not sure if this is the representation in the ATLAS\n");
      fi;
    fi;
    ExportMtx(Concatenation(directory, "/", file, "-A"), y[1]);
    ExportMtx(Concatenation(directory, "/", file, "-B"), y[2]);
    ExportMtx(Concatenation(directory, "/", file, "-R"), y[3]);
  od;

  for dim in Set(List(partitions, p -> p.andim)) do;
    if Size(Filtered(partitions, p -> p.andim = dim)) > 1 then;
      PrintFormatted("WARNING: Multiple (equivalence classes of) reps of A{1} of dimension {2}\n", n, dim);      
      PrintFormatted("  Traces on 3a: {1}\n", traces(dim, "3a"));
      bclass := StringFormatted("{1}a", Order(gens[2]));
      PrintFormatted("  Traces on {1}: {2}\n", bclass, traces(dim, bclass));
    fi;
  od;

  return gens;
end;

