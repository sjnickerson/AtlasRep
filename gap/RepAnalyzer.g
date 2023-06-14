RepAnalyzer := function(dir, files)
  local analyze;

  analyze := function(group)
    local g, n;
    for n in [2..Size(GeneratorsOfGroup(group))] do;
      g := Group(List([1..n], x->GeneratorsOfGroup(group)[x]));
      Print("  ZEntryBound<g1..g", n, "> = ",
          ZEntryBound(group).zentrybound, "\n");
    od;
    return true;
  end;

  ForEach(dir, files, analyze);
end;       
