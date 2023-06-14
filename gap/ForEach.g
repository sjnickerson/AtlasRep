ForEach := function(dir, files, fn)
  local file, g;
  for file in files do;
    g := ImportGroup(Concatenation(dir, "/", file));
    Print(file, "\n");
    CallFuncList(fn, [g]);
    Print("\n");
  od;
end;        
