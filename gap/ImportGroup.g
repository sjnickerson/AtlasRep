ImportGroup := function(file)
  local ext, filename, mat, mats;
  mats := [];
  for ext in ["A", "B", "Z", "R"] do;
    filename := Concatenation(file, "-", ext);
    if IsExistingFile(filename) then
      mat := ImportMtx(filename);
      Add(mats, mat);
    fi;
  od;
  return Group(mats);
end;
