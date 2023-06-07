ExportMtx := function(file, mat)
  local row, elt, fieldorder, stream;
  stream := OutputTextFile(file, false);
  if IsFFE(mat[1][1]) then
    fieldorder := Size(Field(mat[1][1]));
    PrintTo(stream, "1 ", fieldorder, " ", Size(mat), " ", Size(mat[1]), "\n");
    for row in mat do;
      for elt in row do;
        PrintTo(stream, IntFFE(elt));
      od;
      PrintTo(stream, "\n");
    od;
  else;
    PrintTo(stream, "5 1 ", Size(mat), " ", Size(mat[1]), "\n");
    for row in mat do;
      for elt in row do;
        PrintTo(stream, elt, " ");
      od;
      PrintTo(stream, "\n");
    od;
  fi;
  CloseStream(stream);
end;
