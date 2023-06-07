ImportMtx := function(file)
  local text, line, params, one, fieldorder,
    type, rows, cols, row, i, mat;
  text := InputTextFile(file);
  if text = fail then
    return fail;
  fi;
  line := ReadLine(text);
  params := List(SplitString(line, "", " \t\n"), Int);
  type := params[1];
  if type = 1 then;
    fieldorder := params[2];
    if not IsPrime(fieldorder) then;
      return fail;
    fi;
    one := Z(fieldorder)^0;
  elif type = 5 then;
    one := 1;
  else
    return fail;
  fi;
  rows := params[3];
  cols := params[4];
  mat := [];
  for i in [1..rows] do;
    row := ReadLine(text);
    NormalizeWhitespace(row);
    if ' ' in row then;
      row := SplitString(row, " ", " \t\n");
      Add(mat, List(row, x -> Int(x) * one));
    else;
      Add(mat, List(row, x -> Int([x]) * one));
    fi;
  od;
  CloseStream(text);
  return mat;
end;
