ImportReference := function(groupname)
  local rep, path;
  rep := rec();
  rep.groupname := groupname;
  path := Concatenation("../data/", groupname);
  rep.a := ImportMtx(Concatenation(path, "/RefA"));
  rep.b := ImportMtx(Concatenation(path, "/RefB"));
  rep.r := ImportMtx(Concatenation(path, "/RefR"));
  rep.z := ImportMtx(Concatenation(path, "/RefZ"));
  rep.group := Group(rep.a, rep.b, rep.r, rep.z);
  return rep;
end;
