AtlasZEntryBound := function(groupname, maxdim)
  local G, reps, zentry;
  reps := AllAtlasGeneratingSetInfos(groupname, IsMatrixGroup, Ring, Integers, Dimension, [1..maxdim]);
  for rep in reps do;
    Print(rep.repname);
    Print(": ");
    G := Group(AtlasGenerators(groupname, rep.repnr).generators);
    zentry := ZEntryBound(G);
    Print(zentry.zentrybound);
    Print("; ");
    Print(zentry.suborbits);
    Print("\n");
  od;
end;
