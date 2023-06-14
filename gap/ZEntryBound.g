# ZEntryBound
# Computes an absolute bound for finite matrix groups over Z.
ZEntryBound := function(G)
  local i, j, g, gens, n,
        v, v2,
        vecs, allvecs, newvecs,
        unitvec, isunitvec, maxabsvalue,
        printwhen,
        startpoints,
        suborbits,
        zentrybound,
        result;

  # Set up helper functions

  # Returns the unit vector e_i
  unitvec := function(i)
    local v;
    v := List([1..n], x -> 0);
    v[i] := 1;
    return v;
  end;

  # Returns 0 if not a unit vector, otherwise
  # returns j for which v = e_j.
  isunitvec := function(v)
    local j, i;
    j := 0;
    for i in [1..n] do;
      if v[i] = 1 then;
        if j > 0 then;
          return 0;
        fi;
        j := i;
      elif v[i] <> 0 then;
        return 0;
      fi;
    od;
    return j;
  end;

  # Returns the maximum absolute value in a Z-vector
  maxabsvalue := function(v)
    return Maximum(List(v, AbsInt));
  end;

  # Check some preconditions
  if not IsMatrixGroup(G) then;
    return fail;
  fi;
  if not (FieldOfMatrixGroup(G) = Rationals) then;
    return fail;
  fi;

  gens := GeneratorsOfGroup(G);
  n := Size(gens[1]);

  # For each unit vector, this points to the unit
  # vector known to be in the same G orbit, or 0
  # if it's not yet known.
  startpoints := List([1..n], x -> 0);
  suborbits := [];

  # The maximum absolute value, which we'll update
  # as we go through.
  zentrybound := 0;

  for i in [1..n] do;
    if startpoints[i] <> 0 then
      continue;
    fi;
    startpoints[i] := i;

    # All vectors we've seen so far.
    allvecs := NewDictionary(unitvec(i), false);
    # All vectors we've seen on this iteration.
    newvecs := NewDictionary(unitvec(i), false);
    AddDictionary(newvecs, unitvec(i));

    repeat
      # The vectors we're going to iterate through.
      vecs := ShallowCopy(newvecs);
      newvecs := NewDictionary(unitvec(i), false);

      for g in gens do;
        for v in vecs do;
          v2 := v*g;
          j := isunitvec(v2);
          if j > 0 then;
            if startpoints[j] = 0 then;
              startpoints[j] := i;
            fi;
          fi;

          if not KnowsDictionary(vecs, v2) and not KnowsDictionary(allvecs, v2) then;
            AddDictionary(newvecs, v2);
          fi;
        od;
      od;

      for v in newvecs do;
        AddDictionary(allvecs, v);
      od;
    until IsEmpty(newvecs);

    Add(suborbits, Size(allvecs));
    zentrybound := Maximum(zentrybound, Maximum(List(allvecs, maxabsvalue)));
  od;

  result := rec();
  result.zentrybound := zentrybound;
  result.suborbits := suborbits;
  return result;
end;
