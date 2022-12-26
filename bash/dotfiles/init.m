(** User Mathematica initialization file **)

AppendTo[$Path, "/home/jarinf/.Mathematica/SciDraw/packages"]

(* Checks if an axis number is negative, and if so applies a bar over the positive value for use in axis labeling *)
checkNegNum[num_] := Module[{str},
  If[num < 0, str = ToString[
\!\(\*OverscriptBox[\(Abs[num]\), \(_\)]\),
     FormatType -> StandardForm], str = ToString[num]];
  str]

CreateBCC[nx_, ny_, nz_, a_] := Module[{n, m, x, y, z, data},
   n = 2*nx*ny*nz;
   x = Table[0, {i, 1, n}];
   y = Table[0, {i, 1, n}];
   z = Table[0, {i, 1, n}];
   x[[2]] = a/2;
   y[[2]] = a/2;
   z[[2]] = a/2;
   m = 0;
   For[i = 0, i < nx, i++, {
     For[j = 0, j < ny, j++, {
       For[k = 0, k < nz, k++, {
         For[l = 0, l < 2, l++, {
            x[[l + m + 1]] = x[[l + 1]] + a*k;
            y[[l + m + 1]] = y[[l + 1]] + a*j;
            z[[l + m + 1]] = z[[l + 1]] + a*i;
            };
          ];(*l loop*)
         m = m + 2;
         }](*k loop*)
       }](*j loop*)
     }];(*i loop*)

   data = Transpose[{x, y, z}]];
CreateFCC[nx_, ny_, nz_, a_] := Module[{n, m, x, y, z, data},
   n = 4*nx*ny*nz;
   x = Table[0, {i, 1, n}];
   y = Table[0, {i, 1, n}];
   z = Table[0, {i, 1, n}];
   x[[2]] = a/2;
   y[[2]] = a/2;
   y[[3]] = a/2;
   z[[3]] = a/2;
   x[[4]] = a/2;
   z[[4]] = a/2;
   m = 0;
   For[i = 0, i < nx, i++, {
     For[j = 0, j < ny, j++, {
       For[k = 0, k < nz, k++, {
         For[l = 0, l < 4, l++, {
            x[[l + m + 1]] = x[[l + 1]] + a*k;
            y[[l + m + 1]] = y[[l + 1]] + a*j;
            z[[l + m + 1]] = z[[l + 1]] + a*i;
            };
          ];(*l loop*)

         m = m + 4;}](*k loop*)
       }](*j loop*)
     }];(*i loop*)

      data = Transpose[{x, y, z}]];

(* Compares the rotated coordinate axes to the original <100> coordinate axes *)
compareAxes[newX_, newY_, newZ_] := Legended[Graphics3D[{
    {Red, Arrow[{{0, 0, 0}, {1, 0, 0}}]}, {Red, Dashed,
     Arrow[{{0, 0, 0}, Normalize[newX]}]},
    {Blue, Arrow[{{0, 0, 0}, {0, 1, 0}}]}, {Blue, Dashed,
     Arrow[{{0, 0, 0}, Normalize[newY]}]},
    {Green, Arrow[{{0, 0, 0}, {0, 0, 1}}]}, {Green, Dashed,
     Arrow[{{0, 0, 0}, Normalize[newZ]}]},
     Text["100", {1.1, 0, 0}],
    Text[checkNegNum[newX[[1]]] <> checkNegNum[newX[[2]]] <>
      checkNegNum[newX[[3]]], {(newX[[1]] + newX[[1]]/10)/
      Norm[newX], (newX[[2]] + newX[[2]]/10)/Norm[newX], (
      newX[[3]] + newX[[3]]/10)/Norm[newX]}],
    Text["010", {0, 1.1, 0}],
    Text[checkNegNum[newY[[1]]] <> checkNegNum[newY[[2]]] <>
      checkNegNum[newY[[3]]], {(newY[[1]] + newY[[1]]/10)/
      Norm[newY], (newY[[2]] + newY[[2]]/10)/Norm[newY], (
      newY[[3]] + newY[[3]]/10)/Norm[newY]}],
    Text["001", {0, 0, 1.1}],
    Text[checkNegNum[newZ[[1]]] <> checkNegNum[newZ[[2]]] <>
      checkNegNum[newZ[[3]]], {(newZ[[1]] + newZ[[1]]/10)/
      Norm[newZ], (newZ[[2]] + newZ[[2]]/10)/Norm[newZ], (
      newZ[[3]] + newZ[[3]]/10)/Norm[newZ]}]},
    Boxed -> False, ViewVertical -> newZ, ImageSize -> Large],
  SwatchLegend[{Red, Blue, Green}, {"x axis", "y axis", "z axis"}]]

(* Rotates a vector array by an angle theta about the z axis.  This is already implemented within Mathematica however *)
rotateVector[vecArray_, \[Theta]_] := Module[{xtemp, ytemp, newVec},
  newVec = Reap[Do[
       Sow[{(*xtemp=*)
         Cos[\[Theta]]*vecArray[[i]][[1]] -
          Sin[\[Theta]]*vecArray[[i]][[2]],
         (*ytemp=*)
         Sin[\[Theta]]*vecArray[[i]][[1]] +
          Cos[\[Theta]]*vecArray[[i]][[2]],
         vecArray[[i]][[3]]}], {i, 1, Length[vecArray]}]
      ][[2]][[1]]]

(* Calculates the sigma number for the specified hkl *)
calculateSigma[h_, k_, l_] := Module[{sumsq, sigNum},
  If[h == k == l == 0, Print["Invalid: hkl = 000"]; Return[0]];
  sumsq = h^2 + k^2 + l^2;
  sigNum = sumsq;
  While[EvenQ[sigNum], sigNum = sigNum/2];
  sigNum
  ]

(*box is a vector of 3 elements - length in x, length in y, and length in z*)
generateCellLinkedList[atoms_, box_, rCut_] :=
 Module[{rCutSq, ncell, id, lcell, nAtomsPerCell, drijSq, rij, icell,
   pcell, tmp1, tmp2, ia, ja, ka, iatom, n, i, j, k, l, ii, jj, kk, m},
  rCutSq = rCut^2;
  ncell = Floor[box/rCut] + 1;
  lcell = box/ncell;
  n = Length[atoms];
  nAtomsPerCell =
   Max[Floor[Length[atoms] / (ncell[[1]]*ncell[[2]]*ncell[[3]])], 200];
  icell =
   Table[Table[Table[0, {i, 1, ncell[[1]]}], {j, 1, ncell[[2]]}], {k,
     1, ncell[[3]]}];
  pcell =
   Table[Table[
     Table[Table[0, {i, 1, ncell[[1]]}], {j, 1, ncell[[2]]}], {k, 1,
      ncell[[3]]}], {l, 1, nAtomsPerCell}];
  iatom =
   Table[Table[1, {i, 1, Length[atoms]}], {j, 1, nAtomsPerCell}];

  For[i = 1, i <= n, i++,
   {id = Floor[atoms[[i]]/lcell] + 1;
    For[j = 1, j <= 3, j++,
     If[id[[j]] > ncell[[j]],
      id[[j]] =
       ncell[[j]]]] (*make sure we don't go past memory bounds*);
    For[j = 1, j <= 3, j++, If[id[[j]] == 0, id[[j]] = 1]];
    icell[[ id[[1]], id[[2]], id[[3]] ]] =
     icell[[ id[[1]], id[[2]], id[[3]] ]] + 1;
    pcell[[ id[[1]], id[[2]], id[[3]],
      icell[[ id[[1]], id[[2]], id[[3]] ]] ]] = i;
    }];

  For[i = 1, i <= ncell[[1]], i++,
   {For[j = 1, j <= ncell[[2]], j++,
      {For[k = 1, k <= ncell[[3]], k++,
         {For[l = 1, l <= icell[[i, j, k]], l++,
            {tmp1 = pcell[[i, j, k, l]];
             For[ii = 0, ii <= 2, ii++,
              {For[jj = 0, jj <= 2, jj++,
                 {For[kk = 0, kk <= 2, kk++,
                    {ia = i + ii;
                    ja = j + jj;
                    ka = k + kk;
                    If[ia > ncell[[1]], ia = 1];
                    If[ja > ncell[[2]], ja = 1];
                    If[ka > ncell[[3]], ka = 1];
                    If[ia < 1, ia = ncell[[1]]];
                    If[ja < 1, ja = ncell[[2]]];
                    If[ka < 1, ka = ncell[[3]]];

                    For[m = 1, m <= icell[[ia, ja, ka]], m++,
                    {tmp2 = pcell[[ia, ja, ka, m]];
                    If[tmp2 <= tmp1, Continue];
                    rij = atoms[[tmp1]] - atoms[[tmp2]];
                    rij = rij - Round[rij/box]*box;

                    drijSq = Sum[rij[[z]]^2, {z, 1, 3}];

                    If[drijSq > rCutSq, Continue];

                    iatom[[1, tmp1]] += 1;
                    iatom[[iatom[[1, tmp1]], tmp1]] = tmp2;
                    iatom[[1, tmp2]] += 1;
                    iatom[[iatom[[1, tmp2]], tmp2]] = tmp1;
                    }];
                    }];
                  }];
               }];
             }];
          }];
       }];
    }];
  iatom
  ]
