#!/usr/bin/env wolframscript

(* 
Author : Suzie Tallon & Thomas Lenormand 
Date : February 2026
This script is used to detect change points in topology data using likelihood-based models.
It identifies clusters of positions with high likelihood, estimates BIC for 0,1,2,3,4,5 change point models, and exports both model summaries and position-specific segments for downstream analyses.
*)

(* ---------------- Basic helper functions ---------------- *)

T[x_] := Transpose[x]; (* Transpose a matrix/list *)
L[x_] := Length[x]; (* Length of a list *)

minslice = 10; (* Minimum size for a segment *)

(* Compute final list of positions to keep, given clusters and margin *)
ListFinal[clusters_List, marge_?IntegerQ] := Module[
  {ranges, allRanges, result}, 
  allRanges = 
   Range[Min[Flatten[clusters] - (marge - 1)], 
    Max[Flatten[clusters] + (marge - 1)]];
  ranges = 
   Table[Range[Max[clusters[[i]]] + marge, 
     Min[clusters[[i + 1]]] - marge], {i, 1, (clusters // L) - 1}];
  result = Complement[allRanges, Sequence @@ ranges];
  result = Select[result, # > minslice &];
  Select[result, # < maxpos - minslice &]
  ]

(* Retrieve contiguous positions around a maximum that exceed a threshold *)
getContig[lst_, maxPos_, seuil_, dir_] := 
  Module[{i = maxPos, idx = {}}, 
   While[1 <= i <= Length[lst] && lst[[i]] < seuil, AppendTo[idx, i];
    i += dir;];
   idx];

(* Identify positions surrounding the global maximum likelihood in clusters *)
getContigAroundMax[clusters_, mod_, nb_] := 
 Module[{maxPos, lst, seuil},
  maxPos = (Position[clusters, mod[[nb]]] // 
      Flatten)[[1]] \
(*position of the highest peak (the one with max likelihood) in the list of clusters*); 
  lst = Table[
    Max[cluster[[All, nb]]], {cluster, 
     clusters}] (*list of maximum likelihood per peak *); 
  seuil = 
   Min[clusters[[maxPos, All, nb]]] (*Minimum of the highest peak*);
  Join[Reverse@getContig[lst, maxPos - 1, seuil, -1], {maxPos}, 
   getContig[lst, maxPos + 1, seuil, 1]] 
  ]

(* Append a 2D list to a CSV file *)
AppendCSV[file_, data_List] := Module[{strm}, strm = OpenAppend[file];
   WriteString[strm, 
    StringRiffle[ToString /@ #, ","] <> "\n" & /@ data // StringJoin];
   Close[strm];];

(* Identify clusters using DBSCAN *)
cluster[data_] := Module[{constant, variable, clust, result},
  constant = Select[data, Length[Union[#]] == 1 &];
  variable = Select[data, Length[Union[#]] > 1 &];
  clust = 
   If[variable == {}, {}, 
    FindClusters[Flatten[variable], Method -> "DBSCAN"]];
  result = SortBy[Join[constant, clust], Min]
  ]

(* Subsample positions to reduce computational load; select candidate positions, cluster them, and compute final list *)
SubSampling[list_, sousech_, model_] :=
 Module[{listsousech, validPairsFunction, tsousech, modsousech, 
   tabSEsousech, candidateposition, clust, listfinal},
  listsousech = list[[;; ;; sousech]];
  validPairsFunction = ToExpression["validPairs" <> ToString[model]];
  tsousech = validPairsFunction[listsousech, minslice, cutSE] // Quiet;
  modsousech = Last[SortBy[tsousech, #[[model + 1]] &]];
  tabSEsousech = 
   Select[tsousech, #[[model + 1]] > modsousech[[model + 1]] - cutSE &];
  candidateposition = tabSEsousech[[All, 1 ;; model]] // T;
  clust = cluster[candidateposition];
  listfinal = ListFinal[clust, sousech];
  {modsousech, tabSEsousech, listfinal}
  ]

(* Export limits for a chosen model *)
ExportLimits[Nbcp_, TabSE_, Mod_] := 
 Flatten[{{{Mod[[Nbcp + 1]]}}, {datatopo[[Mod[[Range[1, Nbcp]]], 2]]},
    Table[datatopo[[{Min[TabSE[[All, i]]], Max[TabSE[[All, i]]]}, 
      2]], {i, 1, Nbcp}]}, 1]

limits[Nbcp_, TabSE_, Mod_] := Flatten[{Mod[[Range[1, Nbcp]]]}, 1]

eps = 0.01
errvect = {0,0,0,0,0,0,0,0,1,1,1,1,0,0} (* we expect some errors on turnover topologies *) 

(* Compute log-likelihood *)
mll[vect_] := (tot = vect // Total; 
  mlv = Table[
    If[vect[[i]] > 0, vect[[i]] Log[(vect[[i]]/tot + errvect[[i]] eps) (1/(1 + 4 eps))], 0], {i, 1, 14}]; 
  mlv // Total // N)

(* ---------------- Data import and initialization ---------------- *)

If[Length[$ScriptCommandLine] < 3, 
  Echo["Erreur : need a dataset name and a directory in argument !"];
  Exit[]];

directory = $ScriptCommandLine[[2]];

SetDirectory[StringJoin["/Users/suzie/Documents/THESE/GENOMIQUE/Analyses/", directory, "/"]]

dataname = $ScriptCommandLine[[3]];

datatopo = 
  SortBy[Import[StringJoin[dataname, "_topo_limits.csv"], "CSV"] // 
    Rest, #[[2]] &];
datatopo[[All, 2]] = ToExpression[datatopo[[All, 2]]];
datatopo = SortBy[datatopo, #[[2]] &];

cc = datatopo[[All, -14 ;; -1]] // ToExpression;

(* Precompute constants for all models *)

maxpos = cc // Length

listbase = Range[minslice, maxpos - minslice];

cutSE = Log[cc // L // N]  (*penalité du BIC*)


Share[]

(* ---------------- models ---------------- *)


(* The following sections fit 0,1,2 change point models and compute BICs *)
(* model0: single segment, global likelihood *)

mod0 = mll[cc // Total];

BICmod0 = -2 mod0 + cutSE 13 ;

Export[StringJoin["BIC_", dataname, 
  ".csv"], {{"model0", BICmod0}}, "CSV"]

(* model1: single change point *)

MaxLik[lim_?IntegerQ] := (cc1 = Table[cc[[i]], {i, 1, lim}] // Total; 
  cc2 = Table[cc[[i]], {i, lim + 1, maxpos}] // Total; 
  mll[cc1] + mll[cc2])

t1 = Table[{lim, MaxLik[lim]}, {lim, minslice, maxpos - minslice}] // 
   Quiet;

mod1 = Last[SortBy[t1, #[[2]] &]];

BICmod1 = -2 mod1[[2]] + cutSE (2*13 + 1)

tabSE1 = Select[t1, #[[2]] > mod1[[2]] - cutSE &];

t1lim = ExportLimits[1, tabSE1, mod1]

(*Export[StringJoin["topo1_", dataname, ".csv"], t1lim, "CSV"];*)
AppendCSV[StringJoin["BIC_", dataname, ".csv"], {{"model1", BICmod1}}];

Clear[t1];

(* model2: two change points *)

MaxLik2[lim1_?IntegerQ, 
  lim2_?IntegerQ] := (cc1 = Table[cc[[i]], {i, 1, lim1}] // Total; 
  cc2 = Table[cc[[i]], {i, lim1 + 1, lim2}] // Total; 
  cc3 = Table[cc[[i]], {i, lim2 + 1, maxpos}] // Total; 
  mll[cc1] + mll[cc2] + mll[cc3])

AbsoluteTiming[
 t2 = Table[{lim1, lim2, MaxLik2[lim1, lim2]}, {lim1, minslice, 
      maxpos - 2 minslice}, {lim2, lim1 + minslice, 
      maxpos - minslice}] // Quiet;]

mod2 = Last[SortBy[Flatten[t2, 1], #[[3]] &]];

BICmod2 = -2 mod2[[3]] + cutSE (3*13 + 2);

tabSE2 = Select[Flatten[t2, 1], #[[3]] > mod2[[3]] - cutSE &];

t2lim = ExportLimits[2, tabSE2, mod2]

(*Export[StringJoin["topo2_", dataname, ".csv"], t2lim, "CSV"];*)
(*Export[StringJoin["topo_tabSE2_", dataname, ".csv"], tabSE2, "CSV"];*)
AppendCSV[StringJoin["BIC_", dataname, ".csv"], {{"model2", BICmod2}}];

Clear[t2]

(*model3: three change points*)

MaxLik3[lim1_?IntegerQ, lim2_?IntegerQ, 
  lim3_?IntegerQ] := (cc1 = Table[cc[[i]], {i, 1, lim1}] // Total; 
  cc2 = Table[cc[[i]], {i, lim1 + 1, lim2}] // Total; 
  cc3 = Table[cc[[i]], {i, lim2 + 1, lim3}] // Total; 
  cc4 = Table[cc[[i]], {i, lim3 + 1, maxpos}] // Total; 
  mll[cc1] + mll[cc2] + mll[cc3] + mll[cc4])

validPairs3[list_, minslice_, cutSE_] := 
 Module[{vals = Sort[list], n = Length[list], res, 
   currentMaxLik = -Infinity}, 
  res = Reap[
     Do[i2Candidates = 
       Select[Range[i1 + 1, n], vals[[#]] - vals[[i1]] >= minslice &];
      Do[
       i3Candidates = 
        Select[Range[i2 + 1, n], vals[[#]] - vals[[i2]] >= minslice &];
       Do[Module[{v1, v2, v3, lik}, {v1, v2, v3} = vals[[{i1, i2, i3}]];
         lik = MaxLik3[v1, v2, v3];
         (*Filtre mémoire*)
         If[lik >= currentMaxLik - cutSE, Sow[{v1, v2, v3, lik}];
          If[lik > currentMaxLik, currentMaxLik = lik;];];], {i3, 
         i3Candidates}], {i2, i2Candidates}], {i1, 1, n - 2}]][[2, 1]];
  res]

(* Subsampling to reduce computational load *)
sousech3 = 5;

AbsoluteTiming[{mod3sousech, tabSE3sousech, list3sousech} = 
   SubSampling[listbase, sousech3, 3];]

AbsoluteTiming[{mod3, tabSE3, list3} = SubSampling[list3sousech, 1, 3];]

BICmod3 = -2 mod3[[4]] + cutSE (4*13 + 3);

t3lim = ExportLimits[3, tabSE3, mod3];

(*Export[StringJoin["topo_tabSE3_", dataname, ".csv"], tabSE3, "CSV"];*)
(*Export[StringJoin["topo3_", dataname, ".csv"], t3lim, "CSV"];*)
AppendCSV[StringJoin["BIC_", dataname, ".csv"], {{"model3", BICmod3}}];

Share[]

(*model4 : four change points *)

MaxLik4[lim1_?IntegerQ, lim2_?IntegerQ, lim3_?IntegerQ, 
  lim4_?IntegerQ] := (cc1 = Table[cc[[i]], {i, 1, lim1}] // Total; 
  cc2 = Table[cc[[i]], {i, lim1 + 1, lim2}] // Total; 
  cc3 = Table[cc[[i]], {i, lim2 + 1, lim3}] // Total; 
  cc4 = Table[cc[[i]], {i, lim3 + 1, lim4}] // Total;
  cc5 = Table[cc[[i]], {i, lim4 + 1, maxpos}] // Total; 
  mll[cc1] + mll[cc2] + mll[cc3] + mll[cc4] + mll[cc5])

validPairs4[list_, minslice_, cutSE_] := 
 Module[{vals = Sort[list], n = Length[list], res, 
   currentMaxLik = -Infinity}, 
  res = Reap[
     Do[i2Candidates = 
       Select[Range[i1 + 1, n], vals[[#]] - vals[[i1]] >= minslice &];
      Do[
       i3Candidates = 
        Select[Range[i2 + 1, n], vals[[#]] - vals[[i2]] >= minslice &];
       Do[
        i4Candidates = 
         Select[Range[i3 + 1, n], vals[[#]] - vals[[i3]] >= minslice &];
        Do[
         Module[{v1, v2, v3, v4, lik}, {v1, v2, v3, v4} = 
           vals[[{i1, i2, i3, i4}]];
          lik = MaxLik4[v1, v2, v3, v4];
          If[lik >= currentMaxLik - cutSE, Sow[{v1, v2, v3, v4, lik}];
           If[lik > currentMaxLik, currentMaxLik = lik;];];], {i4, 
          i4Candidates}], {i3, i3Candidates}], {i2, 
        i2Candidates}], {i1, 1, n - 3}]][[2, 1]];
  res]

(* Subsampling to reduce computational load *)
sousech4 = 5;

AbsoluteTiming[{mod4sousech, tabSE4sousech, list4sousech} = 
   SubSampling[listbase, sousech4, 4];]

AbsoluteTiming[{mod4, tabSE4, list4} = SubSampling[list4sousech, 1, 4];]

BICmod4 = -2 mod4[[5]] + cutSE (5*13 + 4);

t4lim = ExportLimits[4, tabSE4, mod4];

(*Export[StringJoin["topo_tabSE4_", dataname, ".csv"], tabSE4, "CSV"];*)
(*Export[StringJoin["topo4_", dataname, ".csv"], t4lim, "CSV"];*)
AppendCSV[StringJoin["BIC_", dataname, ".csv"], {{"model4", BICmod4}}];

Share[]

(*model5: five change points *)

MaxLik5[lim1_?IntegerQ, lim2_?IntegerQ, lim3_?IntegerQ, 
  lim4_?IntegerQ, 
  lim5_?IntegerQ] := (cc1 = Table[cc[[i]], {i, 1, lim1}] // Total; 
  cc2 = Table[cc[[i]], {i, lim1 + 1, lim2}] // Total; 
  cc3 = Table[cc[[i]], {i, lim2 + 1, lim3}] // Total; 
  cc4 = Table[cc[[i]], {i, lim3 + 1, lim4}] // Total;
  cc5 = Table[cc[[i]], {i, lim4 + 1, lim5}] // Total; 
  cc6 = Table[cc[[i]], {i, lim5 + 1, maxpos}] // Total; 
  mll[cc1] + mll[cc2] + mll[cc3] + mll[cc4] + mll[cc5] + mll[cc6])

validPairs5[list_, minslice_, cutSE_] := 
 Module[{vals = Sort[list], n = Length[list], res, 
   currentMaxLik = -Infinity}, 
  res = Reap[
     Do[i2Candidates = 
       Select[Range[i1 + 1, n], vals[[#]] - vals[[i1]] >= minslice &];
      Do[
       i3Candidates = 
        Select[Range[i2 + 1, n], vals[[#]] - vals[[i2]] >= minslice &];
       Do[
        i4Candidates = 
         Select[Range[i3 + 1, n], vals[[#]] - vals[[i3]] >= minslice &];
        Do[
         i5Candidates = 
          Select[Range[i4 + 1, n], 
           vals[[#]] - vals[[i4]] >= minslice &];
         
         Do[Module[{v1, v2, v3, v4, v5, lik}, {v1, v2, v3, v4, v5} = 
            vals[[{i1, i2, i3, i4, i5}]];
           lik = MaxLik5[v1, v2, v3, v4, v5];
           (*Seulement si la combinaison est assez bonne*)
           If[lik >= currentMaxLik - cutSE, 
            Sow[{v1, v2, v3, v4, v5, lik}];
            If[lik > currentMaxLik, currentMaxLik = lik];];], {i5, 
           i5Candidates}], {i4, i4Candidates}], {i3, 
         i3Candidates}], {i2, i2Candidates}], {i1, 1, n - 4}]][[2, 1]];
  res]

(* Subsampling (several times) to reduce computational load *)
sousech5 = 10;

AbsoluteTiming[{mod5sousech, tabSE5sousech, list5sousech} = 
   SubSampling[listbase, sousech5, 5];]

AbsoluteTiming[{mod5sousech2, tabSE5sousech2, list5sousech2} = 
   SubSampling[list5sousech, 4, 5];]

AbsoluteTiming[{mod5, tabSE5, list5} = 
   SubSampling[list5sousech2, 1, 5];]

BICmod5 = -2 mod5[[6]] + cutSE (6*13 + 5);

t5lim = ExportLimits[5, tabSE5, mod5];

(*Export[StringJoin["topo_tabSE5_", dataname, ".csv"], tabSE5, "CSV"];*)
(*Export[StringJoin["topo5_", dataname, ".csv"], t5lim, "CSV"];*)
AppendCSV[StringJoin["BIC_", dataname, ".csv"], {{"model5", BICmod5}}];

Share[]

(* Determine the chosen model with the lowest BIC and extract clusters around max likelihood to compute more precise likelihood confidence bounds for the chosen model*)

nbcp = First[
    Ordering[{BICmod0, BICmod1, BICmod2, BICmod3, BICmod4, BICmod5}, 1]] - 1;

chosenmodel = nbcp // ToString;

tabSE = StringJoin["tabSE", chosenmodel] // ToExpression;

mod = Last[SortBy[tabSE, #[[nbcp + 1]] &]];

clusters = 
  SortBy[FindClusters[tabSE[[All, Range[1, nbcp]]], Method -> "DBSCAN"],
    Min];
clustersFull = 
  Table[Select[tabSE, 
    MemberQ[cluster, #[[Range[1, nbcp]]]] &], {cluster, clusters}];

(* Identify positions around maximum likelihood *)
pic = getContigAroundMax[clustersFull, mod, nbcp + 1];

limcutted = ExportLimits[nbcp, Flatten[clustersFull[[pic]], 1], mod];
lims = limits[nbcp, Flatten[clustersFull[[pic]], 1], mod];

Export[StringJoin["topo", chosenmodel, "_cutted_", dataname, 
  ".csv"], limcutted, "CSV"]

(* If 5 change points, export each segment separately and create a master file listing segments on which the analysis will be performed recursively *)

If[nbcp == 5,
  (*Build the 6 segments from lims*)
  boundaries = Join[{1}, lims, {maxpos + 1}];
  segments = Table[datatopo[[boundaries[[i]] ;; boundaries[[i + 1]] - 1]], {i, 1, Length[boundaries] - 1}];
  subsetNames = {};
(*pour mémoriser seulement les dataname des subsets*)
(*Export each segment and record its subset dataname*)
  Do[seg = segments[[i]];
   lim1 = seg[[1, 2]] // ToString;
   lim2 = seg[[-1, 2]] // ToString;
   subsetName = StringJoin[dataname, "_", lim1, "_", lim2];
   fname = StringJoin[subsetName, "_topo_limits.csv"];
   Export[fname, seg, "CSV"];
   AppendTo[subsetNames, subsetName];, {i, 1, Length[segments]}];
  (*Export master file listing subset datanames ONLY*)
  masterFile = StringJoin[dataname, "_subsets.csv"];
  Export[masterFile, subsetNames, "CSV"];];