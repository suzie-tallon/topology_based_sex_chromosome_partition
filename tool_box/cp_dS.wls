#!/usr/bin/env wolframscript

(* 
Author : Suzie Tallon & Thomas Lenormand 
Date : February 2026
This script is used to detect change points in dS data using likelihood-based models.
It identifies clusters of positions with high likelihood, estimates BIC for 0,1,2 change point models, and exports both model summaries and position-specific segments for downstream analyses.
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
 Flatten[{{{Mod[[Nbcp + 1]]}}, {data[[Mod[[Range[Nbcp]]], 1]]}, 
   Table[data[[{Min[TabSE[[All, i]]], Max[TabSE[[All, i]]]}, 1]], {i, 
     1, Nbcp}]}, 1]

limits[Nbcp_, TabSE_, Mod_] := Flatten[{Mod[[Range[1, Nbcp]]]}, 1]

(* Compute log-likelihood *)
ll[data_] := FindMaximum[Sum[If[Subscript[y, i] == 0, Subscript[w, i]*Log[p], Subscript[w, i]*(-(Subscript[y, i]/b) - a*Log[b] + Log[1 - p] + (-1 + a)*Log[Subscript[y, i]] - Log[Gamma[a]])] /. 
       Subscript[y, i_] :> data[[i,2]] /. Subscript[w, i_] :> data[[i,3]], {i, 1, Length[data]}], {{a, 1, 0.0001, 4}, {b, 0.02, 0.0001, 4}, {p, 0.1, 0.00001, 0.9999}}][[1]]


(* ---------------- Data import and initialization ---------------- *)


If[Length[$ScriptCommandLine] < 3, 
  Echo["Erreur : need a dataset name and a directory in argument !"];
  Exit[]];

directory = $ScriptCommandLine[[2]];

SetDirectory[StringJoin["/Users/suzie/Documents/THESE/GENOMIQUE/Analyses/", directory, "/"]]

dataname = $ScriptCommandLine[[3]];

data = Import[StringJoin["dS_", dataname, ".csv"], "CSV"] // Rest // 
   ToExpression;

weights = data[[All, 3]]/(data[[All, 3]] // Mean) // N;

dat = Table[{(data[[i, 1]] - Min[data[[All, 1]]])/(
     Max[data[[All, 1]]] - Min[data[[All, 1]]]) // N, data[[i, 2]], 
    weights[[i]]}, {i, 1, Length[data]}];

(* Precompute constants for all models *)

maxpos = dat // Length

listbase = Range[minslice, maxpos - minslice];

cutSE = Log[dat // L // N]  (* BIC penalty *)

Share[]

(* ---------------- models ---------------- *)


(* The following sections fit 0,1,2 change point models and compute BICs *)
(* model0: single segment, global likelihood *)

mod0 = ll[dat];

BICmod0 = -2 mod0 + cutSE 3 ;

Export[StringJoin["BIC_", dataname, 
  ".csv"], {{"model0", BICmod0}}, "CSV"]

(* model1: single change point *)

MaxLik[lim_?IntegerQ] := (cc1 = Table[dat[[i]], {i, 1, lim}]; 
  cc2 = Table[dat[[i]], {i, lim + 1, maxpos}]; ll[cc1] + ll[cc2])

t1 = Table[{lim, MaxLik[lim]}, {lim, minslice, maxpos - minslice}] // 
   Quiet;

mod1 = Last[SortBy[t1, #[[2]] &]];

BICmod1 = -2 mod1[[2]] + cutSE (2*3 + 1)

tabSE1 = Select[t1, #[[2]] > mod1[[2]] - cutSE &];

t1lim = ExportLimits[1, tabSE1, mod1]

(*Export[StringJoin["topo1_", dataname, ".csv"], t1lim, "CSV"];*)
AppendCSV[StringJoin["BIC_", dataname, ".csv"], {{"model1", BICmod1}}];

Clear[t1];

(* model2: two change points *)

MaxLik2[lim1_?IntegerQ, 
  lim2_?IntegerQ] := (cc1 = Table[dat[[i]], {i, 1, lim1}]; 
  cc2 = Table[dat[[i]], {i, lim1 + 1, lim2}]; 
  cc3 = Table[dat[[i]], {i, lim2 + 1, maxpos}]; 
  ll[cc1] + ll[cc2] + ll[cc3])

validPairs2[list_, minslice_, cutSE_] := 
 Module[{vals = Sort[list], n = Length[list], res, 
   currentMaxLik = -Infinity}, 
  res = Reap[
     Do[i2Candidates = 
       Select[Range[i1 + 1, n], vals[[#]] - vals[[i1]] >= minslice &];
      Do[Module[{v1, v2, lik}, {v1, v2} = vals[[{i1, i2}]];
        lik = MaxLik2[v1, v2];
        (*filtre mémoire*)
        If[lik >= currentMaxLik - cutSE, Sow[{v1, v2, lik}];
         If[lik > currentMaxLik, currentMaxLik = lik];];], {i2, 
        i2Candidates}], {i1, 1, n - 1}]][[2, 1]];
  res]

(* Subsampling to reduce computational load for 2-change-point model *)
sousech2 = 2;

AbsoluteTiming[{mod2sousech, tabSE2sousech, list2sousech} = 
   SubSampling[listbase, sousech2, 2];]

AbsoluteTiming[{mod2, tabSE2, list2} = SubSampling[list2sousech, 1, 2];]


BICmod2 = -2 mod2[[3]] + cutSE (3*3 + 2);


t2lim = ExportLimits[2, tabSE2, mod2]

(*Export[StringJoin["topo2_", dataname, ".csv"], t2lim, "CSV"];*)
(*Export[StringJoin["topo_tabSE2_", dataname, ".csv"], tabSE2, "CSV"];*)
AppendCSV[StringJoin["BIC_", dataname, ".csv"], {{"model2", BICmod2}}];

Share[]

(* Determine the chosen model with the lowest BIC and extract clusters around max likelihood to compute more precise likelihood confidence bounds for the chosen model*)

nbcp = First[
    Ordering[{BICmod0, BICmod1, BICmod2}, 1]] - 1;

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

(* If 2 change points, export each segment separately and create a master file listing segments on which the analysis will be performed recursively *)

If[nbcp == 2,
  (*Build the 6 segments from lims*)
  boundaries = Join[{1}, lims, {maxpos + 1}];
  segments = Table[data[[boundaries[[i]] ;; boundaries[[i + 1]] - 1]], {i, 1, Length[boundaries] - 1}];
  subsetNames = {};
(*pour mémoriser seulement les dataname des subsets*)
(*Export each segment and record its subset dataname*)
  Do[seg = segments[[i]];
   lim1 = seg[[1, 1]] // ToString;
   lim2 = seg[[-1, 1]] // ToString;
   subsetName = StringJoin[dataname, "_", lim1, "_", lim2];
   fname = StringJoin["dS_",subsetName, ".csv"];
   Export[fname, seg, "CSV"];
   AppendTo[subsetNames, subsetName];, {i, 1, Length[segments]}];
  (*Export master file listing subset datanames ONLY*)
  masterFile = StringJoin[dataname, "_subsets.csv"];
  Export[masterFile, subsetNames, "CSV"];];