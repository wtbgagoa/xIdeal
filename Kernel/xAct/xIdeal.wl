(* ::Package:: *)

(*
    xIdeal
    Identification of exact solutions in xAct
    Alfonso Garc\[IAcute]a-Parrado
    agparrado@uco.es
    Universidad de C\[OAcute]rdoba, Spain

    Salvador Mengual Sendra
    Salvador.Mengual@uv.es
    Universidad  de Valencia, Spain

    (c) 2023, under GPL

    http://www.xAct.es/
    http://groups.google.com/group/xAct
    https://github.com/xAct-contrib
    xIdeal is a package for exploiting ideal characterizations in General Relativity.

    xIdeal is distributed under the GNU General Public License, and runs on top of xTensor and xCoba which are free packages for fast
    manipulation of abstract and component tensor expressions. xTensor and xCoba can be downloaded from http://www.xact.es/
*)



(* ::Input::Initialization:: *)
xAct`xIdeal`$xTensorVersionExpected = {"1.1.2", {2015, 8, 23}};

xAct`xIdeal`$Version = {"0.0.1", {2023, 10, 3}};

(******************************************************************************)

(********************* 1. Initialization **************************************)

(******************************************************************************)

(************************ 1.1 GPL *********************************************)

(* xTerior: Identification of exact solutions in xAct *)

(* Copyright (C) 2023 Alfonso Garcia-Parrado Gomez-Lobo and Salvador Mengual Sendra*)

(* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License,or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place-Suite 330, Boston, MA 02111-1307, USA. 
*)

(*********************** 1.2 Info Package ************************************)

(* :Title: x *)

(* :Author: Alfonso Garcia-Parrado Gomez-Lobo and Salvador Mengual Sendra *)

(* :Summary: Identification of exact solutions in xAct *)

(* :Brief Discussion:
   - xIdeal
   
*)

(* :Context: xAct`xIdeal` *)

(* :Package Version: 0.0.1 *)

(* :Copyright: Alfonso Garcia-Parrado Gomez-Lobo and Salvador Mengual Sendra (2023) *)

(* :History: See git log *)

(* :Keywords: *)

(* :Source: xIdeal.wl *)

(* :Warning: *)

(* :Mathematica Version: 9.0 and later *)

(* :Limitations: - ?? *)



(* ::Section:: *)
(* BeginPackage *)


With[
	{
		xAct`xIdeal`Private`xIdealSymbols = 
			DeleteCases[
				Join[Names["xAct`xIdeal`*"], Names["xAct`xIdeal`Private`*"]], 
				"$Version" | "xAct`xIdeal`$Version" | "$xTensorVersionExpected" | 
				"xAct`xIdeal`$xTensorVersionExpected"
			]
	},
	Unprotect /@ xAct`xIdeal`Private`xIdealSymbols;
	Clear /@ xAct`xIdeal`Private`xIdealSymbols;
]

If[Unevaluated[xAct`xCore`Private`$LastPackage] === xAct`xCore`Private`$LastPackage,
	
	xAct`xCore`Private`$LastPackage = "xAct`xIdeal`"
];

(* Explicit (not hidden) import of xTensor, xPerm and xCore: *)

BeginPackage["xAct`xIdeal`", {"xAct`xCoba`", "xAct`xTensor`", "xAct`xPerm`",
	 "xAct`xCore`"}]

If[Not @ OrderedQ @ Map[Last, {xAct`xIdeal`$xTensorVersionExpected, xAct`xTensor`$Version
	}],
	Throw @ Message[General::versions, "xTensor", xAct`xTensor`$Version,
		 xAct`xIdeal`$xTensorVersionExpected]
]

Print[xAct`xCore`Private`bars]

Print["Package xAct`xIdeal`  version ", xAct`xIdeal`$Version[[1]], ", ",
	 xAct`xIdeal`$Version[[2]]];

Print["Copyright (C) 2023, Alfonso Garcia-Parrado Gomez-Lobo and Salvador Mengual Sendra, under the General Public License."
	];

Off[General::shdw]

xAct`xIdeal`Disclaimer[] :=
	Print[                                                   "These are points 11 and 12 of the General Public License:\n\n
BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. 
EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM `AS IS\.b4 
WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. 
SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.\n\n
IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO 
MAY MODIFY AND/OR REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, 
INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF 
DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM 
TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES."
		]

On[General::shdw]

(* If xIdeal is not being called from other package then write this GPL short disclaimer: *)

If[xAct`xCore`Private`$LastPackage === "xAct`xIdeal`",
	Unset[xAct`xCore`Private`$LastPackage];
	Print[xAct`xCore`Private`bars];
	Print["These packages come with ABSOLUTELY NO WARRANTY; for details type Disclaimer[]. This is free software, and you are welcome to redistribute it under certain conditions. See the General Public License for details."
		];
	Print[xAct`xCore`Private`bars]
]

(************************* 1.4. Non-standard setup ***********************************)

(* Screen all dollar indices: *)

$PrePrint = ScreenDollarIndices;



(* ::Section:: *)
(* Usage information *)


PetrovType::usage = " ";

DebeverNullDirections::usage = " ";

TypeDClassify::usage = " ";

PSimplify::usage = " ";

IsometryGroupDimension::usage = " ";

KerrSolutionQ::usage = " ";


(* ::Section:: *)
(* Messages *)


PetrovType::nometric = "Metric `1` has not been registered as a metric";



(* ::Section:: *)
(* BeginPrivate *)


Begin["`Private`"]

(* ::Section:: *)
(* Computation of the metric concomitants *)

metricConcomitant["G2Form"][metric_CTensor, opts___] :=
(metricConcomitant["G2Form"][metric, opts] = 
	Module[{simplf, cart, a1, b1, c1, d1, e1, f1, epsilonmetric},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = 	Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];
		epsilonmetric = epsilon[metric];
		simplf[
			HeadOfTensor[
				1/2 (-I epsilonmetric[-a1, -b1, -c1, -d1] + 
				metric[-a1, -c1] metric[-b1, -d1] - metric[-a1, -d1] metric[-b1, -c1]), 
				{-a1, -b1, -c1, -d1}
			]
		]
	]
)

(* ::Section:: *)
(* Computation of the Weyl concomitants *)


weylConcomitant["Weyl"][metric_CTensor, opts___] :=
(weylConcomitant["Weyl"][metric, opts] = 
	Module[{cart, cd},
		cart = 	Part[metric, 2, 1, -1];
		cd = CovDOfMetric[metric];
		MetricCompute[metric, cart, "Weyl"[-1, -1, -1, -1], Parallelize -> True, Verbose -> True];
		Weyl[cd]
	]
)


weylConcomitant["WeylDual"][metric_CTensor, opts___] :=
(weylConcomitant["WeylDual"][metric, opts] = 
	Module[{simplf, cart, a1, b1, c1, d1, e1, f1, cd, weylcd, epsilonmetric, weyldual},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];
		cd = CovDOfMetric[metric];
		epsilonmetric = epsilon[metric];
		weylcd = weylConcomitant["Weyl"][metric, opts];
		weyldual = simplf[HeadOfTensor[1/2 epsilonmetric[-c1, -d1, -e1, -f1] weylcd[e1, f1, -a1, -b1], {-c1, -d1, -a1, -b1}]]
	]
)


weylConcomitant["WeylSelfDual"][metric_CTensor, opts___] :=
(weylConcomitant["WeylSelfDual"][metric, opts] = 
	Module[{simplf},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		simplf[1/2 (weylConcomitant["Weyl"][metric, opts] - I * weylConcomitant["WeylDual"][metric, opts])]
	]
)

weylConcomitant["WeylSelfDual2"][metric_CTensor, opts___] :=
(weylConcomitant["WeylSelfDual2"][metric, opts] = 
	Module[{simplf, weylselfdual, cart, a1, b1, c1, d1, e1, f1},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		simplf[HeadOfTensor[1/2 weylselfdual[-a1, -b1, e1, f1] weylselfdual[-e1, -f1, -c1, -d1], {-a1, -b1, -c1, -d1}]]
	]
)

weylConcomitant["TraceWeylSelfDual2"][metric_CTensor, opts___] :=
(weylConcomitant["TraceWeylSelfDual2"][metric, opts] = 
	Module[{weylselfdual2, simplf, cart, a1, b1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		weylselfdual2 = weylConcomitant["WeylSelfDual2"][metric, opts];
		simplf[weylselfdual2[-a1, -b1, a1, b1] / 2]
	]
)

weylConcomitant["WeylSelfDual3"][metric_CTensor, opts___] :=
(weylConcomitant["WeylSelfDual3"][metric, opts] = 
	Module[{simplf, weylselfdual, weylselfdual2, cart, a1, b1, c1, d1, e1, f1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		weylselfdual2 = weylConcomitant["WeylSelfDual2"][metric, opts];
		simplf[HeadOfTensor[1/2 weylselfdual2[-a1, -b1, e1, f1] weylselfdual[-e1, -f1, -c1, -d1], {-a1, -b1, -c1, -d1}]]
	]
)

weylConcomitant["TraceWeylSelfDual3"][metric_CTensor, opts___] :=
(weylConcomitant["TraceWeylSelfDual3"][metric, opts] = 
	Module[{weylselfdual3, simplf, cart, a1, b1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		weylselfdual3 = weylConcomitant["WeylSelfDual3"][metric, opts];
		simplf[weylselfdual3[-a1, -b1, a1, b1] / 2]
	]
)

weylConcomitant["ScalarW"][metric_CTensor, opts___] :=
(weylConcomitant["ScalarW"][metric, opts] = 
	Module[{simplf, cart, weylselfdual, g2form, aa, bb, w, cd, a1, b1, c1, d1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{b1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		g2form = metricConcomitant["G2Form"][metric, opts];
		aa = 2 weylConcomitant["TraceWeylSelfDual2"][metric, opts];
		bb = 2 weylConcomitant["TraceWeylSelfDual3"][metric, opts];
		w = simplf[-bb / (2 aa)]
	]
)

weylConcomitant["ScalarZ"][metric_CTensor, opts___] :=
(weylConcomitant["ScalarZ"][metric, opts] = 
	Module[{simplf, cart, weylselfdual, g2form, w, cd, a1, b1, c1, d1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{b1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];		
		cd = CovDOfMetric[metric];
		w = weylConcomitant["ScalarW"][metric, opts];
		simplf[(metric[b1, d1]) cd[-b1][w] cd[-d1][w]]
	]
)

weylConcomitant["TensorXi"][metric_CTensor, opts___] :=
(weylConcomitant["TensorXi"][metric, opts] = 
	Module[{simplf, cart, weylselfdual, g2form, w, cd, a1, b1, c1, d1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 4];
		weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		g2form = metricConcomitant["G2Form"][metric, opts];
		w = weylConcomitant["ScalarW"][metric, opts];
		cd = CovDOfMetric[metric];
		simplf[(weylselfdual[-a1, b1, -c1, d1] - w g2form[-a1, b1, -c1, d1]) cd[-b1][w] cd[-d1][w]]
	]
)


(* ::Section:: *)
(* Computation of the R-frame concomitants *)

(*
NOTE: I put H as an input in the functions that compute its derivative concomitants
instead of calling RframeConcomitant["ConnectionTensor"] in case it is not computed from the R-frame.
*)

(* TODO: Add different options to compute H to RframeConcomitant["ConnectionTensor"] *)

RframeConcomitant["ConnectionTensor"][metric_CTensor, e0_CTensor, e1_CTensor, e2_CTensor, e3_CTensor, opts___] :=
(RframeConcomitant["ConnectionTensor"][metric, e0, e1, e2, e3, opts] = 
	Module[{simplf, cart, cd, a1, b1, c1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 3];
		cd = CovDOfMetric[metric];
		simplf[HeadOfTensor[-(1/2) (-(cd[a1][e0[b1]] e0[c1] - cd[a1][e0[c1]] e0[b1]) + (cd[a1][e1[b1]] e1[c1] - cd[a1][e1[c1]] e1[b1]) + (cd[a1][e2[b1]] e2[c1] - 
           			cd[a1][e2[c1]] e2[b1]) + (cd[a1][e3[b1]] e3[c1] - cd[a1][e3[c1]] e3[b1])), {-a1, -b1, -c1}]]
	]
 )

RframeConcomitant["C1"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C1"][metric, H, opts] = 
	Module[{simplf, cart, cd, a1, b1, c1, d1, i1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, i1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 5];
		cd = CovDOfMetric[metric];
		simplf[HeadOfTensor[cd[-a1][H[-b1, -c1, -d1]] + H[-a1, -b1, i1] H[-i1, -c1, -d1] + H[-a1, -c1, i1] H[-b1, -i1, -d1] + 
				H[-a1, -d1, i1] H[-b1, -c1, -i1], {-a1, -b1, -c1, -d1}]]
  	]
)

RframeConcomitant["C11"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C11"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 10];
  		epsilonmetric = epsilon[metric];
    		ce1 = RframeConcomitant["C1"][metric, H, opts];
		simplf[HeadOfTensor[epsilonmetric[i1, j1, a1, b1] C1[-i1, -c1, -d1, -e1] C1[-j1, -f1, -g1, -h1], {a1, b1, -c1, -d1, -e1, -f1, -g1, -h1}]]
  	]
)

RframeConcomitant["C2"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C2"][metric, H, opts] = 
	Module[{simplf, cart, cd, ce1, a1, b1, c1, d1, e1, i1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, i1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];
  		cd = CovDOfMetric[metric];
    		ce1 = RframeConcomitant["C1"][metric, H, opts];
		simplf[HeadOfTensor[cd[-a1][ce1[-b1, -c1, -d1, -e1]] + H[-a1, -b1, i1] ce1[-i1, -c1, -d1, -e1] +
	   		H[-a1, -c1, i1] ce1[-b1, -i1, -d1, -e1] + H[-a1, -d1, i1] ce1[-b1, -c1, -i1, -e1] + 
         		H[-a1, -e1, i1] ce1[-b1, -c1, -d1, -i1], {-a1, -b1, -c1, -d1, -e1}]]
  	]
)

RframeConcomitant["C12"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C12"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 11];
  		epsilonmetric = epsilon[metric];
    		ce1 = RframeConcomitant["C1"][metric, H, opts];
      		ce2 = RframeConcomitant["C2"][metric, H, opts];
		simplf[HeadOfTensor[epsilonmetric[i1, j1, a1, b1] ce1[-i1, -c1, -d1, -e1] ce2[-j1, -f1, -g1, -h1, -k1], 
	   		{a1, b1, -c1, -d1, -e1, -f1, -g1, -h1, -k1}]]
  	]
)

RframeConcomitant["C122"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C122"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 15];
  		epsilonmetric = epsilon[metric];
    		ce1 = RframeConcomitant["C1"][metric, H, opts];
      		ce2 = RframeConcomitant["C2"][metric, H, opts];
		simplf[HeadOfTensor[epsilonmetric[i1, j1, k1, a1] ce1[-i1, -b1, -c1, -d1] ce2[-j1, -e1, -f1, -g1, -h1] 
	   		ce2[-k1, -l1, -m1, -n1, -o1], {a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -l1, -m1, -n1, -o1}]]
  	]
)

RframeConcomitant["C3"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C3"][metric, H, opts] = 
	Module[{simplf, cart, cd, ce2, a1, b1, c1, d1, e1, f1, i1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, i1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 7];
  		cd = CovDOfMetric[metric];
    		ce2 = RframeConcomitant["C2"][metric, H, opts];
		simplf[HeadOfTensor[cd[-a1][ce2[-b1, -c1, -d1, -e1, -f1]] + H[-a1, -b1, i1] ce2[-i1, -c1, -d1, -e1, -f1] + 
              		H[-a1, -c1, i1] ce2[-b1, -i1, -d1, -e1, -f1] + H[-a1, -d1, i1] ce2[-b1, -c1, -i1, -e1, -f1] + 
          		H[-a1, -e1, i1] ce2[-b1, -c1, -d1, -i1, -f1] + H[-a1, -f1, i1] ce2[-b1, -c1, -d1, -e1, -i1], 
		  	{-a1, -b1, -c1, -d1, -e1, -f1}]]
  	]
)

RframeConcomitant["C123"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C123"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, ce3, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 16];
  		epsilonmetric = epsilon[metric];
    		ce1 = RframeConcomitant["C1"][metric, H, opts];
     		ce2 = RframeConcomitant["C2"][metric, H, opts];
    		ce3 = RframeConcomitant["C3"][metric, H, opts];
		simplf[HeadOfTensor[epsilonmetric[i1, j1, k1, a1] ce1[-i1, -b1, -c1, -d1] ce2[-j1, -e1, -f1, -g1, -h1] 
	      		ce3[-k1, -l1, -m1, -n1, -o1, -p1], {a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -l1, -m1, -n1, -o1, -p1}]]
  	]
)

RframeConcomitant["C1233"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C1233"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, ce3, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, 
 		q1, r1, s1, t1, u1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1, s1, t1, u1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 21];
  		epsilonmetric = epsilon[metric];
    		ce1 = RframeConcomitant["C1"][metric, H, opts];
     		ce2 = RframeConcomitant["C2"][metric, H, opts];
    		ce3 = RframeConcomitant["C3"][metric, H, opts];
		simplf[HeadOfTensor[epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce2[-j1, -d1, -e1, -f1, -g1] 
	      		ce3[-k1, -h1, -m1, -n1, -o1, -p1] ce3[-l1, -q1, -r1, -s1, -t1, -u1], 
	      		{-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1, -q1, -r1, -s1, -t1, -u1}]]
  	]
)

RframeConcomitant["C4"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C4"][metric, H, opts] = 
	Module[{simplf, cart, cd, ce3, a1, b1, c1, d1, e1, f1, g1, i1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, i1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 8];
  		cd = CovDOfMetric[metric];
    		ce3 = RframeConcomitant["C3"][metric, H, opts];
		simplf[HeadOfTensor[cd[-a1][ce3[-b1, -c1, -d1, -e1, -f1, -g1]] + H[-a1, -b1, i1] ce3[-i1, -c1, -d1, -e1, -f1, -g1] + 
          		H[-a1, -c1, i1] ce3[-b1, -i1, -d1, -e1, -f1, -g1] + H[-a1, -d1, i1] ce3[-b1, -c1, -i1, -e1, -f1, -g1] + 
           		H[-a1, -e1, i1] ce3[-b1, -c1, -d1, -i1, -f1, -g1] + H[-a1, -f1, i1] ce3[-b1, -c1, -d1, -e1, -i1, -g1] + 
          		H[-a1, -g1, i1] ce3[-b1, -c1, -d1, -e1, -f1, -i1], {-a1, -b1, -c1, -d1, -e1, -f1, -g1}]]
  	]
)

RframeConcomitant["C1234"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C1234"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, ce3, ce4, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, 
 		q1, r1, s1, t1, u1, v1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1, s1, t1, u1, v1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 22];
  		epsilonmetric = epsilon[metric];
    		ce1 = RframeConcomitant["C1"][metric, H, opts];
     		ce2 = RframeConcomitant["C2"][metric, H, opts];
    		ce3 = RframeConcomitant["C3"][metric, H, opts];
      		ce4 = RframeConcomitant["C4"][metric, H, opts];
		simplf[HeadOfTensor[epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce2[-j1, -d1, -e1, -f1, -g1] 
	       		ce3[-k1, -h1, -m1, -n1, -o1, -p1] ce4[-l1, -q1, -r1, -s1, -t1, -u1, -v1], 
			{-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1, -q1, -r1, -s1, -t1, -u1, -v1}]]
  	]
)

RframeConcomitant["C1222"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C1222"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, 
 		q1, r1, s1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1, s1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 19];
  		epsilonmetric = epsilon[metric];
    		ce1 = RframeConcomitant["C1"][metric, H, opts];
     		ce2 = RframeConcomitant["C2"][metric, H, opts];
		simplf[HeadOfTensor[epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce2[-j1, -d1, -e1, -f1, -g1] 
	   		ce2[-k1, -h1, -m1, -n1, -o1] ce2[-l1, -p1, -q1, -r1, -s1], 
	  		{-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1, -q1, -r1, -s1}]]
  	]
)

RframeConcomitant["C1223"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C1223"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, ce3, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, 
 		q1, r1, s1, t1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1, s1, t1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 20];
  		epsilonmetric = epsilon[metric];
    		ce1 = RframeConcomitant["C1"][metric, H, opts];
     		ce2 = RframeConcomitant["C2"][metric, H, opts];
       		ce3 = RframeConcomitant["C3"][metric, H, opts];
		simplf[HeadOfTensor[epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce2[-j1, -d1, -e1, -f1, -g1] 
	   		ce2[-k1, -h1, -m1, -n1, -o1] ce3[-l1, -p1, -q1, -r1, -s1, -t1], 
	  		{-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1, -q1, -r1, -s1, -t1}]]
  	]
)

RframeConcomitant["C111"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C111"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, m1, n1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, m1, n1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 13];
  		epsilonmetric = epsilon[metric];
    		ce1 = RframeConcomitant["C1"][metric, H, opts];
		simplf[HeadOfTensor[epsilonmetric[i1, j1, k1, a1] ce1[-i1, -b1, -c1, -d1] ce1[-j1, -e1, -f1, -g1] 
	   		ce1[-k1, -h1, -m1, -n1], {a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -l1, -m1, -n1}]]
  	]
)

RframeConcomitant["C112"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C112"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, m1, n1, o1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, m1, n1, o1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 14];
  		epsilonmetric = epsilon[metric];
    		ce1 = RframeConcomitant["C1"][metric, H, opts];
      		ce2 = RframeConcomitant["C2"][metric, H, opts];
		simplf[HeadOfTensor[epsilonmetric[i1, j1, k1, a1] ce1[-i1, -b1, -c1, -d1] ce1[-j1, -e1, -f1, -g1] 
	   		ce2[-k1, -h1, -m1, -n1, -o1], {a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -l1, -m1, -n1}]]
  	]
)

RframeConcomitant["C1122"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C1122"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, 
 		q1, r1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 18];
  		epsilonmetric = epsilon[metric];
    		ce1 = RframeConcomitant["C1"][metric, H, opts];
     		ce2 = RframeConcomitant["C2"][metric, H, opts];
		simplf[HeadOfTensor[epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce1[-j1, -d1, -e1, -f1] 
	   		ce2[-k1, -g1, -h1, -m1, -n1,] ce2[-l1, -o1, -p1, -q1, -r1], 
	  		{-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1, -q1, -r1}]]
  	]
)

RframeConcomitant["C1123"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C1123"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, ce3, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, 
 		q1, r1, s1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1,s1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 19];
  		epsilonmetric = epsilon[metric];
    		ce1 = RframeConcomitant["C1"][metric, H, opts];
     		ce2 = RframeConcomitant["C2"][metric, H, opts];
       		ce3 = RframeConcomitant["C3"][metric, H, opts];
		simplf[HeadOfTensor[epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce1[-j1, -d1, -e1, -f1] 
	   		ce2[-k1, -g1, -h1, -m1, -n1,] ce3[-l1, -o1, -p1, -q1, -r1, -s1], 
	  		{-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1, -q1, -r1, -s1}]]
  	]
)

RframeConcomitant["C1111"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C1111"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 16];
  		epsilonmetric = epsilon[metric];
    		ce1 = RframeConcomitant["C1"][metric, H, opts];
		simplf[HeadOfTensor[epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce1[-j1, -d1, -e1, -f1] ce1[-k1, -g1, -h1, -m1] 
  			ce1[-l1, -n1, -o1, -p1], {-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1}]]
  	]
)

RframeConcomitant["C1112"][metric_CTensor, H_CTensor, opts___] :=
(RframeConcomitant["C1112"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 17];
  		epsilonmetric = epsilon[metric];
    		ce1 = RframeConcomitant["C1"][metric, H, opts];
      		ce2 = RframeConcomitant["C2"][metric, H, opts];
		simplf[HeadOfTensor[epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce1[-j1, -d1, -e1, -f1] ce1[-k1, -g1, -h1, -m1] 
  			ce2[-l1, -n1, -o1, -p1, -q1], {-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1, -q1}]]
  	]
)


(* ::Section:: *)
(* Computation of the Petrov types *)


(*
TODO: there are different algorithms for doing this computation. Add Method option
to be able to choose between them. Add names for each method option.
*)
Options[PetrovType] = {Method -> "Default", PSimplify -> $CVSimplify}
PetrovType[metric_CTensor, opts : OptionsPattern[]] :=
	Catch @
		Module[{cart, cd, weylcd, epsilonmetric, weyldual, weylselfdual, g2form, weylselfdual2, weylselfdual3, aa, bb,
			 rho, simplf},
			If[Not @ MetricQ @ metric,
				Throw[Message[PetrovType::nometric, metric]]
			];
			simplf = OptionValue[PSimplify];
			cart = Part[metric, 2, 1, -1];
			epsilonmetric = epsilon[metric];
			weylcd = weylConcomitant["Weyl"][metric, opts];
			weyldual = weylConcomitant["WeylDual"][metric, opts];
			weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
			g2form = metricConcomitant["G2Form"][metric, opts];
			weylselfdual2 = weylConcomitant["WeylSelfDual2"][metric, opts];
			aa = weylConcomitant["TraceWeylSelfDual2"][metric, opts];
			weylselfdual3 = weylConcomitant["WeylSelfDual3"][metric, opts];
			bb = weylConcomitant["TraceWeylSelfDual3"][metric, opts];
			rho = simplf[- bb / aa];
			Which[
				WeylSelfDual2 === Zero,
					Print["Type N"]
				,
				WeylSelfDual3 === Zero,
					Print["Type III"]
				,
				simplf[aa weylselfdual2 - aa^2 / 3 g2form - bb weylselfdual] === Zero,
					Print["Type D"]
				,
				6 bb^2 - aa^3 === 0,
					Print["Type II"]
				,
				True,
					Print["Type I"]
			]
		]


(* ::Section:: *)
(* Computation of Deveber null directions for each Petrov type *)


(* TODO: we should have private or public functions for the different concommitants *)
(*
TODO: there are different algorithms for doing this computation. Add Method option
to be able to choose between them. Add names for each method option.
*)

DebeverNullDirections[metric_CTensor, u_CTensor, w_CTensor] :=
	Module[{cart, a, b, c, d, e, f, i, j, CD, WeylCD, RiemannCD, RicciCD,
		 RicciScalarCD, epsilonmetric, WeylDual, WeylSelfDual, Q, gamma, Q2, 
		aa, Q3, bb, rho, P, Pdag, scrP, scrP2, dseda, S, Ch2, v0, v1},
		If[Not @ MetricQ @ metric,
			Throw[Message[PetrovType::nometric, metric]]
		];
		cart = Part[metric, 2, 1, -1];
		{a, b, c, d, e, f, i, j} = GetIndicesOfVBundle[Tangent @ ManifoldOfChart
			 @ cart, 8];
		MetricCompute[metric, cart, All, Parallelize -> True, Verbose -> True
			];
		CD = CovDOfMetric[metric];
		WeylCD = Weyl[CD];
		RiemannCD = Riemann[CD];
		RicciCD = Ricci[CD];
		RicciScalarCD = RicciScalar[CD];
		epsilonmetric = epsilon[metric];
		WeylDual = Simplify[HeadOfTensor[1/2 epsilonmetric[-c, -d, -e, -f] 
			WeylCD[e, f, -a, -b], {-c, -d, -a, -b}]];
		WeylSelfDual = Simplify[1/2 (WeylCD - I * WeylDual)];
		Q = HeadOfTensor[2 u[a] u[c] WeylSelfDual[-a, -b, -c, -d], {-b, -d}
			] // FullSimplify;
		gamma = HeadOfTensor[metric[-a, -b] + u[-a] u[-b], {-a, -b}] // FullSimplify
			;
		Q2 = HeadOfTensor[Q[-a, -b] Q[b, -c], {-a, -c}] // FullSimplify;
		aa = Q2[-a, a] // FullSimplify;
		Q3 = HeadOfTensor[Q2[-a, -c] Q[c, -d], {-a, -d}] // Simplify;
		bb = -Q3[-a, a] // FullSimplify;
		Which[
			Q === Zero,
				Print["Type O"]
			,
			Q2 === Zero,
				Print["Type N"];
				HeadOfTensor[Dagger[Q[a, b]] (Q[-b, -a] u[c] + I Q[-b, d] epsilonmetric[
					-a, -d, c, -e] u[e]) // Simplify, {c}]
			,
			Q3 === Zero,
				Print["Type III"];
				HeadOfTensor[Dagger[Q2[a, b]] (Q2[-b, -a] u[c] + I Q2[-b, d] epsilonmetric[
					-a, -d, c, -e] u[e]) // Simplify, {c}]
			,
			Simplify[(aa^2 / 3) gamma - aa Q2 - bb Q] === Zero,
				Print["Type D"];
				rho = -bb / aa // FullSimplify;
				P = 1 / (3 rho) Q // Simplify;
				Pdag = Dagger[P] // Simplify;
				scrP = HeadOfTensor[P[-a, -b] Pdag[b, -c], {-a, -c}] // Simplify;
					
				scrP2 = HeadOfTensor[Pdag[-a, -b] P[b, -c], {-a, -c}] // Simplify
					;
				dseda = scrP[-a, a] + 1/3 // Simplify;
				S = Simplify[1/4 (1 + 2 / (3 Sqrt[dseda])) (P + Pdag) + 1 / (4 Sqrt[
					dseda]) (scrP + scrP2) + 1/6 (1 + 1 / (3 Sqrt[dseda])) gamma];
				Ch2 = Simplify[(1 + Sqrt[dseda]) / 2];
				v0 = HeadOfTensor[PowerExpand[Ch2 u[a] + I / (2 Sqrt[dseda]) scrP2[
					c, b] epsilonmetric[-c, -b, a, -d] u[d]] // Simplify, {a}];
				v1 = HeadOfTensor[PowerExpand[S[-c, a] w[c] / Sqrt[S[-d, -b] w[d]
					 w[b]]] // Simplify, {a}];
				{v0 + v1, v0 - v1}
			,
			Simplify[6 bb^2 - aa^3] === 0,
				Print["Type II"];
				rho = bb / aa;
				P = rho Q + 2 rho^2 gamma - Q2;
				HeadOfTensor[Dagger[P[a, b]] (P[-b, -a] u[c] + I P[-b, d] epsilonmetric[
					-a, -d, c, -e] u[e]) // Simplify, {c}]
			,
			True,
				Print["Type I"]
		]
	]



(* ::Section:: *)
(*  Classification of type D metrics*)


(* Test when a symbolic function is non-negative *)
(* 
Recall that the value of Assumptions option is always logical statement. 
Therefore it should be expressed in terms of the logical syntax  
*)
Options[SymbolicPositiveQ] := {Assumptions -> True};
SymbolicPositiveQ[x_, OptionsPattern[]] :=
	Block[{$Assumptions = $Assumptions && OptionValue[Assumptions]},
		Which[
			Simplify[x] === 0,
				False
			,
			Simplify[x] === Zero,
				False
			,
			Simplify[x + Abs[x]] === 0,
				False
			,
			Simplify[x - Abs[x]] === 0,
				True
			,
			True,
				"Undefined"
		]
	]


Options[TypeDClassify] = {Assumptions -> True}

TypeDClassify[metric_CTensor, w_CTensor, OptionsPattern[]] :=
	Catch @
		Module[{cart, CD, W, RiemannCD, RicciCD, RicciScalarCD, epsilonmetric,
			 W2, W3, TrW3, rho, drho, dlogrho, alpha, S, P, Q, C3, a, b, c, d, e,
			 f, i, j, k, l, C5, assumptionsC5},
			If[Not @ MetricQ @ metric,
				Throw[Message[PetrovType::nometric, metric]]
			];
			assumptionsC5 = OptionValue[Assumptions];
			cart = Part[metric, 2, 1, -1];
			{a, b, c, d, e, f, i, j, k, l} = GetIndicesOfVBundle[Tangent @ ManifoldOfChart
				 @ cart, 10];
			MetricCompute[metric, cart, All, Parallelize -> True, Verbose -> True
				];
			CD = CovDOfMetric[metric];
			W = Weyl[CD];
			RiemannCD = Riemann[CD];
			RicciCD = Ricci[CD];
			RicciScalarCD = RicciScalar[CD];
			epsilonmetric = epsilon[metric];
			W2 = HeadOfTensor[1/2 W[-a, -b, i, j] W[-i, -j, -c, -d], {-a, -b, -c, -d}] // Simplify;
			W3 = HeadOfTensor[1/2 W2[-a, -b, i, j] W[-i, -j, -c, -d], {-a, -b, -c, -d}] // Simplify;
			TrW3 = 1/2 W3[-i, -j, i, j] // Simplify;
			rho = -(1/12 TrW3) ^ (1/3) // FullSimplify;
			drho = CTensor[Grad[rho, ScalarsOfChart @ cart], {-cart}] // Simplify;
			dlogrho = CTensor[Grad[Log[rho], ScalarsOfChart @ cart], {-cart}] // Simplify;
			alpha = 1/9 metric[-i, -j] dlogrho[i] dlogrho[j] - 2 rho // FullSimplify;
			S = HeadOfTensor[1 / (3 rho) (W[-a, -b, -c, -d] - 
			rho (metric[-a, -c] metric[-b, -d] - metric[-a, -d] metric[-b, -c])), {-a, -b, -c, -d}] // Simplify;
			P = HeadOfTensor[epsilonmetric[-a, -b, -i, -j] W[-k, i, -l, j] drho[k] drho[l], {-a, -b}] // Simplify;
			Q = HeadOfTensor[S[-i, -a, -j, -b] drho[i] drho[j], {-a, -b}] // Simplify;
			C3 = 1/2 S[-a, -b, -i, -j] S[i, j, -c, -d] + S[-a, -b, -c, -d] // Simplify;
			C5 = 2 Q[-i, -j] w[i] w[j] + Q[-k, k] // Simplify;
			Which[
				RicciCD =!= Zero,
					Print["No vacuum"]
				,
				rho === Zero || C3 =!= 0,
					Print["Vacuum no Type D"]
				,
				P === Zero,
					Which[
						SymbolicPositiveQ[C5, Assumptions -> assumptionsC5] === "Undefined",
							
							"Undefined sign in C5"
						,
						SymbolicPositiveQ[C5, Assumptions -> assumptionsC5],
							Which[
								SymbolicPositiveQ[alpha, Assumptions -> assumptionsC5] === "Undefined",
									
									"Undefined sign in \[Alpha]"
								,
								SymbolicPositiveQ[alpha, Assumptions -> assumptionsC5],
									"\!\(\*SubscriptBox[\(A\), \(1\)]\)-metric"
								,
								alpha === Zero,
									"\!\(\*SubscriptBox[\(A\), \(3\)]\)-metric"
								,
								Not @ SymbolicPositiveQ[alpha, Assumptions -> assumptionsC5],
									
									"\!\(\*SubscriptBox[\(A\), \(2\)]\)-metric"
								,
								True,
									(*Should be "Undefined" *)SymbolicPositiveQ[alpha, Assumptions
										 -> assumptionsC5]
							]
						,
						Not @ SymbolicPositiveQ[C5, Assumptions -> assumptionsC5],
							Which[
								SymbolicPositiveQ[alpha, Assumptions -> assumptionsC5] === "Undefined",
									
									"Undefined sign in \[Alpha]"
								,
								SymbolicPositiveQ[alpha, Assumptions -> assumptionsC5],
									"\!\(\*SubscriptBox[\(B\), \(1\)]\)-metric"
								,
								alpha === Zero,
									"\!\(\*SubscriptBox[\(B\), \(3\)]\)-metric"
								,
								Not @ SymbolicPositiveQ[alpha, Assumptions -> assumptionsC5],
									
									"\!\(\*SubscriptBox[\(B\), \(2\)]\)-metric"
								,
								True,
									(*Should be "Undefined" *)
									SymbolicPositiveQ[alpha, Assumptions -> assumptionsC5]
							]
						,
						True,
							(*Should be "Undefined" *)
							SymbolicPositiveQ[alpha, Assumptions -> assumptionsC5]
					]
				,
				True,
					Print["C-metric"];
			]
		]

(* ::Section:: *)
(* Identification of the Kerr metric *)

(* Real and imaginary parts of symbolic complex quantities. *)

SymbolicRe[expr_] := ComplexExpand[(expr + Dagger[expr]) / 2]

symbolicIm[expr_] := ComplexExpand[(expr - Dagger[expr]) / (2 I)]

symbolicComplexNorm2[expr_] := ComplexExpand[expr Dagger[expr]] 

Options[KerrSolutionQ] = {PSimplify -> Simplify}
KerrSolutionQ[metric_CTensor, opts : OptionsPattern[]] :=
Catch @
		Module[{cart, cd, weylcd, epsilonmetric, weyldual, weylselfdual, g2form, weylselfdual2, weylselfdual3, aa, bb,
			 rho, a1, b1, c1, d1, e1, f1, simplf, w, z, riccicd, z1, z2, modz},
			If[Not @ MetricQ @ metric,
				Throw[Message[PetrovType::nometric, metric]]
			];
			simplf = OptionValue[PSimplify];
			cart = Part[metric, 2, 1, -1];
			{a1, b1, c1, d1, e1, f1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];
			(* If we compute Weyl here then we speed up the computation of Ricci below *)
			weylcd = weylConcomitant["Weyl"][metric, opts];
			cd = CovDOfMetric[metric];
			riccicd = simplf[Ricci[cd]];
			If[
				riccicd === Zero,
				Print["Vacuum"];
				If[PetrovType[metric, opts] === "Type D",
					Print["Type D"];
					If[simplf@ Antisymmetrize[weylConcomitant["TensorXi"][metric, opts][-a1, -b1] weylConcomitant["TensorXi"][metric, opts][-c1, -d1], {-b1, -c1}] === 0,
						Print["Kerr NUT"];
						If[symbolicIm[z^3 Dagger[w^8]] === 0,
							Print["Complex Kerr"];
							weyldual = weylConcomitant["WeylDual"][metric, opts];
							weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
							g2form = metricConcomitant["G2Form"][metric, opts];
							weylselfdual2 = weylConcomitant["WeylSelfDual2"][metric, opts];
							weylselfdual3 = weylConcomitant["WeylSelfDual3"][metric, opts];
							aa = 2 weylConcomitant["TraceWeylSelfDual2"][metric, opts];
							bb = 2 weylConcomitant["TraceWeylSelfDual3"][metric, opts];
							w = simplf[-bb / (2 aa)];
							z = weylConcomitant["ScalarZ"][metric, opts];
							z1 = simplf@ symbolicRe[z^3 Dagger[w^8]];
							z2 = simplf@ symbolicRe[w^3 Dagger[z]];
							modz = simplf@ symbolicComplexNorm2[z];
							If[SymbolicPositiveQ[-z1 / (18 z2 - modz)^3, opts],
								Print["Kerr"];
								True
							],
							False
						],
						False
					],
					False
				],
				False
			]
		]


(* ::Section:: *)
(*  Determination of the dimension of the isometry group*)

IsometryGroupDimension[metric_CTensor, H_CTensor, opts___] :=
 	Catch@ 
  		Module[{simplf, C1, C2, C3, C4, C11, C12, C122, C123, C1233, 
    			C1234, C1222, C1223, C111, C112, C1122, C1123, C1111, C1112},
   			If[Not@MetricQ@metric, 
    				Throw[Message[IsometryGroupDimension::nometric, metric]]];
			simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
   			C1 = RframeConcomitant["C1"][metric, H, opts];
   			Which[
      
    				C1 === Zero,
    					Print["\!\(\*SubscriptBox[\(G\), \(4\)]\)"],
	 
    				C11 = RframeConcomitant["C11"][metric, H, opts];
    				C11 === Zero,
    					Which[
	 
     						C12 = RframeConcomitant["C12"][metric, H, opts];
     						C12 === Zero,
     							Print["\!\(\*SubscriptBox[\(G\), \(3\)]\)"],
     				
     						C122 = RframeConcomitant["C122"][metric, H, opts];
     						C122 === Zero,
     							Which[
	    
      								C123 = RframeConcomitant["C123"][metric, H, opts];
      								C123 === Zero,
      									Print["\!\(\*SubscriptBox[\(G\), \(2  b\)]\)"],
      						
      								C1233 = RframeConcomitant["C1233"][metric, H, opts];
      								Not[C1233 == Zero],
      									Print["No symmetries"],
      						
      								C1234 = RframeConcomitant["C1234"][metric, H, opts];
      								C1234 === Zero,
      									Print["\!\(\*SubscriptBox[\(G\), \(1  d\)]\)"],
	       
      								True,
      									Print["No symmetries"]
      							],
     				
     						C1222 = RframeConcomitant["C1222"][metric, H, opts];
     						Not[C1222 === Zero],
     							Print["No symmetries"],
     				
     						C1223 = RframeConcomitant["C1223"][metric, H, opts];
     						C1223 === Zero,
     							Print["\!\(\*SubscriptBox[\(G\), \(1  c\)]\)"],
	    
     						True,
     							Print["No symmetries"]
     					],
    				C111 = RframeConcomitant["C111"][metric, H, opts];
    				C111 === Zero,
    					Which[
     				
     						C112 = RframeConcomitant["C112"][metric, H, opts];
     						C112 === Zero,
     							Print["\!\(\*SubscriptBox[\(G\), \(2  a\)]\)"],
     				
     						C1122 = RframeConcomitant["C1122"][metric, H, opts];
     						Not[C1122 === Zero],
     							Print["No symmetries"],
     				
     						C1123 = RframeConcomitant["C1123"][metric, H, opts];
     						C1123 === Zero,
     							Print["\!\(\*SubscriptBox[\(G\), \(1  b\)]\)"],
	    
     						True,
     							Print["No symmetries"]
     					],
    				C1111 = RframeConcomitant["C1111"][metric, H, opts];
    				Not[C1111 === Zero],
    					Print["No symmetries"],
	 
    				C1112 = RframeConcomitant["C1112"][metric, H, opts];
    				C1112 === Zero,
    					Print["\!\(\*SubscriptBox[\(G\), \(1  a\)]\)"],
	 
    				True,
    					Print["No symmetries C1112"]		
    			]
   
   		]



(****************************************************************)

(****************** 5. End private and package ******************)

(****************************************************************)

End[];

EndPackage[];
