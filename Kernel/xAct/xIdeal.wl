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
		{a1, b1, c1, d1, e1, f1} = GetIndicesOfVBundle[Tangent @ ManifoldOfChart@ cart, 6];
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
		{a1, b1, c1, d1, e1, f1} = GetIndicesOfVBundle[Tangent @ ManifoldOfChart@ cart, 6];
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
		{a1, b1, c1, d1, e1, f1} = GetIndicesOfVBundle[Tangent @ ManifoldOfChart@ cart, 6];
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		simplf[HeadOfTensor[1/2 weylselfdual[-a1, -b1, e1, f1] weylselfdual[-e1, -f1, -c1, -d1], {-a1, -b1, -c1, -d1}]]
	]
)

weylConcomitant["WeylSelfDual3"][metric_CTensor, opts___] :=
(weylConcomitant["WeylSelfDual3"][metric, opts] = 
	Module[{simplf, weylselfdual, weylselfdual2, cart, a1, b1, c1, d1, e1, f1},
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1} = GetIndicesOfVBundle[Tangent @ ManifoldOfChart@ cart, 6];
		simplf = (PSimplify /. FilterRules[{opts}, PSimplify]);
		weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		weylselfdual2 = weylConcomitant["WeylSelfDual2"][metric, opts];
		simplf[HeadOfTensor[1/2 weylselfdual2[-a1, -b1, e1, f1] weylselfdual[-e1, -f1, -c1, -d1], {-a1, -b1, -c1, -d1}]]
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
			 rho, a1, b1, c1, d1, e1, f1, simplf},
			If[Not @ MetricQ @ metric,
				Throw[Message[PetrovType::nometric, metric]]
			];
			simplf = OptionValue[PSimplify];
			cart = Part[metric, 2, 1, -1];
			{a1, b1, c1, d1, e1, f1} = GetIndicesOfVBundle[Tangent @ ManifoldOfChart@ cart, 6];
			epsilonmetric = epsilon[metric];
			weylcd = weylConcomitant["Weyl"][metric, opts];
			weyldual = weylConcomitant["WeylDual"][metric, opts];
			weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
			g2form = metricConcomitant["G2Form"][metric, opts];
			weylselfdual2 = weylConcomitant["WeylSelfDual2"][metric, opts];
			aa = simplf[1/2 weylselfdual2[-a1, -b1, a1, b1]];
			weylselfdual3 = weylConcomitant["WeylSelfDual3"][metric, opts];
			bb = simplf[weylselfdual3[-a1, -b1, a1, b1] / 2];
			rho = simplf[-bb / aa];
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
(*  Determination of the dimension of the isometry group*)


IsometryGroupDimension[metric_CTensor, e0_CTensor, e1_CTensor, e2_CTensor, e3_CTensor] :=
 	Catch@ 
  		Module[{cart, CD,
    			epsilonmetric, H, C1, C2, C3, C4, C11, C12, C122, C123, C1233, 
    			C1234, C1222, C1223, C111, C112, C1122, C1123, C1111, C1112, a, b,
     			c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z},
   			If[Not@MetricQ@metric, 
    				Throw[Message[IsometryGroupDimension::nometric, metric]]];
   			cart = Part[metric, 2, 1, -1];
   			{a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, r, s, t, u, v, w,
     			 x, y, z} = GetIndicesOfVBundle[Tangent@ManifoldOfChart@cart, 25];
   			MetricCompute[metric, cart, All, Parallelize -> True, 
    				Verbose -> True];
   			CD = CovDOfMetric[metric];
   			epsilonmetric = epsilon[metric];
   			H = HeadOfTensor[-(1/2) (-(CD[a][e0[b]] e0[c] - CD[a][e0[c]] e0[b]) + (CD[a][e1[b]] e1[c] - CD[a][e1[c]] e1[b]) + (CD[a][e2[b]] e2[c] - 
           			CD[a][e2[c]] e2[b]) + (CD[a][e3[b]] e3[c] - CD[a][e3[c]] e3[b])), {a, b, c}] // Simplify;
   			C1 = Simplify[HeadOfTensor[CD[-a][H[-b, -c, -d]] + H[-a, -b, i] H[-i, -c, -d] + H[-a, -c, j] H[-b, -j, -d] + 
				H[-a, -d, k] H[-b, -c, -k], {-a, -b, -c, -d}]];
   			Which[
      
    				C1 === Zero,
    					Print["\!\(\*SubscriptBox[\(G\), \(4\)]\)"],
	 
    				C11 = Simplify[HeadOfTensor[epsilonmetric[i, j, a, b] C1[-i, -c, -d, -e] C1[-j, -f, -g, -h], {a, b, -c, -d, -e, -f, -g, -h}]];
    				C11 === Zero,
    					Which[
	 
     						C2 = Simplify[HeadOfTensor[CD[-a][C1[-b, -c, -d, -e]] + H[-a, -b, i] C1[-i, -c, -d, -e] +
	   						H[-a, -c, j] C1[-b, -j, -d, -e] + H[-a, -d, k] C1[-b, -c, -k, -e] + 
         						H[-a, -e, l] C1[-b, -c, -d, -l], {-a, -b, -c, -d, -e}]];
     						C12 = Simplify[HeadOfTensor[epsilonmetric[i, j, a, b] C1[-i, -c, -d, -e] C2[-j, -f, -g, -h, -k], 
	   						{a, b, -c, -d, -e, -f, -g, -h, -k}]];
     						C12 === Zero,
     							Print["\!\(\*SubscriptBox[\(G\), \(3\)]\)"],
     				
     						C122 = Simplify[HeadOfTensor[epsilonmetric[i, j, k, a] C1[-i, -b, -c, -d] C2[-j, -e, -f, -g, -h] 
	   						C2[-k, -l, -m, -n, -o], {a, -b, -c, -d, -e, -f, -g, -h, -l, -m, -n, -o}]];
     						C122 === Zero,
     							Which[
	    
      						  		C3 = Simplify[HeadOfTensor[CD[-a][C2[-b, -c, -d, -e, -f]] + H[-a, -b, i] C2[-i, -c, -d, -e, -f] + 
              								H[-a, -c, j] C2[-b, -j, -d, -e, -f] + H[-a, -d, k] C2[-b, -c, -k, -e, -f] + 
          								H[-a, -e, l] C2[-b, -c, -d, -l, -f] + H[-a, -f, m] C2[-b, -c, -d, -e, -m], 
		  							{-a, -b, -c, -d, -e, -f}]];
      								C123 = Simplify[HeadOfTensor[epsilonmetric[i, j, k, a] C1[-i, -b, -c, -d] C2[-j, -e, -f, -g, -h] 
	      								C3[-k, -l, -m, -n, -o, -p], {a, -b, -c, -d, -e, -f, -g, -h, -l, -m, -n, -o, -p}]];
      								C123 === Zero,
      									Print["\!\(\*SubscriptBox[\(G\), \(2  b\)]\)"],
      						
      								C1233 = Simplify[HeadOfTensor[epsilonmetric[i, j, k, l] C1[-i, -a, -b, -c] C2[-j, -d, -e, -f, -g] 
	      								C3[-k, -h, -m, -n, -o, -p] C3[-l, -q, -r, -s, -t, -u], 
	      								{-a, -b, -c, -d, -e, -f, -g, -h, -m, -n, -o, -p, -q, -r, -s, -t, -u}]];
      								Not[C1233 == Zero],
      									Print["No symmetries"],
      						
      								C4 = Simplify[HeadOfTensor[CD[-a][C3[-b, -c, -d, -e, -f, -g]] + H[-a, -b, i] C3[-i, -c, -d, -e, -f, -g] + 
          								H[-a, -c, j] C3[-b, -j, -d, -e, -f, -g] + H[-a, -d, k] C3[-b, -c, -k, -e, -f, -g] + 
           								H[-a, -e, l] C3[-b, -c, -d, -l, -f, -g] + H[-a, -f, m] C3[-b, -c, -d, -e, -m, -g] + 
          								H[-a, -g, n] C3[-b, -c, -d, -e, -f, -n], {-a, -b, -c, -d, -e, -f, -g}]];
      								 C1234 = Simplify[HeadOfTensor[epsilonmetric[i, j, k, l] C1[-i, -a, -b, -c] C2[-j, -d, -e, -f, -g] 
	       								C3[-k, -h, -m, -n, -o, -p] C4[-l, -q, -r, -s, -t, -u, -v], 
									{-a, -b, -c, -d, -e, -f, -g, -h, -m, -n, -o, -p, -q, -r, -s, -t, -u, -v}]];
      								C1234 === Zero,
      									Print["\!\(\*SubscriptBox[\(G\), \(1  d\)]\)"],
	       
      								True,
      									Print["No symmetries"]
      							],
     				
     						C1222 = Simplify[HeadOfTensor[epsilonmetric[i, j, k, l] C1[-i, -a, -b, -c] C2[-j, -d, -e, -f, -g] 
	   						C2[-k, -h, -m, -n, -o] C2[-l, -p, -q, -r, -s], 
	  						{-a, -b, -c, -d, -e, -f, -g, -h, -m, -n, -o, -p, -q, -r, -s}]];
     						Not[C1222 === Zero],
     							Print["No symmetries"],
     				
     						C3 = Simplify[HeadOfTensor[CD[-a][C2[-b, -c, -d, -e, -f]] + H[-a, -b, i] C2[-i, -c, -d, -e, -f] + 
         						H[-a, -c, j] C2[-b, -j, -d, -e, -f] + H[-a, -d, k] C2[-b, -c, -k, -e, -f] + 
         						H[-a, -e, l] C2[-b, -c, -d, -l, -f] + H[-a, -f, m] C2[-b, -c, -d, -e, -m],
	       						{-a, -b, -c, -d, -e, -f}]];
     						C1223 = Simplify[HeadOfTensor[epsilonmetric[i, j, k, l] C1[-i, -a, -b, -c] C2[-j, -d, -e, -f, -g] 
	   						C2[-k, -h, -m, -n, -o] C3[-l, -p, -q, -r, -s, -t], 
	  						{-a, -b, -c, -d, -e, -f, -g, -h, -m, -n, -o, -p, -q, -r, -s, -t}]];
     						C1223 === Zero,
     							Print["\!\(\*SubscriptBox[\(G\), \(1  c\)]\)"],
	    
     						True,
     							Print["No symmetries"]
     					],
    				C111 = Simplify[HeadOfTensor[epsilonmetric[i, j, k, a] C1[-i, -b, -c, -d] C1[-j, -e, -f, -g] C1[-k, -h, -m, -n],
					{a, -b, -c, -d, -e, -f, -g, -h, -m, -n}]];
    				C111 === Zero,
    					Which[
     				
     						C2 = Simplify[HeadOfTensor[CD[-a][C1[-b, -c, -d, -e]] + H[-a, -b, i] C1[-i, -c, -d, -e] +
          						H[-a, -c, j] C1[-b, -j, -d, -e] + H[-a, -d, k] C1[-b, -c, -k, -e] + H[-a, -e, l] C1[-b, -c, -d, -l], 
	       						{-a, -b, -c, -d, -e}]];
     						C112 = Simplify[HeadOfTensor[epsilonmetric[i, j, k, a] C1[-i, -b, -c, -d] C1[-j, -e, -f, -g] C2[-k, -h, -m, -n, -o], 
	   						{a, -b, -c, -d, -e, -f, -g, -h, -m, -n, -o}]];
     						C112 === Zero,
     							Print["\!\(\*SubscriptBox[\(G\), \(2  a\)]\)"],
     				
     						C1122 = Simplify[HeadOfTensor[epsilonmetric[i, j, k, l] C1[-i, -a, -b, -c] C1[-j, -d, -e, -f] C2[-k, -g, -h, -m, -n] 
	   						C2[-l, -o, -p, -q, -r], {-a, -b, -c, -d, -e, -f, -g, -h, -m, -n, -o, -p, -q, -r}]];
     						Not[C1122 === Zero],
     							Print["No symmetries"],
     				
     						C3 = Simplify[HeadOfTensor[CD[-a][C2[-b, -c, -d, -e, -f]] + H[-a, -b, i] C2[-i, -c, -d, -e, -f] + 
         						H[-a, -c, j] C2[-b, -j, -d, -e, -f] + H[-a, -d, k] C2[-b, -c, -k, -e, -f] + 
         						H[-a, -e, l] C2[-b, -c, -d, -l, -f] + H[-a, -f, m] C2[-b, -c, -d, -e, -m], 
	       							{-a, -b, -c, -d, -e, -f}]];
     						C1123 = Simplify[HeadOfTensor[epsilonmetric[i, j, k, l] C1[-i, -a, -b, -c] C1[-j, -d, -e, -f] C2[-k, -g, -h, -m, -n] 
	   						C3[-l, -o, -p, -q, -r, -s], {-a, -b, -c, -d, -e, -f, -g, -h, -m, -n, -o, -p, -q, -r, -s}]];
     						C1123 === Zero,
     							Print["\!\(\*SubscriptBox[\(G\), \(1  b\)]\)"],
	    
     						True,
     							Print["No symmetries"]
     					],
    				C1111 = Simplify[HeadOfTensor[epsilonmetric[i, j, k, l] C1[-i, -a, -b, -c] C1[-j, -d, -e, -f] C1[-k, -g, -h, -m] C1[-l, -n, -o, -p], 
					{-a, -b, -c, -d, -e, -f, -g, -h, -m, -n, -o, -p}]];
    				Not[C1111 === Zero],
    					Print["No symmetries"],
	 
    				C2 = Simplify[HeadOfTensor[CD[-a][C1[-b, -c, -d, -e]] + H[-a, -b, i] C1[-i, -c, -d, -e] + H[-a, -d, j] C1[-b, -c, -j, -e] + 
        				H[-a, -d, k] C1[-b, -c, -k, -e] + H[-a, -e, l] C1[-b, -c, -d, -l], {-a, -b, -c, -d, -e}]];
    				C1112 = Simplify[HeadOfTensor[epsilonmetric[i, j, k, l] C1[-i, -a, -b, -c] C1[-j, -d, -e, -f] C1[-k, -g, -h, -m] C2[-l, -n, -o, -p, -q], 
					{-a, -b, -c, -d, -e, -f, -g, -h, -m, -n, -o, -p, -q}]];
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
