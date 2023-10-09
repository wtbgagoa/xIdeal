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

(* ::Section:: *)
(* Messages *)

PetrovType::nometric = "Metric `1` has not been registered as a metric";

(* ::Section:: *)
(* BeginPrivate *)

Begin["`Private`"]

(******************************************************************************)

(********************* 2. Basic structures ************************************)

(******************************************************************************)

(* ::Section:: *)
(* Computation of the Petrov types *)

(*
Explanation
*)

PetrovType[metric_CTensor] :=
	Catch @
		Module[{cart, CD, WeylCD, RiemannCD, RicciCD, RicciScalarCD, epsilonmetric,
			 WeylDual, WeylSelfDual, G2Form, WeylSelfDual2, WeylSelfDual3, aa, bb,
			 rho, a, b, c, d, e, f, i, j},
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
			G2Form = Simplify[HeadOfTensor[1/2 (-I epsilonmetric[-a, -b, -c, -
				d] + metric[-a, -c] metric[-b, -d] - metric[-a, -d] metric[-b, -c]), 
				{-a, -b, -c, -d}]];
			WeylSelfDual2 = Simplify[HeadOfTensor[1/2 WeylSelfDual[-a, -b, i, 
				j] WeylSelfDual[-i, -j, -c, -d], {-a, -b, -c, -d}]];
			aa = Simplify[1/2 WeylSelfDual2[-a, -b, a, b]];
			WeylSelfDual3 = Simplify[HeadOfTensor[1/2 WeylSelfDual2[-a, -b, i,
				 j] WeylSelfDual[-i, -j, -c, -d], {-a, -b, -c, -d}]];
			bb = Simplify[WeylSelfDual3[-a, -b, a, b] / 2];
			rho = Simplify[-bb / aa];
			Which[
				WeylSelfDual2 === Zero,
					Print["Type N"]
				,
				WeylSelfDual3 === Zero,
					Print["Type III"]
				,
				Simplify[aa WeylSelfDual2 - aa^2 / 3 G2Form - bb WeylSelfDual] ===
					 Zero,
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

(* TODO: we should have private or public functions for the different 
concommitants *)

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
(*To do: fix the fact that assumptions added through options don't work properly.*)
Options[SymbolicPositiveQ] := {Assumptions -> True};
SymbolicPositiveQ[x_, opts : OptionsPattern[]] :=
	Block[{$Assumptions = $Assumptions},	
		$Assumptions = Join[{OptionValue[Assumptions]}, $Assumptions];
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

(****************************************************************)

(****************** 5. End private and package ******************)

(****************************************************************)

End[];

EndPackage[];
