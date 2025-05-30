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
	 "xAct`xCore`", "xAct`ExactSolsData`"}]

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


PetrovType::usage = "PetrovType[metric] returns the Petrov Type of metric.";

DebeverNullDirections::usage = "DebeverNullDirections[metric] returns the multiple Debever null directions of metric.";

TypeDClassify::usage = "TypeDClassify[metric, w] returns the subfamily of vacuum Type D solutions to which metric belongs. To do so, it needs an arbitrary unitary time-like vector w.";

PSimplify::usage = " ";

ConnectionTensor::usage = "ConnectionTensor[metric, opts] returns the connection tensor of metric.";

IsometryGroupDimension::usage = "IsometryGroupDimension[metric, H] returns the dimension of the maximal isometry group admitted by metric.";

KerrSolutionQ::usage = "KerrSolutionQ[metric] returns True if metric is Kerr metric and False otherwise.";

ClearxIdealCache::usage = " ";

SaveExactSolution::usage = " ";

PerfectFluidQ::usage = "PerfectFluidQ[metric, w] returns True if metric is of the perfect fluid type. To do so, it needs an arbitrary unitary time-like vector w.";

PerfectFluidVariables::usage = "PerfectFluidVariables[metric, w] returns a list with the energy density, the pressure and the fluid flow of metric if it is of the perfect fluid type. To do so, it needs an arbitrary unitary time-like vector w.";

ThermodynamicPerfectFluidQ::usage = "ThermodynamicPerfectFluidQ[metric, w] returns True if metric is of the thermodynamic perfect fluid type. To do so, it needs an arbitrary unitary time-like vector w.";

GenericIdealGasQ::usage = "GenericIdealGasQ[metric, w] returns True if metric is of the generic ideal gas type. To do so, it needs an arbitrary unitary time-like vector w.";

Rframe::usage = " ";

(* ::Section:: *)
(* Messages *)


PetrovType::nometric = "Metric `1` has not been registered as a metric";

PetrovType::noobserver = "Value `1` for \"Observer\" is invalid";

PetrovType::nospatialmetric = "Invalid spatial metric for spacetime metric `1` or \" Observer\" `2`";

PetrovType::nopsimplify = "Value `1` for \"PSimplify\" is invalid\" ";

DebeverNullDirections::nometric = "Metric `1` has not been registered as a metric";

TypeDClassify::nometric = "Metric `1` has not been registered as a metric";

KerrSolutionQ::nometric = "Metric `1` has not been registered as a metric";

PerfectFluidQ::nometric = "Metric `1` has not been registered as a metric";

PerfectFluidVariables::nometric = "Metric `1` has not been registered as a metric";

PerfectFluidVariables::noperfectfluid = "Metric `1` is not of the perfect fluid type";

ThermodynamicPerfectFluidQ::nometric = "Metric `1` has not been registered as a metric";

ThermodynamicPerfectFluidQ::noperfectfluid = "Metric `1` is not of the perfect fluid type";

GenericIdealGasQ::nometric = "Metric `1` has not been registered as a metric";

GenericIdealGasQ::noperfectfluid = "Metric `1` is not of the perfect fluid type";

GenericIdealGasQ::nothermodynamicperfectfluid = "Metric `1` is not of the thermodynamic perfect fluid type";

ConnectionTensor::nometric = "Metric `1` has not been registered as a metric";

IsometryGroupDimension::nometric = "Metric `1` has not been registered as a metric";

(* ::Section:: *)
(* BeginPrivate *)


Begin["`Private`"]

(* ::Section:: *)
(* Utilities *)
applyfunc[expr_, func_, head_, errortag_String] := 
	With[{res = func[expr]},
		If[Head[res] === func, Throw[Message[MessageName[head, errortag], func]], res]
	]

(* ::Section:: *)
(* Computation of the metric concomitants *)
(* TODO: we should avoid defining default options for private functions *)
Options[metricConcomitant] = {PSimplify -> $CVSimplify, Parallelize -> True, Verbose -> True, "Observer" -> Null, "Vector" -> Null, "Bivector" -> Null, "Method" -> Default}

metricConcomitant["G"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["G"][metric, opts] = 
	Module[{simplf, cart, a1, b1, c1, d1, gmetric, time, vb},
		{simplf, vb} = OptionValue[metricConcomitant, {opts}, {PSimplify, Verbose}];
		cart = 	Part[metric, 2, 1, -1];
		{a1, b1, c1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 4];
		time = AbsoluteTime[];
		gmetric = metric[-a1, -b1] metric[-c1, -d1] - metric[-a1, -d1] metric[-c1, -b1];
		If[vb,
			Print["** ReportCompute: computing metric concomitant \"G\" in", AbsoluteTime[] - time, " seconds:"];
		];
		gmetric = HeadOfTensor[gmetric, {-a1, -b1, -c1, -d1}];
		time = AbsoluteTime[];
		simplf[gmetric];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"G\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		gmetric
	]
)

metricConcomitant["G2Form"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["G2Form"][metric, opts] = 
	Module[{simplf, cart, a1, b1, c1, d1, epsilonmetric, gmetric, time, vb},
		{simplf, vb} = OptionValue[metricConcomitant, {opts}, {PSimplify, Verbose}];
		cart = 	Part[metric, 2, 1, -1];
		{a1, b1, c1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 4];
		epsilonmetric = epsilon[metric];
		time = AbsoluteTime[];
		gmetric = 1/2 (-I epsilonmetric[-a1, -b1, -c1, -d1] + 
			metric[-a1, -c1] metric[-b1, -d1] - metric[-a1, -d1] metric[-b1, -c1]);
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"G2Form\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		gmetric = HeadOfTensor[gmetric, {-a1, -b1, -c1, -d1}];
		time = AbsoluteTime[];
		gmetric = simplf[gmetric];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"G2Form\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		gmetric
	]
)

metricConcomitant["SpatialMetric"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["SpatialMetric"][metric, opts] = 
	Module[{simplf, cart, obs, a1, b1, vb, time, smetric},
		{obs, simplf, vb} = OptionValue[metricConcomitant, {opts}, {"Observer", PSimplify, Verbose}];
		If[Head[obs] =!= CTensor || TensorRank[Part[obs, 1]] =!= 1,
			Throw[Message[PetrovType::noobserver, obs]]
		];
		cart = 	Part[metric, 2, 1, -1];
		{a1, b1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		Which[
			SignatureOfMetric[metric] === {3, 1, 0},
				If[Simplify[obs[a1]metric[-a1, -b1]obs[b1]] =!= -1,
					Throw[Message[PetrovType::noobserver, obs]]
				];
				time = AbsoluteTime[];
				smetric = HeadOfTensor[metric[-a1, -b1] + obs[-a1] obs[-b1], {-a1, -b1}];
				If[vb, 
					Print["** ReportCompute: computing metric concomitant \"SpatialMetric\" in ", AbsoluteTime[] - time, " seconds:"]
				];
				time = AbsoluteTime[];
				smetric = applyfunc[smetric, simplf, PetrovType, "nopsimplify"];
				If[vb,
					Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"SpatialMetric\" in ", AbsoluteTime[] - time, " seconds:"]
				];
				smetric,
			SignatureOfMetric[metric] === {1, 3, 0},
				If[Simplify[obs[a1]metric[-a1, -b1]obs[b1]]  =!= 1,
					Throw[Message[PetrovType::noobserver, obs]]
				];
				time = AbsoluteTime[];
				smetric = HeadOfTensor[metric[-a1, -b1] - obs[-a1] obs[-b1], {-a1, -b1}];
				If[vb, 
					Print["** ReportCompute: computing metric concomitant \"SpatialMetric\" in ", AbsoluteTime[] - time, " seconds:"]
				];
				time = AbsoluteTime[];
				smetric = applyfunc[smetric, simplf, PetrovType, "nopsimplify"];
				If[vb,
					Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"SpatialMetric\" in ", AbsoluteTime[] - time, " seconds:"]
				];
				smetric,
			True,
				Throw[Message[PetrovType::nospatialmetric, metric, obs]]
		]
	]
)

metricConcomitant["Ricci"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["Ricci"][metric, opts] = 
	Module[{cd},
		cd = CovDOfMetric[metric];
		weylConcomitant["Weyl"][metric, opts];
		Ricci[cd]
	]
)

metricConcomitant["RicciScalar"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["RicciScalar"][metric, opts] = 
	Module[{cd},
		cd = CovDOfMetric[metric];
		weylConcomitant["Weyl"][metric, opts];
		RicciScalar[cd]
	]
)

metricConcomitant["SchoutenTensor"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["SchoutenTensor"][metric, opts] = 
	Module[{simplf, cart, cd, ricciscalarcd, ricci, a1, b1, schouten, vb, time},
		{simplf, vb} = OptionValue[metricConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		cd = CovDOfMetric[metric];
		ricciscalarcd = metricConcomitant["RicciScalar"][metric, opts];
		ricci = Ricci[cd];
		time = AbsoluteTime[];
		schouten = 1/2 (ricci[-a1, -b1] - 1/6 ricciscalarcd[] metric[-a1, -b1]);
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"SchoutenTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		schouten = HeadOfTensor[schouten, {-a1, -b1}];
		time = AbsoluteTime[];
		schouten = simplf[schouten];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"SchoutenTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		schouten
	]
)

metricConcomitant["STensor"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["STensor"][metric, opts] = 
	Module[{simplf, cart, cd, ricciscalarcd, ricci, a1, b1, stensor, vb, time},
		{simplf, vb} = OptionValue[metricConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		cd = CovDOfMetric[metric];
		ricciscalarcd = metricConcomitant["RicciScalar"][metric, opts];
		ricci = Ricci[cd];
		time = AbsoluteTime[];
		stensor = ricci[-a1, -b1] - 1/4 ricciscalarcd[] metric[-a1, -b1];
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"STensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		stensor = HeadOfTensor[stensor, {-a1, -b1}];
		time = AbsoluteTime[];
		stensor = simplf[stensor];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"STensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		stensor
	]
)

metricConcomitant["qScalar"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["qScalar"][metric, opts] = 
	Module[{simplf, cart, st, s2, s3, trS2,trS3, a1, b1, c1, q, vb, time},
		{simplf, vb} = OptionValue[metricConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 3];
		st = metricConcomitant["STensor"][metric, opts];
		time = AbsoluteTime[];
		s2 = st[-a1, -c1] st[c1, -b1];
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"STensor2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		s2 = HeadOfTensor[s2, {-a1, -b1}];
		time = AbsoluteTime[];
		s2 = simplf[s2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"STensor2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		s3 = s2[-a1, -c1] st[c1, -b1];
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"STensor3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		s3 = HeadOfTensor[s3, {-a1, -b1}];
		time = AbsoluteTime[];
		s3 = simplf[s3];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"STensor3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		trs2 = s2[-a1, a1];
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"TrSTensor2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		trs2 = simplf[trs2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"TrSTensor2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		trs3 = s3[-a1, a1];
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"TrSTensor3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		trs3 = simplf[trs3];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"TrSTensor3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		q = -2 trs3 / trs2;
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"qScalar\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		q = simplf[q];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"qScalar\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		q
	]
)

metricConcomitant["QTensor"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["QTensor"][metric, opts] = 
	Module[{simplf, cart, st, q, a1, b1, qtensor, vb, time},
		{simplf, vb} = OptionValue[metricConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		st = metricConcomitant["STensor"][metric, opts];
		q = metricConcomitant["qScalar"][metric, opts];
		time = AbsoluteTime[];
		qtensor = st - 1/4 q metric;
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"QTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		qtensor = simplf[qtensor];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"QTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		qtensor
	]
)

metricConcomitant["QTensor2"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["QTensor2"][metric, opts] = 
	Module[{simplf, cart, qt, q2, a1, b1, c1, vb, time},
		{simplf, vb} = OptionValue[metricConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 3];
		qt = metricConcomitant["QTensor"][metric, opts];
		time = AbsoluteTime[];
		q2 = qt[-a1, -b1] qt[b1, -c1];
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"QTensor2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		q2 = HeadOfTensor[q2, {-a1, -c1}];
		time = AbsoluteTime[];
		q2 = simplf[q2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"QTensor2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		q2
	]
)

metricConcomitant["FlowProjector"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["FlowProjector"][metric, opts] = 
	Module[{simplf, cart, qt, q, a1, b1, fp, vb, time},
		{simplf, vb} = OptionValue[metricConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		qt = metricConcomitant["QTensor"][metric, opts];
		q = metricConcomitant["qScalar"][metric, opts];
		time = AbsoluteTime[];
		fp = qt / q;
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"FlowProjector\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		fp = simplf[fp];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"FlowProjector\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		fp
	]
)

metricConcomitant["FluPerCond1"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["FluPerCond1"][metric, opts] = 
	Module[{simplf, cart, qt, q, q2, a1, b1, cond1, vb, time},
		{simplf, vb} = OptionValue[metricConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		qt = metricConcomitant["QTensor"][metric, opts];
		q = metricConcomitant["qScalar"][metric, opts];
		q2 = metricConcomitant["QTensor2"][metric, opts];
		time = AbsoluteTime[];
		cond1 = q2 + q qt;
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"FluPerCond1\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		cond1 = simplf[cond1];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"FluPerCond1\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		cond1
	]
)

metricConcomitant["FluPerCond2"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["FluPerCond2"][metric, opts] = 
	Module[{simplf, cart, qt, q, v, a1, b1, cond1, vb, time},
		{v, simplf, vb} = OptionValue[metricConcomitant, {opts}, {"Vector", PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		qt = metricConcomitant["QTensor"][metric, opts];
		q = metricConcomitant["qScalar"][metric, opts];
		time = AbsoluteTime[];
		cond2 = q qt[-a1, -b1] v[a1] v[b1];
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"FluPerCond2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		cond2 = simplf[cond2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"FluPerCond2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		cond2
	]
)

metricConcomitant["rScalar"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["rScalar"][metric, opts] = 
	Module[{simplf, cart, cd, riccicd, r, a1, vb, time},
		{simplf, vb} = OptionValue[metricConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 1];
		cd = CovDOfMetric[metric];
		riccicd = metricConcomitant["Ricci"][metric, opts];
		time = AbsoluteTime[];
		r = riccicd[-a1, a1];
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"rScalar\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		r = simplf[r];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"rScalar\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		r
	]
)

metricConcomitant["EnergyDensity"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["EnergyDensity"][metric, opts] = 
	Module[{simplf, q, r, edens, vb, time},
		{simplf, vb} = OptionValue[metricConcomitant, {opts}, {PSimplify, Verbose}];
		q = metricConcomitant["qScalar"][metric, opts];
		r = metricConcomitant["rScalar"][metric, opts];
		time = AbsoluteTime[];
		edens = (3 q + r) / 4;
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"EnergyDensity\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		edens = simplf[edens];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"EnergyDensity\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		edens
	]
)

metricConcomitant["Pressure"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["Pressure"][metric, opts] = 
	Module[{simplf, q, r, press, vb, time},
		{simplf, vb} = OptionValue[metricConcomitant, {opts}, {PSimplify, Verbose}];
		q = metricConcomitant["qScalar"][metric, opts];
		r = metricConcomitant["rScalar"][metric, opts];
		time = AbsoluteTime[];
		press = (q - r) / 4;
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"Pressure\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		press = simplf[press];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"Pressure\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		press
	]
)

metricConcomitant["FluPerFlow"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["FluPerFlow"][metric, opts] = 
	Module[{simplf, cart, proj, v, a1, b1, flow, vb, time},
		{v, simplf, vb} = OptionValue[metricConcomitant, {opts}, {"Vector", PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		proj = metricConcomitant["FlowProjector"][metric, opts];
		time = AbsoluteTime[];
		flow = proj[-a1, -b1] v[b1] / Sqrt[proj[-a1, -b1] v[b1] v[a1]];
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"FluPerFlow\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		flow = HeadOfTensor[flow, {-a1}];
		time = AbsoluteTime[];
		flow = simplf[flow];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"FluPerFlow\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		flow
	]
)

metricConcomitant["IndicatrixFunction"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["IndicatrixFunction"][metric, opts] = 
	Module[{simplf, cart, cd, flow, press, edens, indicatrix, a1, b1, vb, time},
		{simplf, vb} = OptionValue[metricConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		cd = CovDOfMetric[metric];
		flow = metricConcomitant["FluPerFlow"][metric, opts];
		press = metricConcomitant["Pressure"][metric, opts];
		edens = metricConcomitant["EnergyDensity"][metric, opts];
		time = AbsoluteTime[];
		indicatrix = (flow[a1] TensorDerivative[CTensor[press, {}], cd][-a1]) / (flow[b1] TensorDerivative[CTensor[edens, {}], cd][-b1]);
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"IndicatrixFunction\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		indicatrix = simplf[indicatrix];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"IndicatrixFunction\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		indicatrix
	]
)

metricConcomitant["ThermoFluPerCond"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["ThermoFluPerCond"][metric, opts] = 
	Module[{simplf, cart, cd, indicatrix, press, edens, dindicatrix, dpress, dedens, a1, b1, c1, cond, vb, time},
		{simplf, vb} = OptionValue[metricConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 3];
		cd = CovDOfMetric[metric];
		press = metricConcomitant["Pressure"][metric, opts];
		edens = metricConcomitant["EnergyDensity"][metric, opts];
		indicatrix = metricConcomitant["IndicatrixFunction"][metric, opts];
		time = AbsoluteTime[];
		dpress = TensorDerivative[CTensor[press, {}], cd];
		If[vb, 
			Print["** ReportCompute: computing the exterior derivative of metric concomitant \"Pressure\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dpress = simplf[dpress];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to the exterior derivative of metric concomitant \"Pressure\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dedens = TensorDerivative[CTensor[edens, {}], cd];
		If[vb, 
			Print["** ReportCompute: computing the exterior derivative of metric concomitant \"EnergyDensity\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dedens = simplf[dedens];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to the exterior derivative of metric concomitant \"EnergyDensity\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dindicatrix = TensorDerivative[CTensor[indicatrix, {}], cd];
		If[vb, 
			Print["** ReportCompute: computing the exterior derivative of metric concomitant \"IndicatrixFunction\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dindicatrix = simplf[dindicatrix];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to the exterior derivative of metric concomitant \"IndicatrixFunction\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		cond = Antisymmetrize[dindicatrix[-a1] dpress[-b1] dedens[-c1], {-a1, -b1, -c1}];
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"ThermoFluPerCond\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		cond = simplf[cond];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"ThermoFluPerCond\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		cond
	]
)

metricConcomitant["GenericIdealGasCond"][metric_CTensor, opts : OptionsPattern[]] :=
(metricConcomitant["ThermoFluPerCond"][metric, opts] = 
	Module[{simplf, cart, cd, indicatrix, press, edens, pi, dindicatrix, dpi, a1, b1, cond, vb, time},
		{simplf, vb} = OptionValue[metricConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		cd = CovDOfMetric[metric];
		press = metricConcomitant["Pressure"][metric, opts];
		edens = metricConcomitant["EnergyDensity"][metric, opts];
		indicatrix = metricConcomitant["IndicatrixFunction"][metric, opts];
		time = AbsoluteTime[];
		pi = press / edens;
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"Pi\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		pi = simplf[pi];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"Pressure\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dpi = TensorDerivative[CTensor[pi, {}], cd];
		If[vb, 
			Print["** ReportCompute: computing the exterior derivative of metric concomitant \"Pi\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dpi = simplf[dpi];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to the exterior derivative of metric concomitant \"Pi\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dindicatrix = TensorDerivative[CTensor[indicatrix, {}], cd];
		If[vb, 
			Print["** ReportCompute: computing the exterior derivative of metric concomitant \"IndicatrixFunction\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dindicatrix = simplf[dindicatrix];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to the exterior derivative of metric concomitant \"IndicatrixFunction\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		cond = Antisymmetrize[dindicatrix[-a1] dpi[-b1], {-a1, -b1}];
		If[vb, 
			Print["** ReportCompute: computing metric concomitant \"GenericIdealGasCond\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		cond = simplf[cond];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to metric concomitant \"GenericIdealGasCond\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		cond
	]
)

(* This deletes metric concomitants for all metrics  *)

ClearxIdealCache["MetricConcomitants"] := 
	Module[{},
		SubValues[metricConcomitant] = DeleteCases[SubValues[metricConcomitant], _?(FreeQ[First[#], Pattern] &)];
	]
(* ::Section:: *)
(* Computation of the Weyl concomitants *)

(* TODO: we should avoid defining default options for private functions *)
Options[weylConcomitant] = {PSimplify -> $CVSimplify, Parallelize -> True, Assumptions -> True, Verbose -> True, "Observer" -> Null, "Vector" -> Null, "Bivector" -> Null, "Bivector2" -> Null, Method -> "Default"}

weylConcomitant["Weyl"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["Weyl"][metric, opts] = 
	Module[{cart, cd, parallel, vb},
		cart = 	Part[metric, 2, 1, -1];
		cd = CovDOfMetric[metric];
		{parallel, vb} = OptionValue[weylConcomitant, {opts}, {Parallelize, Verbose}];
		MetricCompute[metric, cart, "Weyl"[-1, -1, -1, -1], Parallelize -> parallel, Verbose -> vb];
		Weyl[cd]
	]
)

weylConcomitant["Weyl2"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["Weyl2"][metric, opts] = 
	Module[{simplf, cart, a1, b1, c1, d1, i1, j1, weylcd, weyl2cd, vb, time},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, i1, j1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];
		weylcd = weylConcomitant["Weyl"][metric, opts];
		time = AbsoluteTime[];
		weyl2cd = 1/2 weylcd[-a1, -b1, i1, j1] weylcd[-i1, -j1, -c1, -d1]; 
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"Weyl2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		weyl2cd = HeadOfTensor[weyl2cd, {-a1, -b1, -c1, -d1}]; 
		time = AbsoluteTime[];
		simplf[weyl2cd];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"Weyl2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		weyl2cd
	]
)

weylConcomitant["Weyl3"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["Weyl3"][metric, opts] = 
	Module[{simplf, cart, a1, b1, c1, d1, i1, j1, weylcd, weyl2cd, weyl3cd, vb, time},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, i1, j1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];
		weylcd = weylConcomitant["Weyl"][metric, opts];
  		weyl2cd = weylConcomitant["Weyl2"][metric, opts];
		time = AbsoluteTime[];
		weyl3cd = 1/2 weyl2cd[-a1, -b1, i1, j1] weylcd[-i1, -j1, -c1, -d1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"Weyl3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		weyl3cd = HeadOfTensor[weyl3cd, {-a1, -b1, -c1, -d1}];
		time = AbsoluteTime[];
		weyl3cd = simplf[weyl3cd];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"Weyl3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		weyl3cd
	]
)

weylConcomitant["TraceWeyl3"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["TraceWeyl3"][metric, opts] = 
	Module[{simplf, cart, a1, b1, weyl3cd, trweyl3, vb, time},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
  		weyl3cd = weylConcomitant["Weyl3"][metric, opts];
		time = AbsoluteTime[];
		trweyl3 = 1/2 weyl3cd[-a1, -b1, a1, b1]; 
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TraceWeyl3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		trweyl3 = simplf[trweyl3];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TraceWeyl3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		trweyl3
	]
)

weylConcomitant["WeylDual"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["WeylDual"][metric, opts] = 
	Module[{simplf, cart, a1, b1, c1, d1, e1, f1, cd, weylcd, epsilonmetric, weyldual, vb, time},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];
		cd = CovDOfMetric[metric];
		epsilonmetric = epsilon[metric];
		weylcd = weylConcomitant["Weyl"][metric, opts];
		time = AbsoluteTime[];
		weyldual = 1/2 epsilonmetric[-c1, -d1, -e1, -f1] weylcd[e1, f1, -a1, -b1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"WeylDual\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		weyldual = ToCCanonical[weyldual];
		weyldual = HeadOfTensor[weyldual, {-c1, -d1, -a1, -b1}];
		time = AbsoluteTime[];
		weyldual = simplf[weyldual];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"WeylDual\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		weyldual
	]
)


weylConcomitant["WeylSelfDual"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["WeylSelfDual"][metric, opts] = 
	Module[{simplf, weylcd, weyldual, weylsd, vb, time},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		weylcd = weylConcomitant["Weyl"][metric, opts];
		weyldual = weylConcomitant["WeylDual"][metric, opts];
		time = AbsoluteTime[];
		weylsd = 1/2 (weylcd - I * weyldual);
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"WeylSelfDual\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		weylsd = simplf[weylsd];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"WeylSelfDual\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		weylsd
	]
)

weylConcomitant["WeylMatrixQ"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["WeylMatrixQ"][metric, opts] = 
	Module[{cart, simplf, obs, mq, a1, b1, c1, d1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 4];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		(* In this particular case we need to input the opts arg. Why? *)
		obs = OptionValue[weylConcomitant, {opts}, "Observer"];
		time = AbsoluteTime[];
		mq = 2 obs[a1] obs[c1] weylConcomitant["WeylSelfDual"][metric, opts][-a1, -b1, -c1, -d1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"WeylMatrixQ\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		mq = HeadOfTensor[mq, {-b1, -d1}];
		time = AbsoluteTime[];
		mq = simplf[mq];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"WeylMatrixQ\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		mq
	]
)

weylConcomitant["TraceWeylMatrixQ"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["TraceWeylMatrixQ"][metric, opts] = 
	Module[{cart, simplf, obs, mq, trmq, a1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 1];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		mq = weylConcomitant["WeylMatrixQ"][metric, opts];
		time = AbsoluteTime[];
		trmq = mq[a1, -a1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TraceWeylMatrixQ\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		trmq = simplf[trmq];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TraceWeylMatrixQ\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		trmq
	]
)

weylConcomitant["WeylMatrixQ2"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["WeylMatrixQ2"][metric, opts] = 
	Module[{cart, simplf, mq, mq2, a1, b1, c1, d1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 3];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		mq = weylConcomitant["WeylMatrixQ"][metric, opts];
		time = AbsoluteTime[];
		mq2 = mq[-a1, -b1] mq[b1, -c1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"WeylMatrixQ2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		mq2 = HeadOfTensor[mq2, {-a1, -c1}];
		time = AbsoluteTime[];
		mq2 = simplf[mq2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"WeylMatrixQ2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		mq2
	]
)

weylConcomitant["TraceWeylMatrixQ2"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["TraceWeylMatrixQ2"][metric, opts] = 
	Module[{cart, simplf, obs, mq2, trmq2, a1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 1];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		mq2 = weylConcomitant["WeylMatrixQ2"][metric, opts];
		time = AbsoluteTime[];
		trmq2 = mq2[a1, -a1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TraceWeylMatrixQ2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		trmq2 = simplf[trmq2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TraceWeylMatrixQ2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		trmq2
	]
)

weylConcomitant["WeylMatrixQ3"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["WeylMatrixQ3"][metric, opts] = 
	Module[{cart, simplf, mq, mq2, mq3, a1, b1, c1, d1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 4];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		mq = weylConcomitant["WeylMatrixQ"][metric, opts];
		mq2 = weylConcomitant["WeylMatrixQ2"][metric, opts];
		time = AbsoluteTime[];
		mq3 = mq2[-a1, -b1] mq[b1, -c1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"WeylMatrixQ3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		mq3 = HeadOfTensor[mq3, {-a1, -c1}];
		time = AbsoluteTime[];
		mq3 = simplf[mq3];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"WeylMatrixQ3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		mq3
	]
)

weylConcomitant["TraceWeylMatrixQ3"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["TraceWeylMatrixQ3"][metric, opts] = 
	Module[{cart, simplf, mq3, trmq3, a1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 1];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		mq3 = weylConcomitant["WeylMatrixQ3"][metric, opts];
		time = AbsoluteTime[];
		trmq3 = mq3[a1, -a1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TraceWeylMatrixQ3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		trmq3 = simplf[trmq3];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TraceWeylMatrixQ3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		trmq3
	]
)

weylConcomitant["WeylSelfDual2"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["WeylSelfDual2"][metric, opts] = 
	Module[{simplf, weylselfdual, weylselfdual2, cart, a1, b1, c1, d1, e1, f1, time, vb},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		time = AbsoluteTime[];
		weylselfdual2 = 1/2 weylselfdual[-a1, -b1, e1, f1] weylselfdual[-e1, -f1, -c1, -d1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"WeylSelfDual2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		weylselfdual2 = HeadOfTensor[weylselfdual2, {-a1, -b1, -c1, -d1}];
		time = AbsoluteTime[];
		weylselfdual2 = simplf[weylselfdual2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"WeylSelfDual2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		weylselfdual2
	]
)

weylConcomitant["TraceWeylSelfDual2"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["TraceWeylSelfDual2"][metric, opts] = 
	Module[{weylselfdual2, trweylselfdual2, simplf, cart, a1, b1, vb, time},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		weylselfdual2 = weylConcomitant["WeylSelfDual2"][metric, opts];
		time = AbsoluteTime[];
		trweylselfdual2 = weylselfdual2[-a1, -b1, a1, b1] / 2;
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TraceWeylSelfDual2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		trweylselfdual2 = simplf[trweylselfdual2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TraceWeylSelfDual2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		trweylselfdual2
	]
)

weylConcomitant["WeylSelfDual3"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["WeylSelfDual3"][metric, opts] = 
	Module[{simplf, weylselfdual, weylselfdual2, weylselfdual3, cart, a1, b1, c1, d1, e1, f1, vb, time},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];
		weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		weylselfdual2 = weylConcomitant["WeylSelfDual2"][metric, opts];
		time = AbsoluteTime[];
		weylselfdual3 = 1/2 weylselfdual2[-a1, -b1, e1, f1] weylselfdual[-e1, -f1, -c1, -d1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"WeylSelfDual3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		weylselfdual3 = HeadOfTensor[weylselfdual3, {-a1, -b1, -c1, -d1}];
		time = AbsoluteTime[];
		weylselfdual3 = simplf[weylselfdual3];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"WeylSelfDual3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		weylselfdual3
	]
)

weylConcomitant["TraceWeylSelfDual3"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["TraceWeylSelfDual3"][metric, opts] = 
	Module[{weylselfdual3, trweylselfdual3, simplf, cart, a1, b1, vb, time},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		weylselfdual3 = weylConcomitant["WeylSelfDual3"][metric, opts];
		time = AbsoluteTime[];
		trweylselfdual3 = weylselfdual3[-a1, -b1, a1, b1] / 2;
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TraceWeylSelfDual3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		trweylselfdual3 = simplf[trweylselfdual3];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TraceWeylSelfDual3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		trweylselfdual3
	]
)

weylConcomitant["ScalarW"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["ScalarW"][metric, opts] = 
	Module[{simplf, cart, weylselfdual, g2form, aa, bb, w, cd, a1, b1, c1, d1},
		simplf = OptionValue[weylConcomitant, PSimplify];
		cart = Part[metric, 2, 1, -1];
		{b1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];
		weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		g2form = metricConcomitant["G2Form"][metric, opts];
		aa = 4 weylConcomitant["TraceWeylSelfDual2"][metric, opts];
		bb = 8 weylConcomitant["TraceWeylSelfDual3"][metric, opts];
		w = CTensor[simplf[-bb / (2 aa)], {}]
	]
)

weylConcomitant["ScalarZ"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["ScalarZ"][metric, opts] = 
	Module[{simplf, cart, weylselfdual, g2form, w, cd, cdw, a1, b1, c1, d1},
		simplf = OptionValue[weylConcomitant, PSimplify];
		cart = Part[metric, 2, 1, -1];
		{b1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 2];		
		cd = CovDOfMetric[metric];
		w = weylConcomitant["ScalarW"][metric, opts];
		cdw = TensorDerivative[w, cd];
		simplf[(metric[b1, d1]) cdw[-b1] cdw[-d1]]
	]
)

weylConcomitant["TensorXi"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["TensorXi"][metric, opts] = 
	Module[{simplf, cart, weylselfdual, g2form, w, cd, a1, b1, c1, d1, cdw, tensorxi},
		simplf = OptionValue[weylConcomitant, PSimplify];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 4];
		weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		g2form = metricConcomitant["G2Form"][metric, opts];
		w = weylConcomitant["ScalarW"][metric, opts];
		cd = CovDOfMetric[metric];
		cdw = TensorDerivative[w, cd];
		tensorxi = (weylselfdual[-a1, -b1, -c1, -d1] - w[] g2form[-a1, -b1, -c1, -d1]) cdw[b1] cdw[d1];
		tensorxi = HeadOfTensor[tensorxi, {-a1, -c1}];
		simplf[tensorxi]
	]
)

weylConcomitant["ConformalTensorXi"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["ConformalTensorXi"][metric, opts] = 
	Module[{simplf, cart, weylselfdual, g2form, w, cd, a1, b1, c1, d1, cdw, schouten, cdschouten, weylcd, weylcdweylcd, lambda, tensorxi},
		simplf = OptionValue[weylConcomitant, PSimplify];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 4];
		weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		g2form = metricConcomitant["G2Form"][metric, opts];
		w = weylConcomitant["ScalarW"][metric, opts];
		cd = CovDOfMetric[metric];
		cdw = TensorDerivative[w, cd];
		schouten = metricConcomitant["SchoutenTensor"][metric, opts];
		cdschouten = TensorDerivative[schouten, cd];
		weylcd = Weyl[cd];
		weylcdweylcd = weylcd[-a1, -b1, -c1, -d1]weylcd[a1, b1, c1, d1];
		lambda = (2 / weylcdweylcd) weylcd[-a1, b1, c1, d1]Antisymmetrize[cdschouten[-c1, -d1, -b1], {-d1, -b1}];
		lambda = simplf[HeadOfTensor[lambda, {-a1}]];
		lambda = simplf[cdw + 2 w lambda];
		tensorxi = (weylselfdual[-a1, -b1, -c1, -d1] - w[] g2form[-a1, -b1, -c1, -d1]) lambda[b1] lambda[d1];
		tensorxi = HeadOfTensor[tensorxi, {-a1, -c1}];
		simplf[tensorxi]
	]
)

weylConcomitant["ConformalScalarZ"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["ConformalScalarZ"][metric, opts] = 
	Module[{simplf, cart, weylselfdual, g2form, w, cd, cdw, a1, b1, c1, d1, schouten, cdschouten, weylcd, weylcdweylcd, lambda},
		simplf = OptionValue[weylConcomitant, PSimplify];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 4];		
		cd = CovDOfMetric[metric];
		w = weylConcomitant["ScalarW"][metric, opts];
		cdw = TensorDerivative[w, cd];
		schouten = metricConcomitant["SchoutenTensor"][metric, opts];
		cdschouten = TensorDerivative[schouten, cd];
		weylcd = Weyl[cd];
		weylcdweylcd = weylcd[-a1, -b1, -c1, -d1]weylcd[a1, b1, c1, d1];
		lambda = (2 / weylcdweylcd) weylcd[-a1, b1, c1, d1]Antisymmetrize[cdschouten[-c1, -d1, -b1], {-d1, -b1}];
		lambda = simplf[HeadOfTensor[lambda, {-a1}]];
		lambda = simplf[cdw + 2 w lambda];
		simplf[(metric[b1, d1]) lambda[-b1] lambda[-d1]]
	]
)

weylConcomitant["ConformalLambda"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["ConformalLambda"][metric, opts] = 
	Module[{simplf, cart, weylselfdual, g2form, w, cd, cdw, a1, b1, c1, d1, schouten, cdschouten, weylcd, weylcdweylcd, lambda},
		simplf = OptionValue[weylConcomitant, PSimplify];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 4];		
		cd = CovDOfMetric[metric];
		w = weylConcomitant["ScalarW"][metric, opts];
		cdw = TensorDerivative[w, cd];
		schouten = metricConcomitant["SchoutenTensor"][metric, opts];
		BreakPoint[];
		cdschouten = TensorDerivative[schouten, cd];
		weylcd = Weyl[cd];
		weylcdweylcd = simplf[weylcd[-a1, -b1, -c1, -d1]weylcd[a1, b1, c1, d1]];
		lambda = (2 / weylcdweylcd) weylcd[-a1, b1, c1, d1]Antisymmetrize[cdschouten[-b1, -c1, -d1], {-d1, -b1}];
		lambda = simplf[HeadOfTensor[lambda, {-a1}]];
		simplf[cdw + 2 w lambda]
	]
)

weylConcomitant["ScalarRho"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["ScalarRho"][metric, opts] = 
	Module[{simplf, cart, trw3, rho, vb, time},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];	
		trw3 = weylConcomitant["TraceWeyl3"][metric, opts];
		time = AbsoluteTime[];
		rho = -(trw3 / 12)^(1/3);
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"ScalarRho\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		rho = simplf[rho];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"ScalarRho\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		rho
	]
)

weylConcomitant["TensorS"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["TensorS"][metric, opts] = 
	Module[{simplf, cart, w, G, rho, s, vb, time},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 4];		
		w = weylConcomitant["Weyl"][metric, opts];
		G = metricConcomitant["G"][metric, opts];
		rho = weylConcomitant["ScalarRho"][metric, opts];
		time = AbsoluteTime[];
		s = (w - rho G) / (3 rho);
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorS\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		s = simplf[s];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorS\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		s
	]
)

weylConcomitant["VectorPhi"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["TensorPhi"][metric, opts] = 
	Module[{simplf, cart, S, cd, cds, a1, b1, c1, d1, phi, vb, time},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 4];		
		S = weylConcomitant["TensorS"][metric, opts];
		cd = CovDOfMetric[metric];
		time = AbsoluteTime[];
		cds = cd[-a1][S[a1, -b1, -c1, -d1]];
		If[vb, 
			Print["** ReportCompute: computing covariant derivative of Weyl concomitant \"TensorS\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		cds = HeadOfTensor[cds, {-b1, -c1, -d1}];
		time = AbsoluteTime[];
		cds = simplf[cds];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to covariant derivative of Weyl concomitant \"TensorS\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		phi = S[a1, -b1, c1, d1]cds[-a1, -c1, -d1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"VectorPhi\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		phi = HeadOfTensor[phi, {-b1}];
		time = AbsoluteTime[];
		phi = simplf[phi];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"VectorPhi\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		phi
	]
)

weylConcomitant["TensorB"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["TensorB"][metric, opts] = 
	Module[{simplf, cart, S, R, r, a1, b1, c1, d1, B, vb, time},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 4];		
		S = weylConcomitant["TensorS"][metric, opts];
		R = weylConcomitant["Ricci"][metric, opts];
		r = weylConcomitant["RicciScalar"][metric, opts];
		time = AbsoluteTime[];
		B = R[-a1, -b1] - 1/2 S[-a1, -b1, -c1, -d1] R[c1, d1] - 1/2 r metric[-a1, -b1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorB\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		B = HeadOfTensor[R[-a1, -b1] - 1/2 S[-a1, -b1, -c1, -d1] R[c1, d1] - 1/2 r metric[-a1, -b1], {-a1, -b1}];
		time = AbsoluteTime[];
		B = simplf[B];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorB\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		B
	]
)

(* Canonical bivectors as Weyl concomitants *)

weylConcomitant["PTNCanonicalBivector"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["PTNCanonicalBivector"][metric, opts] = 
	Module[{simplf, cart, obs, X, weylselfdual, cbv, a1, b1, i1, j1, k1, l1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, i1, j1, k1, l1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];		
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
  		obs = OptionValue[weylConcomitant, {opts}, "Observer"];
    	X = OptionValue[weylConcomitant, {opts}, "Bivector"];
      	weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		time = AbsoluteTime[];
		cbv = weylselfdual[-a1, -b1, -i1, -j1] X[i1, j1] / Sqrt[weylselfdual[-i1, -j1, -k1, -l1] X[i1, j1] X[k1, l1]];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"PTNCanonicalBivector\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		cbv = HeadOfTensor[cbv, {-a1, -b1}];
		time = AbsoluteTime[];
    	cbv = simplf[cbv];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"PTNCanonicalBivector\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		cbv
	]
)

weylConcomitant["PTIIICanonicalBivector1"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["PTIIICanonicalBivector1"][metric, opts] = 
	Module[{simplf, cart, obs, X, weylselfdual2, cbv, a1, b1, i1, j1, k1, l1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, i1, j1, k1, l1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];		
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
  		obs = OptionValue[weylConcomitant, {opts}, "Observer"];
    	X = OptionValue[weylConcomitant, {opts}, "Bivector"];
      	weylselfdual2 = weylConcomitant["WeylSelfDual2"][metric, opts];
		time = AbsoluteTime[];
		cbv = weylselfdual2[-a1, -b1, -i1, -j1] X[i1, j1] / Sqrt[-weylselfdual2[-i1, -j1, -k1, -l1] X[i1, j1] X[k1, l1]];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"PTIIICanonicalBivector1\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		cbv = HeadOfTensor[cbv, {-a1, -b1}];
		time = AbsoluteTime[];
    	cbv = simplf[cbv];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"PTIIICanonicalBivector1\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		cbv
	]
)

weylConcomitant["PTIIICanonicalBivector2"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["PTIIICanonicalBivector2"][metric, opts] = 
	Module[{simplf, cart, obs, X, weylselfdual, g2form, scrh, cbv, a1, b1, i1, j1, k1, l1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, i1, j1, k1, l1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];		
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
  		obs = OptionValue[weylConcomitant, {opts}, "Observer"];
    	X = OptionValue[weylConcomitant, {opts}, "Bivector"];
      	weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		g2form = metricConcomitant["G2Form"][metric, opts];
  		scrh = weylConcomitant["PTIIICanonicalBivector1"][metric, opts];
		time = AbsoluteTime[];
		cbv = (2 1/2 scrh[i1, j1] X[-i1, -j1] 1/2 weylselfdual[-a1, -b1, -i1, -j1] X[i1, j1] - 1/4 weylselfdual[-i1, -j1, -k1, -l1] X[i1, j1] X[k1, l1] scrh[-a1, -b1]) 
       			/ (2 (1/4 scrg[-i1, -j1, -k1, -l1] scrh[i1, j1] X[k1, l1])^2);
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"PTIIICanonicalBivector2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		cbv = HeadOfTensor[cbv, {-a1, -b1}];
		time = AbsoluteTime[];
    	cbv = simplf[cbv];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"PTIIICanonicalBivector2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		cbv
	]
)

weylConcomitant["PTDCanonicalBivector"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["PTDCanonicalBivector"][metric, opts] = 
	Module[{simplf, cart, obs, X, aa, bb, rho, weylselfdual, g2form, scrp, scrp2, cbv, a1, b1, i1, j1, k1, l1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, i1, j1, k1, l1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];		
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
  		obs = OptionValue[weylConcomitant, {opts}, "Observer"];
    	X = OptionValue[weylConcomitant, {opts}, "Bivector"];
      	bb = -weylConcomitant["TraceWeylSelfDual3"][metric, opts];
		aa = weylConcomitant["TraceWeylSelfDual2"][metric, opts];
		rho = bb / aa;
  		g2form = metricConcomitant["G2Form"][metric, opts];
    	weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		time = AbsoluteTime[];
		scrp = weylselfdual - rho g2form;
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorScP\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		scrp = simplf[scrp];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorScP\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		scrp2 = 1/2 scrp[-a1, -b1, -i1, -j1] scrp[i1, j1, -k1, -l1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorScP2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
  		scrp2 = HeadOfTensor[scrp2, {-a1, -b1, -k1, -l1}];
		time = AbsoluteTime[];
		scrp2 = simplf[scrp2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorScP2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		cbv = scrp[-a1, -b1, -i1, -j1] X[i1, j1] / Sqrt[-scrp2[-i1, -j1, -k1, -l1] X[i1, j1] X[k1, l1]];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"PTDCanonicalBivector\" in ", AbsoluteTime[] - time, " seconds:"]
		];
  		cbv = HeadOfTensor[cbv, {-a1, -b1}];
		time = AbsoluteTime[];
    	cbv = simplf[cbv];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"PTDCanonicalBivector\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		cbv
	]
)

weylConcomitant["PTIICanonicalBivector1"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["PTIICanonicalBivector1"][metric, opts] = 
	Module[{simplf, cart, obs, X, aa, bb, rho, weylselfdual, weylselfdual2, g2form, scrq, cbv, a1, b1, i1, j1, k1, l1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, i1, j1, k1, l1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];		
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
  		obs = OptionValue[weylConcomitant, {opts}, "Observer"];
    	X = OptionValue[weylConcomitant, {opts}, "Bivector"];
      	bb = -weylConcomitant["TraceWeylSelfDual3"][metric, opts];
		aa = weylConcomitant["TraceWeylSelfDual2"][metric, opts];
		rho = bb / aa;
  		g2form = metricConcomitant["G2Form"][metric, opts];
    	weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		weylselfdual2 = weylConcomitant["WeylSelfDual2"][metric, opts];
		time = AbsoluteTime[];
		scrq = (weylselfdual2 + rho weylselfdual - 2 rho^2 g2form) / (3 rho);
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorScQ\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		scrq = simplf[scrq];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorScQ\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		cbv = scrq[-a1, -b1, -i1, -j1] X[i1, j1] / Sqrt[scrq[-i1, -j1, -k1, -l1] X[i1, j1] X[k1, l1]];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"PTIICanonicalBivector1\" in ", AbsoluteTime[] - time, " seconds:"]
		];
  		cbv = HeadOfTensor[cbv, {-a1, -b1}];
		time = AbsoluteTime[];
    	cbv = simplf[cbv];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"PTIICanonicalBivector1\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		cbv
	]
)

weylConcomitant["PTIICanonicalBivector2"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["PTIICanonicalBivector2"][metric, opts] = 
	Module[{simplf, cart, obs, X, aa, bb, rho, weylselfdual, g2form, p, scrp, scrp2, cbv, a1, b1, i1, j1, k1, l1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, i1, j1, k1, l1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];		
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
  		obs = OptionValue[weylConcomitant, {opts}, "Observer"];
    	X = OptionValue[weylConcomitant, {opts}, "Bivector"];
      	bb = -weylConcomitant["TraceWeylSelfDual3"][metric, opts];
		aa = weylConcomitant["TraceWeylSelfDual2"][metric, opts];
		rho = bb / aa;
  		g2form = metricConcomitant["G2Form"][metric, opts];
    	weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
		time = AbsoluteTime[];
		p = weylselfdual - rho g2form;
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorP\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		p = simplf[];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorP\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		scrp = 1/2 p[-a1, -b1, -i1, -j1] p[i1, j1, -k1, -l1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorScP\" in ", AbsoluteTime[] - time, " seconds:"]
		];
  		scrp = HeadOfTensor[scrp, {-a1, -b1, -k1, -l1}];
		time = AbsoluteTime[];
		scrp = simplf[scrp];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorScP\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		scrp2 = 1/2 scrp[-a1, -b1, -i1, -j1] scrp[i1, j1, -k1, -l1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorScP2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
  		scrp2 = HeadOfTensor[scrp2, {-a1, -b1, -k1, -l1}];
		time = AbsoluteTime[];
		scrp2 = simplf[scrp2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorScP2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		cbv = scrp[-a1, -b1, -i1, -j1] X[i1, j1] / Sqrt[-scrp2[-i1, -j1, -k1, -l1] X[i1, j1] X[k1, l1]];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"PTIICanonicalBivector2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
  		cbv = HeadOfTensor[cbv, {-a1, -b1}];
		time = AbsoluteTime[];
		cbv = simplf[cbv];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"PTIICanonicalBivector2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
    	cbv
	]
)

(* Connection tensor for different Petrov Types as Weyl concomitants *)

weylConcomitant["PTIConnectionTensor"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["PTIConnectionTensor"][metric, opts] = 
	Module[{cart, simplf, aa, bb, selfdualW, selfdualW2, G2form, cd, scX, scY, selfdualH, H, a1, b1, c1, d1, e1, f1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		selfdualW = weylConcomitant["WeylSelfDual"][metric, opts];
		selfdualW2 = weylConcomitant["WeylSelfDual2"][metric, opts];
		aa = weylConcomitant["TraceWeylSelfDual2"][metric, opts];
		bb = weylConcomitant["TraceWeylSelfDual3"][metric, opts];
		G2form = metricConcomitant["G2Form"];
		cd = CovDOfMetric[metric];
		time = AbsoluteTime[];
		scX = cd[-a1][selfdualW[-b1, -c1, -d1, -e1]] selfdualW[d1, e1, c1, -f1] / 2;
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"scXtensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		scX = HeadOfTensor[scX, {-a1, -b1, -f1}];
		time = AbsoluteTime[];
		scX = simplf[scX];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"scXtensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		scY = (3 aa selfdualW2 + 6 bb selfdualW + aa^2 G2form / 2) / (aa^3 - 6 bb^2);
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"scYtensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		scY = simplf[scY];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"scYtensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		selfdualH = scX[-a1, -b1, -c1] scY[b1, c1, -d1, -e1] / Sqrt[2];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"SelfDualPTIConnectionTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		selfdualH = HeadOfTensor[selfdualH, {-a1, -d1, -e1}];
		time = AbsoluteTime[];
		selfdualH = simplf[selfdualH];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"SelfDualPTIConnectionTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		H = (selfdualH + Dagger[selfdualH]) / Sqrt[2];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"PTIConnectionTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		H = simplf[H];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"PTIConnectionTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		H
	]
)

weylConcomitant["PTIIConnectionTensor"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["PTIIConnectionTensor"][metric, opts] = 
	Module[{cart, simplf, cd, scU, scLp, scS, scLm, selfdualH, H, a1, b1, c1, d1, e1, f1, g1, h1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 8];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		scU = weylConcomitant["PTIICanonicalBivector2"];
		scLp = weylConcomitant["PTIICanonicalBivector1"];
		G2form = metricConcomitant["G2Form"];
		cd = CovDOfMetric[metric];
		X = OptionValue[weylConcomitant, {opts}, "Bivector"];
		time = AbsoluteTime[];
		scS = G2form[-a1, -b1, -c1, -d1] + scU[-a1, -b1] scU[-c1, -d1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"scStensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		scS = HeadOfTensor[scS, {-a1, -b1, -c1, -d1}];
		time = AbsoluteTime[];
		scS = simplf[scS];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"scStensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		scLm = (scS[-a1, -b1, -c1, -d1] X[a1, b1] X[c1, d1] scLp[-e1, -f1, -g1, -h1] - 2 scLp[-e1, -f1, -a1, -b1] X[a1, b1] scS[-g1, -h1, -c1, -d1] X[c1, d1]) / (2 scLp[-a1, -b1, -c1, -d1] X[c1, d1] scLp[a1, b1, -e1, -f1] X[e1, f1]);
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"PTIINullEigenBivectorm\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		scLm = HeadOfTensor[scLm, {-e1, -f1, -g1, -h1}];
		time = AbsoluteTime[];
		scLm = simplf[scLm];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"PTIINullEigenBivectorm\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		selfdualH = - (cd[-a1][scU[-b1, -c1]] scU[c1, -d1] + cd[-a1][scLp[-b1, -c1]] scLm[c1, -d1] + cd[-a1][scLm[-b1, -c1]] scLp[c1, -d1]) / Sqrt[2];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"SelfDualPTIIConnectionTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		selfdualH = HeadOfTensor[selfdualH, {-a1, -b1, -d1}];
		time = AbsoluteTime[];
		selfdualH = simplf[selfdualH];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"SelfDualPTIIConnectionTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		H = (selfdualH + Dagger[selfdualH]) / Sqrt[2];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"PTIIConnectionTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		H = simplf[H];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"PTIIConnectionTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		H
	]
)

weylConcomitant["PTIIIConnectionTensor"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["PTIIIConnectionTensor"][metric, opts] = 
	Module[{cart, simplf, cd, scU, scLp, scS, scLm, selfdualH, H, a1, b1, c1, d1, e1, f1, g1, h1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 8];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		scU = weylConcomitant["PTIIICanonicalBivector2"];
		scLp = weylConcomitant["PTIIICanonicalBivector1"];
		G2form = metricConcomitant["G2Form"];
		cd = CovDOfMetric[metric];
		X = OptionValue[weylConcomitant, {opts}, "Bivector"];
		time = AbsoluteTime[];
		scS = G2form[-a1, -b1, -c1, -d1] + scU[-a1, -b1] scU[-c1, -d1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"scStensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		scS = HeadOfTensor[scS, {-a1, -b1, -c1, -d1}];
		time = AbsoluteTime[];
		scS = simplf[scS];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"scStensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		scLm = (scS[-a1, -b1, -c1, -d1] X[a1, b1] X[c1, d1] scLp[-e1, -f1, -g1, -h1] - 2 scLp[-e1, -f1, -a1, -b1] X[a1, b1] scS[-g1, -h1, -c1, -d1] X[c1, d1]) / (2 scLp[-a1, -b1, -c1, -d1] X[c1, d1] scLp[a1, b1, -e1, -f1] X[e1, f1]);
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"PTIIINullEigenBivectorm\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		scLm = HeadOfTensor[scLm, {-e1, -f1, -g1, -h1}];
		time = AbsoluteTime[];
		scLm = simplf[scLm];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"PTIIINullEigenBivectorm\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		selfdualH = - (cd[-a1][scU[-b1, -c1]] scU[c1, -d1] + cd[-a1][scLp[-b1, -c1]] scLm[c1, -d1] + cd[-a1][scLm[-b1, -c1]] scLp[c1, -d1]) / Sqrt[2];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"SelfDualPTIIIConnectionTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		selfdualH = HeadOfTensor[selfdualH, {-a1, -b1, -d1}];
		time = AbsoluteTime[];
		selfdualH = simplf[selfdualH];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"SelfDualPTIIIConnectionTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		H = (selfdualH + Dagger[selfdualH]) / Sqrt[2];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"PTIIIConnectionTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		H = simplf[H];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"PTIIIConnectionTensor\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		H
	]
)

(* Debever directions as Weyl concomitants *)

weylConcomitant["NullDirectionTypeN"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["NullDirectionTypeN"][metric, opts] = 
	Module[{cart, simplf, obs, mq, dir, a1, b1, c1, d1, e1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 5];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		obs = OptionValue[weylConcomitant, {opts}, "Observer"];
		mq = weylConcomitant["WeylMatrixQ"][metric, opts];
		time = AbsoluteTime[];
		dir = Dagger[mq[a1, b1]] (mq[-b1, -a1] obs[c1] + I mq[-b1, d1] epsilon[metric][-a1, -d1, c1, -e1] obs[e1]);
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"NullDirectionTypeN\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		dir = HeadOfTensor[dir, {c1}];
		time = AbsoluteTime[];
		dir = simplf[dir];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"NullDirectionTypeN\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		dir
	]
)

weylConcomitant["WNullDirectionTypeN"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["WNullDirectionTypeN"][metric, opts] = 
	Module[{cart, simplf, obs, weylselfdual, canonicalbivector, mh, mh2, dir, a1, b1, i1, j1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, i1, j1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 4];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
  		obs = OptionValue[weylConcomitant, {opts}, "Observer"];
  		canonicalbivector = weylConcomitant["PTNCanonicalBivector"][metric, opts];
		time = AbsoluteTime[];
	 	mh = ComplexExpand[(canonicalbivector + Dagger[canonicalbivector])/Sqrt[2]];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorH\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		mh2 = mh[-a1, -i1] mh[i1, -b1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorH2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
   		mh2 = HeadOfTensor[mh2, {-a1, -b1}];
		time = AbsoluteTime[];
		mh2 = simplf[mh2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorH2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dir = mh2[-a1, -i1] obs[i1] / Sqrt[-mh2[i1, -j1] obs[i1] obs[j1]];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"WNullDirectionTypeN\" in ", AbsoluteTime[] - time, " seconds:"]
		];
     	dir = HeadOfTensor[dir, {-a1}];
		time = AbsoluteTime[];
		dir = simplf[dir];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"WNullDirectionTypeN\" in ", AbsoluteTime[] - time, " seconds:"]
		];
       	dir
	]
)

weylConcomitant["NullDirectionTypeIII"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["NullDirectionTypeIII"][metric, opts] = 
	Module[{cart, simplf, obs, mq, a1, b1, c1, d1, e1, dir, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 5];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		obs = OptionValue[weylConcomitant, {opts}, "Observer"];
		mq = weylConcomitant["WeylMatrixQ2"][metric, opts];
		time = AbsoluteTime[];
		dir = Dagger[mq[a1, b1]] (mq[-b1, -a1] obs[c1] + I mq[-b1, d1] epsilon[metric][-a1, -d1, c1, -e1] obs[e1]);
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"NullDirectionTypeIII\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		dir = HeadOfTensor[dir, {c1}];
		time = AbsoluteTime[];
		dir = simplf[dir];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"NullDirectionTypeIII\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		dir
	]
)

weylConcomitant["WNullDirectionTypeIII"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["WNullDirectionTypeIII"][metric, opts] = 
	Module[{cart, simplf, obs, weylselfdual, canonicalbivector, mh, mh2, dir, a1, b1, i1, j1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, i1, j1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 4];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
  		obs = OptionValue[weylConcomitant, {opts}, "Observer"];
  		canonicalbivector = weylConcomitant["PTIIICanonicalBivector1"][metric, opts];
		time = AbsoluteTime[];
	 	mh = ComplexExpand[(canonicalbivector + Dagger[canonicalbivector])/Sqrt[2]];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorH\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		mh2 = mh[-a1, -i1] mh[i1, -b1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorH2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
   		mh2 = HeadOfTensor[mh2, {-a1, -b1}];
		time = AbsoluteTime[];
		mh2 = simplf[mh2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorH2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dir = mh2[-a1, -i1] obs[i1] / Sqrt[-mh2[i1, -j1] obs[i1] obs[j1]];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"WNullDirectionTypeIII\" in ", AbsoluteTime[] - time, " seconds:"]
		];
     	dir = HeadOfTensor[dir, {-a1}];
		time = AbsoluteTime[];
       	dir = simplf[dir];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"WNullDirectionTypeIII\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		dir
	]
)

(* TODO: when doing simplifications we are not able to handle outputs with Piecewise *)
weylConcomitant["NullDirectionTypeD"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["NullDirectionTypeD"][metric, opts] = 
	Module[{cart, simplf, obs, w, a1, b1, c1, d1, e1, bb, aa, mq, rho, gamma, P, Pdag, scrP, scrP2, S, dseda, Ch2, v0, v1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 5];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		obs = OptionValue[weylConcomitant, {opts}, "Observer"];
		w = OptionValue[weylConcomitant, {opts}, "Vector"];
		bb = -weylConcomitant["TraceWeylMatrixQ3"][metric, opts];
		aa = weylConcomitant["TraceWeylMatrixQ2"][metric, opts];
		mq = weylConcomitant["WeylMatrixQ"][metric, opts];
		gamma = metricConcomitant["SpatialMetric"][metric, opts];
		time = AbsoluteTime[];
		rho = -bb / aa;
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"ScalarRho\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		rho = simplf[rho];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"ScalarRho\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		P = 1 / (3 rho) mq;
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorP\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		P = simplf[P];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorP\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		Pdag = Dagger[P];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorPdag\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		Pdag = simplf[Pdag];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorPdag\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		scrP = P[-a1, -b1] Pdag[b1, -c1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorScrP\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		scrP = HeadOfTensor[scrP, {-a1, -c1}];
		time = AbsoluteTime[];
		scrP = simplf[scrP];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorScrP\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		scrP2 = Pdag[-a1, -b1] P[b1, -c1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorScrP2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		scrP2 = HeadOfTensor[scrP2, {-a1, -c1}];
		time = AbsoluteTime[];
		scrP2 = simplf[scrP2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorScrP2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dseda = scrP[-a1, a1] + 1/3;
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"ScalarDseda\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dseda = simplf[dseda];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"ScalarDseda\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		S = 1/4 (1 + 2 / (3 Sqrt[dseda])) (P + Pdag) + 1 / (4 Sqrt[dseda]) (scrP + scrP2) + 1/6 (1 + 1 / (3 Sqrt[dseda])) gamma;
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorS\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		S = simplf[S];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorS\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		Ch2 = (1 + Sqrt[dseda]) / 2;
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"Cosh2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		Ch2 = simplf[Ch2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"Cosh2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		v0 = Ch2 obs[a1] + I / (2 Sqrt[dseda]) scrP2[c1, b1] epsilon[metric][-c1, -b1, a1, -d1] obs[d1];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"v0\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		v0 = PowerExpand[v0];
		v0 = HeadOfTensor[v0, {a1}];
		time = AbsoluteTime[];
		v0 = simplf[v0];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"v0\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		v1 = S[-c1, a1] w[c1] / Sqrt[S[-d1, -b1] w[d1] w[b1]];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"v1\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		v1 = PowerExpand[v1];
		v1 = HeadOfTensor[v1, {a1}];
		time = AbsoluteTime[];
		v1 = simplf[v1];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"v1\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dir1 = v0 + v1;
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"NullDirectionTypeDPlus\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dir1 = simplf[dir1];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"NullDirectionTypeDPlus\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		dir2 = v0 - v1;
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"NullDirectionTypeDMinus\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dir2 = simplf[dir2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"NullDirectionTypeDMinus\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		{dir1, dir2}
	]
)

weylConcomitant["WNullDirectionTypeD"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["WNullDirectionTypeD"][metric, opts] = 
	Module[{cart, simplf, obs, canonicalbivector, mu, lp, lm, a1, b1, i1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, i1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 3];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
  		obs = OptionValue[weylConcomitant, {opts}, "Observer"];
  		canonicalbivector = weylConcomitant["PTDCanonicalBivector"][metric, opts];
	 	mu = ComplexExpand[(canonicalbivector + Dagger[canonicalbivector])/Sqrt[2]];
		time = AbsoluteTime[];
     	lp = HeadOfTensor[(mu[-b1, i1] mu[-i1, a1] + mu[-b1, a1]) obs[b1], {a1}];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"WNullDirectionTypeDPlus\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		lp = simplf[lp];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"WNullDirectionTypeDPlus\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
       	lm = HeadOfTensor[(mu[-b1, i1] mu[-i1, a1] - mu[-b1, a1]) obs[b1], {a1}];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"WNullDirectionTypeDMinus\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		lm = simplf[lm];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"WNullDirectionTypeDMinus\" in ", AbsoluteTime[] - time, " seconds:"]
		];
	 	{lp, lm}
	]
)


weylConcomitant["NullDirectionTypeII"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["NullDirectionTypeII"][metric, opts] = 
	Module[{cart, simplf, obs, mq, mq2, a1, b1, c1, d1, e1, bb, aa, gamma, rho, P, dvb, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 5];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		obs = OptionValue[weylConcomitant, {opts}, "Observer"];
		(* In the following calls the "Observer" option is supossed to be non-Null *)
		mq = weylConcomitant["WeylMatrixQ"][metric, opts];
		mq2 = weylConcomitant["WeylMatrixQ2"][metric, opts];
		bb = weylConcomitant["TraceWeylMatrixQ3"][metric, opts];
		aa = weylConcomitant["TraceWeylMatrixQ2"][metric, opts];
		gamma = metricConcomitant["SpatialMetric"][metric, opts];
		time = AbsoluteTime[];
		rho = bb / aa;
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"ScalarRho\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		rho = simplf[rho];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"ScalarRho\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		P = rho mq + 2 rho^2 gamma - mq2;
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorP\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		P = simplf[P];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorP\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dvb = Dagger[P[a1, b1]] (P[-b1, -a1] obs[c1] + I P[-b1, d1] epsilon[metric][-a1, -d1, c1, -e1] obs[e1]);
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"NullDirectionTypeII\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		dvb = HeadOfTensor[dvb, {c1}];
		time = AbsoluteTime[];
		dvb = simplf[dvb];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"NullDirectionTypeII\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		dvb
	]
)

weylConcomitant["WNullDirectionTypeII"][metric_CTensor, opts : OptionsPattern[]] :=
(weylConcomitant["WNullDirectionTypeII"][metric, opts] = 
	Module[{cart, simplf, obs, canonicalbivector, mh, mh2, dir, a1, b1, i1, j1, k1, l1, vb, time},
		cart = Part[metric, 2, 1, -1];
		{a1, b1, i1, j1, k1, l1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
  		obs = OptionValue[weylConcomitant, {opts}, "Observer"];
  		canonicalbivector = weylConcomitant["PTIICanonicalBivector1"][metric, opts];
		time = AbsoluteTime[];
	 	mh = ComplexExpand[(canonicalbivector + Dagger[canonicalbivector])/Sqrt[2]];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorH\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		mh2 = HeadOfTensor[mh[-a1, -i1] mh[i1, -b1], {-a1, -b1}];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"TensorH2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
   		mh2 = simplf[mh2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"TensorH2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		time = AbsoluteTime[];
		dir = mh2[-a1, -i1] obs[i1] / Sqrt[-mh2[-i1, -j1] obs[i1] obs[j1]];
		If[vb, 
			Print["** ReportCompute: computing Weyl concomitant \"WNullDirectionTypeII\" in ", AbsoluteTime[] - time, " seconds:"]
		];
     	dir = HeadOfTensor[dir, {-a1}];
		time = AbsoluteTime[];
       	dir = simplf[dir];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Weyl concomitant \"WNullDirectionTypeII\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		dir
	]
)

(* This deletes the computed Weyl concomitants for all metrics  *)

ClearxIdealCache["WeylConcomitants"] := 
	Module[{},
		SubValues[weylConcomitant] = DeleteCases[SubValues[weylConcomitant], _?(FreeQ[First[#], Pattern] &)];
	]


(* ::Section:: *)
(* Computation of the R-frame concomitants *)

(*
NOTE: I put H as an input in the functions that compute its derivative concomitants
instead of calling ConnectionTensorConcomitant["ConnectionTensor"] in case it is not computed from the R-frame.
*)

(* TODO: Add different options to compute H to ConnectionTensorConcomitant["ConnectionTensor"] *)

Options[ConnectionTensorConcomitant] = {PSimplify -> $CVSimplify, Parallelize -> True, Verbose -> True}

ConnectionTensorConcomitant["C1"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C1"][metric, H, opts] = 
	Module[{simplf, cart, cd, a1, b1, c1, d1, i1, vb, time, ce1},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, i1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 5];
		cd = CovDOfMetric[metric];
		time = AbsoluteTime[];
		ce1 = cd[-a1][H[-b1, -c1, -d1]] + H[-a1, -b1, i1] H[-i1, -c1, -d1] + H[-a1, -c1, i1] H[-b1, -i1, -d1] + 
			H[-a1, -d1, i1] H[-b1, -c1, -i1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C1\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1 = HeadOfTensor[ce1, {-a1, -b1, -c1, -d1}];
		time = AbsoluteTime[];
		ce1 = simplf[ce1];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C1\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1
  	]
)

ConnectionTensorConcomitant["C11"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C11"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, vb, time, ce11},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 10];
  		epsilonmetric = epsilon[metric];
    	ce1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
		time = AbsoluteTime[];
		ce11 = epsilonmetric[i1, j1, a1, b1] ce1[-i1, -c1, -d1, -e1] ce1[-j1, -f1, -g1, -h1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C11\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce11 = HeadOfTensor[ce11, {a1, b1, -c1, -d1, -e1, -f1, -g1, -h1}];
		time = AbsoluteTime[];
		ce11 = simplf[ce11];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C11\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce11
    ]
)

ConnectionTensorConcomitant["C2"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C2"][metric, H, opts] = 
	Module[{simplf, cart, cd, ce1, a1, b1, c1, d1, e1, i1, vb, time, ce2},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, i1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 6];
  		cd = CovDOfMetric[metric];
    	ce1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
		time = AbsoluteTime[];
		ce2 = cd[-a1][ce1[-b1, -c1, -d1, -e1]] + H[-a1, -b1, i1] ce1[-i1, -c1, -d1, -e1] +
	   		H[-a1, -c1, i1] ce1[-b1, -i1, -d1, -e1] + H[-a1, -d1, i1] ce1[-b1, -c1, -i1, -e1] + 
         	H[-a1, -e1, i1] ce1[-b1, -c1, -d1, -i1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce2 = HeadOfTensor[ce2, {-a1, -b1, -c1, -d1, -e1}];
		time = AbsoluteTime[];
		ce2 = simplf[ce2];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C2\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce2
  	]
)

ConnectionTensorConcomitant["C12"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C12"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, vb, time, ce12},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 11];
  		epsilonmetric = epsilon[metric];
    	ce1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
      	ce2 = ConnectionTensorConcomitant["C2"][metric, H, opts];
		time = AbsoluteTime[];
		ce12 = epsilonmetric[i1, j1, a1, b1] ce1[-i1, -c1, -d1, -e1] ce2[-j1, -f1, -g1, -h1, -k1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C12\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce12 = HeadOfTensor[ce12, {a1, b1, -c1, -d1, -e1, -f1, -g1, -h1, -k1}];
		time = AbsoluteTime[];
		ce12 = simplf[ce12];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C12\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce12
  	]
)

ConnectionTensorConcomitant["C122"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C122"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, vb, time, ce122},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 15];
  		epsilonmetric = epsilon[metric];
    	ce1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
      	ce2 = ConnectionTensorConcomitant["C2"][metric, H, opts];
		ce122 = epsilonmetric[i1, j1, k1, a1] ce1[-i1, -b1, -c1, -d1] ce2[-j1, -e1, -f1, -g1, -h1] ce2[-k1, -l1, -m1, -n1, -o1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C122\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce122 = HeadOfTensor[ce122, {a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -l1, -m1, -n1, -o1}];
		time = AbsoluteTime[];
		ce122 = simplf[ce122];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C122\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce122
  	]
)

ConnectionTensorConcomitant["C3"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C3"][metric, H, opts] = 
	Module[{simplf, cart, cd, ce2, a1, b1, c1, d1, e1, f1, i1, vb, time, ce3},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, i1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 7];
  		cd = CovDOfMetric[metric];
    	ce2 = ConnectionTensorConcomitant["C2"][metric, H, opts];
		time = AbsoluteTime[];
		ce3 = cd[-a1][ce2[-b1, -c1, -d1, -e1, -f1]] + H[-a1, -b1, i1] ce2[-i1, -c1, -d1, -e1, -f1] + 
            H[-a1, -c1, i1] ce2[-b1, -i1, -d1, -e1, -f1] + H[-a1, -d1, i1] ce2[-b1, -c1, -i1, -e1, -f1] + 
          	H[-a1, -e1, i1] ce2[-b1, -c1, -d1, -i1, -f1] + H[-a1, -f1, i1] ce2[-b1, -c1, -d1, -e1, -i1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce3 = HeadOfTensor[ce3, {-a1, -b1, -c1, -d1, -e1, -f1}];
		time = AbsoluteTime[];
		ce3 = simplf[ce3];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C3\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce3
  	]
)

ConnectionTensorConcomitant["C123"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C123"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, ce3, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, vb, time, ce123},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 16];
  		epsilonmetric = epsilon[metric];
    	ce1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
     	ce2 = ConnectionTensorConcomitant["C2"][metric, H, opts];
    	ce3 = ConnectionTensorConcomitant["C3"][metric, H, opts];
		time = AbsoluteTime[];
		ce123 = epsilonmetric[i1, j1, k1, a1] ce1[-i1, -b1, -c1, -d1] ce2[-j1, -e1, -f1, -g1, -h1] ce3[-k1, -l1, -m1, -n1, -o1, -p1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C123\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce123 = HeadOfTensor[ce123, {a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -l1, -m1, -n1, -o1, -p1}];
		time = AbsoluteTime[];
		ce123 = simplf[ce123];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C123\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce123
  	]
)

ConnectionTensorConcomitant["C1233"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C1233"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, ce3, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, 
 		q1, r1, s1, t1, u1, vb, time, ce1233},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1, s1, t1, u1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 21];
  		epsilonmetric = epsilon[metric];
    	ce1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
     	ce2 = ConnectionTensorConcomitant["C2"][metric, H, opts];
    	ce3 = ConnectionTensorConcomitant["C3"][metric, H, opts];
		time = AbsoluteTime[];
		ce1233 = 
			epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce2[-j1, -d1, -e1, -f1, -g1] ce3[-k1, -h1, -m1, -n1, -o1, -p1] ce3[-l1, -q1, -r1, -s1, -t1, -u1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C1233\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1233 = HeadOfTensor[ce1233, {-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1, -q1, -r1, -s1, -t1, -u1}];
		time = AbsoluteTime[];
		ce1233 = simplf[ce1233];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C1233\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1233
  	]
)

ConnectionTensorConcomitant["C4"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C4"][metric, H, opts] = 
	Module[{simplf, cart, cd, ce3, a1, b1, c1, d1, e1, f1, g1, i1, vb, time, ce4},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, i1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 8];
  		cd = CovDOfMetric[metric];
    	ce3 = ConnectionTensorConcomitant["C3"][metric, H, opts];
		time = AbsoluteTime[];
		ce4 = cd[-a1][ce3[-b1, -c1, -d1, -e1, -f1, -g1]] + H[-a1, -b1, i1] ce3[-i1, -c1, -d1, -e1, -f1, -g1] + 
          	H[-a1, -c1, i1] ce3[-b1, -i1, -d1, -e1, -f1, -g1] + H[-a1, -d1, i1] ce3[-b1, -c1, -i1, -e1, -f1, -g1] + 
           	H[-a1, -e1, i1] ce3[-b1, -c1, -d1, -i1, -f1, -g1] + H[-a1, -f1, i1] ce3[-b1, -c1, -d1, -e1, -i1, -g1] + 
          	H[-a1, -g1, i1] ce3[-b1, -c1, -d1, -e1, -f1, -i1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C4\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce4 = HeadOfTensor[ce4, {-a1, -b1, -c1, -d1, -e1, -f1, -g1}];
		time = AbsoluteTime[];
		ce4 = simplf[ce4];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C4\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce4
  	]
)

ConnectionTensorConcomitant["C1234"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C1234"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, ce3, ce4, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, 
 		q1, r1, s1, t1, u1, v1, vb, time, ce1234},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1, s1, t1, u1, v1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 22];
  		epsilonmetric = epsilon[metric];
    	ce1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
     	ce2 = ConnectionTensorConcomitant["C2"][metric, H, opts];
    	ce3 = ConnectionTensorConcomitant["C3"][metric, H, opts];
      	ce4 = ConnectionTensorConcomitant["C4"][metric, H, opts];
		time = AbsoluteTime[];
		ce1234 = epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce2[-j1, -d1, -e1, -f1, -g1] ce3[-k1, -h1, -m1, -n1, -o1, -p1] ce4[-l1, -q1, -r1, -s1, -t1, -u1, -v1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C1234\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1234 = HeadOfTensor[ce1234, {-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1, -q1, -r1, -s1, -t1, -u1, -v1}];
		time = AbsoluteTime[];
		ce1234 = simplf[ce1234];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C1234\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1234
  	]
)

ConnectionTensorConcomitant["C1222"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C1222"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, 
 		q1, r1, s1, vb, time, ce1222},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1, s1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 19];
  		epsilonmetric = epsilon[metric];
    	ce1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
     	ce2 = ConnectionTensorConcomitant["C2"][metric, H, opts];
		time = AbsoluteTime[];
		ce1222 = epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce2[-j1, -d1, -e1, -f1, -g1] ce2[-k1, -h1, -m1, -n1, -o1] ce2[-l1, -p1, -q1, -r1, -s1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C1222\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1222 = HeadOfTensor[ce1222, {-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1, -q1, -r1, -s1}];
		time = AbsoluteTime[];
		ce1222 = simplf[ce1222];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C1222\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1222
  	]
)

ConnectionTensorConcomitant["C1223"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C1223"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, ce3, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, 
 		q1, r1, s1, t1, vb, time, ce1223},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1, s1, t1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 20];
  		epsilonmetric = epsilon[metric];
    	ce1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
     	ce2 = ConnectionTensorConcomitant["C2"][metric, H, opts];
       	ce3 = ConnectionTensorConcomitant["C3"][metric, H, opts];
		time = AbsoluteTime[];
		ce1223 = epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce2[-j1, -d1, -e1, -f1, -g1] ce2[-k1, -h1, -m1, -n1, -o1] ce3[-l1, -p1, -q1, -r1, -s1, -t1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C1223\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1223 = HeadOfTensor[ce1223, {-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1, -q1, -r1, -s1, -t1}];
		time = AbsoluteTime[];
		ce1223 = simplf[];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C1223\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1223
  	]
)

ConnectionTensorConcomitant["C111"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C111"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, m1, n1, vb, time, ce111},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, m1, n1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 13];
  		epsilonmetric = epsilon[metric];
    	ce1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
		time = AbsoluteTime[];
		ce111 = epsilonmetric[i1, j1, k1, a1] ce1[-i1, -b1, -c1, -d1] ce1[-j1, -e1, -f1, -g1] ce1[-k1, -h1, -m1, -n1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C111\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce111 = HeadOfTensor[ce111, {a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1}];
		time = AbsoluteTime[];
		ce111 = simplf[ce111];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C111\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce111
  	]
)

ConnectionTensorConcomitant["C112"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C112"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, m1, n1, o1, vb, time, ce112},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, m1, n1, o1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 14];
  		epsilonmetric = epsilon[metric];
    	ce1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
      	ce2 = ConnectionTensorConcomitant["C2"][metric, H, opts];
		time = AbsoluteTime[];
		ce112 = epsilonmetric[i1, j1, k1, a1] ce1[-i1, -b1, -c1, -d1] ce1[-j1, -e1, -f1, -g1] ce2[-k1, -h1, -m1, -n1, -o1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C112\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce112 = HeadOfTensor[ce112, {a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1}];
		time = AbsoluteTime[];
		ce112 = simplf[ce112];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C112\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce112
  	]
)

ConnectionTensorConcomitant["C1122"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C1122"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, 
 		q1, r1, vb, time, ce1122},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 18];
  		epsilonmetric = epsilon[metric];
    	ce1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
     	ce2 = ConnectionTensorConcomitant["C2"][metric, H, opts];
		time = AbsoluteTime[];
		ce1122 = epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce1[-j1, -d1, -e1, -f1] ce2[-k1, -g1, -h1, -m1, -n1] ce2[-l1, -o1, -p1, -q1, -r1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C1122\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1122 = HeadOfTensor[ce1122, {-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1, -q1, -r1}];
		time = AbsoluteTime[];
		ce1122 = simplf[ce1122];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C1122\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1122
  	]
)

ConnectionTensorConcomitant["C1123"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C1123"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, ce3, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, 
 		q1, r1, s1, vb, time, ce1123},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1,s1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 19];
  		epsilonmetric = epsilon[metric];
    	ce1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
     	ce2 = ConnectionTensorConcomitant["C2"][metric, H, opts];
       	ce3 = ConnectionTensorConcomitant["C3"][metric, H, opts];
		time = AbsoluteTime[];
		ce1123 = epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce1[-j1, -d1, -e1, -f1] ce2[-k1, -g1, -h1, -m1, -n1] ce3[-l1, -o1, -p1, -q1, -r1, -s1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C1123\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1123 = HeadOfTensor[ce1123, {-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1, -q1, -r1, -s1}];
		time = AbsoluteTime[];
		ce1123 = simplf[ce1123];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C1123\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1123
  	]
)

ConnectionTensorConcomitant["C1111"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C1111"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, vb, time, ce1111},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 16];
  		epsilonmetric = epsilon[metric];
    	ce1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
		time = AbsoluteTime[];
		ce1111 = epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce1[-j1, -d1, -e1, -f1] ce1[-k1, -g1, -h1, -m1] ce1[-l1, -n1, -o1, -p1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C1111\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1111 = HeadOfTensor[ce1111, {-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1}];
		time = AbsoluteTime[];
		ce1111 = simplf[ce1111];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C1111\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1111
  	]
)

ConnectionTensorConcomitant["C1112"][metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
(ConnectionTensorConcomitant["C1112"][metric, H, opts] = 
	Module[{simplf, cart, epsilonmetric, ce1, ce2, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, vb, time, ce1112},
		{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
		cart = Part[metric, 2, 1, -1];
		{a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1} = 
  			GetIndicesOfVBundle[VBundleOfBasis @ cart, 17];
  		epsilonmetric = epsilon[metric];
    	ce1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
      	ce2 = ConnectionTensorConcomitant["C2"][metric, H, opts];
		time = AbsoluteTime[];
		ce1112 = epsilonmetric[i1, j1, k1, l1] ce1[-i1, -a1, -b1, -c1] ce1[-j1, -d1, -e1, -f1] ce1[-k1, -g1, -h1, -m1] ce2[-l1, -n1, -o1, -p1, -q1];
		If[vb, 
			Print["** ReportCompute: computing Connection Tensor concomitant \"C1112\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1112 = HeadOfTensor[ce1112, {-a1, -b1, -c1, -d1, -e1, -f1, -g1, -h1, -m1, -n1, -o1, -p1, -q1}];
		time = AbsoluteTime[];
		ce1112 = simplf[ce1112];
		If[vb,
			Print["** ReportCompute: applying  ", simplf, " to Connection Tensor concomitant \"C1112\" in ", AbsoluteTime[] - time, " seconds:"]
		];
		ce1112
  	]
)

(* This deletes ConnectionTensor concomitants for all metrics  *)

ClearxIdealCache["ConnectionTensorConcomitants"] := 
	Module[{},
		SubValues[ConnectionTensorConcomitant] = DeleteCases[SubValues[ConnectionTensorConcomitant], _?(FreeQ[First[#], Pattern] &)];
	]

(* ::Section:: *)
(* Clear cache *)

ClearxIdealCache[] :=
    (
        ClearxIdealCache["MetricConcomitants"];
        ClearxIdealCache["WeylConcomitants"];
        ClearxIdealCache["ConnectionTensorConcomitants"]
    )

(* ::Section:: *)
(* Computation of the Petrov types *)


(*
TODO: there are two algorithms for doing this computation. Merge them in the same function.
*)
Options[PetrovType] = {Method -> "Default", PSimplify -> $CVSimplify, Parallelize -> True, Verbose -> True, "Observer" -> Null, "Vector" -> Null, "Bivector" -> Null}

PetrovType[metric_CTensor, opts : OptionsPattern[]] :=
    Module[{method},
        method = OptionValue[Method];
        Switch[method,
            "Default" || "WeylSelfDual",
                petrovType1[metric, opts]
            ,
            "PetrovMatrix",
                petrovType2[metric, opts]
            ,
            _,
                petrovType1[metric, opts]
        ]
    ]

(* Method 1: "WeylSelfdual" *)
petrovType1[metric_CTensor, opts : OptionsPattern[]] :=
	Catch @
		Module[{cart, cd, weylcd, epsilonmetric, weyldual, weylselfdual, g2form, weylselfdual2, weylselfdual3, aa, bb,
			 rho, simplf},
			If[Not @ MetricQ @ metric,
				Throw[Message[PetrovType::nometric, metric]]
			];
			simplf = OptionValue[PetrovType, {opts}, PSimplify];
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
			Which[
   				weylcd === Zero,
       					"Type O"
	    			,
				weylselfdual2 === Zero,
					"Type N"
				,
				weylselfdual3 === Zero,
					"Type III"
				,
				simplf[aa weylselfdual2 - aa^2 / 3 g2form - bb weylselfdual] === Zero,
					"Type D"
				,
				6 bb^2 - aa^3 === 0,
					"Type II"
				,
				True,
					"Type I"
			]
		]

(* Method 2: "PetrovMatrix" *)
petrovType2[metric_CTensor, opts : OptionsPattern[]] :=
	Catch @
		Module[{Q, gamma, Q2, aa, Q3, bb, simplf},
			If[Not @ MetricQ @ metric,
				Throw[Message[PetrovType::nometric, metric]]
			];
			(*TODO: make sure that an observer is included among the options*)
			simplf = OptionValue[PetrovType, {opts}, PSimplify];
			Q = weylConcomitant["WeylMatrixQ"][metric, opts];
			gamma = metricConcomitant["SpatialMetric"][metric, opts];
			Q2 = weylConcomitant["WeylMatrixQ2"][metric, opts];
			aa = weylConcomitant["TraceWeylMatrixQ2"][metric, opts];
			Q3 = weylConcomitant["WeylMatrixQ3"][metric, opts];
			bb = -weylConcomitant["TraceWeylMatrixQ3"][metric, opts];
			Which[
				Q === Zero,
       					"Type O"
	    		,
				Q2 === Zero,
					"Type N"
				,
				Q3 === Zero,
					"Type III"
				,
				simplf[(aa^2 / 3) gamma - aa Q2 - bb Q] === Zero,
					"Type D"
				,
				simplf[6 bb^2 - aa^3] === 0,
					"Type II"
				,
				True,
					"Type I"
			]
		]


(* ::Section:: *)
(* Computation of multiple Deveber null directions for each Petrov type *)


Options[DebeverNullDirections] = {Method -> "Default", PSimplify -> $CVSimplify, Verbose -> True, Parallelize -> True, "Observer" -> Null, "Vector" -> Null, "Bivector" -> Null}

DebeverNullDirections[metric_CTensor, opts : OptionsPattern[]] :=
    Module[{method},
        method = OptionValue[Method];
        Switch[method,
            "Default" || "WeylSelfDual",
                debeverNullDirections1[metric, opts]
            ,
            "PetrovMatrix",
                debeverNullDirections2[metric, opts]
            ,
            _,
                debeverNullDirections1[metric, opts]
        ]
    ]


(* Method 1: "WeylSelfdual" *)
debeverNullDirections1[metric_CTensor, opts : OptionsPattern[]] :=
Catch@ Module[{ptype},
		If[Not @ MetricQ @ metric,
			Throw[Message[DebeverNullDirections::nometric, metric]]
		];
  		ptype = PetrovType[metric, opts];
		Which[
			ptype === "Type O",
				Print["Type O"]
			,
			ptype === "Type N",
				Print["Type N"];
				weylConcomitant["WNullDirectionTypeN"][metric, opts]
			,
			ptype === "Type III",
				Print["Type III"];
				weylConcomitant["WNullDirectionTypeIII"][metric, opts]
			,
			ptype === "Type D",
				Print["Type D"];
				weylConcomitant["WNullDirectionTypeD"][metric, opts]
			,
			ptype === "Type II",
				Print["Type II"];
				weylConcomitant["WNullDirectionTypeII"][metric, opts]
			,
			True,
				Print["Type I"]
		]
	]

(* Method 2: "PetrovMatrix" *)
(* TODO: Make sure that, if Type D, "Observer" is a string with two arguments *)
debeverNullDirections2[metric_CTensor, opts : OptionsPattern[]] :=
Catch@ Module[{ptype},
		If[Not @ MetricQ @ metric,
			Throw[Message[DebeverNullDirections::nometric, metric]]
		];
  		ptype = PetrovType[metric, opts];
		Which[
			ptype === "Type O",
				Print["Type O"]
			,
			ptype === "Type N",
				Print["Type N"];
				weylConcomitant["NullDirectionTypeN"][metric, opts]
			,
			ptype === "Type III",
				Print["Type III"];
				weylConcomitant["NullDirectionTypeIII"][metric, opts]
			,
			ptype === "Type D",
				Print["Type D"];
				weylConcomitant["NullDirectionTypeD"][metric, opts]
			,
			ptype === "Type II",
				Print["Type II"];
				weylConcomitant["NullDirectionTypeII"][metric, opts]
			,
			True,
				Print["Type I"]
		]
	]



(* ::Section:: *)
(*  Classification of type D metrics*)
(* TODO: Add the new Weyl concomitants that appear here *)

(* Test when a symbolic function is non-negative *)
(* 
Recall that the value of Assumptions option is always logical statement. 
Therefore it should be expressed in terms of the logical syntax  
*)
Options[SymbolicPositiveQ] = {Assumptions -> True, PSimplify -> Simplify, Verbose -> True};

SymbolicPositiveQ[x_, OptionsPattern[]] :=
	Block[{$Assumptions = $Assumptions && OptionValue[Assumptions]},
		Which[
			Simplify[x] === 0,
				False
			,
			Simplify[x] === Zero,
				False
			,
			Simplify[ComplexExpand[x + Abs[x]]] === 0,
				False
			,
			Simplify[ComplexExpand[x - Abs[x]]] === 0,
				True
			,
			True,
				"Unknown"
		]
	]


Options[TypeDClassify] = {Assumptions -> True, Method -> "Default", PSimplify -> $CVSimplify}

TypeDClassify[metric_CTensor, w_CTensor, opts : OptionsPattern[]] :=
	Catch @
		Module[{cart, cd, W, RicciCD, epsilonmetric, logrho, TrW3, rho, drho, dlogrho, alpha, S, P, Q, C3, a, b, c, d, e,
			 f, i, j, k, l, C5, assumptions, simplf},
			If[Not @ MetricQ @ metric,
				Throw[Message[TypeDClassify::nometric, metric]]
			];
			assumptions = OptionValue[Assumptions];
   			simplf = OptionValue[PSimplify];
			cart = Part[metric, 2, 1, -1];
			{a, b, c, d, e, f, i, j, k, l} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 10];
			W = weylConcomitant["Weyl"][metric, opts];
			cd = CovDOfMetric[metric];
			RicciCD = Ricci[cd];
			epsilonmetric = epsilon[metric];
			TrW3 = weylConcomitant["TraceWeyl3"][metric, opts];
			rho = simplf[-(1/12 TrW3) ^ (1/3)];
			logrho = CTensor[Log[rho], {}];
			rho = CTensor[rho, {}];
			drho = simplf[TensorDerivative[rho, cd]];
			dlogrho = simplf[TensorDerivative[logrho, cd]];
			alpha = simplf[1/9 metric[-i, -j] dlogrho[i] dlogrho[j] - 2 rho[]];
			S = HeadOfTensor[
       				1 / (3 rho[]) (W[-a, -b, -c, -d] - rho[] (metric[-a, -c] metric[-b, -d] - 
					metric[-a, -d] metric[-b, -c])), {-a, -b, -c, -d}
     			];
			S = simplf[S];
			P = HeadOfTensor[epsilonmetric[-a, -b, -i, -j] W[-k, i, -l, j] drho[k] drho[l], {-a, -b}];
			P = simplf[P];
			Q = HeadOfTensor[S[-i, -a, -j, -b] drho[i] drho[j], {-a, -b}];
			Q = simplf[Q];
			C3 = 1/2 S[-a, -b, -i, -j] S[i, j, -c, -d] + S[-a, -b, -c, -d];
			C3 = simplf[C3];
			C5 = 2 Q[-i, -j] w[i] w[j] + Q[-k, k];
			C5 = simplf[C5];
			Which[
				RicciCD =!= Zero,
					"No vacuum"
				,
				rho[] === Zero || C3 =!= 0,
					"Vacuum no Type D"
				,
				P === Zero,
					Which[
						SymbolicPositiveQ[C5, Assumptions -> assumptions] === "Undefined",
							"Undefined sign in C5"
						,
						SymbolicPositiveQ[C5, Assumptions -> assumptions],
							Which[
								SymbolicPositiveQ[alpha, Assumptions -> assumptions] === "Undefined",
									
									"Undefined sign in \[Alpha]"
								,
								SymbolicPositiveQ[alpha, Assumptions -> assumptions],
									"\!\(\*SubscriptBox[\(A\), \(1\)]\)-metric"
								,
								alpha === Zero,
									"\!\(\*SubscriptBox[\(A\), \(3\)]\)-metric"
								,
								Not @ SymbolicPositiveQ[alpha, Assumptions -> assumptions],
									
									"\!\(\*SubscriptBox[\(A\), \(2\)]\)-metric"
								,
								True,
									(*Should be "Undefined" *)SymbolicPositiveQ[alpha, Assumptions-> assumptions]
							]
						,
						Not @ SymbolicPositiveQ[C5, Assumptions -> assumptions],
							Which[
								SymbolicPositiveQ[alpha, Assumptions -> assumptions] === "Undefined",
									
									"Undefined sign in \[Alpha]"
								,
								SymbolicPositiveQ[alpha, Assumptions -> assumptions],
									"\!\(\*SubscriptBox[\(B\), \(1\)]\)-metric"
								,
								alpha === Zero,
									"\!\(\*SubscriptBox[\(B\), \(3\)]\)-metric"
								,
								Not @ SymbolicPositiveQ[alpha, Assumptions -> assumptions],
									
									"\!\(\*SubscriptBox[\(B\), \(2\)]\)-metric"
								,
								True,
									(*Should be "Undefined" *)
									SymbolicPositiveQ[alpha, Assumptions -> assumptions]
							]
						,
						True,
							(*Should be "Undefined" *)
							SymbolicPositiveQ[alpha, Assumptions -> assumptions]
					]
				,
				True,
					"C-metric";
			]
		]

(* ::Section:: *)
(* Identification of the Kerr metric *)

(* Real and imaginary parts of symbolic complex quantities. *)

symbolicRe[expr_] := ComplexExpand[(expr + Dagger[expr]) / 2]

symbolicIm[expr_] := ComplexExpand[(expr - Dagger[expr]) / (2 I)]

symbolicComplexNorm2[expr_] := ComplexExpand[expr Dagger[expr]] 

Options[KerrSolutionQ] = {PSimplify -> $CVSimplify, Parallelize -> True, Verbose -> True}
KerrSolutionQ[metric_CTensor, opts : OptionsPattern[]] :=
Catch @
		Module[{cart, cd, weylcd, epsilonmetric, weyldual, weylselfdual, g2form, weylselfdual2, weylselfdual3, aa, bb,
			 rho, a1, b1, c1, d1, e1, f1, simplf, w, z, xi, riccicd, z1, z2, modz},
			If[Not @ MetricQ @ metric,
				Throw[Message[KerrSolutionQ::nometric, metric]]
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
					xi = weylConcomitant["TensorXi"][metric, opts];
					xi = Antisymmetrize[xi[-a1, -b1] xi[-c1, -d1], {-b1, -c1}];
					xi = HeadOfTensor[xi, {-a1, -b1, -c1, -d1}];
					If[simplf[xi] === Zero,
						Print["Kerr NUT"];
						weyldual = weylConcomitant["WeylDual"][metric, opts];
						weylselfdual = weylConcomitant["WeylSelfDual"][metric, opts];
						g2form = metricConcomitant["G2Form"][metric, opts];
						weylselfdual2 = weylConcomitant["WeylSelfDual2"][metric, opts];
						weylselfdual3 = weylConcomitant["WeylSelfDual3"][metric, opts];
						aa = 4 weylConcomitant["TraceWeylSelfDual2"][metric, opts];
						bb = 8 weylConcomitant["TraceWeylSelfDual3"][metric, opts];
						w = simplf[- bb / (2 aa)];
						z = weylConcomitant["ScalarZ"][metric, opts];
						If[symbolicIm[z^3 Dagger[w^8]] === 0,
							Print["Complex Kerr"];
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
(* Perfect fluid characterization *)

(* TODO: Check that an arbitrary time-like vector is given *)
Options[PerfectFluidQ] = {Assumptions -> True, PSimplify -> $CVSimplify, Verbose -> True, Parallelize -> True, "Vector" -> Null}

PerfectFluidQ[metric_CTensor, opts : OptionsPattern[]] :=
	Catch@ 
		Module[{cond1, cond2},
			If[Not @ MetricQ @ metric, 
    					Throw[Message[PerfectFluidQ::nometric, metric]]];
			Block[{$Assumptions = $Assumptions && OptionValue[Assumptions]},
				cond1 = metricConcomitant["FluPerCond1"][metric, opts];
				cond2 = metricConcomitant["FluPerCond2"][metric, opts];
				Which[
					SymbolicPositiveQ[cond2] === "Unknown",
						"Unknown"
					,
					cond1 === Zero && SymbolicPositiveQ[cond2],
						True
					,
					Not[cond1 === Zero] && SymbolicPositiveQ[cond2],
						False
					,
					cond1 === Zero && Not[SymbolicPositiveQ[cond2]],
						False
					,
					Not[cond1 === Zero] && Not[SymbolicPositiveQ[cond2]],
						False
					,
					True,
						"Unknown"
				]
			]
		]

(* TODO: Check that an arbitrary time-like vector is given *)
Options[PerfectFluidVariables] = {Assumptions -> True, PSimplify -> $CVSimplify, Verbose -> True, Parallelize -> True, "Vector" -> Null}

PerfectFluidVariables[metric_CTensor, opts : OptionsPattern[]] :=
	Catch@ 
		Module[{edens, press, flow},
			If[Not @ MetricQ @ metric, 
    					Throw[Message[PerfectFluidVariables::nometric, metric]]
			];
			If[Not @ PerfectFluidQ[metric, opts],
						Throw[Message[PerfectFluidVariables::noperfectfluid, metric]]
			];
			edens = metricConcomitant["EnergyDensity"][metric, opts];
			press = metricConcomitant["Pressure"][metric, opts];
			flow = metricConcomitant["FluPerFlow"][metric, opts];
			{edens, press, flow}
		]

(* TODO: Check that an arbitrary time-like vector is given *)
Options[ThermodynamicPerfectFluidQ] = {Assumptions -> True, PSimplify -> $CVSimplify, Verbose -> True, Parallelize -> True, "Vector" -> Null}

ThermodynamicPerfectFluidQ[metric_CTensor, opts : OptionsPattern[]] :=
	Catch@
		Module[{cond},
			If[Not @ MetricQ @ metric, 
    					Throw[Message[ThermodynamicPerfectFluidQ::nometric, metric]]
			];
			If[Not @ PerfectFluidQ[metric, opts],
						Throw[Message[ThermodynamicPerfectFluidQ::noperfectfluid, metric]]
			];
			Block[{$Assumptions = $Assumptions && OptionValue[Assumptions]},
				cond = metricConcomitant["ThermoFluPerCond"][metric, opts];
				If[cond === 0, True, False, False]
			]
		]

(* TODO: Check that an arbitrary time-like vector is given *)
Options[GenericIdealGasQ] = {Assumptions -> True, PSimplify -> $CVSimplify, Verbose -> True, Parallelize -> True, "Vector" -> Null}

GenericIdealGasQ[metric_CTensor, opts : OptionsPattern[]] :=
	Catch@
		Module[{cond},
			If[Not @ MetricQ @ metric, 
    					Throw[Message[GenericIdealGasQ::nometric, metric]]
			];
			If[Not @ PerfectFluidQ[metric, opts],
						Throw[Message[GenericIdealGasQ::noperfectfluid, metric]]
			];
			If[Not @ ThermodynamicPerfectFluidQ[metric, opts],
						Throw[Message[GenericIdealGasQ::nothermodynamicperfectfluid, metric]]
			];
			Block[{$Assumptions = $Assumptions && OptionValue[Assumptions]},
				cond = metricConcomitant["GenericIdealGasCond"][metric, opts];
				If[cond === 0, True, False, False]
			]
		]

(* ::Section:: *)
(*  Determination of the Connection Tensor*)

(*TODO: Add the documentation of this function*)
(*TODO: If no R-frame is given, depending on the Petrov Type more options will be necessary*)
(*TODO: Fix the fact that weylConcomitants do not recognize the "Rframe" option*)
Options[ConnectionTensor] = {Rframe -> {Null, Null, Null, Null}, PSimplify -> $CVSimplify, Parallelize -> True, Verbose -> True}

ConnectionTensor[metric_CTensor, opts : OptionsPattern[]] :=
	Catch@ 
		Module[{simplf, e0, e1, e2, e3, connectionTens, cart, cd, a1, b1, c1, vb, time},
			If[Not@MetricQ@metric, 
    				Throw[Message[ConnectionTensor::nometric, metric]]];
			{simplf, vb} = OptionValue[weylConcomitant, {opts}, {PSimplify, Verbose}];
			{e0, e1, e2, e3} = OptionValue[Rframe];
			cart = Part[metric, 2, 1, -1];
			{a1, b1, c1} = GetIndicesOfVBundle[VBundleOfBasis @ cart, 3];
			cd = CovDOfMetric[metric];
			If[
				e0 =!= Null && e1 =!= Null && e2 =!= Null && e3 =!= Null,
					time = AbsoluteTime[];
					connectionTens = -(1/2) (-(cd[-a1][e0[-b1]] e0[-c1] - cd[-a1][e0[-c1]] e0[-b1]) + (cd[-a1][e1[-b1]] e1[-c1] - cd[-a1][e1[-c1]] e1[-b1]) + 
						(cd[-a1][e2[-b1]] e2[-c1] - cd[-a1][e2[-c1]] e2[-b1]) + (cd[-a1][e3[-b1]] e3[-c1] - cd[-a1][e3[-c1]] e3[-b1]));
					If[vb, 
						Print["** ReportCompute: computing \"ConnectionTensor\" in ", AbsoluteTime[] - time, " seconds:"]
					];
					connectionTens = HeadOfTensor[connectionTens, {-a1, -b1, -c1}];
					time = AbsoluteTime[];
					connectionTens = simplf[connectionTens];
					If[vb,
						Print["** ReportCompute: applying  ", simplf, " to \"ConnectionTensor\" in ", AbsoluteTime[] - time, " seconds:"]
					],
				ptype = PetrovType[metric, opts];
				Which[
					ptype === "Type O",
						Print["Type O, an R-frame is needed"];
						connectionTens = Null;
					,
					ptype === "Type N",
						Print["Type N, an R-frame is needed"];
						connectionTens = Null;
					,
					ptype === "Type III",
						connectionTens = weylConcomitant["PTIIIConnectionTensor"][metric, opts];
					,
					ptype === "Type D",
						Print["Type D, an R-frame is needed"];
						connectionTens = Null;
					,
					ptype === "Type II",
						connectionTens = weylConcomitant["PTIIConnectionTensor"][metric, opts];
					,
					True,
						connectionTens = weylConcomitant["PTIConnectionTensor"][metric, opts];
				]
			];
			connectionTens
		]

(* ::Section:: *)
(*  Determination of the dimension of the isometry group*)

Options[IsometryGroupDimension] = {PSimplify -> $CVSimplify, Parallelize -> True, Verbose -> True}
IsometryGroupDimension[metric_CTensor, H_CTensor, opts : OptionsPattern[]] :=
 	Catch@ 
  		Module[{simplf, C1, C2, C3, C4, C11, C12, C122, C123, C1233, 
    			C1234, C1222, C1223, C111, C112, C1122, C1123, C1111, C1112},
   			If[Not@MetricQ@metric, 
    				Throw[Message[IsometryGroupDimension::nometric, metric]]];
			simplf = OptionValue[PSimplify];
   			C1 = ConnectionTensorConcomitant["C1"][metric, H, opts];
   			Which[
      
    				C1 === Zero,
    					Print["\!\(\*SubscriptBox[\(G\), \(4\)]\)"],
	 
    				C11 = ConnectionTensorConcomitant["C11"][metric, H, opts];
    				C11 === Zero,
    					Which[
	 
     						C12 = ConnectionTensorConcomitant["C12"][metric, H, opts];
     						C12 === Zero,
     							Print["\!\(\*SubscriptBox[\(G\), \(3\)]\)"],
     				
     						C122 = ConnectionTensorConcomitant["C122"][metric, H, opts];
     						C122 === Zero,
     							Which[
	    
      								C123 = ConnectionTensorConcomitant["C123"][metric, H, opts];
      								C123 === Zero,
      									Print["\!\(\*SubscriptBox[\(G\), \(2  b\)]\)"],
      						
      								C1233 = ConnectionTensorConcomitant["C1233"][metric, H, opts];
      								Not[C1233 == Zero],
      									Print["No symmetries"],
      						
      								C1234 = ConnectionTensorConcomitant["C1234"][metric, H, opts];
      								C1234 === Zero,
      									Print["\!\(\*SubscriptBox[\(G\), \(1  d\)]\)"],
	       
      								True,
      									Print["No symmetries"]
      							],
     				
     						C1222 = ConnectionTensorConcomitant["C1222"][metric, H, opts];
     						Not[C1222 === Zero],
     							Print["No symmetries"],
     				
     						C1223 = ConnectionTensorConcomitant["C1223"][metric, H, opts];
     						C1223 === Zero,
     							Print["\!\(\*SubscriptBox[\(G\), \(1  c\)]\)"],
	    
     						True,
     							Print["No symmetries"]
     					],
    				C111 = ConnectionTensorConcomitant["C111"][metric, H, opts];
    				C111 === Zero,
    					Which[
     				
     						C112 = ConnectionTensorConcomitant["C112"][metric, H, opts];
     						C112 === Zero,
     							Print["\!\(\*SubscriptBox[\(G\), \(2  a\)]\)"],
     				
     						C1122 = ConnectionTensorConcomitant["C1122"][metric, H, opts];
     						Not[C1122 === Zero],
     							Print["No symmetries"],
     				
     						C1123 = ConnectionTensorConcomitant["C1123"][metric, H, opts];
     						C1123 === Zero,
     							Print["\!\(\*SubscriptBox[\(G\), \(1  b\)]\)"],
	    
     						True,
     							Print["No symmetries"]
     					],
    				C1111 = ConnectionTensorConcomitant["C1111"][metric, H, opts];
    				Not[C1111 === Zero],
    					Print["No symmetries"],
	 
    				C1112 = ConnectionTensorConcomitant["C1112"][metric, H, opts];
    				C1112 === Zero,
    					Print["\!\(\*SubscriptBox[\(G\), \(1  a\)]\)"],
	 
    				True,
    					Print["No symmetries"]		
    			]
   
   		]


(* ::Section:: *)
(* SaveExactSolution *)



Options[SaveExactSolution] = {
	"ParameterNames" -> {}, 
	"ParameterAssumptions" -> Null , 
	"IsIdeal" -> True, 
	"CoordinateSystemName" -> " ",
	"CoordinateNames" -> {},
	"CoordinateAssumptions" -> Null,
	"ScalarFunctionNames" -> {},
	"Classes" -> {" "}
}

SaveExactSolution[metric_Function, exactsolname_String, opts : OptionsPattern[]] :=
	Module[{params, paramassmp, isideal, coords, sclrs, coordsassmp, funcs, cls},

		{params, paramassmp, isideal, sclrs, coords, coordsassmp, funcs, cls} = OptionValue[{"ParameterNames", "ParameterAssumptions", "IsIdeal", "CoordinateNames", "CoordinateSystemName", "CoordinateAssumptions", "ScalarFunctionNames", "Classes"}];

		exactSolsData[exactsolname, "ParameterNames"] = params;

		exactSolsData[exactsolname, "ParameterAssumptions"] = paramassmp;

		exactSolsData[exactsolname, "IsIDEAL"] = isideal;

		exactSolsData[exactsolname, "Classes"] = cls;

		exactSolsData[exactsolname, {coords, "CoordinateNames"}] = sclrs;

		exactSolsData[exactsolname, {coords, "CoordinateAssumptions"}] = coordsassmp;

		exactSolsData[exactsolname, {coords, "ParameterNames"}] = exactSolsData[exactsolname, "ParameterNames"];

		exactSolsData[exactsolname, {coords, "ParameterAssumptions"}] = exactSolsData[exactsolname, "ParameterAssumptions"];

		exactSolsData[exactsolname, {coords, "ScalarFunctionNames"}] = funcs;

		exactSolsData[exactsolname, {coords, "Metric"}] = metric;

	]

(****************************************************************)

(****************** 5. End private and package ******************)

(****************************************************************)


End[];

EndPackage[];
