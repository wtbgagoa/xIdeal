BeginPackage["xAct`ExactSolsData`"]

GenRelExactSolsData::usage = " ";

Begin["xAct`xIdeal`Private`"]

(* ::Section:: *)
(* Data *)

allmetrics = {
	"BertottiRobinsonSolution",
 	"Friedmann",
	"FarnsworthKerrI",
	"FarnsworthKerrII",
	"FarnsworthKerrIII",
	"GeneralSphericalSymmetry",
 	"GeneralSzekeresSzafron",
 	"KantowskiSachs",
	"KasnerI",
	"KasnerII",
	"KasnerIII",
	"Kerr",
 	"KerrNUT",
	"LemaitreTolman",
	"TaubI",
	"TaubII",
	"ThermodynamicStephani",
	"OsvathKoutrasI",
	"OsvathKoutrasII",
	"OsvathKoutrasIII",
	"PetrovSolution",
	"PPWave",
	"ReissnerNordstrom",
	"Schwarzschild",
   	"Stephani",
	"StephaniThermodynamic",
	"StephaniThermodynamicSpherical",
    "SzekeresSzafronI",
    "SzekeresSzafronII",
    "Wills",
	"WindmillI"
}

allclasses = {
	"AbelianG3onT3",
	"AxialSymmetry",
	"BarotropicPerfectFluid",
 	"ConformallyFlat",
	"ConformallyStatic", 
	"ConformallyStationary",
	"DMetrics",
 	"Dust",
	"EinsteinMaxwellSolution",
	"G1",
	"G2",
	"G3",
	"G3onS2",
	"G3IonS3", 
	"G3IIonS3", 
	"G3IIIonS3", 
	"G3IVonS3", 
	"G3VonS3", 
	"G3VI0onS3", 
	"G3VIhonS3", 
	"G3VII0onS3", 
	"G3VIIhonS3", 
	"G3VIIIonS3", 
	"G3IXonS3",
	"G3IonT3", 
	"G3IIonT3", 
	"G3IIIonT3", 
	"G3IVonT3", 
	"G3VonT3", 
	"G3VIonT3", 
	"G3VI0onT3", 
	"G3VIhonT3", 
	"G3VII0onT3", 
	"G3VIIhonT3", 
	"G3VIIIonT3", 
	"G3IXonT3",  
	"G4",
	"G5",
	"G6",
	"G7",
	"G8",
	"G9",
	"G10",
 	"Geodesic",
	"PerfectFluid",
	"PetrovTypeI",
	"PetrovTypeII",
	"PetrovTypeIII",
	"PetrovTypeN",
	"PetrovTypeD",
	"SpatialG1", 
	"SpatialG2", 
	"SpatialG3", 
	"SpatialG4", 
	"SpatialG5",
	"SpatialG6", 
	"SphericalSymmetry",
	"Static", 
	"Stationary",
 	"ThermodynamicPerfectFluid",
	"Vacuum",
	"VacuumTypeD",
	"Warped22"
}

allmetricproperties = {
	"IsIDEAL",
	"CoordinateAssumptions",
	"CoordinateSystemName",
	"ParameterAssumptions",
	"CoordinateNames",
	"ParameterNames",
	"ScalarFunctionNames",
	"ExactSolutionName",
	"Metric"
}

allcoordinatesystems = {
	"AdaptedCoordinates",
	"BoyerLindquistCoordinates",
	"CanonicalCoordinates",
	"ComplexCoordinates",
	"ExpansionGradientAdaptedCoordinates",
	"GroupGeneratorsAdaptedCoordinates",
 	"IsotropicCoordinates",
  	"ReducedCircumferencePolarCoordinates",
 	"SchwarzschildCoordinates",
	"SphericalCoordinates"
}

exactSolsData["Solutions"] = allmetrics

exactSolsData["Classes"] = allclasses

exactSolsData["Properties"] = allmetricproperties

exactSolsData["CoordinateSystems"] = allcoordinatesystems

(* Valid metrics *)
Set[metricQ[#], True]& /@ allmetrics;
metricQ[_] := False;

(* Valid classes of exact solutions/metrics *)
Set[exactsolclassQ[#], True]& /@ allclasses;
exactsolclassQ[_] := False;

(* Valid properties of exactsolutions/metrics *)
Set[metricpropertyQ[#], True]& /@ allmetricproperties;
metricpropertyQ[_] := False;

(* Valid coordinate systems *)
Set[coordinatesystemQ[#], True]& /@ allcoordinatesystems;
coordinatesystemQ[_] := False;


(* ::Section:: *)
(* Classification of exact solutions *)



exactSolsData["SphericalSymmetry"] = {
	"GeneralSphericalSymmetry",
	"Schwarzschild",
	"StephaniSpherical"
}

exactSolsData["AxialSymmetry"] = {
	"Kerr"
}

exactSolsData["PetrovTypeD"] = {
	"Kerr",
	"Schwarzschild"
}

exactSolsData["Vacuum"] = {
	"Kerr",
	"Schwarzschild"
}

(* ::Section:: *)
(* Exact solutions database *)


(* ::Subsection:: *)
(* Bertotti-Robinson Solution *)

exactSolsData["BertottiRobinsonSolution", "Classes"] = {"EinsteinMaxwellSolution", "PetrovTypeN", "ConformallyFlat"}

exactSolsData["BertottiRobinsonSolution", "CoordinateSystems"] = {"ComplexCoordinates"}

exactSolsData["BertottiRobinsonSolution", "DefaultCoordinates"] = "ComplexCoordinates"

exactSolsData["BertottiRobinsonSolution", "IsIDEAL"] = False
 
exactSolsData["BertottiRobinsonSolution", "ParameterAssumptions"] = Null
 
exactSolsData["BertottiRobinsonSolution", "ParameterNames"] = {"k"}
 
exactSolsData["BertottiRobinsonSolution", {"ComplexCoordinates", "CoordinateAssumptions"}] = Null
 
exactSolsData["BertottiRobinsonSolution", {"ComplexCoordinates", "CoordinateNames"}] = {"u", "v", "x1", "x2"}
 
exactSolsData["BertottiRobinsonSolution", {"ComplexCoordinates", "Metric"}] = 
    Function[{coords, params, scfuncs}, 
    	With[{u = coords[[1]], v = coords[[2]], x1 = coords[[3]], x2 = coords[[4]], k = params[[1]]}, 
    		{
				{(-k)*v^2, -1, 0, 0}, 
				{-1, 0, 0, 0}, 
				{0, 0, 0, (1 + (k/2)*x1*x2)^(-2)}, 
				{0, 0, (1 + (k/2)*x1*x2)^(-2), 0}
			}
		]
	]
 
exactSolsData["BertottiRobinsonSolution", {"ComplexCoordinates", "ParameterAssumptions"}] = exactSolsData["BertottiRobinsonSolution", "ParameterAssumptions"]
 
exactSolsData["BertottiRobinsonSolution", {"ComplexCoordinates", "ParameterNames"}] = exactSolsData["BertottiRobinsonSolution", "ParameterNames"]
 
exactSolsData["BertottiRobinsonSolution", {"ComplexCoordinates", "ScalarFunctionNames"}] = {}


(* ::Subsection:: *)
(* Farnsworth-Kerr I *)

exactSolsData["FarnsworthKerrI", "Classes"] = {"PerfectFluid", "Homogeneous", "G4", "G3IXonS3"}
 
exactSolsData["FarnsworthKerrI", "CoordinateSystems"] = {"CanonicalCoordinates"}

exactSolsData["FarnsworthKerrI", "DefaultCoordinates"] = "CanonicalCoordinates"

exactSolsData["FarnsworthKerrI", "IsIDEAL"] = False
 
exactSolsData["FarnsworthKerrI", "ParameterAssumptions"] = 
    Function[{coords, params, scfuncs}, 
    	With[{k = params[[1]], a = params[[2]]}, 
    		Element[k, Reals] && Abs[k] < 1/2 && Element[a, Reals];
		]
	]
 
exactSolsData["FarnsworthKerrI", "ParameterNames"] = {"k", "a"}
 
exactSolsData["FarnsworthKerrI", {"CanonicalCoordinates", "CoordinateAssumptions"}] = Null
 
exactSolsData["FarnsworthKerrI", {"CanonicalCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}
 
exactSolsData["FarnsworthKerrI", {"CanonicalCoordinates", "Metric"}] = 
    Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]], k = params[[1]], a = params[[2]]}, 
    		a^2*{
				{-1, 0, (-Sqrt[1 - 2*k^2])*Sin[x], (-Sqrt[1 - 2*k^2])*Cos[x]*Cos[y]}, 
				{0, 1 - k, 0, (1 - k)*Sin[y]}, 
        		{(-Sqrt[1 - 2*k^2])*Sin[x], 0, (1 + k)*Cos[x]^2 + (1 + 2*k^2)*Sin[x]^2, 
					k*(2*k - 1)*Cos[x]*Sin[x]*Cos[y]}, 
				{(-Sqrt[1 - 2*k^2])*Cos[x]*Cos[y], (1 - k)*Sin[y], k*(2*k - 1)*Cos[x]*Sin[x]*Cos[y], 
					(1 - k)*Sin[y]^2 + (1 + k)*Sin[x]^2*Cos[y]^2 + (1 + 2*k^2)*Cos[x]^2*Cos[y]^2}
			}
		]
	]
 
exactSolsData["FarnsworthKerrI", {"CanonicalCoordinates", "ParameterAssumptions"}] = exactSolsData["FarnsworthKerrI", "ParameterAssumptions"]
 
exactSolsData["FarnsworthKerrI", {"CanonicalCoordinates", "ParameterNames"}] = exactSolsData["FarnsworthKerrI", "ParameterNames"]
 
exactSolsData["FarnsworthKerrI", {"CanonicalCoordinates", "ScalarFunctionNames"}] = {}


(* ::Subsection:: *)
(* Farnsworth-Kerr II *)

exactSolsData["FarnsworthKerrII", "Classes"] = {"PerfectFluid", "Homogeneous", "G4", "G3VIIIonT3", "G3IIIonS3"}
 
exactSolsData["FarnsworthKerrII", "CoordinateSystems"] = {"CanonicalCoordinates"}

exactSolsData["FarnsworthKerrII", "DefaultCoordinates"] = "CanonicalCoordinates"

exactSolsData["FarnsworthKerrII", "IsIDEAL"] = False
 
exactSolsData["FarnsworthKerrII", "ParameterAssumptions"] = 
    Function[{coords, params, scfuncs}, 
    	With[{k = params[[1]], a = params[[2]]}, 
      		Element[k, Reals] && Inequality[1, Less, 4*k^2, LessEqual, 2] && Element[a, Reals];
		]
	]
 
exactSolsData["FarnsworthKerrII", "ParameterNames"] = {"k", "a"}
 
exactSolsData["FarnsworthKerrII", {"CanonicalCoordinates", "CoordinateAssumptions"}] = Null
 
exactSolsData["FarnsworthKerrII", {"CanonicalCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}
 
exactSolsData["FarnsworthKerrII", {"CanonicalCoordinates", "Metric"}] = 
    Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]], k = params[[1]], a = params[[2]]}, 
    		(-a^2)*{
				{1, 0, Sqrt[1 - 2*k^2]*Sin[x], Sqrt[1 - 2*k^2]*Cos[x]*Cosh[y]}, 
        		{0, 1 - k, 0, (1 - k)*Sinh[y]}, 
				{Sqrt[1 - 2*k^2]*Sin[x], 0, (1 + k)*Cos[x]^2 - (1 + 2*k^2)*Sin[x]^2, 
					(-(2 + k + 2*k^2))*Cos[x]*Sin[x]*Cosh[y]}, 
				{Sqrt[1 - 2*k^2]*Cos[x]*Cosh[y], (1 - k)*Sinh[y], (-(2 + k + 2*k^2))*Cos[x]*Sin[x]*Cosh[y], 
					(1 - k)*Sinh[y]^2 + (1 + k)*Sin[x]^2*Cosh[y]^2 - (1 + 2*k^2)*Cos[x]^2*Cosh[y]^2}
			}
		]
	]
 
exactSolsData["FarnsworthKerrII", {"CanonicalCoordinates", "ParameterAssumptions"}] = exactSolsData["FarnsworthKerrII", "ParameterAssumptions"]
 
exactSolsData["FarnsworthKerrII", {"CanonicalCoordinates", "ParameterNames"}] = exactSolsData["FarnsworthKerrII", "ParameterNames"]
 
exactSolsData["FarnsworthKerrII", {"CanonicalCoordinates", "ScalarFunctionNames"}] = {}


(* ::Subsection:: *)
(* Farnsworth-Kerr III *)

exactSolsData["FarnsworthKerrIII", "Classes"] = {"DMetrics", "PetrovTypeD", "PerfectFluid", "Homogeneous", "G4", "G3IXonS3"}
 
exactSolsData["FarnsworthKerrIII", "CoordinateSystems"] = {"CanonicalCoordinates"}

exactSolsData["FarnsworthKerrIII", "DefaultCoordinates"] = "CanonicalCoordinates"

exactSolsData["FarnsworthKerrIII", "IsIDEAL"] = False
 
exactSolsData["FarnsworthKerrIII", "ParameterAssumptions"] = 
    Function[{coords, params, scfuncs}, 
    	With[{s = params[[1]], a = params[[2]]}, 
    		Element[s, Reals] && Abs[s] < 1 && Element[a, Reals];
		]
	]
 
exactSolsData["FarnsworthKerrIII", "ParameterNames"] = {"s", "a"}
 
exactSolsData["FarnsworthKerrIII", {"CanonicalCoordinates", "CoordinateAssumptions"}] = Null
 
exactSolsData["FarnsworthKerrIII", {"CanonicalCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}
 
exactSolsData["FarnsworthKerrIII", {"CanonicalCoordinates", "Metric"}] = 
    Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]], s = params[[1]], a = params[[2]]}, 
    		(-a^2)*{
				{1, 0, 0, 0}, 
				{0, 1 - s, 0, (-(1 - s))*Sinh[y]}, 
        		{0, 0, (1 + s)*Cos[x]^2 - 2*Sin[x]^2, (-(3 + s))*Cos[x]*Sin[x]*Cosh[y]}, 
				{0, (-(1 - s))*Sinh[y], (-(3 + s))*Cos[x]*Sin[x]*Cosh[y], 
					(1 - s)*Sinh[y]^2 + (1 + s)*Sin[x]^2*Cosh[y]^2 - 2*Cos[x]^2*Cosh[y]^2}
			}
		]
	]
 
exactSolsData["FarnsworthKerrIII", {"CanonicalCoordinates", "ParameterAssumptions"}] = exactSolsData["FarnsworthKerrIII", "ParameterAssumptions"]
 
exactSolsData["FarnsworthKerrIII", {"CanonicalCoordinates", "ParameterNames"}] = exactSolsData["FarnsworthKerrIII", "ParameterNames"]
 
exactSolsData["FarnsworthKerrIII", {"CanonicalCoordinates", "ScalarFunctionNames"}] = {}


(* ::Subsection:: *)
(* Friedmann in reduced circumference polar coordinates *)

exactSolsData["Friedmann", "Classes"] = {"PerfectFluid", "ThermodynamicPerfectFluid",
     "ConformallyFlat", "SpatiallyHomogeneous", "SpatialG6", "BarotropicPerfectFluid",
     "ConformallyStatic"}

exactSolsData["Friedmann", "CoordinateSystems"] = {"ReducedCircumferencePolarCoordinates"}

exactSolsData["Friedmann", "DefaultCoordinates"] = "ReducedCircumferencePolarCoordinates"

exactSolsData["Friedmann", "IsIDEAL"] = True

exactSolsData["Friedmann", "ParameterNames"] = {"k"}

exactSolsData["Friedmann", "ParameterAssumptions"] = Null

exactSolsData["Friedmann", {"ReducedCircumferencePolarCoordinates", "CoordinateNames"}] = {"t", "r", "\[Theta]", "\[Phi]"}

exactSolsData["Friedmann", {"ReducedCircumferencePolarCoordinates", "CoordinateAssumptions"}] = Null

exactSolsData["Friedmann", {"ReducedCircumferencePolarCoordinates", "ParameterNames"}] = {"k"}

exactSolsData["Friedmann", {"ReducedCircumferencePolarCoordinates", "ParameterAssumptions"}] = Null

exactSolsData["Friedmann", {"ReducedCircumferencePolarCoordinates", "ScalarFunctionNames"}] = {"R"}


exactSolsData["Friedmann", {"ReducedCircumferencePolarCoordinates", "Metric"}] =
    Function[{coords, params, scfuncs},
        With[{t = coords[[1]], r = coords[[2]], theta = coords[[3]], 
            phi = coords[[4]], k = params[[1]], R = scfuncs[[1]]},
            DiagonalMatrix[
				{
					-1, 
					R[t] ^ 2 / (1 - k * r^2), 
					R[t] ^ 2 * r^2, 
					R[t] ^ 2 * r^2 * Sin[theta] ^ 2
				}
			]
        ]
    ]




(* ::Subsection:: *)
(* GeneralSpherical metric in spherical coordinates *)

exactSolsData["GeneralSphericalSymmetry", "Classes"] = {"DMetrics", "PerfectFluid", "PetrovTypeD", 
	"ThermodynamicPerfectFluid", "SphericalSymmetry", "Warped22"}

exactSolsData["GeneralSphericalSymmetry", "CoordinateSystems"] = {"SphericalCoordinates"} 

exactSolsData["GeneralSphericalSymmetry", "DefaultCoordinates"] = "SphericalCoordinates"

exactSolsData["GeneralSphericalSymmetry", "ParameterNames"] = {}

exactSolsData["GeneralSphericalSymmetry", "ParameterAssumptions"] = Null

exactSolsData["GeneralSphericalSymmetry", "IsIDEAL"] = True

exactSolsData["GeneralSphericalSymmetry", {"SphericalCoordinates", "CoordinateNames"}] = {"t", "r", "\[Theta]", "\[Phi]"}

exactSolsData["GeneralSphericalSymmetry", {"SphericalCoordinates", "CoordinateAssumptions"}] = #[[2]] > 0 && Pi > #[[3]] > 0 &

exactSolsData["GeneralSphericalSymmetry", {"SphericalCoordinates", "ParameterNames"}] = exactSolsData["GeneralSpherical", "ParameterNames"]

exactSolsData["GeneralSphericalSymmetry", {"SphericalCoordinates", "ParameterAssumptions"}] = exactSolsData["GeneralSpherical", "ParameterAssumptions"]

exactSolsData["GeneralSphericalSymmetry", {"SphericalCoordinates", "ScalarFunctionNames"}] = {"\[Lambda]", "\[Mu]", "\[Nu]"}

defaultcoordinates["GeneralSphericalSymmetry"] = "SphericalCoordinates"

(* The syntax is exactSolsData[args__][{coords_List, parameters_List, functions_List}] *)

exactSolsData["GeneralSphericalSymmetry", {"SphericalCoordinates", "Metric"}] =
	Function[{coords, params, funcs},
		With[{nu = funcs[[3]], lambda = funcs[[1]], mu = funcs[[2]], t = coords[[1]], r = coords[[2]], theta = coords[[3]], phi = coords[[4]]},
			DiagonalMatrix[
				{
					-E^(2 nu[t, r]), 
					E^(2 lambda[t, r]), 
					E^(2 mu[t, r]), 
					E^(2 mu[t, r]) Sin[theta]^2
				}
			]
		] 
	]


(* ::Subsection:: *)
(* Kerr *)

exactSolsData["Kerr", "Classes"] = {"DMetrics", "PetrovTypeD", 
	"AxialSymmetry", "Vacuum", "Stationary"}
 
exactSolsData["Kerr", "CoordinateSystems"] = {"BoyerLindquistCoordinates"}

exactSolsData["Kerr", "DefaultCoordinates"] = "BoyerLindquistCoordinates"

exactSolsData["Kerr", "IsIDEAL"] = True
 
exactSolsData["Kerr", "ParameterAssumptions"] = 
    Function[{coords, params, scfuncs}, 
    	With[{m = params[[1]], a = params[[2]]}, 
			m > 0 && a >= 0
		]
	]
 
exactSolsData["Kerr", "ParameterNames"] = {"m", "a"}
 
exactSolsData["Kerr", {"BoyerLindquistCoordinates", "CoordinateAssumptions"}] = 
	Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], r = coords[[2]], theta = coords[[3]], 
    		phi = coords[[4]], m = params[[1]], a = params[[2]], Sigma = scfuncs[[1]]}, 
				Sigma = 
					Function[{r, a, theta}, 
						r^2 + a^2*Cos[theta]^2
					]; 
				r > 0 && Pi > theta > 0 && Element[Cos[theta], Reals] && 
					Sin[theta] > 0 && Sigma > 0
		]
	]
 
exactSolsData["Kerr", {"BoyerLindquistCoordinates", "CoordinateNames"}] = 
    {"t", "r", "\[Theta]", "\[Phi]"}
 
exactSolsData["Kerr", {"BoyerLindquistCoordinates", "Metric"}] = 
    Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], r = coords[[2]], theta = coords[[3]], phi = coords[[4]], 
			m = params[[1]], a = params[[2]], Sigma = scfuncs[[1]], Delta = scfuncs[[2]]}, 
    			Sigma = 
					Function[{r, a, theta}, 
						r^2 + a^2*Cos[theta]^2
					]; 
       			Delta = 
					Function[{r, m, a}, 
						r^2 - 2*m*r + a^2
					]; 
        		{
					{-(1 - (2*m*r)/Sigma[r, a, theta]), 0, 0, -((2*m*r*a*Sin[theta]^2)/Sigma[r, a, theta])}, 
        			{0, Sigma[r, a, theta]/Delta[r, m, a], 0, 0}, 
        			{0, 0, Sigma[r, a, theta], 0}, 
        			{-((2*m*r*a*Sin[theta]^2)/Sigma[r, a, theta]), 0, 0, 
						(r^2 + a^2 + ((2*m*r*a^2)/Sigma[r, a, theta])*Sin[theta]^2)*Sin[theta]^2}
				}
		]
	]
 
exactSolsData["Kerr", {"BoyerLindquistCoordinates", "ParameterAssumptions"}] = exactSolsData["Kerr", "ParameterAssumptions"]
 
exactSolsData["Kerr", {"BoyerLindquistCoordinates", "ParameterNames"}] = exactSolsData["Kerr", "ParameterNames"]
 
exactSolsData["Kerr", {"BoyerLindquistCoordinates", "ScalarFunctionNames"}] = 
    {"\[CapitalSigma]", "\[CapitalDelta]"}


(* ::Subsection:: *)
(* Kerr-NUT *)

exactSolsData["KerrNUT", "Classes"] = {"DMetrics", "PetrovTypeD", "Vacuum"}

exactSolsData["KerrNUT", "CoordinateSystems"] = {"ExpansionGradientAdaptedCoordinates"}

exactSolsData["KerrNUT", "DefaultCoordinates"] = "ExpansionGradientAdaptedCoordinates"
 
exactSolsData["KerrNUT", "IsIDEAL"] = True
 
exactSolsData["KerrNUT", "ParameterAssumptions"] = Null
 
exactSolsData["KerrNUT", "ParameterNames"] = {"p", "k", "s"}
 
exactSolsData["KerrNUT", {"ExpansionGradientAdaptedCoordinates", "CoordinateAssumptions"}] = 
	Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]], p = params[[1]], 
			k = params[[2]], s = params[[3]], alpha = scfuncs[[1]], beta = scfuncs[[2]]}, 
      			alpha = 
					Function[{p, x, k, s}, 
						p*x^2 + k*(3 - k^2)*(x/(1 + k^2)^3) + s
					]; 
				beta = 
					Function[{p, y, k, s}, 
						(-p)*y^2 + (3*k^2 - 1)*(y/(1 + k^2)^3) + s
					]; 
				x^2 + y^2 > 0 && alpha[p, x, k, s] > 0 && beta[p, y, k, s] > 0
		]
	]
 
exactSolsData["KerrNUT", {"ExpansionGradientAdaptedCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}
 
exactSolsData["KerrNUT", {"ExpansionGradientAdaptedCoordinates", "Metric"}] = 
    Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]], p = params[[1]], 
			k = params[[2]], s = params[[3]], alpha = scfuncs[[1]], beta = scfuncs[[2]]}, 
      			alpha = 
					Function[{p, x, k, s}, 
						p*x^2 + k*(3 - k^2)*(x/(1 + k^2)^3) + s
					]; 
				beta = 
					Function[{p, y, k, s}, 
						(-p)*y^2 + (3*k^2 - 1)*(y/(1 + k^2)^3) + s
					]; 
       			{
					{-((alpha[p, x, k, s]*y^4 - beta[p, y, k, s]*x^4)/(x^2 + y^2)), 0, 0, 
        				-((alpha[p, x, k, s]*y^2 + beta[p, y, k, s]*x^2)/(x^2 + y^2))}, 
        			{0, (x^2 + y^2)/alpha[p, x, k, s], 0, 0}, 
       				{0, 0, (x^2 + y^2)/beta[p, y, k, s], 0}, 
        			{-((alpha[p, x, k, s]*y^2 + beta[p, y, k, s]*x^2)/(x^2 + y^2)), 0, 0, 
        				(beta[p, y, k, s] - alpha[p, x, k, s])/(x^2 + y^2)}
				}
		]
	]
 
exactSolsData["KerrNUT", {"ExpansionGradientAdaptedCoordinates", "ParameterAssumptions"}] = exactSolsData["KerrNUT", "ParameterAssumptions"]
 
exactSolsData["KerrNUT", {"ExpansionGradientAdaptedCoordinates", "ParameterNames"}] = exactSolsData["KerrNUT", "ParameterNames"]
 
exactSolsData["KerrNUT", {"ExpansionGradientAdaptedCoordinates", "ScalarFunctionNames"}] = {"\[Alpha]", "\[Beta]"}


(* ::Subsection:: *)
(* Lemaitre-Tolman *)

exactSolsData["LemaitreTolman", "Classes"] = {"DMetrics", "PerfectFluid", 
    "PetrovTypeD", "ThermodynamicPerfectFluid", "SphericalSymmetry", "Warped22"}

exactSolsData["LemaitreTolman", "CoordinateSystems"] = {"SphericalCoordinates"}

exactSolsData["LemaitreTolman", "DefaultCoordinates"] = "SphericalCoordinates"

exactSolsData["LemaitreTolman", "IsIDEAL"] = True
 
exactSolsData["LemaitreTolman", "ParameterAssumptions"] = Null
 
exactSolsData["LemaitreTolman", "ParameterNames"] = {}
 
exactSolsData["LemaitreTolman", {"SphericalCoordinates", "CoordinateAssumptions"}] = Null
 
exactSolsData["LemaitreTolman", {"SphericalCoordinates", "CoordinateNames"}] = {"t", "r", "\[Theta]", "\[Phi]"}
 
exactSolsData["LemaitreTolman", {"SphericalCoordinates", "Metric"}] = 
    Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], r = coords[[2]], theta = coords[[3]], phi = coords[[4]], R = scfuncs[[1]]}, 
			DiagonalMatrix[
       			{
					-1, 
					D[R[t, r], r]^2, 
					R[t, r]^2, 
					R[t, r]^2*Sin[theta]^2
				}
			]
		]
	]
 
exactSolsData["LemaitreTolman", {"SphericalCoordinates", "ParameterAssumptions"}] = exactSolsData["LemaitreTolman", "ParameterAssumptions"]
 
exactSolsData["LemaitreTolman", {"SphericalCoordinates", "ParameterNames"}] = exactSolsData["LemaitreTolman", "ParameterNames"]
 
exactSolsData["LemaitreTolman", {"SphericalCoordinates", "ScalarFunctionNames"}] = {"R"}


(* ::Subsection:: *)
(* Osvath-Koutras I *)

exactSolsData["OsvathKoutrasI", "Classes"] = {"PetrovTypeI", "PerfectFluid", "Homogeneous", "G4", 
	"AbelianG3onT3", "G3VIhonS3", "G3VIonT3"}

exactSolsData["OsvathKoutrasI", "CoordinateSystems"] = {"GroupGeneratorsAdaptedCoordinates"}

exactSolsData["OsvathKoutrasI", "DefaultCoordinates"] = "GroupGeneratorsAdaptedCoordinates"
 
exactSolsData["OsvathKoutrasI", "IsIDEAL"] = False
 
exactSolsData["OsvathKoutrasI", "ParameterAssumptions"] = 
    Function[{coords, params, scfuncs}, 
    	With[{s = params[[1]], a = params[[2]]}, 
    		Inequality[1/2, LessEqual, s^2, Less, 1.2296814706969093] && Element[a, Reals]
		]
	]
 
exactSolsData["OsvathKoutrasI", "ParameterNames"] = {"s", "a", "\[Beta]", "A", "B", "F", "b"}
 
exactSolsData["OsvathKoutrasI", {"GroupGeneratorsAdaptedCoordinates", "CoordinateAssumptions"}] = Null
 
exactSolsData["OsvathKoutrasI", {"GroupGeneratorsAdaptedCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}
 
exactSolsData["OsvathKoutrasI", {"GroupGeneratorsAdaptedCoordinates", "Metric"}] = 
	Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]], s = params[[1]], a = params[[2]], 
			beta = params[[3]], A = params[[4]], B = params[[5]], F = params[[6]], b = params[[7]]}, 
      			beta = 
					Function[{s}, 
						Sqrt[1 + 2*s^2*(1 - s^2)*(3 - s^2)]
					]; 
       			A = Function[{beta, s}, 
						(1 - beta[s])/2
					]; 
       			B = 
					Function[{beta, s}, 
						(1 + beta[s])/2]; F = Function[{s}, 1 - s^2
					]; 
       			b = 
					Function[{s}, 
						Sqrt[2]*s*(3 - s^2)
					]; 
    			a^2*{
					{((4/b[s]^2)*A[beta, s]^2 - 1)*Exp[2*A[beta, s]*z], 
						((4/b[s]^2)*A[beta, s]*B[beta, s] - 1)*Exp[(A[beta, s] + B[beta, s])*z], 0, 0}, 
        			{((4/b[s]^2)*A[beta, s]*B[beta, s] - 1)*Exp[(A[beta, s] + B[beta, s])*z], 
						((4/b[s]^2)*B[beta, s]^2 - 1)*Exp[2*B[beta, s]*z], 0, 0}, 
					{0, 0, Exp[2*F[s]*z], 0}, 
        			{0, 0, 0, 1}
				}
		]
	]
 
exactSolsData["OsvathKoutrasI", {"GroupGeneratorsAdaptedCoordinates", "ParameterAssumptions"}] = exactSolsData["OsvathKoutrasI", "ParameterAssumptions"]
 
exactSolsData["OsvathKoutrasI", {"GroupGeneratorsAdaptedCoordinates", "ParameterNames"}] = exactSolsData["OsvathKoutrasI", "ParameterNames"]
 
exactSolsData["OsvathKoutrasI", {"GroupGeneratorsAdaptedCoordinates", "ScalarFunctionNames"}] = {}


(* ::Subsection:: *)
(* Osvath-Koutras II *)

exactSolsData["OsvathKoutrasII", "Classes"] = {"PetrovTypeI", "PerfectFluid", "Homogeneous", "G4", 
	"AbelianG3onT3", "G3VIhonS3", "G3IVonT3"}

exactSolsData["OsvathKoutrasII", "CoordinateSystems"] = {"GroupGeneratorsAdaptedCoordinates"}

exactSolsData["OsvathKoutrasII", "DefaultCoordinates"] = "GroupGeneratorsAdaptedCoordinates"
 
exactSolsData["OsvathKoutrasII", "IsIDEAL"] = False
 
exactSolsData["OsvathKoutrasII", "ParameterAssumptions"] = 
    Function[{coords, params, scfuncs}, 
		With[{a = params[[2]]}, 
    		Element[a, Reals]
		]
	]
 
exactSolsData["OsvathKoutrasII", "ParameterNames"] = {"s", "a", "F", "b"}
 
exactSolsData["OsvathKoutrasII", {"GroupGeneratorsAdaptedCoordinates", "CoordinateAssumptions"}] = Null
 
exactSolsData["OsvathKoutrasII", {"GroupGeneratorsAdaptedCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}
 
exactSolsData["OsvathKoutrasII", {"GroupGeneratorsAdaptedCoordinates", "Metric"}] = 
	Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]], s = params[[1]], a = params[[2]], 
			F = params[[3]], b = params[[4]]}, 
				s = Sqrt[1.2296814706969093]; 
       			F = Function[{s}, 1 - s^2]; 
				b = Function[{s}, Sqrt[2]*s*(3 - s^2)]; 
       			a^2*{
					{(-4^(-1))*(b[s] - 1/b[s])^2*E^z, (-4^(-1))*(b[s] - 1/b[s])^2*E^z*z, 0, 0}, 
					{(-4^(-1))*(b[s] - 1/b[s])^2*E^z*z, E^z*(1 - (z^2/4)*(b[s] - 1/b[s])^2), 0, 0}, 
					{0, 0, E^(2*F[s]*z), 0}, 
					{0, 0, 0, 1}
				}
		]
	]
 
exactSolsData["OsvathKoutrasII", {"GroupGeneratorsAdaptedCoordinates", "ParameterAssumptions"}] = exactSolsData["OsvathKoutrasII", "ParameterAssumptions"]
 
exactSolsData["OsvathKoutrasII", {"GroupGeneratorsAdaptedCoordinates", "ParameterNames"}] = exactSolsData["OsvathKoutrasII", "ParameterNames"]
 
exactSolsData["OsvathKoutrasII", {"GroupGeneratorsAdaptedCoordinates", "ScalarFunctionNames"}] = {}


(* ::Subsection:: *)
(* Osvath-Koutras III *)

exactSolsData["OsvathKoutrasIII", "Classes"] = {"PetrovTypeI", "PerfectFluid", "Homogeneous", "G4", 
	"AbelianG3onT3", "G3VIIhonS3"}

exactSolsData["OsvathKoutrasIII", "CoordinateSystems"] = {"GroupGeneratorsAdaptedCoordinates"}

exactSolsData["OsvathKoutrasIII", "DefaultCoordinates"] = "GroupGeneratorsAdaptedCoordinates"
 
exactSolsData["OsvathKoutrasIII", "IsIDEAL"] = False
 
exactSolsData["OsvathKoutrasIII", "ParameterAssumptions"] = 
    Function[{coords, params, scfuncs}, 
    	With[{s = params[[1]], a = params[[2]]}, 
    		Inequality[1.2296814706969093, Less, s^2, LessEqual, 2] && Element[a, Reals]
		]
	]
 
exactSolsData["OsvathKoutrasIII", "ParameterNames"] = {"s", "a", "\[Beta]2", "F", "b", "k"}
 
exactSolsData["OsvathKoutrasIII", {"GroupGeneratorsAdaptedCoordinates", "CoordinateAssumptions"}] = Null
 
exactSolsData["OsvathKoutrasIII", {"GroupGeneratorsAdaptedCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}
 
exactSolsData["OsvathKoutrasIII", {"GroupGeneratorsAdaptedCoordinates", "Metric"}] = 
	Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]], s = params[[1]], a = params[[2]], 
       		beta2 = params[[3]], F = params[[4]], b = params[[5]], k = params[[6]]}, 
	   			beta2 = Function[{s}, 1 + 2*s^2*(1 - s^2)*(3 - s^2)]; 
				F = Function[{s}, 1 - s^2]; 
       			b = Function[{s}, Sqrt[2]*s*(3 - s^2)]; 
       			k = Function[{beta2, s}, Sqrt[-beta2]/2]; 
       			a^2*{
					{
						E^z*((1/b[s]^2)*(Cos[k[beta2, s]*z] - 2*k[beta2, s]*Sin[k[beta2, s]*z])^2 - Cos[k[beta2, s]*z]^2), 
          				E^z*((1/b[s]^2)*(Cos[k[beta2, s]*z] - 2*k[beta2, s]*Sin[k[beta2, s]*z])*(2*k[beta2, s]*Cos[k[beta2, s]*z] + 
							Sin[k[beta2, s]*z]) - Cos[k[beta2, s]*z]*Sin[k[beta2, s]*z]), 
          				0, 
						0
					}, 
		  			{
						E^z*((1/b[s]^2)*(Cos[k[beta2, s]*z] - 2*k[beta2, s]*Sin[k[beta2, s]*z])*(2*k[beta2, s]*Cos[k[beta2, s]*z] + 
              				Sin[k[beta2, s]*z]) - Cos[k[beta2, s]*z]*Sin[k[beta2, s]*z]), 
          				E^z*((1/b[s]^2)*(2*k[beta2, s]*Cos[k[beta2, s]*z] + Sin[k[beta2, s]*z])^2 - Sin[k[beta2, s]*z]^2), 
						0, 
						0
					}, 
         			{0, 0, Exp[2*F[s]*z], 0}, 
		 			{0, 0, 0, 1}
				}
		]
	]
 
exactSolsData["OsvathKoutrasIII", {"GroupGeneratorsAdaptedCoordinates", "ParameterAssumptions"}] = exactSolsData["OsvathKoutrasIII", "ParameterAssumptions"]
 
exactSolsData["OsvathKoutrasIII", {"GroupGeneratorsAdaptedCoordinates", "ParameterNames"}] = exactSolsData["OsvathKoutrasIII", "ParameterNames"]
 
exactSolsData["OsvathKoutrasIII", {"GroupGeneratorsAdaptedCoordinates", "ScalarFunctionNames"}] = {}


(* ::Subsection:: *)
(* Petrov Solution *)

exactSolsData["PetrovSolution", "Classes"] = {"Vacuum", "G4", "G3IonT3", "G3VIIhonT3", "PetrovTypeI"}

exactSolsData["PetrovSolution", "CoordinateSystems"] = {"CanonicalCoordinates"}

exactSolsData["PetrovSolution", "DefaultCoordinates"] = "CanonicalCoordinates"
 
exactSolsData["PetrovSolution", "IsIDEAL"] = False
 
exactSolsData["PetrovSolution", "ParameterAssumptions"] = Null
 
exactSolsData["PetrovSolution", "ParameterNames"] = {"k"}
 
exactSolsData["PetrovSolution", {"CanonicalCoordinates", "CoordinateAssumptions"}] = Null
 
exactSolsData["PetrovSolution", {"CanonicalCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}
 
exactSolsData["PetrovSolution", {"CanonicalCoordinates", "Metric"}] = 
    Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]], k = params[[1]]}, 
    		{
				{(-E^x)*Cos[Sqrt[3]*x], 0, 0, (-E^x)*Sin[Sqrt[3]*x]}, 
				{0, 1, 0, 0}, 
        		{0, 0, E^(-2*x), 0}, 
				{(-E^x)*Sin[Sqrt[3]*x], 0, 0, E^x*Cos[Sqrt[3]*x]}
			}/k^2
		]
	]
 
exactSolsData["PetrovSolution", {"CanonicalCoordinates", "ParameterAssumptions"}] = Null
 
exactSolsData["PetrovSolution", {"CanonicalCoordinates", "ParameterNames"}] = {"k"}
 
exactSolsData["PetrovSolution", {"CanonicalCoordinates", "ScalarFunctionNames"}] = {}


(* ::Subsection:: *)
(* Reissner-Nordström *)

exactSolsData["ReissnerNordstrom", "CoordinateSystems"] = {"SchwarzschildCoordinates"}

exactSolsData["ReissnerNordstrom", "DefaultCoordinates"] = "SchwarzschildCoordinates"

exactSolsData["ReissnerNordstrom", "ParameterNames"] = {"m", "rQ"}

exactSolsData["ReissnerNordstrom", "ParameterAssumptions"] = 
    Function[{coords, params, scfuncs}, 
    	With[{m = params[[1]], rQ = params[[2]]}, 
			m > 0 && rQ > 0
		]
	]

exactSolsData["ReissnerNordstrom", "IsIDEAL"] = True

exactSolsData["ReissnerNordstrom", "Classes"] = {"DMetrics", "PetrovTypeD", 
	"SphericalSymmetry", "Warped22", "Static", "EinsteinMaxwellSolution"}

exactSolsData["ReissnerNordstrom", {"SchwarzschildCoordinates", "CoordinateNames"}] = {"t", "r", "\[Theta]", "\[Phi]"}

exactSolsData["ReissnerNordstrom", {"SchwarzschildCoordinates", "CoordinateAssumptions"}] = 
	Function[{coords, params, scfuncs}, 
     	With[{t = coords[[1]], r = coords[[2]], theta = coords[[3]], phi = coords[[4]]}, 
			r > 0 && Pi > theta > 0
		]
	]

exactSolsData["ReissnerNordstrom", {"SchwarzschildCoordinates", "ScalarFunctionNames"}] = {}

exactSolsData["ReissnerNordstrom", {"SchwarzschildCoordinates", "ParameterNames"}] = exactSolsData["ReissnerNordstrom", "ParameterNames"]

exactSolsData["ReissnerNordstrom", {"SchwarzschildCoordinates", "ParameterAssumptions"}] = exactSolsData["ReissnerNordstrom", "ParameterAssumptions"]

exactSolsData["ReissnerNordstrom", {"SchwarzschildCoordinates", "Metric"}] = 
    Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], r = coords[[2]], theta = coords[[3]], phi = coords[[4]], m = params[[1]], rQ = params[[2]]}, 
    		DiagonalMatrix[
				{
					-(1 - (2*m)/r + rQ^2/r^2), 
					(1 - (2*m)/r + rQ^2/r^2)^(-1), 
					r^2, 
					r^2*Sin[theta]^2
				}
			]
		]
	]

(* ::Subsection:: *)
(* Schwarzschild in Schwarzschild coordinates *)
exactSolsData["Schwarzschild", "Classes"] = {"PetrovTypeD", "Static", "SphericalSymmetry", "Vacuum", "VacuumTypeD"}

exactSolsData["Schwarzschild", "CoordinateSystems"] = {"SchwarzschildCoordinates", "IsotropicCoordinates"}

exactSolsData["Schwarzschild", "DefaultCoordinates"] = "SchwarzschildCoordinates"

exactSolsData["Schwarzschild", "ParameterNames"] = {"m"}

exactSolsData["Schwarzschild", "ParameterAssumptions"] = #[[1]] > 0 &

exactSolsData["Schwarzschild", "IsIDEAL"] = True

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "CoordinateNames"}] = {"t", "r", "\[Theta]", "\[Phi]"}

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "CoordinateAssumptions"}] = #[[2]] > 0 && Pi > #[[3]] > 0 &

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "ParameterNames"}] = exactSolsData["Schwarzschild", "ParameterNames"]

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "ParameterAssumptions"}] = exactSolsData["Schwarzschild", "ParameterAssumptions"]

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "ScalarFunctionNames"}] = {}

defaultcoordinates["Schwarzschild"] = "SchwarzschildCoordinates"

(* The syntax is exactSolsData[args__][{coords_List, parameters_List, functions_List}] *)

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "Metric"}] = 
	Function[{coords, params, funcs},
		With[{t = coords[[1]], r = coords[[2]], theta = coords[[3]], phi = coords[[4]], m = params[[1]]},
			DiagonalMatrix[{-(1 - (2 m) / r), (1 - (2 m) / r)^-1, r^2, r^2 Sin[theta]^2}] 
		]
	]

(* ::Subsection:: *)
(* Schwarzschild in Isotropic coordinates *)

exactSolsData["Schwarzschild", {"IsotropicCoordinates", "ParameterNames"}] = {"m"}

exactSolsData["Schwarzschild", {"IsotropicCoordinates", "ParameterAssumptions"}] =
    Function[{coords, params, scfuncs},
        With[{m = params[[1]]},
            m > 0
        ]
    ]

exactSolsData["Schwarzschild", {"IsotropicCoordinates", "ScalarFunctionNames"}] = {"R"}

exactSolsData["Schwarzschild", {"IsotropicCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}

exactSolsData["Schwarzschild", {"IsotropicCoordinates", "CoordinateAssumptions"}] =
    Function[{coords, params, scfuncs},
        With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = 
            coords[[4]], m = params[[1]], R = scfuncs[[1]]},
            R =
                Function[{x, y, z},
                    Sqrt[x^2 + y^2 + z^2]
                ];
            m >= 0 && R[x, y, z] > 2 * m
        ]
    ]

exactSolsData["Schwarzschild", {"IsotropicCoordinates", "Metric"}] =
    Function[{coords, params, scfuncs},
        With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = 
            coords[[4]], m = params[[1]], R = scfuncs[[1]]},
            R =
                Function[{x, y, z},
                    Sqrt[x^2 + y^2 + z^2]
                ];
            DiagonalMatrix[
				{
					-((1 - m / (2 * R[x, y, z])) ^ 2 / (1 + m / (2 * R[x, y, z])) ^ 2), 
					(1 + m / (2 * R[x, y, z])) ^ 4, 
					(1 + m / (2 * R[x, y, z])) ^ 4, 
					(1 + m / (2 * R[x, y, z])) ^ 4
				}
			]
        ]
    ]


(* ::Subsection:: *)
(* Stephani in adapted coordinates *)

exactSolsData["Stephani", "Classes"] = {"PerfectFluid", "SpatialG6", "ConformallyFlat"}

exactSolsData["Stephani", "CoordinateSystems"] = {"AdaptedCoordinates"}

exactSolsData["Stephani", "DefaultCoordinates"] = "AdaptedCoordinates"
 
exactSolsData["Stephani", "IsIDEAL"] = True
 
exactSolsData["Stephani", "ParameterAssumptions"] = exactSolsData["Stephani", "ParameterAssumptions"]
 
exactSolsData["Stephani", "ParameterNames"] = exactSolsData["Stephani", "ParameterNames"]
 
exactSolsData["Stephani", {"AdaptedCoordinates", "CoordinateAssumptions"}] = Null
 
exactSolsData["Stephani", {"AdaptedCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}
 
exactSolsData["Stephani", {"AdaptedCoordinates", "Metric"}] = 
    Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]], Omega = scfuncs[[1]], alpha = scfuncs[[2]], 
    		R = scfuncs[[3]], b1 = scfuncs[[4]], b2 = scfuncs[[5]], b3 = scfuncs[[6]], k = scfuncs[[7]]}, 
      			Omega = 
					Function[{R, t, b1, b2, b3, x, y, z, k}, 
         				R[t]/(1 + 2*(b1[t]*x + b2[t]*y + b3[t]*z) + k[t]*((x^2 + y^2 + z^2)/4))
					]; 
       			alpha = 
					Function[{R, t, b1, b2, b3, x, y, z, k}, 
        				R[t]*(D[Omega[R, t, b1, b2, b3, x, y, z, k], t]/(Omega[R, t, b1, b2, b3, x, y, z, k]*D[R[t], t]))
					]; 
       			DiagonalMatrix[
					{
						-alpha[R, t, b1, b2, b3, x, y, z, k]^2, 
        				Omega[R, t, b1, b2, b3, x, y, z, k]^2, 
						Omega[R, t, b1, b2, b3, x, y, z, k]^2, 
						Omega[R, t, b1, b2, b3, x, y, z, k]^2
					}
				]
		]
	]
 
exactSolsData["Stephani", {"AdaptedCoordinates", "ParameterAssumptions"}] = Null
 
exactSolsData["Stephani", {"AdaptedCoordinates", "ParameterNames"}] = {}
 
exactSolsData["Stephani", {"AdaptedCoordinates", "ScalarFunctionNames"}] = {"\[CapitalOmega]", "\[Alpha]", "R", "b1", "b2", "b3", "K"}


(* ::Subsection:: *)
(* Stephani Thermodynamic in adapted coordinates *)

exactSolsData["StephaniThermodynamic", "Classes"] = {"PerfectFluid", "ThermodynamicPerfectFluid",
     "G3S2", "SpatialG6", "ConformallyFlat"}

exactSolsData["StephaniThermodynamic", "CoordinateSystems"] = {"AdaptedCoordinates"}

exactSolsData["StephaniThermodynamic", "DefaultCoordinates"] = "AdaptedCoordinates"

exactSolsData["StephaniThermodynamic", "IsIDEAL"] = True

exactSolsData["StephaniThermodynamic", "ParameterAssumptions"] = Null

exactSolsData["StephaniThermodynamic", "ParameterNames"] = {"\[CurlyEpsilon]"}

exactSolsData["StephaniThermodynamic", {"AdaptedCoordinates", "CoordinateAssumptions"}] = Null

exactSolsData["StephaniThermodynamic", {"AdaptedCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}

exactSolsData["StephaniThermodynamic", {"AdaptedCoordinates", "Metric"}] =
    Function[{coords, params, scfuncs},
        With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = 
            coords[[4]], epsilon = params[[1]], w = scfuncs[[1]], L = scfuncs[[2]],
            Omega = scfuncs[[3]], alpha = scfuncs[[4]], R = scfuncs[[5]], b = scfuncs[[6]]},
            w =
                Function[{x, y, z},
                    2 * (z / (1 + epsilon * ((x^2 + y^2 + z^2) / 4)))            
                ];
            L =
                Function[{R, b, t, w, x, y, z},
                    R[t] / (1 + b[t] * w[x, y, z])
                ];
            Omega =
                Function[{w, x, y, z, L, R, b, t},
                    w[x, y, z] * (L[R, b, t, w, x, y, z] / (2 * z))
                ];
            alpha =
                Function[{R, t, L, b, w, x, y, z},
                    R[t] * (D[L[R, b, t, w, x, y, z], t] / (L[R, b, t, w, x, y, z] * D[R[t], t]))
                ];
            DiagonalMatrix[
				{
					-alpha[R, t, L, b, w, x, y, z] ^ 2, 
					Omega[w, x, y, z, L, R, b, t] ^ 2, 
					Omega[w, x, y, z, L, R, b, t] ^ 2, 
					Omega[w, x, y, z, L, R, b, t] ^ 2
				}
			]
        ]
    ]

exactSolsData["StephaniThermodynamic", {"AdaptedCoordinates", "ParameterAssumptions"}] = Null

exactSolsData["StephaniThermodynamic", {"AdaptedCoordinates", "ParameterNames"}] = {"\[CurlyEpsilon]"}

exactSolsData["StephaniThermodynamic", {"AdaptedCoordinates", "ScalarFunctionNames"}] = {"w", "L", "\[CapitalOmega]", "\[Alpha]", "R", "b"}


(* ::Subsection:: *)
(* Stephani Thermodynamic and Spherically symmetric in adapted coordinates *)

exactSolsData["StephaniThermodynamicSpherical", "Classes"] = {"PerfectFluid", "ThermodynamicPerfectFluid", "G3onS2", "SpatialG6", 
    "ConformallyFlat", "SphericalSymmetry", "Warped22"}

exactSolsData["StephaniThermodynamicSpherical", "CoordinateSystems"] = {"SphericalCoordinates"}

exactSolsData["StephaniThermodynamicSpherical", "DefaultCoordinates"] = "SphericalCoordinates"
 
exactSolsData["StephaniThermodynamicSpherical", "IsIDEAL"] = True
 
exactSolsData["StephaniThermodynamicSpherical", "ParameterAssumptions"] = Null
 
exactSolsData["StephaniThermodynamicSpherical", "ParameterNames"] = {}
 
exactSolsData["StephaniThermodynamicSpherical", {"SphericalCoordinates", "CoordinateAssumptions"}] = Null
 
exactSolsData["StephaniThermodynamicSpherical", {"SphericalCoordinates", "CoordinateNames"}] = {"t", "r", "\[Theta]", "\[Phi]"}
 
exactSolsData["StephaniThermodynamicSpherical", {"SphericalCoordinates", "Metric"}] = 
	Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], r = coords[[2]], theta = coords[[3]], phi = coords[[4]], 
			Omega = scfuncs[[1]], alpha = scfuncs[[2]], R = scfuncs[[3]], k = scfuncs[[4]]}, 
    			Omega = 
					Function[{R, t, k, r}, 
						R[t]/(1 + k[t]*(r^2/4))
					]; 
    			alpha = 
					Function[{R, t, Omega, k, r}, 
						R[t]*(D[Omega[R, t, k, r], t]/(Omega[R, t, k, r]*D[R[t], t]))
					]; 
				DiagonalMatrix[
       				{
						-alpha[R, t, Omega, k, r]^2, 
						Omega[R, t, k, r]^2, 
        				Omega[R, t, k, r]^2*r^2, 
						Omega[R, t, k, r]^2*r^2*Sin[theta]^2
					}
				]; 
		]
	]
 
exactSolsData["StephaniThermodynamicSpherical", {"SphericalCoordinates", "ParameterAssumptions"}] = Null
 
exactSolsData["StephaniThermodynamicSpherical", {"SphericalCoordinates", "ParameterNames"}] = {}
 
exactSolsData["StephaniThermodynamicSpherical", {"SphericalCoordinates", "ScalarFunctionNames"}] = {"\[CapitalOmega]", "\[Alpha]", "R", "k"}

(* ::Section:: *)
(* GenRelExactSolsData *)

(* ::Subsection:: *)
(* Metrics: default coordinates *)

iGenRelExactSolsData[metric_?metricQ] := exactSolsData[metric, "CoordinateSystems"]

iGenRelExactSolsData[metric_?metricQ, "CoordinateSystems"] := exactSolsData[metric, "CoordinateSystems"]

iGenRelExactSolsData[metric_?metricQ, coords_?coordinatesystemQ] := exactSolsData[metric, {coords, "Metric"}]

(* The syntax is GRData[args__String, {coords_List, parameters_List, functions_List}] *)

iGenRelExactSolsData[metric_?metricQ, coordname_?coordinatesystemQ, {coords_List, parameters_List, functions_List}] := exactSolsData[metric, {coordname, "Metric"}][coords, parameters, functions]

iGenRelExactSolsData[metric_?metricQ, "Properties"] := allmetricproperties

iGenRelExactSolsData[metric_?metricQ, "IsIDEAL"] := exactSolsData[metric, "IsIDEAL"]

iGenRelExactSolsData[metric_?metricQ, "CoordinateAssumptions"] := Module[{coords},
	coords = exactSolsData[metric, "DefaultCoordinates"];
	exactSolsData[metric, {coords, "CoordinateAssumptions"}]
]

iGenRelExactSolsData[metric_?metricQ, "CoordinateSystemName"] := Module[{coords},
	coords = exactSolsData[metric, "DefaultCoordinates"]
]

iGenRelExactSolsData[metric_?metricQ, "CoordinateNames"] := Module[{coords},
	coords = exactSolsData[metric, "DefaultCoordinates"];
	exactSolsData[metric, {coords, "CoordinateNames"}]
]

iGenRelExactSolsData[metric_?metricQ, "ParameterAssumptions"] := exactSolsData[metric, "ParameterAssumptions"]

iGenRelExactSolsData[metric_?metricQ, "ParameterNames"] := exactSolsData[metric, "ParameterNames"]

iGenRelExactSolsData[metric_?metricQ, "ScalarFunctionNames"] := Module[{coords},
	coords = exactSolsData[metric, "DefaultCoordinates"];
	exactSolsData[metric, {coords, "ScalarFunctionNames"}]
]

iGenRelExactSolsData[metric_?metricQ, "Metric"] := Module[{coords},
	coords = exactSolsData[metric, "DefaultCoordinates"];
	exactSolsData[metric, {coords, "Metric"}]
]

(* ::Subsection:: *)
(* Metrics: user given coordinates *)

iGenRelExactSolsData[metric_?metricQ, {coordname_?coordinatesystemQ, "CoordinateAssumptions"}] := 
	If[
		FreeQ[exactSolsData[metric, "CoordinateSystems"], coordname],
		(* TODO: error control *)
		Throw[$Failed],
		exactSolsData[metric, {coordname, "CoordinateAssumptions"}]
	]
	
 
iGenRelExactSolsData[metric_?metricQ, {coordname_?coordinatesystemQ, "CoordinateNames"}] := 
	If[
		FreeQ[exactSolsData[metric, "CoordinateSystems"], coordname],
		(* TODO: error control *)
		Throw[$Failed],
		exactSolsData[metric, {coordname, "CoordinateNames"}]
	]

iGenRelExactSolsData[metric_?metricQ, {coordname_?coordinatesystemQ, "Metric"}] :=
	If[
		FreeQ[exactSolsData[metric, "CoordinateSystems"], coordname],
		(* TODO: error control *)
		Throw[$Failed],
		exactSolsData[metric, {coordname, "Metric"}]
	]

iGenRelExactSolsData[metric_?metricQ, {coordname_?coordinatesystemQ, "ParameterNames"}] :=
	If[
		FreeQ[exactSolsData[metric, "CoordinateSystems"], coordname],
		(* TODO: error control *)
		Throw[$Failed],
		exactSolsData[metric, {coordname, "ParameterNames"}]
	]

iGenRelExactSolsData[metric_?metricQ, {coordname_?coordinatesystemQ, "ParameterAssumptions"}] :=
	If[
		FreeQ[exactSolsData[metric, "CoordinateSystems"], coordname],
		(* TODO: error control *)
		Throw[$Failed],
		exactSolsData[metric, {coordname, "ParameterAssumptions"}]
	]

iGenRelExactSolsData[metric_?metricQ, {coordname_?coordinatesystemQ, "ScalarFunctionNames"}] :=
	If[
		FreeQ[exactSolsData[metric, "CoordinateSystems"], coordname],
		(* TODO: error control *)
		Throw[$Failed],
		exactSolsData[metric, {coordname, "ScalarFunctionNames"}]
	]


(* ::Subsection:: *)
(* General definitions *)
iGenRelExactSolsData[] = allmetrics;
iGenRelExactSolsData["ExactSolutions"] = allmetrics;
iGenRelExactSolsData[All] = allmetrics;
iGenRelExactSolsData["Classes"] = allclasses;
iGenRelExactSolsData[All, "Classes"] = allclasses;
iGenRelExactSolsData["Properties"] = allmetricproperties;
iGenRelExactSolsData[All, "Properties"] = allmetricproperties;
iGenRelExactSolsData["CoordinateSystems"] = allcoordinatesystems;
iGenRelExactSolsData[All, "CoordinateSystems"] = allcoordinatesystems;

iGenRelExactSolsData[___] := $Failed;

(* Entry point of GenRelExactSolsData *)
GenRelExactSolsData[args___]:= Module[{res = iGenRelExactSolsData[args]}, (res/; UnsameQ[res, $Failed])]

(****************************************************************)

(****************** 5. End private and package ******************)

(****************************************************************)


End[]

EndPackage[]