BeginPackage["xAct`ExactSolsData`"]

GenRelExactSolsData::usage = " ";

Begin["xAct`xIdeal`Private`"]

(* ::Section:: *)
(* Data *)

allmetrics = {
	"BertottiRobinson",
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
 	"KerrNut",
	"LemaitreTolman",
	"TaubI",
	"TaubII",
	"ThermodynamicStephani",
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
	"FunctionNames",
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


(* The syntax is GRData[args__String, {coords_List, parameters_List, functions_List}] *)

iGenRelExactSolsData[solname_String, coordname_String, "Metric", {coords_List, parameters_List, functions_List}] := exactSolsData[solname, {coordname, "Metric"}][coords, parameters, functions]

(* ::Subsection:: *)
(* Bertotti-Robinson Solution *)

exactSolsData["BertottiRobinsonSolution", "Classes"] = {"EinsteinMaxwellSolution", "PetrovTypeN", "ConformallyFlat"}

exactSolsData["BertottiRobinsonSolution", "CoordinateSystems"] = {"ComplexCoordinates"}

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
			};
		]
	]
 
exactSolsData["BertottiRobinsonSolution", {"ComplexCoordinates", "ParameterAssumptions"}] = exactSolsData["BertottiRobinsonSolution", "ParameterAssumptions"]
 
exactSolsData["BertottiRobinsonSolution", {"ComplexCoordinates", "ParameterNames"}] = exactSolsData["BertottiRobinsonSolution", "ParameterNames"]
 
exactSolsData["BertottiRobinsonSolution", {"ComplexCoordinates", "ScalarFunctionNames"}] = {}


(* ::Subsection:: *)
(* Farnsworth-Kerr I *)

exactSolsData["FarnsworthKerrI", "Classes"] = {"PerfectFluid", "Homogeneous", "G4", "G3IXonS3"}
 
exactSolsData["FarnsworthKerrI", "CoordinateSystems"] = {"CanonicalCoordinates"}

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
			};
		]
	]
 
exactSolsData["FarnsworthKerrI", {"CanonicalCoordinates", "ParameterAssumptions"}] = exactSolsData["FarnsworthKerrI", "ParameterAssumptions"]
 
exactSolsData["FarnsworthKerrI", {"CanonicalCoordinates", "ParameterNames"}] = exactSolsData["FarnsworthKerrI", "ParameterNames"]
 
exactSolsData["FarnsworthKerrI", {"CanonicalCoordinates", "ScalarFunctionNames"}] = {}


(* ::Subsection:: *)
(* Farnsworth-Kerr II *)

exactSolsData["FarnsworthKerrII", "Classes"] = {"PerfectFluid", "Homogeneous", "G4", "G3VIIIonT3", "G3IIIonS3"}
 
exactSolsData["FarnsworthKerrII", "CoordinateSystems"] = {"CanonicalCoordinates"}

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
			};
		]
	]
 
exactSolsData["FarnsworthKerrII", {"CanonicalCoordinates", "ParameterAssumptions"}] = exactSolsData["FarnsworthKerrII", "ParameterAssumptions"]
 
exactSolsData["FarnsworthKerrII", {"CanonicalCoordinates", "ParameterNames"}] = exactSolsData["FarnsworthKerrII", "ParameterNames"]
 
exactSolsData["FarnsworthKerrII", {"CanonicalCoordinates", "ScalarFunctionNames"}] = {}


(* ::Subsection:: *)
(* Farnsworth-Kerr III *)

exactSolsData["FarnsworthKerrIII", "Classes"] = {"DMetrics", "PetrovTypeD", "PerfectFluid", "Homogeneous", "G4", "G3IXonS3"}
 
exactSolsData["FarnsworthKerrIII", "CoordinateSystems"] = {"CanonicalCoordinates"}

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
			};
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

exactSolsData["KerrNut", "Classes"] = {"DMetrics", "PetrovTypeD", "Vacuum"}

exactSolsData["KerrNut", "CoordinateSystems"] = {"ExpansionGradientAdaptedCoordinates"}
 
exactSolsData["KerrNut", "IsIDEAL"] = True
 
exactSolsData["KerrNut", "ParameterAssumptions"] = Null
 
exactSolsData["KerrNut", "ParameterNames"] = {"p", "k", "s"}
 
exactSolsData["KerrNut", {"ExpansionGradientAdaptedCoordinates", "CoordinateAssumptions"}] = 
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
				x^2 + y^2 > 0 && alpha > 0 && beta > 0
		]
	]
 
exactSolsData["KerrNut", {"ExpansionGradientAdaptedCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}
 
exactSolsData["KerrNut", {"ExpansionGradientAdaptedCoordinates", "Metric"}] = 
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
				};
		]
	]
 
exactSolsData["KerrNut", {"ExpansionGradientAdaptedCoordinates", "ParameterAssumptions"}] = exactSolsData["KerrNut", "ParameterAssumptions"]
 
exactSolsData["KerrNut", {"ExpansionGradientAdaptedCoordinates", "ParameterNames"}] = exactSolsData["KerrNut", "ParameterNames"]
 
exactSolsData["KerrNut", {"ExpansionGradientAdaptedCoordinates", "ScalarFunctionNames"}] = {"\[Alpha]", "\[Beta]"}


(* ::Subsection:: *)
(* Lemaitre-Tolman *)

exactSolsData["LemaitreTolman", "Classes"] = {"DMetrics", "PerfectFluid", 
    "PetrovTypeD", "ThermodynamicPerfectFluid", "SphericalSymmetry", "Warped22"}
 
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
(* Petrov Solution *)

exactSolsData["PetrovSolution", "Classes"] = {"Vacuum", "G4", "G3IonT3", "G3VIIhonT3", "PetrovTypeI"}
 
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
			}/k^2; 
		]
	]
 
exactSolsData["PetrovSolution", {"CanonicalCoordinates", "ParameterAssumptions"}] = Null
 
exactSolsData["PetrovSolution", {"CanonicalCoordinates", "ParameterNames"}] = {"k"}
 
exactSolsData["PetrovSolution", {"CanonicalCoordinates", "ScalarFunctionNames"}] = {}


(* ::Subsection:: *)
(* Reissner-Nordström *)

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

exactSolsData["Schwarzschild", "ParameterNames"] = {"m"}

exactSolsData["Schwarzschild", "ParameterAssumptions"] = #[[1]] > 0 &

exactSolsData["Schwarzschild", "IsIDEAL"] = True

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "CoordinateNames"}] = {"t", "r", "\[Theta]", "\[Phi]"}

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "CoordinateAssumptions"}] = #[[2]] > 0 && Pi > #[[3]] > 0 &

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "ParameterNames"}] = exactSolsData["Schwarzschild", "ParameterNames"]

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "ParameterAssumptions"}] = exactSolsData["Schwarzschild", "ParameterAssumptions"]

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

(* Metrics *)


(* Special definitions *)
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