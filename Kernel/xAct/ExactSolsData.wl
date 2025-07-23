BeginPackage["xAct`ExactSolsData`"]

(* ::Section:: *)
(* Usage information *)

GRExactSolsData::usage = " ";

(* ::Section:: *)
(* Messages *)

GRExactSolsData::noprop = "Unknown argument or property";

(* ::Section:: *)
(* BeginPrivate *)

Begin["xAct`xIdeal`Private`"]

(* ::Section:: *)
(* all metrics, classes, properties and coordinate systems *)
allmetrics = {
	"BertottiRobinsonSolution",
	"ElectroVacTypeD",
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
	"Classes",
	"CoordinateAssumptions",
	"CoordinateNames",
	"CoordinateSystemName",
	"CoordinateSystems",
	"IsIDEAL",
	"Metric",
	"ParameterAssumptions",
	"ParameterNames",
	"ScalarFunctionValues",
	"ScalarFunctionNames"
}

allcoordinatesystems = {
	"AdaptedCoordinates",
	"BoyerLindquistCoordinates",
	"CanonicalCoordinates",
	"ComplexCoordinates",
	"ExpansionGradientAdaptedCoordinates",
	"GroupGeneratorsAdaptedCoordinates",
	"HarmonicCoordinates",
 	"IsotropicCoordinates",
	"PlanarCoordinates",
  	"ReducedCircumferencePolarCoordinates",
 	"SchwarzschildCoordinates",
	"SphericalCoordinates",
	"TypeDCoordinates"
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

(* TODO: finish this section *)

exactSolsData["AbelianG3onT3"] = {
	"OsvathKoutrasI", 
	"OsvathKoutrasII", 
	"OsvathKoutrasIII"
}

exactSolsData["AxialSymmetry"] = {
	"Kerr"
}

exactSolsData["BarotropicPerfectFluid"] = {
	"Friedmann"
}

exactSolsData["ConformallyFlat"] = {
	"BertottiRobinsonSolution", 
	"Friedmann", 
	"Stephani", 
	"StephaniThermodynamic", 
	"StephaniThermodynamicSpherical"
}

exactSolsData["ConformallyStatic"] = {
	"Friedmann"
}

exactSolsData["PetrovTypeD"] = {
	"Kerr",
	"Schwarzschild"
}

exactSolsData["DMetrics"] = {
	"FarnsworthKerrIII", 
	"GeneralSphericalSymmetry", 
	"Kerr", 
	"KerrNUT", 
	"LemaitreTolman", 
	"ReissnerNordstrom"
}

exactSolsData["EinsteinMaxwellSolution"] = {
	"BertottiRobinsonSolution", 
	"ReissnerNordstrom"
}

exactSolsData["G3IIIonS3"] = {
	"FarnsworthKerrII"
}

exactSolsData["G3IonT3"] = {
	"PetrovSolution"
} 

exactSolsData["G3IVonT3"] = {
	"OsvathKoutrasII"
}

exactSolsData["G3IXonS3"] = {
	"FarnsworthKerrI", 
	"FarnsworthKerrIII"
} 

exactSolsData["G3onS2"] = {
	"StephaniThermodynamicSpherical"
}

exactSolsData["G3S2"] = {
	"StephaniThermodynamic"
}

exactSolsData["G3VIhonS3"] = {
	"OsvathKoutrasI", 
	"OsvathKoutrasII"
}

exactSolsData["G3VIIhonS3"] = {
	"OsvathKoutrasIII"
}

exactSolsData["G3VIIhonT3"] = {
	"PetrovSolution"
}

exactSolsData["G3VIIIonT3"] = {
	"FarnsworthKerrII"
}

exactSolsData["G3VIonT3"] = {
	"OsvathKoutrasI"
}

exactSolsData["G4"] = {
	"FarnsworthKerrI", 
	"FarnsworthKerrII", 
	"FarnsworthKerrIII", 
	"OsvathKoutrasI", 
	"OsvathKoutrasII", 
	"OsvathKoutrasIII", 
	"PetrovSolution"
}

exactSolsData["Homogeneous"] = {
	"FarnsworthKerrI", 
	"FarnsworthKerrII", 
	"FarnsworthKerrIII", 
	"OsvathKoutrasI", 
	"OsvathKoutrasII", 
	"OsvathKoutrasIII"
}

exactSolsData["PerfectFluid"] = {
	"Friedmann", 
	"FarnsworthKerrI", 
	"FarnsworthKerrII", 
	"FarnsworthKerrIII", 
	"GeneralSphericalSymmetry",
	"LemaitreTolman",
	"OsvathKoutrasI", 
	"OsvathKoutrasII", 
	"OsvathKoutrasIII", 
	"Stephani",
	"StephaniThermodynamic",
	"StephaniThermodynamicSpherical"
}

exactSolsData["PetrovTypeD"] = {
	"FarnsworthKerrIII", 
	"GeneralSphericalSymmetry", 
	"Kerr", 
	"KerrNUT", 
	"LemaitreTolman", 
	"ReissnerNordstrom",
	"Schwarzschild"
}

exactSolsData["PetrovTypeI"] = {
	"OsvathKoutrasI", 
	"OsvathKoutrasII",
	"OsvathKoutrasIII",
	"PetrovSolution"
}

exactSolsData["PetrovTypeN"] = {
	"BertottiRobinsonSolution"
}

exactSolsData["SpatialG6"] = {
	"Friedmann", 
	"Stephani", 
	"StephaniThermodynamic", 
	"StephaniThermodynamicSpherical"
}

exactSolsData["SpatiallyHomogeneous"] = {
	"Friedmann"
}

exactSolsData["SphericalSymmetry"] = {
	"GeneralSphericalSymmetry",
	"LemaitreTolman", 
	"ReissnerNordstrom", 
	"Schwarzschild", 
	"StephaniThermodynamicSpherical"
}

exactSolsData["Static"] = {
	"ReissnerNordstrom", 
	"Schwarzschild"
}

exactSolsData["Stationary"] = {
	"Kerr"
}

exactSolsData["ThermodynamicPerfectFluid"] = {
	"Friedmann", 
	"GeneralSphericalSymmetry", 
	"LemaitreTolman", 
	"StephaniThermodynamic", 
	"StephaniThermodynamicSpherical"
}

exactSolsData["Vacuum"] = {
	"Kerr", 
	"KerrNUT", 
	"PetrovSolution",
	"Schwarzschild"
}

exactSolsData["VacuumTypeD"] = {
	"Kerr", 
	"KerrNUT",
	"PetrovSolution",
	"Schwarzschild"
}

exactSolsData["Warped22"] = {
	"GeneralSphericalSymmetry", 
	"LemaitreTolman", 
	"ReissnerNordstrom", 
	"StephaniThermodynamicSpherical"
}

exactSolsData["SphericalSymmetry"] = {
	"GeneralSphericalSymmetry",
	"Schwarzschild",
	"StephaniSpherical"
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

exactSolsData["BertottiRobinsonSolution", {"ComplexCoordinates", "ScalarFunctionValues"}] = {}

(* ::Subsection:: *)
(* General electrovac type D solution. See https://arxiv.org/abs/2407.14863v1 *)

exactSolsData["ElectroVacTypeD", "Classes"] = exactSolsData["Kerr", "Classes"]

exactSolsData["ElectroVacTypeD", "CoordinateSystems"] = {"TypeDCoordinates"}

exactSolsData["ElectroVacTypeD", "DefaultCoordinates"] = "TypeDCoordinates"
 
exactSolsData["ElectroVacTypeD", "IsIDEAL"] = True
 
exactSolsData["ElectroVacTypeD", "ParameterAssumptions"] = Null
 
exactSolsData["ElectroVacTypeD", "ParameterNames"] = {"m", "a", "l", "e", "g", "\[Alpha]", "\[Lambda]"}
 
exactSolsData["ElectroVacTypeD", {"TypeDCoordinates", "CoordinateAssumptions"}] = -Infinity < #[[1]] < Infinity && -Infinity < #[[2]] < Infinity && 0 < #[[3]] < Pi &

exactSolsData["ElectroVacTypeD", {"TypeDCoordinates", "CoordinateNames"}] = {"t", "q", "\[Theta]", "\[Phi]"}
 
exactSolsData["ElectroVacTypeD", {"TypeDCoordinates", "Metric"}] = 
	Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], q = coords[[2]], theta = coords[[3]], phi = coords[[4]],
				m = params[[1]], a = params[[2]], l = params[[3]], e = params[[4]], g = params[[5]], alpha = params[[6]], lambda = params[[7]], 
				P = scfuncs[[1]], Q = scfuncs[[2]], rho = scfuncs[[3]], Omega = scfuncs[[4]]}, 
    			Quiet[P = 
					Function[{theta}, 
						1 - 2 ((alpha a / (a^2 + l^2)) m - (lambda/3) l) (l + a Cos[theta]) + 
						((alpha^2 a^2 / (a^2 + l^2)^2) (a^2 - l^2 + e^2 + g^2) + lambda/3) (l + a Cos[theta])^2
					]
				]; 
    			Quiet[Q = 
					Function[{q}, 
						(1 - 2 m q + (a^2 - l^2 + e^2 + g^2) q^2) (q + alpha a (a - l) / (a^2 + l^2)) (q - alpha a (a + l) / (a^2 + l^2))
    					- (lambda / 3) (1 + 2 alpha a l (a^2 - l^2) / (a^2 + l^2) q + (a^2 + 3 l^2) q^2)
					]
				];
				Quiet[Omega = 
					Function[{q, theta}, 
						q - alpha a / (a^2 + l^2) (l + a Cos[theta])
					]
				];
				Quiet[rho = 
					Function[{q, theta}, 
						 1 + q^2 (l + a Cos[theta])^2
					]
				];
				{
					{
    					-Q[q] / (rho[q, theta]^2 Omega[q, theta]^2), 0, 0, 
    					-Q[q] / (rho[q, theta]^2 Omega[q, theta]^2) (a Sin[theta]^2 + 4 l Sin[theta/2]^2)
  					},
  					{
    					0, rho[q, theta]^2 / (Q[q] Omega[q, theta]^2), 0, 0
  					},
					{
    					0, 0, rho[q, theta]^2 / (P[theta] Omega[q, theta]^2), 0
  					},
  					{
    					-Q[q] / (rho[q, theta]^2 Omega[q, theta]^2) (a Sin[theta]^2 + 4 l Sin[theta/2]^2), 0, 0,  P[theta] Sin[theta]^2 / (rho[q, theta]^2 Omega[q, theta]^2) (1 + (a + l)^2 q^2)
  					}
				}
		]
	]
 
exactSolsData["ElectroVacTypeD", {"TypeDCoordinates", "ParameterAssumptions"}] = Null
 
exactSolsData["ElectroVacTypeD", {"TypeDCoordinates", "ParameterNames"}] = {"m", "a", "l", "e", "g", "\[Alpha]", "\[Lambda]"}
 
exactSolsData["ElectroVacTypeD", {"TypeDCoordinates", "ScalarFunctionNames"}] = {"P", "Q", "\[Rho]", "\[CapitalOmega]"}

exactSolsData["ElectroVacTypeD", {"TypeDCoordinates", "ScalarFunctionValues"}] =
Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], q = coords[[2]], theta = coords[[3]], phi = coords[[4]],
				m = params[[1]], a = params[[2]], l = params[[3]], e = params[[4]], g = params[[5]], alpha = params[[6]], lambda = params[[7]], 
				P = scfuncs[[1]], Q = scfuncs[[2]], rho = scfuncs[[3]], Omega = scfuncs[[4]]}, 
    			Quiet[P = 
					Function[{theta}, 
						1 - 2 ((alpha a / (a^2 + l^2)) m - (lambda/3) l) (l + a Cos[theta]) + 
						((alpha^2 a^2 / (a^2 + l^2)^2) (a^2 - l^2 + e^2 + g^2) + lambda/3) (l + a Cos[theta])^2
					]
				]; 
    			Quiet[Q = 
					Function[{q}, 
						(1 - 2 m q + (a^2 - l^2 + e^2 + g^2) q^2) (q + alpha a (a - l) / (a^2 + l^2)) (q - alpha a (a + l) / (a^2 + l^2))
    					- (lambda / 3) (1 + 2 alpha a l (a^2 - l^2) / (a^2 + l^2) q + (a^2 + 3 l^2) q^2)
					]
				];
				Quiet[Omega = 
					Function[{q, theta}, 
						q - alpha a / (a^2 + l^2) (l + a Cos[theta])
					]
				];
				Quiet[rho = 
					Function[{q, theta}, 
						 1 + q^2 (l + a Cos[theta])^2
					]
				];
				{
					P[theta],
					Q[q],
					Omega[q, theta],
					rho[q, theta]
				}
		]
]
	

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

exactSolsData["FarnsworthKerrI", {"CanonicalCoordinates", "ScalarFunctionValues"}] = {}

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

exactSolsData["FarnsworthKerrII", {"CanonicalCoordinates", "ScalarFunctionValues"}] = {}

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

exactSolsData["FarnsworthKerrII", {"CanonicalCoordinates", "ScalarFunctionValues"}] = {}

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

exactSolsData["Friedmann", {"ReducedCircumferencePolarCoordinates", "ScalarFunctionValues"}] =
    Function[{coords, params, scfuncs},
        With[{t = coords[[1]], r = coords[[2]], theta = coords[[3]], 
            phi = coords[[4]], k = params[[1]], R = scfuncs[[1]]},
				{
					R[t]
				}
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

exactSolsData["GeneralSphericalSymmetry", {"SphericalCoordinates", "ScalarFunctionValues"}] =
	Function[{coords, params, funcs},
		With[{nu = funcs[[3]], lambda = funcs[[1]], mu = funcs[[2]], t = coords[[1]], r = coords[[2]], theta = coords[[3]], phi = coords[[4]]},
				{
					lambda[t, r], 
					mu[t, r], 
					nu[t, r]
				}
			]
		]

(* ::Subsection:: *)
(* GeneralSzekeresSzafron metric in planar coordinates *)

exactSolsData["GeneralSzekeresSzafron", "Classes"] = {"PerfectFluid"}

exactSolsData["GeneralSzekeresSzafron", "CoordinateSystems"] = {"PlanarCoordinates"} 

exactSolsData["GeneralSzekeresSzafron", "DefaultCoordinates"] = "PlanarCoordinates"

exactSolsData["GeneralSzekeresSzafron", "ParameterNames"] = {}

exactSolsData["GeneralSzekeresSzafron", "ParameterAssumptions"] = Null

exactSolsData["GeneralSzekeresSzafron", "IsIDEAL"] = False

exactSolsData["GeneralSzekeresSzafron", {"PlanarCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}

exactSolsData["GeneralSzekeresSzafron", {"PlanarCoordinates", "CoordinateAssumptions"}] = 
	-Infinity < #[[1]] < Infinity && -Infinity < #[[2]] < Infinity && -Infinity < #[[3]] < Infinity && -Infinity < #[[4]] < Infinity &

exactSolsData["GeneralSzekeresSzafron", {"PlanarCoordinates", "ParameterNames"}] = exactSolsData["GeneralSzekeresSzafron", "ParameterNames"]

exactSolsData["GeneralSzekeresSzafron", {"PlanarCoordinates", "ParameterAssumptions"}] = exactSolsData["GeneralSzekeresSzafron", "ParameterAssumptions"]

exactSolsData["GeneralSzekeresSzafron", {"PlanarCoordinates", "ScalarFunctionNames"}] = {"\[Alpha]", "\[Beta]"}

defaultcoordinates["GeneralSzekeresSzafron"] = "PlanarCoordinates"

(* The syntax is exactSolsData[args__][{coords_List, parameters_List, functions_List}] *)

exactSolsData["GeneralSzekeresSzafron", {"PlanarCoordinates", "Metric"}] =
	Function[{coords, params, funcs},
		With[{alpha = funcs[[1]], beta = funcs[[2]], t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]]},
			DiagonalMatrix[
				{
					-1,  
					E^(2 beta[t, x, y, z]), 
					E^(2 beta[t, x, y, z]),
					E^(2 alpha[t, x, y, z])
				}
			]
		] 
	]

exactSolsData["GeneralSzekeresSzafron", {"PlanarCoordinates", "ScalarFunctionValues"}] =
	Function[{coords, params, funcs},
		With[{alpha = funcs[[1]], beta = funcs[[2]], t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]]},
				{
					alpha[t, x, y, z], 
					beta[t, x, y, z]
				}
			]
		]

(* ::Subsection:: *)
(* KantowskiSachs metric in spherical coordinates *)

exactSolsData["KantowskiSachs", "Classes"] = {"SphericalSymmetry", "Warped22"}

exactSolsData["KantowskiSachs", "CoordinateSystems"] = {"SphericalCoordinates"} 

exactSolsData["KantowskiSachs", "DefaultCoordinates"] = "SphericalCoordinates"

exactSolsData["KantowskiSachs", "ParameterNames"] = {}

exactSolsData["KantowskiSachs", "ParameterAssumptions"] = Null

exactSolsData["KantowskiSachs", "IsIDEAL"] = False

exactSolsData["KantowskiSachs", {"SphericalCoordinates", "CoordinateNames"}] = {"t", "r", "\[Theta]", "\[Phi]"}

exactSolsData["KantowskiSachs", {"SphericalCoordinates", "CoordinateAssumptions"}] = #[[2]] > 0 && Pi > #[[3]] > 0 &

exactSolsData["KantowskiSachs", {"SphericalCoordinates", "ParameterNames"}] = exactSolsData["KantowskiSachs", "ParameterNames"]

exactSolsData["KantowskiSachs", {"SphericalCoordinates", "ParameterAssumptions"}] = exactSolsData["KantowskiSachs", "ParameterAssumptions"]

exactSolsData["KantowskiSachs", {"SphericalCoordinates", "ScalarFunctionNames"}] = {"R", "S"}

defaultcoordinates["KantowskiSachs"] = "SphericalCoordinates"

(* The syntax is exactSolsData[args__][{coords_List, parameters_List, functions_List}] *)

exactSolsData["KantowskiSachs", {"SphericalCoordinates", "Metric"}] =
	Function[{coords, params, funcs},
		With[{R = funcs[[1]], S = funcs[[2]], t = coords[[1]], r = coords[[2]], theta = coords[[3]], phi = coords[[4]]},
			DiagonalMatrix[
				{
					-1, 
					R[t], 
					S[t], 
					S[t] Sin[theta]^2
				}
			]
		] 
	]

exactSolsData["KantowskiSachs", {"SphericalCoordinates", "ScalarFunctionValues"}] =
	Function[{coords, params, funcs},
		With[{R = funcs[[1]], S = funcs[[2]], t = coords[[1]], r = coords[[2]], theta = coords[[3]], phi = coords[[4]]},
				{
					R[t], 
					S[t]
				}
			]
		]

(* ::Subsection:: *)
(* KasnerI *)

exactSolsData["KasnerI", "Classes"] = {"Vacuum", "G3IXonS3"}

exactSolsData["KasnerI", "CoordinateSystems"] = {"CanonicalCoordinates"} 

exactSolsData["KasnerI", "DefaultCoordinates"] = "CanonicalCoordinates"

exactSolsData["KasnerI", "ParameterNames"] = {"p1", "p2", "p3"}

exactSolsData["KasnerI", "ParameterAssumptions"] =  
    Function[{coords, params, scfuncs}, 
    	With[{p1 = params[[1]], p2 = params[[2]], p3 = params[[3]]}, 
    		Element[p1, Reals] && Element[p2, Reals] && Element[p3, Reals] && 
			p1 + p2 + p3 == 1 && p1^2 + p2^2 + p3^2 == 1
		]
	]

exactSolsData["KasnerI", "IsIDEAL"] = False

exactSolsData["KasnerI", {"CanonicalCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}

exactSolsData["KasnerI", {"CanonicalCoordinates", "CoordinateAssumptions"}] = 
	-Infinity < #[[1]] < Infinity && -Infinity < #[[2]] < Infinity && -Infinity < #[[3]] < Infinity && -Infinity < #[[4]] < Infinity &

exactSolsData["KasnerI", {"CanonicalCoordinates", "ParameterNames"}] = exactSolsData["KasnerI", "ParameterNames"]

exactSolsData["KasnerI", {"CanonicalCoordinates", "ParameterAssumptions"}] = exactSolsData["KasnerI", "ParameterAssumptions"]

exactSolsData["KasnerI", {"CanonicalCoordinates", "ScalarFunctionNames"}] = {}

defaultcoordinates["KasnerI"] = "CanonicalCoordinates"

(* The syntax is exactSolsData[args__][{coords_List, parameters_List, functions_List}] *)

exactSolsData["KasnerI", {"CanonicalCoordinates", "Metric"}] =
	Function[{coords, params, funcs},
		With[{p1 = params[[1]], p2 = params[[2]], p3 = params[[3]], 
			t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]]},
			DiagonalMatrix[
				{
					-1,  
					t^(2p1), 
					t^(2p2),
					t^(2p3)
				}
			]
		] 
	]

exactSolsData["KasnerI", {"CanonicalCoordinates", "ScalarFunctionValues"}] = {}

(* ::Subsection:: *)
(* KasnerII *)

exactSolsData["KasnerII", "Classes"] = {"Vacuum", "G3IXonS3"}

exactSolsData["KasnerII", "CoordinateSystems"] = {"CanonicalCoordinates"} 

exactSolsData["KasnerII", "DefaultCoordinates"] = "CanonicalCoordinates"

exactSolsData["KasnerII", "ParameterNames"] = {"p1", "p2", "p3", "b", "epsilon"}

exactSolsData["KasnerII", "ParameterAssumptions"] =  
    Function[{coords, params, scfuncs}, 
    	With[{p1 = params[[1]], p2 = params[[2]], p3 = params[[3]]}, 
    		Element[p1, Reals] && Element[p2, Reals] && Element[p3, Reals] && 
			p1 + p2 + p3 == 1 && p1^2 + p2^2 + p3^2 == 1
		]
	]

exactSolsData["KasnerII", "IsIDEAL"] = False

exactSolsData["KasnerII", {"CanonicalCoordinates", "CoordinateNames"}] = {"x", "y", "z", "t"}

exactSolsData["KasnerII", {"CanonicalCoordinates", "CoordinateAssumptions"}] = 
	-Infinity < #[[1]] < Infinity && -Infinity < #[[2]] < Infinity && -Infinity < #[[3]] < Infinity && -Infinity < #[[4]] < Infinity &

exactSolsData["KasnerII", {"CanonicalCoordinates", "ParameterNames"}] = exactSolsData["KasnerII", "ParameterNames"]

exactSolsData["KasnerII", {"CanonicalCoordinates", "ParameterAssumptions"}] = exactSolsData["KasnerII", "ParameterAssumptions"]

exactSolsData["KasnerII", {"CanonicalCoordinates", "ScalarFunctionNames"}] = {"G"}

defaultcoordinates["KasnerII"] = "CanonicalCoordinates"

(* The syntax is exactSolsData[args__][{coords_List, parameters_List, functions_List}] *)

exactSolsData["KasnerII", {"CanonicalCoordinates", "Metric"}] =
	Function[{coords, params, funcs},
		With[{p1 = params[[1]], p2 = params[[2]], p3 = params[[3]], b = params[[4]], epsilon = params[[5]],
			t = coords[[4]], x = coords[[1]], y = coords[[2]], z = coords[[3]], G2 = funcs[[1]]},
			(* Define the auxiliary function G^2 *)
			Quiet[
				G2 = Function[{t, b, p1}, 1 + b^2 t^(4 p1)]
			];

			(* Metric in the order {x,y,z,t} *)
			{
  				{
    				epsilon t^(2 p1)/G2[t, b, p1],
    				4 epsilon p1 b z t^(2 p1)/G2[t, b, p1],
    				0,
    				0
  				},
  				{
    				4 epsilon p1 b z t^(2 p1)/G2[t, b, p1],
    				16 epsilon p1^2 b^2 z^2 t^(2 p1)/G2[t, b, p1] + G2[t, b, p1] t^(2 p2),
    				0,
    				0
  				},
  				{
    				0,
    				0,
    				G2[t, b, p1] t^(2 p3),
    				0
  				},
  				{
    				0,
    				0,
    				0,
    				-epsilon G2[t, b, p1]
  				}
			}

		] 
	]

exactSolsData["KasnerII", {"CanonicalCoordinates", "ScalarFunctionValues"}] =
	Function[{coords, params}, 
    	With[{t = coords[[1]], p1 = params[[1]], b = params[[4]]}, 
        		{
					1 + b^2*t^(4*p1)
				}
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
				Quiet[Sigma = 
					Function[{r, a, theta}, 
						r^2 + a^2*Cos[theta]^2
					]
				]; 
				r > 0 && Pi > theta > 0 && Element[Cos[theta], Reals] && 
					Sin[theta] > 0 && Sigma[r, a, theta] > 0
		]
	]
 
exactSolsData["Kerr", {"BoyerLindquistCoordinates", "CoordinateNames"}] = 
    {"t", "r", "\[Theta]", "\[Phi]"}
 
exactSolsData["Kerr", {"BoyerLindquistCoordinates", "Metric"}] = 
    Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], r = coords[[2]], theta = coords[[3]], phi = coords[[4]], 
			m = params[[1]], a = params[[2]], Sigma = scfuncs[[1]], Delta = scfuncs[[2]]}, 
    			Quiet[Sigma = 
					Function[{r, a, theta}, 
						r^2 + a^2*Cos[theta]^2
					]
				];
       			Quiet[Delta = 
					Function[{r, m, a}, 
						r^2 - 2*m*r + a^2
					]
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

exactSolsData["Kerr", {"BoyerLindquistCoordinates", "ScalarFunctionValues"}] =
	Function[{coords, params}, 
    	With[{t = coords[[1]], r = coords[[2]], theta = coords[[3]], phi = coords[[4]], 
			m = params[[1]], a = params[[2]]}, 
        		{
					r^2 + a^2*Cos[theta]^2,
					r^2 - 2*m*r + a^2
				}
		]
	]

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

				Quiet[alpha = 
					Function[{p, x, k, s}, 
						p*x^2 + k*(3 - k^2)*(x/(1 + k^2)^3) + s
					]
				];

				Quiet[beta = 
					Function[{p, y, k, s}, 
						(-p)*y^2 + (3*k^2 - 1)*(y/(1 + k^2)^3) + s
					]
				]; 
				x^2 + y^2 > 0 && alpha[p, x, k, s] > 0 && beta[p, y, k, s] > 0
		]
	]
 
exactSolsData["KerrNUT", {"ExpansionGradientAdaptedCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}
 
exactSolsData["KerrNUT", {"ExpansionGradientAdaptedCoordinates", "Metric"}] = 
    Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]], p = params[[1]], 
			k = params[[2]], s = params[[3]], alpha = scfuncs[[1]], beta = scfuncs[[2]]}, 
      			Quiet[alpha = 
					Function[{p, x, k, s}, 
						p*x^2 + k*(3 - k^2)*(x/(1 + k^2)^3) + s
					]
				]; 
				Quiet[beta = 
					Function[{p, y, k, s}, 
						(-p)*y^2 + (3*k^2 - 1)*(y/(1 + k^2)^3) + s
					]
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

exactSolsData["KerrNUT", {"ExpansionGradientAdaptedCoordinates", "ScalarFunctionValues"}] =
	Function[{coords, params}, 
    	With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]], p = params[[1]], 
			k = params[[2]], s = params[[3]]}, 
       			{
					p*x^2 + k*(3 - k^2)*(x/(1 + k^2)^3) + s,
					(-p)*y^2 + (3*k^2 - 1)*(y/(1 + k^2)^3) + s
				}
		]
	]

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

exactSolsData["LemaitreTolman", {"SphericalCoordinates", "ScalarFunctionValues"}] =
    Function[{coords, params, scfuncs},
        With[{t = coords[[1]], r = coords[[2]], theta = coords[[3]], 
            phi = coords[[4]], k = params[[1]], R = scfuncs[[1]]},
				{
					R[t, r]
				}
			]
        ]



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
      			Quiet[beta = 
					Function[{s}, 
						Sqrt[1 + 2*s^2*(1 - s^2)*(3 - s^2)]
					]
				]; 
       			Quiet[A = Function[{beta, s}, 
						(1 - beta[s])/2
					]
				]; 
       			Quiet[B = 
					Function[{beta, s}, 
						(1 + beta[s])/2]; F = Function[{s}, 1 - s^2
					]
				]; 
       			Quiet[b = 
					Function[{s}, 
						Sqrt[2]*s*(3 - s^2)
					]
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

exactSolsData["OsvathKoutrasI", {"GroupGeneratorsAdaptedCoordinates", "ScalarFunctionValues"}] = {}

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
				Quiet[s = Sqrt[1.2296814706969093]]; 
       			Quiet[F = Function[{s}, 1 - s^2]]; 
				Quiet[b = Function[{s}, Sqrt[2]*s*(3 - s^2)]]; 
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

exactSolsData["OsvathKoutrasII", {"GroupGeneratorsAdaptedCoordinates", "ScalarFunctionValues"}] = {}

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
				Quiet[F = Function[{s}, 1 - s^2]]; 
       			Quiet[b = Function[{s}, Sqrt[2]*s*(3 - s^2)]]; 
       			Quiet[k = Function[{beta2, s}, Sqrt[-beta2]/2]]; 
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

exactSolsData["OsvathKoutrasIII", {"GroupGeneratorsAdaptedCoordinates", "ScalarFunctionValues"}] = {}

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

exactSolsData["PetrovSolution", {"CanonicalCoordinates", "ScalarFunctionValues"}] = {}

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

exactSolsData["ReissnerNordstrom", {"SchwarzschildCoordinates", "ScalarFunctionValues"}] = {}

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

exactSolsData["Schwarzschild", "CoordinateSystems"] = {"SchwarzschildCoordinates", "IsotropicCoordinates", "HarmonicCoordinates"}

exactSolsData["Schwarzschild", "DefaultCoordinates"] = "SchwarzschildCoordinates"

exactSolsData["Schwarzschild", "ParameterNames"] = {"m"}

exactSolsData["Schwarzschild", "ParameterAssumptions"] = #[[1]] > 0 &

exactSolsData["Schwarzschild", "IsIDEAL"] = True

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "CoordinateNames"}] = {"t", "r", "\[Theta]", "\[Phi]"}

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "CoordinateAssumptions"}] = #[[2]] > 0 && Pi > #[[3]] > 0 &

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "ParameterNames"}] = exactSolsData["Schwarzschild", "ParameterNames"]

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "ParameterAssumptions"}] = exactSolsData["Schwarzschild", "ParameterAssumptions"]

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "ScalarFunctionNames"}] = {}

exactSolsData["Schwarzschild", {"SchwarzschildCoordinates", "ScalarFunctionValues"}] = {}

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
            Quiet[R =
                Function[{x, y, z},
                    Sqrt[x^2 + y^2 + z^2]
                ]
			];
            m >= 0 && R[x, y, z] > 2 * m
        ]
    ]

exactSolsData["Schwarzschild", {"IsotropicCoordinates", "Metric"}] :=
  Function[
      DiagonalMatrix[{
        -((1 - Slot[2][[1]]/(2 Slot[3][[1]][Slot[1][[2]], Slot[1][[3]], Slot[1][[4]]]))^2/(1 + Slot[2][[1]]/(2 Slot[3][[1]][Slot[1][[2]], Slot[1][[3]], Slot[1][[4]]]))^2),
         (1 + Slot[2][[1]]/(2 Slot[3][[1]][Slot[1][[2]], Slot[1][[3]], Slot[1][[4]]]))^4,
         (1 + Slot[2][[1]]/(2 Slot[3][[1]][Slot[1][[2]], Slot[1][[3]], Slot[1][[4]]]))^4,
         (1 + Slot[2][[1]]/(2 Slot[3][[1]][Slot[1][[2]], Slot[1][[3]], Slot[1][[4]]]))^4
      }]
    ];

exactSolsData["Schwarzschild", {"IsotropicCoordinates", "ScalarFunctionValues"}] =
    Function[
			{Sqrt[Slot[1][[2]]^2 + Slot[1][[3]]^2 + Slot[1][[4]]^2]}
    ]

	(* ::Subsection:: *)
(* Schwarzschild in harmonic coordinates *)

exactSolsData["Schwarzschild", {"HarmonicCoordinates", "ParameterNames"}] = {"m"}

exactSolsData["Schwarzschild", {"HarmonicCoordinates", "ParameterAssumptions"}] =
    Function[{coords, params, scfuncs},
        With[{m = params[[1]]},
            m > 0
        ]
    ]

exactSolsData["Schwarzschild", {"HarmonicCoordinates", "ScalarFunctionNames"}] = {"R"}

exactSolsData["Schwarzschild", {"HarmonicCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}

exactSolsData["Schwarzschild", {"HarmonicCoordinates", "CoordinateAssumptions"}] =
    Function[{coords, params, scfuncs},
        With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = 
            coords[[4]], m = params[[1]], R = scfuncs[[1]]},
            Quiet[R =
                Function[{x, y, z},
                    Sqrt[x^2 + y^2 + z^2]
                ]
			];
            m >= 0 && R[x, y, z] >  m
        ]
    ]

exactSolsData["Schwarzschild", {"HarmonicCoordinates", "Metric"}] =
    Function[{coords, params, scfuncs},
        With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = 
            coords[[4]], m = params[[1]], r = scfuncs[[1]]},
            Quiet[r =
                Function[{x, y, z},
                    Sqrt[x^2 + y^2 + z^2]
                ]
			];
			{
				{(1 - m/r[x, y, z])/(1 + m/r[x, y, z]), 0, 0, 0}, 
				{	
					0, 
					-(1 + m/r[x, y, z])^2 - (m^2 x^2 (1 + m / r[x, y, z]))/((1 - m / r[x, y, z]) r[x, y, z]^4), 
					-((m^2 x y (1 + m/r[x, y, z]))/((1 - m/r[x, y, z]) r[x, y, z]^4)), 
					-((m^2 x z (1 + m/r[x, y, z]))/((1 - m/r[x, y, z]) r[x, y, z]^4))
				}, 
				{	
					0, 
					-((m^2 x y (1 + m/r[x, y, z]))/((1 - m/r[x, y, z]) r[x, y, z]^4)), 
					-(1 + m/r[x, y, z])^2 - (m^2 y^2 (1 + m/r[x, y, z]))/((1 - m/r[x, y, z]) r[x, y, z]^4), 
					-((m^2 y z (1 + m/r[x, y, z]))/((1 - m/r[x, y, z]) r[x, y, z]^4))
				}, 
				{	
					0, 
					-((m^2 x z (1 + m/r[x, y, z]))/((1 - m/r[x, y, z]) r[x, y, z]^4)), 
					-((m^2 y z (1 + m/r[x, y, z]))/((1 - m/r[x, y, z]) r[x, y, z]^4)), 
					-(1 + m/r[x, y, z])^2 - (m^2 z^2 (1 + m/r[x, y, z]))/((1 - m/r[x, y, z]) r[x, y, z]^4)
				}
   			}
        ]
    ]

exactSolsData["Schwarzschild", {"HarmonicCoordinates", "ScalarFunctionValues"}] =
    Function[{coords, params},
        With[{t = coords[[1]], x = coords[[2]], y = coords[[c3]], z = coords[[4]], m = params[[1]]},
			{Sqrt[x^2 + y^2 + z^2]}
        ]
    ]

(* ::Subsection:: *)
(* Stephani in adapted coordinates *)

exactSolsData["Stephani", "Classes"] = {"PerfectFluid", "SpatialG6", "ConformallyFlat"}

exactSolsData["Stephani", "CoordinateSystems"] = {"AdaptedCoordinates"}

exactSolsData["Stephani", "DefaultCoordinates"] = "AdaptedCoordinates"
 
exactSolsData["Stephani", "IsIDEAL"] = True
 
exactSolsData["Stephani", "ParameterAssumptions"] = exactSolsData["Stephani", {"AdaptedCoordinates", "ParameterAssumptions"}]
 
exactSolsData["Stephani", "ParameterNames"] = exactSolsData["Stephani", {"AdaptedCoordinates", "ParameterNames"}]
 
exactSolsData["Stephani", {"AdaptedCoordinates", "CoordinateAssumptions"}] = Null
 
exactSolsData["Stephani", {"AdaptedCoordinates", "CoordinateNames"}] = {"t", "x", "y", "z"}
 
exactSolsData["Stephani", {"AdaptedCoordinates", "Metric"}] = 
    Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]], Omega = scfuncs[[1]], alpha = scfuncs[[2]], 
    		R = scfuncs[[3]], b1 = scfuncs[[4]], b2 = scfuncs[[5]], b3 = scfuncs[[6]], k = scfuncs[[7]]}, 
      			Quiet[Omega = 
					Function[{R, t, b1, b2, b3, x, y, z, k}, 
         				R[t]/(1 + 2*(b1[t]*x + b2[t]*y + b3[t]*z) + k[t]*((x^2 + y^2 + z^2)/4))
					]
				]; 
       			Quiet[alpha = 
					Function[{R, t, b1, b2, b3, x, y, z, k}, 
        				R[t]*(D[Omega[R, t, b1, b2, b3, x, y, z, k], t]/(Omega[R, t, b1, b2, b3, x, y, z, k]*D[R[t], t]))
					]
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

exactSolsData["Stephani", {"AdaptedCoordinates", "ScalarFunctionValues"}] =
    Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = coords[[4]], Omega = scfuncs[[1]], alpha = scfuncs[[2]], 
    		R = scfuncs[[3]], b1 = scfuncs[[4]], b2 = scfuncs[[5]], b3 = scfuncs[[6]], k = scfuncs[[7]]},
			Quiet[Omega = 
					Function[{R, t, b1, b2, b3, x, y, z, k}, 
         				R[t]/(1 + 2*(b1[t]*x + b2[t]*y + b3[t]*z) + k[t]*((x^2 + y^2 + z^2)/4))
					]
				]; 
       			Quiet[alpha = 
					Function[{R, t, b1, b2, b3, x, y, z, k}, 
        				R[t]*(D[Omega[R, t, b1, b2, b3, x, y, z, k], t]/(Omega[R, t, b1, b2, b3, x, y, z, k]*D[R[t], t]))
					]
				];
				{
					Omega[R, t, b1, b2, b3, x, y, z, k],
					alpha[R, t, b1, b2, b3, x, y, z, k],
					R[t],
					b1[t], 
					b2[t], 
					b3[t], 
					k[t]	
				}
		]
	]

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
            Quiet[w =
                Function[{x, y, z},
                    2 * (z / (1 + epsilon * ((x^2 + y^2 + z^2) / 4)))            
                ]
			];
            Quiet[L =
                Function[{R, b, t, w, x, y, z},
                    R[t] / (1 + b[t] * w[x, y, z])
                ]
			];
            Quiet[Omega =
                Function[{w, x, y, z, L, R, b, t},
                    w[x, y, z] * (L[R, b, t, w, x, y, z] / (2 * z))
                ]
			];
            Quiet[alpha =
                Function[{R, t, L, b, w, x, y, z},
                    R[t] * (D[L[R, b, t, w, x, y, z], t] / (L[R, b, t, w, x, y, z] * D[R[t], t]))
                ]
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

exactSolsData["StephaniThermodynamic", {"AdaptedCoordinates", "ScalarFunctionValues"}] =
    Function[{coords, params, scfuncs},
        With[{t = coords[[1]], x = coords[[2]], y = coords[[3]], z = 
            coords[[4]], epsilon = params[[1]], w = scfuncs[[1]], L = scfuncs[[2]],
            Omega = scfuncs[[3]], alpha = scfuncs[[4]], R = scfuncs[[5]], b = scfuncs[[6]]},
            Quiet[w =
                Function[{x, y, z},
                    2 * (z / (1 + epsilon * ((x^2 + y^2 + z^2) / 4)))            
                ]
			];
            Quiet[L =
                Function[{R, b, t, w, x, y, z},
                    R[t] / (1 + b[t] * w[x, y, z])
                ]
			];
            Quiet[Omega =
                Function[{w, x, y, z, L, R, b, t},
                    w[x, y, z] * (L[R, b, t, w, x, y, z] / (2 * z))
                ]
			];
            Quiet[alpha =
                Function[{R, t, L, b, w, x, y, z},
                    R[t] * (D[L[R, b, t, w, x, y, z], t] / (L[R, b, t, w, x, y, z] * D[R[t], t]))
                ]
			];
           {
			w[x, y, z],
			L[R, b, t, w, x, y, z],
			Omega[w, x, y, z, L, R, b, t],
			alpha[R, t, L, b, w, x, y, z],
			R[t],
			b[r]
		   }
        ]
    ]

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
    			Quiet[Omega = 
					Function[{R, t, k, r}, 
						R[t]/(1 + k[t]*(r^2/4))
					]
				]; 
    			Quiet[alpha = 
					Function[{R, t, Omega, k, r}, 
						R[t]*(D[Omega[R, t, k, r], t]/(Omega[R, t, k, r]*D[R[t], t]))
					]
				];
				DiagonalMatrix[
       				{
						-alpha[R, t, Omega, k, r]^2, 
						Omega[R, t, k, r]^2, 
        				Omega[R, t, k, r]^2*r^2, 
						Omega[R, t, k, r]^2*r^2*Sin[theta]^2
					}
				] 
		]
	]
 
exactSolsData["StephaniThermodynamicSpherical", {"SphericalCoordinates", "ParameterAssumptions"}] = Null
 
exactSolsData["StephaniThermodynamicSpherical", {"SphericalCoordinates", "ParameterNames"}] = {}
 
exactSolsData["StephaniThermodynamicSpherical", {"SphericalCoordinates", "ScalarFunctionNames"}] = {"\[CapitalOmega]", "\[Alpha]", "R", "k"}

exactSolsData["StephaniThermodynamicSpherical", {"SphericalCoordinates", "ScalarFunctionValues"}] =
	Function[{coords, params, scfuncs}, 
    	With[{t = coords[[1]], r = coords[[2]], theta = coords[[3]], phi = coords[[4]], 
			Omega = scfuncs[[1]], alpha = scfuncs[[2]], R = scfuncs[[3]], k = scfuncs[[4]]}, 
    			Quiet[Omega = 
					Function[{R, t, k, r}, 
						R[t]/(1 + k[t]*(r^2/4))
					]
				]; 
    			Quiet[alpha = 
					Function[{R, t, Omega, k, r}, 
						R[t]*(D[Omega[R, t, k, r], t]/(Omega[R, t, k, r]*D[R[t], t]))
					]
				];
				{
					Omega[R, t, k, r],
					alpha[R, t, Omega, k, r],
					R[t],
					k[t]
				}
		]
	]

(* ::Section:: *)
(* GRExactSolsData *)

(* ::Subsection:: *)
(* Metrics: default coordinates *)

iGRExactSolsData[metric_?metricQ] := exactSolsData[metric, "CoordinateSystems"]

iGRExactSolsData[metric_?metricQ, "CoordinateSystems"] := exactSolsData[metric, "CoordinateSystems"]

iGRExactSolsData[metric_?metricQ, "Classes"] := exactSolsData[metric, "Classes"]

iGRExactSolsData[metric_?metricQ, coords_?coordinatesystemQ] := exactSolsData[metric, {coords, "Metric"}]

(* The syntax is GRData[args__String, {coords_List, parameters_List, functions_List}] *)

iGRExactSolsData[metric_?metricQ, coordname_?coordinatesystemQ, {coords_List, parameters_List, functions_List}] := exactSolsData[metric, {coordname, "Metric"}][coords, parameters, functions]

iGRExactSolsData[metric_?metricQ, "Properties"] := allmetricproperties

iGRExactSolsData[metric_?metricQ, "IsIDEAL"] := exactSolsData[metric, "IsIDEAL"]

iGRExactSolsData[metric_?metricQ, "CoordinateAssumptions"] := Module[{coords},
	coords = exactSolsData[metric, "DefaultCoordinates"];
	exactSolsData[metric, {coords, "CoordinateAssumptions"}]
]

iGRExactSolsData[metric_?metricQ, "CoordinateSystemName"] := Module[{coords},
	coords = exactSolsData[metric, "DefaultCoordinates"]
]

iGRExactSolsData[metric_?metricQ, "CoordinateNames"] := Module[{coords},
	coords = exactSolsData[metric, "DefaultCoordinates"];
	exactSolsData[metric, {coords, "CoordinateNames"}]
]

iGRExactSolsData[metric_?metricQ, "ParameterAssumptions"] := exactSolsData[metric, "ParameterAssumptions"]

iGRExactSolsData[metric_?metricQ, "ParameterNames"] := exactSolsData[metric, "ParameterNames"]

iGRExactSolsData[metric_?metricQ, "ScalarFunctionNames"] := Module[{coords},
	coords = exactSolsData[metric, "DefaultCoordinates"];
	exactSolsData[metric, {coords, "ScalarFunctionNames"}]
]

iGRExactSolsData[metric_?metricQ, "ScalarFunctionValues"] := Module[{coords},
	coords = exactSolsData[metric, "DefaultCoordinates"];
	exactSolsData[metric, {coords, "ScalarFunctionValues"}]
]

iGRExactSolsData[metric_?metricQ, "Metric"] := Module[{coords},
	coords = exactSolsData[metric, "DefaultCoordinates"];
	exactSolsData[metric, {coords, "Metric"}]
]

(* ::Subsection:: *)
(* Metrics: user given coordinates *)

(* Available properties for a given coordinate system *)
allcoordinateproperties = {
	"CoordinateAssumptions", 
	"CoordinateNames",
	(* TODO: it should be "MetricTensor" *)
	"Metric",
	"ParameterNames",
	"ParameterAssumptions",
	"ScalarFunctionNames",
	"ScalarFunctionValues"
	}

(* Valid coordinate systems *)
Set[coordinatepropertyQ[#], True]& /@ allcoordinateproperties;
coordinatesystemQ[_] := False;


iGRExactSolsData[metric_?metricQ, {coordname_?coordinatesystemQ, "Properties"}] := allcoordinateproperties

iGRExactSolsData[metric_?metricQ, {coordname_?coordinatesystemQ, "CoordinateAssumptions"}] := exactSolsData[metric, {coordname, "CoordinateAssumptions"}]
	 
iGRExactSolsData[metric_?metricQ, {coordname_?coordinatesystemQ, "CoordinateNames"}] := exactSolsData[metric, {coordname, "CoordinateNames"}]

iGRExactSolsData[metric_?metricQ, {coordname_?coordinatesystemQ, "Metric"}] := exactSolsData[metric, {coordname, "Metric"}]

iGRExactSolsData[metric_?metricQ, {coordname_?coordinatesystemQ, "ParameterNames"}] := exactSolsData[metric, {coordname, "ParameterNames"}]

iGRExactSolsData[metric_?metricQ, {coordname_?coordinatesystemQ, "ParameterAssumptions"}] := exactSolsData[metric, {coordname, "ParameterAssumptions"}]

iGRExactSolsData[metric_?metricQ, {coordname_?coordinatesystemQ, "ScalarFunctionNames"}] := exactSolsData[metric, {coordname, "ScalarFunctionNames"}]

iGRExactSolsData[metric_?metricQ, {coordname_?coordinatesystemQ, "ScalarFunctionValues"}] := exactSolsData[metric, {coordname, "ScalarFunctionValues"}]

(* ::Subsection:: *)
(* General definitions *)
iGRExactSolsData[] = allmetrics;
iGRExactSolsData["ExactSolutions"] = allmetrics;
iGRExactSolsData[All] = allmetrics;
iGRExactSolsData["Classes"] = allclasses;
iGRExactSolsData[All, "Classes"] = allclasses;
iGRExactSolsData["Properties"] = allmetricproperties;
iGRExactSolsData[All, "Properties"] = allmetricproperties;
iGRExactSolsData["CoordinateSystems"] = allcoordinatesystems;
iGRExactSolsData[All, "CoordinateSystems"] = allcoordinatesystems;

(* Metrics contained in a given class *)
iGRExactSolsData[class_?exactsolclassQ] := exactSolsData[class]

iGRExactSolsData[___] := $Failed;

(* Entry point of GRExactSolsData *)
GRExactSolsData[args___]:= Module[{res = iGRExactSolsData[args]}, If[UnsameQ[res, $Failed], res, Message[GRExactSolsData::noprop]]]

(****************************************************************)

(****************** 5. End private and package ******************)

(****************************************************************)


End[]

EndPackage[]