BeginPackage["xAct`xIdeal`ExactSolsData`"]

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
	"G3S2",
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
	"ExactSolutionName"
}

allcoordinatesystems = {
	"AdaptedCoordinates",
	"BoyerLindquistCoordinates",
	"CanonicalCoordinates",
	"ExpansionGradientAdaptedCoordinates",
 	"IsotropicCoordinates",
  	"ReducedCircumferencePolarCoordinates",
 	"SchwarzschildCoordinates",
	"SphericalCoordinates"
}

exactSolsData["Solutions"] = allmetrics

exactSolsData["Classes"] = allclasses

exactSolsData["MetricProperties"] = allmetricproperties

exactSolsData["CoordinateSystems"] = allcoordinatesystems


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

xAct`xIdeal`GRData[solname_String, coordname_String, "Metric", {coords_List, parameters_List, functions_List}] := exactSolsData[solname, {coordname, "Metric"}][coords, parameters, functions]


(* ::Subsection:: *)
(* GeneralSpherical metric in spherical coordinates *)

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
(* Schwarzschild in Schwarzschild coordinates *)

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


(****************************************************************)

(****************** 5. End private and package ******************)

(****************************************************************)


End[]

EndPackage[]