#!/usr/bin/env python3
"""
setup_asm_case.py -- OpenFOAM 13 case for a 3-D assembly (half domain) meshed
with snappyHexMesh from the STL patches written by build_fluid.py.

Same physics as the 2-D Stage 1 cases (water 25 C, k-omega SST, gravity along
the Venturi flow via buoyancyForce, fixed flow rate in, atmospheric total
pressure out) except near-wall treatment: wall functions (nutkWallFunction,
omegaWallFunction), because a wall-resolved 3-D mesh does not fit in memory.
y+ is written and reported.

Usage: python3 setup_asm_case.py <case dir> <fluid dir with fluid.stl + assembly_info.json> <Q L/min> [level]
level: coarse | medium | fine  (background cell 2.2 / 1.6 / 1.15 mm)
"""
import json, math, os, shutil, sys

NU, RHO, G = 8.9e-7, 997.0, 9.81
TI, MIX_FRAC = 0.05, 0.07
END_ITER = int(os.environ.get("ASM_ITERS", "3000"))   # iterations (env ASM_ITERS overrides)
WRITE_INT = 500 if END_ITER % 500 == 0 else END_ITER   # always write the last iteration
LEVELS = {"coarse": 2.2e-3, "medium": 1.6e-3, "fine": 1.15e-3}

case, fluid_dir, q_lpm = sys.argv[1], sys.argv[2], float(sys.argv[3])
level = sys.argv[4] if len(sys.argv) > 4 else "medium"
info = json.load(open(os.path.join(fluid_dir, "assembly_info.json")))
if os.path.isdir(case):
    shutil.rmtree(case)
for d in ("0", "constant/geometry", "system"):
    os.makedirs(os.path.join(case, d))
shutil.copy(os.path.join(fluid_dir, "fluid.stl"), os.path.join(case, "constant/geometry/fluid.stl"))
shutil.copy(os.path.join(fluid_dir, "assembly_info.json"), os.path.join(case, "assembly_info.json"))


def hdr(cls, loc, obj):
    return ("FoamFile\n{\n    format      ascii;\n    class       %s;\n    location    \"%s\";\n"
            "    object      %s;\n}\n\n" % (cls, loc, obj))


def write(rel, text):
    with open(os.path.join(case, rel), "w", newline="\n") as fh:
        fh.write(text)


Q_half = q_lpm / 60000.0 / 2.0
R_in = info["inlet_radius"] * 1e-3
U_in = (q_lpm / 60000.0) / (math.pi * R_in ** 2)
k_in = 1.5 * (TI * U_in) ** 2
L_mix = MIX_FRAC * 2 * R_in
om_in = math.sqrt(k_in) / (0.09 ** 0.25 * L_mix)

walls = ["wall_esig", "wall_elbow", "wall_venturi"]


def field(name, cls, dims, internal, bcs):
    body = hdr(cls, "0", name) + f"dimensions      {dims};\n\ninternalField   {internal};\n\nboundaryField\n{{\n"
    for patch, lines in bcs.items():
        body += f"    {patch}\n    {{\n" + "".join(f"        {l}\n" for l in lines) + "    }\n"
    body += "}\n"
    write(f"0/{name}", body)


def with_walls(d, wall_lines):
    d = dict(d)
    d['"wall_.*"'] = wall_lines
    d["symmetry"] = ["type            symmetry;"]
    return d


field("U", "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)", with_walls({
    "inlet": ["type            flowRateInletVelocity;",
              f"volumetricFlowRate {Q_half:.8g};   // {q_lpm:g} L/min, half domain",
              "value           uniform (0 0 0);"],
    "outlet": ["type            pressureInletOutletVelocity;", "value           uniform (0 0 0);"]},
    ["type            noSlip;"]))
field("p", "volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0", with_walls({
    "inlet": ["type            fixedFluxPressure;", "value           uniform 0;"],
    "outlet": ["type            totalPressure;", "p0              uniform 0;", "value           uniform 0;"]},
    ["type            fixedFluxPressure;", "value           uniform 0;"]))
field("k", "volScalarField", "[0 2 -2 0 0 0 0]", f"uniform {k_in:.6g}", with_walls({
    "inlet": ["type            turbulentIntensityKineticEnergyInlet;", f"intensity       {TI};",
              f"value           uniform {k_in:.6g};"],
    "outlet": ["type            inletOutlet;", f"inletValue      uniform {k_in:.6g};", f"value           uniform {k_in:.6g};"]},
    ["type            kqRWallFunction;", f"value           uniform {k_in:.6g};"]))
field("omega", "volScalarField", "[0 0 -1 0 0 0 0]", f"uniform {om_in:.6g}", with_walls({
    "inlet": ["type            turbulentMixingLengthFrequencyInlet;", f"mixingLength    {L_mix:.6g};",
              f"value           uniform {om_in:.6g};"],
    "outlet": ["type            inletOutlet;", f"inletValue      uniform {om_in:.6g};", f"value           uniform {om_in:.6g};"]},
    ["type            omegaWallFunction;", f"value           uniform {om_in:.6g};"]))
field("nut", "volScalarField", "[0 2 -1 0 0 0 0]", "uniform 0", with_walls({
    "inlet": ["type            calculated;", "value           uniform 0;"],
    "outlet": ["type            calculated;", "value           uniform 0;"]},
    ["type            nutkWallFunction;", "value           uniform 0;"]))

write("constant/physicalProperties", hdr("dictionary", "constant", "physicalProperties") +
      f"viscosityModel  constant;\n\nnu              {NU};   // water 25 C (rho = {RHO:g})\n")
write("constant/momentumTransport", hdr("dictionary", "constant", "momentumTransport") +
      "simulationType  RAS;\n\nRAS\n{\n    model           kOmegaSST;\n    turbulence      on;\n    viscosityModel  Newtonian;\n}\n")
write("constant/g", hdr("uniformDimensionedVectorField", "constant", "g") +
      f"dimensions      [0 1 -2 0 0 0 0];\nvalue           (-{G} 0 0);   // -x is down (Venturi flow direction)\n")
write("constant/fvModels", hdr("dictionary", "constant", "fvModels") +
      "buoyancyForce\n{\n    type            buoyancyForce;\n    U               U;\n}\n")

# ---------------------------------------------------------------- mesh dictionaries
h = LEVELS[level]
xmin, xmax = info["outlet_x"] * 1e-3 - 3e-3, 40e-3   # E-Sig spigot reaches x = +36 mm when flat side is down
ymin, ymax = -16e-3, info["inlet_y"] * 1e-3 + 3e-3
zmin, zmax = -1.0e-3, 25e-3
n = [max(1, round((xmax - xmin) / h)), max(1, round((ymax - ymin) / h)), max(1, round((zmax - zmin) / h))]
write("system/blockMeshDict", hdr("dictionary", "system", "blockMeshDict") + f"""\
convertToMeters 1;
vertices
(
    ({xmin} {ymin} {zmin}) ({xmax} {ymin} {zmin}) ({xmax} {ymax} {zmin}) ({xmin} {ymax} {zmin})
    ({xmin} {ymin} {zmax}) ({xmax} {ymin} {zmax}) ({xmax} {ymax} {zmax}) ({xmin} {ymax} {zmax})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({n[0]} {n[1]} {n[2]}) simpleGrading (1 1 1) );
boundary ( background {{ type patch; faces ( (0 3 2 1) (4 5 6 7) (0 1 5 4) (2 3 7 6) (0 4 7 3) (1 2 6 5) ); }} );
""")
write("system/surfaceFeaturesDict", hdr("dictionary", "system", "surfaceFeaturesDict") +
      "surfaces (\"fluid.stl\");\nincludedAngle   150;\nsubsetFeatures { nonManifoldEdges yes; openEdges yes; }\n")
vx = info["venturi_inlet_x"] * 1e-3
loc = (0.0, 0.10, 0.005)   # inside the E-Sig small bore
write("system/snappyHexMeshDict", hdr("dictionary", "system", "snappyHexMeshDict") + f"""\
castellatedMesh true;
snap            true;
addLayers       true;

geometry
{{
    fluid
    {{
        type triSurfaceMesh;
        file "fluid.stl";
        regions
        {{
            inlet        {{ name inlet; }}
            outlet       {{ name outlet; }}
            symmetry     {{ name symmetry; }}
            wall_esig    {{ name wall_esig; }}
            wall_elbow   {{ name wall_elbow; }}
            wall_venturi {{ name wall_venturi; }}
        }}
    }}
    fineZone
    {{
        type box;
        min ({info["outlet_x"] * 1e-3 - 1e-3} -0.016 -0.001);
        max (0.014 {info["esig_end_face_y"] * 1e-3 + 0.012} 0.025);
    }}
}}

castellatedMeshControls
{{
    maxLocalCells 4000000;
    maxGlobalCells 8000000;
    minRefinementCells 10;
    maxLoadUnbalance 0.10;
    nCellsBetweenLevels 3;
    features ( {{ file "fluid.eMesh"; level 2; }} );
    refinementSurfaces
    {{
        fluid
        {{
            level (1 1);
            regions
            {{
                inlet        {{ level (1 1); patchInfo {{ type patch; }} }}
                outlet       {{ level (2 2); patchInfo {{ type patch; }} }}
                symmetry     {{ level (0 0); patchInfo {{ type symmetry; }} }}
                wall_esig    {{ level (1 1); patchInfo {{ type wall; }} }}
                wall_elbow   {{ level (2 2); patchInfo {{ type wall; }} }}
                wall_venturi {{ level (2 2); patchInfo {{ type wall; }} }}
            }}
        }}
    }}
    resolveFeatureAngle 30;
    refinementRegions
    {{
        fineZone {{ mode inside; level 2; }}
    }}
    insidePoint ({loc[0]} {loc[1]} {loc[2]});
    allowFreeStandingZoneFaces true;
}}

snapControls
{{
    nSmoothPatch 5;
    tolerance 2.0;
    nSolveIter 100;
    nRelaxIter 8;
    nFeatureSnapIter 15;
    implicitFeatureSnap false;
    explicitFeatureSnap true;
    multiRegionFeatureSnap true;
}}

addLayersControls
{{
    relativeSizes true;
    layers
    {{
        "wall_.*" {{ nSurfaceLayers 3; }}
    }}
    expansionRatio 1.2;
    finalLayerThickness 0.45;
    minThickness 0.1;
    nGrow 0;
    featureAngle 120;
    slipFeatureAngle 30;
    nRelaxIter 5;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedialAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
}}

meshQualityControls
{{
    maxNonOrtho 65;
    maxBoundarySkewness 4;
    maxInternalSkewness 4;
    maxConcave 80;
    minVol 1e-16;
    minTetQuality 1e-15;
    minArea -1;
    minTwist 0.02;
    minDeterminant 0.001;
    minFaceWeight 0.05;
    minVolRatio 0.01;
    minTriangleTwist -1;
    nSmoothScale 4;
    errorReduction 0.75;
    relaxed {{ maxNonOrtho 70; maxInternalSkewness 4; maxBoundarySkewness 4; }}
}}

mergeTolerance 1e-6;
""")

# ---------------------------------------------------------------- solver dictionaries
mm = 1e-3
port_x = info["port_x"] * mm
# total pressure sections: inlet face of the E-Sig spigot, E-Sig end, Venturi inlet, throat, outlet
write("system/controlDict", hdr("dictionary", "system", "controlDict") + f"""\
solver          incompressibleFluid;

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {END_ITER};
deltaT          1;

writeControl    timeStep;
writeInterval   {WRITE_INT};
purgeWrite      1;
writeFormat     binary;
writePrecision  10;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

functions
{{
    #includeFunc residuals(p, U, k, omega)
    #includeFunc patchAverage(name=inletAvg, patch=inlet, fields=(p U))
    #includeFunc patchFlowRate(name=inletFlow, patch=inlet)
    #includeFunc patchFlowRate(name=outletFlow, patch=outlet)
}}
""")
write("system/fvSchemes", hdr("dictionary", "system", "fvSchemes") + """\
ddtSchemes      { default steadyState; }
gradSchemes
{
    default         Gauss linear;
    grad(U)         cellLimited Gauss linear 1;
    grad(k)         cellLimited Gauss linear 1;
    grad(omega)     cellLimited Gauss linear 1;
}
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss limitedLinear 1;
    div(phi,omega)  bounded Gauss limitedLinear 1;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear limited corrected 0.5; }
interpolationSchemes { default linear; }
snGradSchemes   { default limited corrected 0.5; }
wallDist        { method meshWave; }
""")
write("system/fvSolution", hdr("dictionary", "system", "fvSolution") + """\
solvers
{
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-08;
        relTol          0.05;
    }
    "(U|k|omega)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-10;
        relTol          0.1;
    }
}
SIMPLE
{
    nNonOrthogonalCorrectors 1;
    consistent      no;
    residualControl
    {
        p               1e-5;
        U               1e-6;
        "(k|omega)"     1e-6;
    }
}
relaxationFactors
{
    fields          { p 0.3; }
    equations
    {
        U               0.7;
        "(k|omega)"     0.6;
    }
}
""")
with open(os.path.join(case, "case_info.txt"), "w", newline="\n") as fh:
    fh.write(f"assembly {os.path.basename(os.path.dirname(os.path.abspath(case)))}\nQ_lpm {q_lpm:g}\nlevel {level}\n"
             f"background_cell_m {h}\nU_inlet_m_s {U_in:.6g}\nk_inlet {k_in:.6g}\nomega_inlet {om_in:.6g}\n"
             f"wall_treatment wall functions (nutkWallFunction, omegaWallFunction)\nhalf_domain yes (symmetry z=0)\n")
print(f"case {case}: background {n} cells of {h*1e3:.2f} mm, U_in {U_in:.4f} m/s")
