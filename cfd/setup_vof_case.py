#!/usr/bin/env python3
"""
setup_vof_case.py -- create one Stage 4 two-phase (VoF) case, OpenFOAM 13
foamRun -solver incompressibleVoF, on the wedge mesh of make_vof_mesh.py
(bore + open air box below the outlet).

Physics
  * transient, incompressible water/air VoF (MULES, interface compression)
  * water 25 C: rho 997, nu 8.9e-7; air 25 C: rho 1.18, nu 1.56e-5;
    surface tension 0.072 N/m; contact angle 90 deg (zeroGradient alpha)
  * gravity +9.81 along +x (flow direction, downward)
  * k-omega SST, wall-resolved, same wall treatment as Stage 1
  * inlet: fixed flow rate of water (alpha = 1) -- steady piped install
  * atmosphere (box side, bottom, annulus level with the outlet face):
    p_rgh prghTotalPressure p0 = 0, U pressureInletOutletVelocity,
    alpha inletOutlet (air flows in)
  * initial state: bore full of water with plug flow Q/A(z), air box empty
    (i.e. the flow is reduced from a running, full device; whether air then
    invades the outlet and reaches the throat is the question)

Monitors (postProcessing): wall probes of p and alpha.water at the port and
along the diffuser/tail, bore water fraction, inlet p, flow rates.

Usage
  python3 setup_vof_case.py --geom current --q 6 [--level medium] [--end 3]
Default case dir: cfd/<geom>_vof/Q<q>
"""
import argparse
import math
import os
import shutil

import make_vof_mesh as mvm
from setup_case import hdr, write, field, NU, RHO, G, BETA1, TURB_INTENSITY, MIXING_LENGTH_FRAC

HERE = os.path.dirname(os.path.abspath(__file__))
RHO_AIR = 1.18
NU_AIR = 1.56e-5
SIGMA = 0.072
END_TIME = 2.0
WRITE_INTERVAL = 0.5


def wall_r(profile, z):
    for (za, ra), (zb, rb) in zip(profile, profile[1:]):
        if za - 1e-9 <= z <= zb + 1e-9:
            return ra + (rb - ra) * (z - za) / (zb - za)
    raise ValueError(z)


def setup(geom, q_lpm, level="vof", case=None, end=END_TIME):
    if case is None:
        case = os.path.join(HERE, f"{geom}_vof", f"Q{q_lpm:g}" + ("" if level == "vof" else f"_{level}"))
    if os.path.isdir(case):
        shutil.rmtree(case)
    os.makedirs(case)
    info = mvm.build_mesh(geom, case, level)
    NR, NXd, n = info["NR"], info["NXd"], info["n_cells"]
    zn, Rn = info["z_nodes_mm"], info["R_nodes_mm"]

    Rin = info["inlet_radius_mm"] * 1e-3
    Q = q_lpm / 60000.0
    Uin = Q / (math.pi * Rin ** 2)
    k_in = 1.5 * (TURB_INTENSITY * Uin) ** 2
    L_mix = MIXING_LENGTH_FRAC * 2 * Rin
    om_in = math.sqrt(k_in) / (0.09 ** 0.25 * L_mix)
    k_air, om_air = 1e-8, 1.0
    om_wall = [10 * 6 * NU / (BETA1 * y ** 2) for y in info["wall_y1_m"]]
    om_list = "nonuniform List<scalar> %d\n(\n%s\n)" % (len(om_wall), "\n".join(f"{w:.6g}" for w in om_wall))
    wedge = ["type            wedge;"]

    # initial fields: bore full with plug flow, box air at rest
    alpha0, U0, k0, om0 = [], [], [], []
    for i in range(NXd):
        Rm = 0.5 * (Rn[i] + Rn[i + 1]) * 1e-3
        u = Q / (math.pi * Rm ** 2)
        for _ in range(NR):
            alpha0.append("1"); U0.append(f"({u:.6g} 0 0)"); k0.append(f"{k_in:.6g}"); om0.append(f"{om_in:.6g}")
    for _ in range(n - NXd * NR):
        alpha0.append("0"); U0.append("(0 0 0)"); k0.append(f"{k_air:.6g}"); om0.append(f"{om_air:.6g}")

    def nl(kind, vals):
        return f"nonuniform List<{kind}> {len(vals)}\n(\n" + "\n".join(vals) + "\n)"

    field(case, "alpha.water", "volScalarField", "[0 0 0 0 0 0 0]", nl("scalar", alpha0), {
        "inlet": ["type            fixedValue;", "value           uniform 1;"],
        "wall": ["type            zeroGradient;"],
        "atmosphere": ["type            inletOutlet;", "inletValue      uniform 0;", "value           uniform 0;"],
        "front": wedge, "back": wedge})
    field(case, "U", "volVectorField", "[0 1 -1 0 0 0 0]", nl("vector", U0), {
        "inlet": ["type            flowRateInletVelocity;",
                  f"meanVelocity    {Uin:.8g};   // Q = {q_lpm:g} L/min over the full inlet area",
                  "value           uniform (0 0 0);"],
        "wall": ["type            noSlip;"],
        "atmosphere": ["type            pressureInletOutletVelocity;", "value           uniform (0 0 0);"],
        "front": wedge, "back": wedge})
    field(case, "p_rgh", "volScalarField", "[1 -1 -2 0 0 0 0]", "uniform 0", {
        "inlet": ["type            fixedFluxPressure;", "value           uniform 0;"],
        "wall": ["type            fixedFluxPressure;", "value           uniform 0;"],
        "atmosphere": ["type            prghTotalPressure;", "p0              uniform 0;", "value           uniform 0;"],
        "front": wedge, "back": wedge})
    field(case, "p", "volScalarField", "[1 -1 -2 0 0 0 0]", "uniform 0", {
        "inlet": ["type            calculated;", "value           uniform 0;"],
        "wall": ["type            calculated;", "value           uniform 0;"],
        "atmosphere": ["type            calculated;", "value           uniform 0;"],
        "front": wedge, "back": wedge})
    field(case, "k", "volScalarField", "[0 2 -2 0 0 0 0]", nl("scalar", k0), {
        "inlet": ["type            turbulentIntensityKineticEnergyInlet;",
                  f"intensity       {TURB_INTENSITY};", f"value           uniform {k_in:.6g};"],
        "wall": ["type            fixedValue;", "value           uniform 0;"],
        "atmosphere": ["type            inletOutlet;", f"inletValue      uniform {k_air:g};", f"value           uniform {k_air:g};"],
        "front": wedge, "back": wedge})
    field(case, "omega", "volScalarField", "[0 0 -1 0 0 0 0]", nl("scalar", om0), {
        "inlet": ["type            turbulentMixingLengthFrequencyInlet;",
                  f"mixingLength    {L_mix:.6g};", f"value           uniform {om_in:.6g};"],
        "wall": ["type            fixedValue;   // Menter: 10*6*nu/(beta1*y1^2)", f"value           {om_list};"],
        "atmosphere": ["type            inletOutlet;", f"inletValue      uniform {om_air:g};", f"value           uniform {om_air:g};"],
        "front": wedge, "back": wedge})
    field(case, "nut", "volScalarField", "[0 2 -1 0 0 0 0]", "uniform 0", {
        "inlet": ["type            calculated;", "value           uniform 0;"],
        "wall": ["type            fixedValue;", "value           uniform 0;"],
        "atmosphere": ["type            calculated;", "value           uniform 0;"],
        "front": wedge, "back": wedge})

    write(case, "constant/phaseProperties", hdr("dictionary", "constant", "phaseProperties") +
          f"phases          (water air);\n\nsigma           {SIGMA};\n")
    write(case, "constant/physicalProperties.water", hdr("dictionary", "constant", "physicalProperties.water") +
          f"viscosityModel  constant;\n\nnu              {NU};\n\nrho             {RHO:g};\n")
    write(case, "constant/physicalProperties.air", hdr("dictionary", "constant", "physicalProperties.air") +
          f"viscosityModel  constant;\n\nnu              {NU_AIR};\n\nrho             {RHO_AIR};\n")
    write(case, "constant/momentumTransport", hdr("dictionary", "constant", "momentumTransport") +
          "simulationType  RAS;\n\nRAS\n{\n    model           kOmegaSST;\n    turbulence      on;\n}\n")
    write(case, "constant/g", hdr("uniformDimensionedVectorField", "constant", "g") +
          f"dimensions      [0 1 -2 0 0 0 0];\nvalue           ({G} 0 0);   // along the flow (downward)\n")

    # wall probes (face centres on the wedge mid-plane)
    prof = info["wall_profile_mm"]
    stn = info["stations"]
    zs = [stn["port"], stn["diffuser_start"]]
    zs += [z for z in (65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 114.0)
           if stn["diffuser_start"] < z < stn["outlet"] - 0.5]
    c = math.cos(math.radians(info["wedge_half_deg"]))
    locs = " ".join(f"({z*1e-3:.6g} {wall_r(prof, z)*c*1e-3:.6g} 0)" for z in zs)
    # axis probes: port level, outlet face, and in the box
    axl = " ".join(f"({z*1e-3:.6g} 1e-5 0)" for z in (stn["port"], stn["outlet"], stn["outlet"] + 10, stn["outlet"] + 30))

    write(case, "system/controlDict", hdr("dictionary", "system", "controlDict") + f"""\
solver          incompressibleVoF;

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {end:g};
deltaT          1e-6;

writeControl    adjustableRunTime;
writeInterval   {WRITE_INTERVAL:g};
purgeWrite      0;
writeFormat     binary;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

adjustTimeStep  yes;
maxCo           1;
maxAlphaCo      1;
maxDeltaT       1e-3;

functions
{{
    wallProbes   // z (mm): {', '.join(f'{z:g}' for z in zs)}
    {{
        type            patchProbes;
        libs            ("libsampling.so");
        patchName       wall;
        fields          (p alpha.water);
        probeLocations  ({locs});
        writeControl    timeStep;
        writeInterval   20;
    }}
    axisProbes   // port, outlet face, outlet+10, outlet+30 mm
    {{
        type            probes;
        libs            ("libsampling.so");
        fields          (p alpha.water U);
        probeLocations  ({axl});
        writeControl    timeStep;
        writeInterval   20;
    }}
    boreWater
    {{
        type            volFieldValue;
        libs            ("libfieldFunctionObjects.so");
        select          cellZone;
        cellZone        bore;
        operation       volAverage;
        fields          (alpha.water);
        writeFields     false;
        writeControl    timeStep;
        writeInterval   20;
    }}
    inletAvg
    {{
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        select          patch;
        patch           inlet;
        operation       areaAverage;
        fields          (p);
        writeFields     false;
        writeControl    timeStep;
        writeInterval   20;
    }}
    #includeFunc patchFlowRate(name=inletFlow, patch=inlet)
    #includeFunc patchFlowRate(name=atmFlow, patch=atmosphere)
    #includeFunc yPlus
}}
""")
    write(case, "system/fvSchemes", hdr("dictionary", "system", "fvSchemes") + """\
ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         Gauss linear;
    grad(k)         cellLimited Gauss linear 1;
    grad(omega)     cellLimited Gauss linear 1;
}

divSchemes
{
    default         none;
    div(rhoPhi,U)   Gauss linearUpwind grad(U);
    div(phi,alpha)  Gauss interfaceCompression vanLeer 1;
    div(phi,k)      Gauss limitedLinear 1;
    div(phi,omega)  Gauss limitedLinear 1;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

wallDist
{
    method          meshWave;
}
""")
    write(case, "system/fvSolution", hdr("dictionary", "system", "fvSolution") + """\
solvers
{
    "alpha.water.*"
    {
        nCorrectors     2;
        nSubCycles      1;
        MULESCorr       yes;
        MULES
        {
            nIter           3;
        }
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0;
    }

    "pcorr.*"
    {
        solver          PCG;
        preconditioner
        {
            preconditioner  GAMG;
            tolerance       1e-5;
            relTol          0;
            smoother        GaussSeidel;
        }
        tolerance       1e-5;
        relTol          0;
        maxIter         50;
    }

    p_rgh
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.05;
        maxIter         50;
    }

    p_rghFinal
    {
        $p_rgh;
        relTol          0;
    }

    "(U|k|omega).*"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0.1;
    }
}

PIMPLE
{
    momentumPredictor no;
    nOuterCorrectors 1;
    nCorrectors     2;
    nNonOrthogonalCorrectors 0;
}

relaxationFactors
{
    equations
    {
        ".*"            1;
    }
}
""")
    with open(os.path.join(case, "case_info.txt"), "w", newline="\n") as fh:
        fh.write(f"stage 4 VoF\ngeometry {geom}\nQ_lpm {q_lpm:g}\nlevel {level}\nU_inlet_m_s {Uin:.6g}\n"
                 f"end_time_s {end:g}\ninitial bore full of water, air box empty\n"
                 f"wall probe z_mm {' '.join(f'{z:g}' for z in zs)}\n")
    print(f"case {case}: {n} cells ({NXd*NR} bore), U_in = {Uin:.4f} m/s")
    return case


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--geom", required=True)
    ap.add_argument("--q", type=float, required=True)
    ap.add_argument("--level", default="vof")
    ap.add_argument("--case", default=None)
    ap.add_argument("--end", type=float, default=END_TIME)
    a = ap.parse_args()
    setup(a.geom, a.q, a.level, a.case, a.end)


if __name__ == "__main__":
    main()
