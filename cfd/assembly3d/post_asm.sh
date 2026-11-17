#!/bin/bash
# post_asm.sh "<case dir under cfd/>"  -- run inside WSL after run_asm.sh.
# Samples the final time in the WSL work dir: cut planes at the part joints, the
# inlet/outlet/wall/symmetry patches (VTK, ASCII) and inlet-seeded streamlines,
# then copies postProcessing back to the case folder.
source /opt/openfoam13/etc/bashrc >/dev/null 2>&1
set -u
SRC="$(cd "$1" && pwd)" || exit 2
REL="${SRC#*/cfd/}"
WORK="$HOME/turitap_cfd/$REL"
cd "$WORK" || exit 2
read -r SPIG ESEND VIN PORT DIFE OUTX < <(python3 - <<'PY'
import json; i=json.load(open("assembly_info.json")); m=1e-3
print(i["spigot_entry_y"]*m, (i["esig_end_face_y"]+0.5)*m, (i["venturi_inlet_x"]-0.3)*m,
      i["port_x"]*m, (i["venturi_inlet_x"]-i.get("diffuser_end_mm",79.4))*m, i["outlet_x"]*m)
PY
)
plane() { # name px py pz nx ny nz
cat <<EOD
        $1 { type cutPlane; planeType pointAndNormal; point ($2 $3 $4); normal ($5 $6 $7); interpolate false; }
EOD
}
cat > system/postFunctions <<EOD
FoamFile { format ascii; class dictionary; object postFunctions; }
    yPlus { type yPlus; libs ("libfieldFunctionObjects.so"); writeControl writeTime; }
    wallShearStress { type wallShearStress; libs ("libfieldFunctionObjects.so"); writeControl writeTime; }
    sections
    {
        type surfaces; libs ("libsampling.so"); writeControl writeTime;
        surfaceFormat vtk; formatOptions { vtk { format ascii; } }; interpolationScheme cellPoint;
        fields (p U); interpolate false;
        surfaces
        {
            inlet  { type patch; patches (inlet); interpolate false; }
            outlet { type patch; patches (outlet); interpolate false; }
$(plane spigot 0 $SPIG 0 0 1 0)
$(plane esigEnd 0 $ESEND 0 0 1 0)
$(plane ventIn $VIN 0 0 1 0 0)
$(plane port $PORT 0 0 1 0 0)
$(plane diffEnd $DIFE 0 0 1 0 0)
        }
    }
    walls
    {
        type surfaces; libs ("libsampling.so"); writeControl writeTime;
        surfaceFormat vtk; formatOptions { vtk { format ascii; } }; interpolationScheme cellPoint;
        fields (p wallShearStress yPlus); interpolate false;
        surfaces
        {
            wall_esig    { type patch; patches (wall_esig); interpolate false; }
            wall_elbow   { type patch; patches (wall_elbow); interpolate false; }
            wall_venturi { type patch; patches (wall_venturi); interpolate false; }
        }
    }
    viz
    {
        type surfaces; libs ("libsampling.so"); writeControl writeTime;
        surfaceFormat vtk; formatOptions { vtk { format ascii; } }; interpolationScheme cellPoint;
        fields (p U wallShearStress); interpolate true;
        surfaces
        {
            centre { type patch; patches (symmetry); interpolate true; }
            wall   { type patch; patches ("wall_.*"); interpolate true; }
        }
    }
    streamlines
    {
        type streamlines; libs ("libfieldFunctionObjects.so");
        executeControl writeTime; writeControl writeTime;
        setFormat vtk; formatOptions { vtk { format ascii; } }; interpolationScheme cellPoint;
        lifeTime 20000; nSubCycle 4; cloudName particleTracks;
        direction forward;
        seedSampleSet { type boundaryRandom; patches (inlet); nPoints 400; }
        fields (p U);
    }
EOD
foamPostProcess -solver incompressibleFluid -latestTime -dict system/postFunctions > log.post 2>&1
echo "post exit $? ; $(grep -c . log.post) log lines; last: $(tail -1 log.post)"
rm -rf "$SRC/postProcessing" && cp -r postProcessing "$SRC/"
cp log.post "$SRC/"
