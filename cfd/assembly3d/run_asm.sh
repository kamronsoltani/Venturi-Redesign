#!/bin/bash
# run_asm.sh "<case dir under cfd/>" [nproc]  -- run inside WSL.
# Stages the case to ~/turitap_cfd/<rel> (OpenFOAM 13 rejects paths with spaces),
# meshes it (blockMesh, surfaceFeatures, snappyHexMesh), runs checkMesh (stop on
# failure or non-orthogonality > 65 deg), solves in parallel, and copies logs,
# postProcessing, mesh_info and the final time back. The work dir keeps the mesh.
source /opt/openfoam13/etc/bashrc >/dev/null 2>&1
set -u
SRC="$(cd "$1" && pwd)" || { echo "no case dir $1"; exit 2; }
NPROC="${2:-8}"
REL="${SRC#*/cfd/}"
WORK="$HOME/turitap_cfd/$REL"
rm -rf "$WORK" && mkdir -p "$WORK"
cp -r "$SRC/0" "$SRC/constant" "$SRC/system" "$SRC/assembly_info.json" "$WORK/"
cd "$WORK"
stage() { echo "$1" > STATUS; date +"%H:%M:%S $1" >> STATUS.log; }
copy_back() {
    cp log.* STATUS.log "$SRC/" 2>/dev/null
    rm -rf "$SRC/postProcessing"; [ -d postProcessing ] && cp -r postProcessing "$SRC/"
}
stage "meshing: blockMesh"
blockMesh > log.blockMesh 2>&1 || { stage "FAILED blockMesh"; copy_back; exit 3; }
stage "meshing: surfaceFeatures"
surfaceFeatures > log.surfaceFeatures 2>&1
cat > system/decomposeParDict <<EOD
FoamFile { format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains $NPROC;
decomposer      scotch;
distributor     scotch;
EOD
stage "meshing: snappyHexMesh ($NPROC cores)"
decomposePar -force > log.decomposePar.mesh 2>&1
mpirun -np "$NPROC" snappyHexMesh -overwrite -parallel > log.snappyHexMesh 2>&1 || { stage "FAILED snappyHexMesh"; copy_back; exit 3; }
reconstructPar -constant > log.reconstructPar.mesh 2>&1
rm -rf processor*
stage "checkMesh"
checkMesh > log.checkMesh 2>&1
NONORTH=$(grep -oP 'Mesh non-orthogonality Max: \K[0-9.eE+-]+' log.checkMesh | head -1)
CELLS=$(grep -oP '^\s*cells:\s*\K[0-9]+' log.checkMesh | head -1)
echo "cells ${CELLS:-?}  max non-orthogonality ${NONORTH:-?} deg  $(grep -E 'Mesh OK|Failed [0-9]+ mesh checks' log.checkMesh | head -1 | sed 's/^ *//')" | tee mesh_summary.txt
cp mesh_summary.txt "$SRC/"
if ! grep -q "Mesh OK" log.checkMesh; then stage "STOP: checkMesh failed"; copy_back; exit 4; fi
if awk -v n="${NONORTH:-999}" 'BEGIN{exit !(n > 65)}'; then stage "STOP: non-orthogonality > 65"; copy_back; exit 4; fi
stage "solving ($NPROC cores)"
decomposePar -force > log.decomposePar 2>&1
mpirun -np "$NPROC" foamRun -parallel > log.foamRun 2>&1
STATUS=$?
stage "reconstructing"
reconstructPar -latestTime > log.reconstructPar 2>&1
rm -rf processor*
copy_back
LATEST=$(foamListTimes -latestTime -noZero 2>/dev/null | tail -1)
if grep -q "^End" log.foamRun; then stage "done (iteration $LATEST)"; else stage "FAILED solver (status $STATUS)"; fi
cp STATUS STATUS.log "$SRC/"
exit $STATUS
