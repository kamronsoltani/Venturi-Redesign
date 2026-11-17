#!/bin/bash
# run_case.sh "<case dir under the project's cfd/>"   -- run inside WSL.
#
# OpenFOAM 13 refuses paths containing spaces ("Venturi Scaling"), so the
# case is staged to ~/turitap_cfd/<same relative path>, run there, and the
# results (logs, postProcessing, final time) are copied back.
#
# 1. checkMesh -> log.checkMesh; stop if it fails or max non-orthogonality > 65 deg
# 2. foamRun   -> log.foamRun   (never printed in full)
source /opt/openfoam13/etc/bashrc >/dev/null 2>&1   # before set -u: bashrc uses unset vars
set -u
SRC="$(cd "$1" && pwd)" || { echo "no case dir $1"; exit 2; }
REL="${SRC#*/cfd/}"
WORK="$HOME/turitap_cfd/$REL"
rm -rf "$WORK" && mkdir -p "$WORK"
cp -r "$SRC/0" "$SRC/constant" "$SRC/system" "$WORK/"
cd "$WORK"

copy_back() {
    rm -rf "$SRC/postProcessing"
    cp log.* "$SRC/" 2>/dev/null
    [ -d postProcessing ] && cp -r postProcessing "$SRC/"
    LATEST=$(foamListTimes -latestTime -noZero 2>/dev/null | tail -1)
    [ -n "$LATEST" ] && rm -rf "$SRC/$LATEST" && cp -r "$LATEST" "$SRC/"
}

checkMesh > log.checkMesh 2>&1
NONORTH=$(grep -oP 'Mesh non-orthogonality Max: \K[0-9.eE+-]+' log.checkMesh | head -1)
FAILED=$(grep -oP 'Failed \K[0-9]+(?= mesh checks)' log.checkMesh | head -1)
echo "checkMesh: max non-orthogonality = ${NONORTH:-?} deg; $(grep -E 'Mesh OK|Failed [0-9]+ mesh checks' log.checkMesh | head -1 | sed 's/^ *//')"
if ! grep -q "Mesh OK" log.checkMesh || [ -n "$FAILED" ]; then
    copy_back; echo "STOP: checkMesh reported failure (see log.checkMesh)"; exit 3
fi
if awk -v n="${NONORTH:-999}" 'BEGIN{exit !(n > 65)}'; then
    copy_back; echo "STOP: non-orthogonality above 65 deg"; exit 4
fi

NPROC="${NPROC:-1}"          # NPROC=n bash run_case.sh <case>  -> run on n cores
if [ "$NPROC" -gt 1 ]; then
    cat > system/decomposeParDict <<EOF
FoamFile { format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains $NPROC;
decomposer      scotch;
method          scotch;
EOF
    decomposePar -force > log.decomposePar 2>&1 || { copy_back; echo "STOP: decomposePar failed"; exit 5; }
    mpirun -np "$NPROC" foamRun -parallel > log.foamRun 2>&1
    STATUS=$?
    reconstructPar -latestTime > log.reconstructPar 2>&1
    rm -rf processor*
else
    foamRun > log.foamRun 2>&1
    STATUS=$?
fi
copy_back
echo "foamRun exit status: $STATUS (work dir $WORK)"
grep -E "solution converged|^End" log.foamRun | tail -2
exit $STATUS
