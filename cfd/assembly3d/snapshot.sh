#!/bin/bash
# snapshot.sh <case rel under cfd/> <time> <snapshot dir (Linux path, under a */cfd/ folder)>
# Copies one written time out of a RUNNING parallel case, reconstructs it in a
# separate folder and runs post_asm.sh on it, so the live solve is not touched.
source /opt/openfoam13/etc/bashrc >/dev/null 2>&1
set -eu
REL=$1; T=$2; SNAP=$3
RUN="$HOME/turitap_cfd/$REL"
SREL="${SNAP#*/cfd/}"; W="$HOME/turitap_cfd/$SREL"
rm -rf "$W"; mkdir -p "$W" "$SNAP"
cp -r "$RUN/constant" "$RUN/system" "$RUN/assembly_info.json" "$W/"
for p in "$RUN"/processor*; do
  b=$(basename "$p"); mkdir -p "$W/$b"; cp -r "$p/constant" "$p/$T" "$W/$b/"
done
cd "$W"; reconstructPar -time "$T" > log.reconstruct 2>&1; rm -rf processor*
cp assembly_info.json "$SNAP/"
bash "$(dirname "$0")/post_asm.sh" "$SNAP"
