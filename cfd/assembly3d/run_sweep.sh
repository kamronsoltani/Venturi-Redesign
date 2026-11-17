#!/bin/bash
# run_sweep.sh <assembly dir under cfd/> <flows> [level]   -- run inside WSL.
#
# Runs one assembly at several flow rates, one case at a time, unattended:
#   bash run_sweep.sh "<repo>/cfd/asm34_current" "20 10 8 6 40 50 60"
#   bash run_sweep.sh "<repo>/cfd/asm34_current" "30" coarse      (mesh study)
#
# For each flow Q (in the order given) it writes the case <asm>/Q<Q>[_<level>],
# meshes it (or reuses the mesh of an earlier case of the same assembly and
# level), starts from the final fields of the nearest finished flow, runs
# checkMesh + the solver (run_asm.sh), samples the result (post_asm.sh),
# computes the reported quantities (post_asm.py -> results.json), builds the
# 3-D viewer file and refreshes cfd/results/assembly_3d_stage5b.csv and
# docs/3d/index.json (collect_asm.py).
#
# A case that already has results.json is skipped, so re-running the same
# command resumes an interrupted sweep. The water volume must already exist in
# <asm>/fluid (fluid.stl + assembly_info.json, from build_fluid.py).
#
# Environment: ITERS (default 2000), NPROC (default 8),
#   LABEL  (viewer label, e.g. "3/4in E-Sig + 3/4in elbow + current Venturi"),
#   PUSH_STATUS=1  push a progress page to the GitHub branch "sweep-status"
#                  every 10 min and after each case (needs git.exe on Windows
#                  with push rights; nothing is pushed to main).
source /opt/openfoam13/etc/bashrc >/dev/null 2>&1
set -u
ASM="$(cd "$1" && pwd)" || { echo "no assembly dir $1"; exit 2; }
FLOWS="$2"
LEVEL="${3:-medium}"
export ASM_ITERS="${ITERS:-2000}"
NPROC="${NPROC:-8}"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
NAME="$(basename "$ASM")"
WROOT="$HOME/turitap_cfd/${ASM#*/cfd/}"
LABEL="${LABEL:-$NAME}"
SUFFIX=""; [ "$LEVEL" != "medium" ] && SUFFIX="_$LEVEL"
LOG="$ASM/sweep.log"
mkdir -p "$WROOT"
echo "$ASM|$FLOWS|$LEVEL|$ASM_ITERS" > "$HOME/turitap_cfd/SWEEP"
note() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a "$LOG"; }

status_push() {
    python3 "$HERE/sweep_status.py" "$ASM" "$FLOWS" "$LEVEL" > "$ASM/SWEEP_STATUS.md" 2>/dev/null || return
    [ "${PUSH_STATUS:-0}" = "1" ] || return
    local RW blob tree commit
    RW="$(wslpath -w "$REPO")"
    blob=$(git.exe -C "$RW" hash-object -w --stdin < "$ASM/SWEEP_STATUS.md" | tr -d '\r') || return
    tree=$(printf '100644 blob %s\tSTATUS.md\n' "$blob" | git.exe -C "$RW" mktree | tr -d '\r') || return
    commit=$(echo "sweep status $(date '+%Y-%m-%d %H:%M')" | git.exe -C "$RW" commit-tree "$tree" | tr -d '\r') || return
    git.exe -C "$RW" push -q -f origin "$commit:refs/heads/sweep-status" >/dev/null 2>&1 || note "status push failed"
}
heartbeat() { while true; do sleep 600; status_push; done; }
heartbeat & HB=$!
trap 'kill $HB 2>/dev/null' EXIT

# a finished work dir (fields at a time > 0) of this assembly and level -> its flow
finished_flows() {
    for d in "$WROOT"/Q*; do
        [ -f "$d/STATUS" ] && grep -q "^done" "$d/STATUS" || continue
        grep -q "^level $LEVEL$" "$d/case_info.txt" 2>/dev/null || continue
        q=$(awk '/^Q_lpm/{print $2}' "$d/case_info.txt"); echo "$q $d"
    done
}
mesh_source() {
    for d in "$WROOT"/Q*; do
        [ -d "$d/constant/polyMesh" ] && [ -f "$d/mesh_summary.txt" ] && grep -q "Mesh OK" "$d/mesh_summary.txt" || continue
        grep -q "^level $LEVEL$" "$d/case_info.txt" 2>/dev/null && { echo "$d"; return; }
    done
}

note "sweep start: $NAME flows [$FLOWS] level $LEVEL, $ASM_ITERS iterations, $NPROC cores"
status_push
for Q in $FLOWS; do
    CASE="$ASM/Q${Q}${SUFFIX}"
    if [ -f "$CASE/results.json" ]; then note "Q$Q$SUFFIX: results exist, skipped"; continue; fi
    rm -rf "$CASE"
    python3 "$HERE/setup_asm_case.py" "$CASE" "$ASM/fluid" "$Q" "$LEVEL" > /dev/null || { note "Q$Q: setup failed"; continue; }
    MESH_FROM="$(mesh_source)"
    INIT_FROM="$(finished_flows | awk -v q="$Q" '{d=$1-q; if(d<0)d=-d; print d, $2}' | sort -n | head -1 | cut -d' ' -f2-)"
    note "Q$Q$SUFFIX: start (mesh from: ${MESH_FROM:-new}; init from: ${INIT_FROM:-uniform})"
    echo "$WROOT/Q${Q}${SUFFIX}" > "$HOME/turitap_cfd/CURRENT_CASE"
    MESH_FROM="$MESH_FROM" INIT_FROM="$INIT_FROM" bash "$HERE/run_asm.sh" "$CASE" "$NPROC"
    RC=$?
    if [ $RC -ne 0 ] && ! grep -q "^done" "$CASE/STATUS" 2>/dev/null; then
        note "Q$Q$SUFFIX: FAILED ($(cat "$CASE/STATUS" 2>/dev/null)); see $CASE/log.*"; status_push; continue
    fi
    bash "$HERE/post_asm.sh" "$CASE" > /dev/null 2>&1
    (cd "$HERE" && python3 post_asm.py "$CASE" > "$CASE/post_summary.txt" 2>&1) || note "Q$Q$SUFFIX: post_asm.py failed"
    (cd "$HERE" && python3 build_viz3d.py "$CASE" "$REPO/docs/3d/cases/${NAME}_Q${Q}${SUFFIX}.json" \
        "$LABEL, $Q L/min$([ -n "$SUFFIX" ] && echo " ($LEVEL mesh)")" > /dev/null 2>&1) || note "Q$Q$SUFFIX: viewer build failed"
    (cd "$HERE" && python3 collect_asm.py "$REPO" > /dev/null 2>&1) || note "collect_asm.py failed"
    note "Q$Q$SUFFIX: done. $(grep -E '^level_inlet_cm|^static_inlet_cm' "$CASE/post_summary.txt" | tr -s ' ' | tr '\n' ' ')"
    status_push
done
note "sweep finished"
status_push
