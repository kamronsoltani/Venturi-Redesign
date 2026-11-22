#!/bin/bash
# run_campaign.sh  -- run inside WSL. Stage 5b ranking campaign, unattended.
#
# Every assembly at 30 L/min first (the ranking run), then every assembly at
# 6 L/min, then every assembly at 8 L/min (dosing cut-off and port suction).
# Each (assembly, flow) case is one run_sweep.sh call, so all of its behaviour
# applies: checkMesh before every solve (stop on failure or non-orthogonality
# > 65 deg), the mesh of an assembly is made once and reused, each flow starts
# from that assembly's nearest finished flow, cases with results.json are
# skipped (re-running the same command resumes after a crash or reboot), and a
# failed case is logged and the campaign moves on to the next one.
#
#   bash run_campaign.sh                       # all 11 assemblies, 30 then 6 then 8
#   FLOWS_ORDER="30" bash run_campaign.sh      # only the ranking run
#   ASSEMBLIES="asm1lr_V1 asm1lr_V4" bash run_campaign.sh
#
# Environment: FLOWS_ORDER (default "30 6 8"), ASSEMBLIES (default all 11),
# plus run_sweep.sh's ITERS, NPROC, PUSH_STATUS.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
CFD="$(cd "$HERE/.." && pwd)"
FLOWS_ORDER="${FLOWS_ORDER:-30 6 8}"
ASSEMBLIES="${ASSEMBLIES:-asm34_current asm34_A asm34_Bvane asm1e_V1 asm1e_V2 asm1e_V3 asm1e_V4 asm1lr_V1 asm1lr_V2 asm1lr_V3 asm1lr_V4}"
LOG="$CFD/results/campaign_stage5b.log"
mkdir -p "$HOME/turitap_cfd" "$CFD/results"
printf 'flows %s\nassemblies %s\n' "$FLOWS_ORDER" "$ASSEMBLIES" > "$HOME/turitap_cfd/CAMPAIGN"

label() {
    local up v
    case "$1" in
        asm34_*)  up="3/4in E-Sig + 3/4in threaded elbow" ;;
        asm1e_*)  up="1in E-Sig + McMaster 4880K464 1in elbow" ;;
        asm1lr_*) up="1in long-radius reducing elbow" ;;
        asm34lr_*) up="3/4in E-Sig + 3/4in long-radius elbow (R 42 mm, no reduction)" ;;
    esac
    case "${1#*_}" in
        current) v="current Venturi" ;;
        A)       v="Venturi A (14 mm throat, 19 mm exit)" ;;
        Bvane)   v="Venturi B (14 mm throat, 26 mm exit) + cross-vane" ;;
        V1)      v="V1 (14 mm throat, 26 mm exit) + vane" ;;
        V2)      v="V2 (16 mm throat, 26 mm exit) + vane" ;;
        V3)      v="V3 (18 mm throat, 26 mm exit) + vane" ;;
        V4)      v="V4 (14 mm throat, 30 mm exit) + 6-blade vane" ;;
        *)       v="${1#*_}" ;;
    esac
    echo "$up + $v"
}

for a in $ASSEMBLIES; do
    [ -f "$CFD/$a/fluid/fluid.stl" ] || { echo "missing water volume $CFD/$a/fluid/fluid.stl"; exit 2; }
done
echo "$(date '+%Y-%m-%d %H:%M:%S') campaign start: flows [$FLOWS_ORDER], assemblies [$ASSEMBLIES]" | tee -a "$LOG"
for Q in $FLOWS_ORDER; do
    for a in $ASSEMBLIES; do
        LABEL="$(label "$a")" bash "$HERE/run_sweep.sh" "$CFD/$a" "$Q"
        echo "$(date '+%Y-%m-%d %H:%M:%S') $a Q$Q: $(cat "$CFD/$a/Q$Q/STATUS" 2>/dev/null || echo 'no status')" | tee -a "$LOG"
    done
done
echo "$(date '+%Y-%m-%d %H:%M:%S') campaign finished" | tee -a "$LOG"
rm -f "$HOME/turitap_cfd/CAMPAIGN"
