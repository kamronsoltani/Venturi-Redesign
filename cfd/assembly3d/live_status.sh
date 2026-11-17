#!/bin/bash
# live_status.sh <work dir in WSL, e.g. ~/turitap_cfd/asm34_current/Q30>  -- refreshing status view
W="${1:-$HOME/turitap_cfd/asm34_current/Q30}"
while true; do
  clear
  echo "3-D CFD  |  $(basename "$(dirname "$W")")/$(basename "$W")  |  $(date +%H:%M:%S)"
  echo "-----------------------------------------------------------"
  [ -f "$W/STATUS.log" ] && tail -n 8 "$W/STATUS.log" || echo "waiting to start..."
  [ -f "$W/mesh_summary.txt" ] && { echo; echo "mesh: $(cat "$W/mesh_summary.txt")"; }
  if [ -f "$W/log.snappyHexMesh" ] && grep -q "meshing: snappy" "$W/STATUS" 2>/dev/null; then
    echo; echo "snappyHexMesh: $(grep -E '^(Surface refinement|Shell refinement|Feature refinement|Snapping|Adding layers|Layer addition|Finished|Morphing|Smoothing)' "$W/log.snappyHexMesh" | tail -1)"
  fi
  R=$(ls -d "$W"/postProcessing/residuals*/*/residuals.dat 2>/dev/null | tail -1)
  if [ -n "$R" ]; then
    echo; tail -1 "$R" | awk '{printf "iteration %-6s of 3000   residuals  p %.1e  Ux %.1e  Uy %.1e  k %.1e  omega %.1e\n",$1,$2,$3,$4,$6,$7}'
    A=$(ls -d "$W"/postProcessing/inletAvg/*/*.dat 2>/dev/null | tail -1)
    [ -n "$A" ] && tail -1 "$A" | awk '{printf "inlet static pressure %.3f kPa  (%.1f cm of water above the spigot entry, static)\n",$2*0.997,$2/9.81*100}'
    FI=$(ls -d "$W"/postProcessing/inletFlow/*/*.dat 2>/dev/null | tail -1); FO=$(ls -d "$W"/postProcessing/outletFlow/*/*.dat 2>/dev/null | tail -1)
    [ -n "$FI" ] && [ -n "$FO" ] && paste <(tail -1 "$FI") <(tail -1 "$FO") | awk '{printf "flow in %.2f L/min, out %.2f L/min (full domain)\n",-$2*120000,$4*120000}'
    echo; echo "residual p, last 10 samples:"; tail -n 1000 "$R" | awk 'NR%100==0{printf "  it %-5s p %.1e\n",$1,$2}' | tail -10
  fi
  sleep 15
done
