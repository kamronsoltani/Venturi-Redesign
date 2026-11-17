# Next steps for the lab computer (updated 2026-11-17)

Read this first, then `LAB_SETUP.md` and `PROJECT_NOTES.md`. Follow the rules at
the end of `LAB_SETUP.md`.

## Where things stand

- Water volumes for all 11 assemblies are in `cfd/<assembly>/fluid`; pictures
  in `cfd/assembly_previews/`. The geometry is final and approved (including
  the off-centre McMaster elbow joint, which matches the real fitting). Do not
  rebuild or change it.
- `asm34_current` (3/4in E-Sig + 3/4in elbow + current Venturi) at 30 L/min is
  done. This machine may already be running it at 8 and 6 L/min (started from
  `run_sweep.sh ../asm34_current '8 6'`).
- Medium-mesh checks (`cfd/assembly_previews/meshcheck_medium.txt`): asm34_A,
  asm34_Bvane and asm1e_V1-V4 pass checkMesh (non-orthogonality 63.4-65.0 deg,
  0.67-1.0 M cells). The four long-radius assemblies (asm1lr_V1-V4) have NOT
  been meshed yet with their current, smaller refinement box; the campaign run
  is their first mesh.

## What to do

1. **Find out what is already running, and do not interrupt it.** In WSL:
   `pgrep -fa "run_sweep|run_campaign|foamRun|snappyHexMesh"` and
   `cat ~/turitap_cfd/CURRENT_CASE ~/turitap_cfd/SWEEP`. If another terminal or
   another agent session on this machine is running a sweep, leave it alone:
   one sweep at a time on this machine (each uses 8 cores).
2. **Do not `git pull` while a sweep is running.** The pull replaces scripts the
   running sweep is still reading. Wait until it has finished (`pgrep` above
   shows nothing, the sweep log says `sweep finished`).
3. Report the finished runs (table as in `LAB_SETUP.md` B5 step 1, and say
   plainly which converged). Do not commit or push until the user asks.
4. Then `git pull`, and **ask the user before** starting the campaign from
   `LAB_SETUP.md` Part B2:
   `wsl -d Ubuntu-22.04 -- bash -lc "cd /mnt/c/cfd/Venturi-Redesign/cfd/assembly3d && PUSH_STATUS=1 bash run_campaign.sh"`
   It runs in its own terminal window, not inside the agent session. It runs
   all 11 assemblies at 30 L/min, then 6, then 8 (about 28 h), and skips any
   case that already has `results.json` (so the 8 and 6 L/min runs of the
   current Venturi are not repeated). Open the live view from B3 and give the
   user the phone status link from B4.

## Known issue to watch: long-radius meshing stall

On the first attempt, the long-radius meshes sat in snappyHexMesh "Shell
refinement" for 90+ minutes, adding ~35 cells per pass (the refinement box
was too big; it is now smaller, but that has not been confirmed to fix it).
If an `asm1lr_*` case shows `meshing: snappyHexMesh` for more than about
30-40 minutes (the others take 15-20 min), or its `log.snappyHexMesh` shows
"Shell refinement iteration" above ~20, tell the user, show the last 20 lines
of that log, and ask before stopping it. Do not change mesh settings without
asking. If the user agrees to stop it, kill only that case's snappyHexMesh;
the campaign logs the failure and moves on, and the other assemblies still run.

## Rules (short version)

- No physics, mesh, geometry or solver changes without asking the user.
- Never print full solver logs; read the last ~50 lines or postProcessing.
- checkMesh before every solve (the scripts do this); stop on failure or
  non-orthogonality above 65 deg.
- Ask before installing anything, and before committing or pushing.
