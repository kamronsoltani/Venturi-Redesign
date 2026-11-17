# Running the 3-D CFD on another computer

Step-by-step setup for a second Windows machine (e.g. the lab PC), then how to
run the overnight sweep, watch it, check it from a phone, and bring the
results back. Everything needed is in this repository; read
`PROJECT_NOTES.md` for the physics, the results so far and the plan.

Time needed: about 1 hour of setup (mostly downloads), then the runs.

---

## Part A: one-time setup

### A1. Windows prerequisites

1. Install **Git for Windows**: https://git-scm.com/download/win (defaults are fine).
2. Optional, for releases and PRs: **GitHub CLI**, https://cli.github.com, then run
   `gh auth login` in PowerShell and follow the prompts.
3. Power settings (Settings > System > Power): when plugged in, set **Sleep: Never**
   and **Turn off screen: whatever you like**. A sleep in the middle of a run stops
   WSL and kills the solve. If you can, pause Windows Update for a week so it does
   not reboot overnight.

### A2. Install WSL2 with Ubuntu 22.04

Open **PowerShell as Administrator** and run:

```powershell
wsl --install -d Ubuntu-22.04
```

Reboot when asked. Ubuntu opens and asks for a Linux user name and password
(any; remember the password, `sudo` asks for it).

### A3. Give WSL enough memory and cores

WSL gets only half the RAM by default. Create the file
`C:\Users\<your Windows user>\.wslconfig` with:

```ini
[wsl2]
memory=48GB
processors=16
```

Then in PowerShell run `wsl --shutdown` (WSL restarts with the new limits the
next time it is opened). On a 64 GB machine 48 GB leaves room for Windows.

### A4. Install OpenFOAM 13 and the Python tools (inside Ubuntu)

Open **Ubuntu** from the Start menu and run these one at a time:

```bash
sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc"
sudo add-apt-repository -y http://dl.openfoam.org/ubuntu
sudo apt update
sudo apt -y install openfoam13 python3-pip python3-numpy python3-scipy python3-matplotlib python3-gmsh
python3 -m pip install --user gmsh==4.15.2
```

(`python3-gmsh` from Ubuntu is 4.8 and too old for the geometry booleans, but
installing it pulls in the graphics libraries the pip `gmsh` 4.15.2 needs.)

Check:

```bash
source /opt/openfoam13/etc/bashrc && foamRun -help | head -3
python3 -c "import gmsh, numpy, scipy; print(gmsh.__version__, numpy.__version__, scipy.__version__)"
```

The second line should print `4.15.2 1.21.5 1.8.0`. The reference machine runs
OpenFOAM Foundation 13 (Debian package `openfoam13` 20260624); a newer 13
package is fine.

### A5. Get the repository

In **PowerShell** (not inside OneDrive and with no spaces in the path):

```powershell
mkdir C:\cfd
cd C:\cfd
git clone https://github.com/kamronsoltani/Venturi-Redesign.git
cd Venturi-Redesign
git config user.name "Kamron Soltani"
git config user.email "kamronsoltani@users.noreply.github.com"
```

Inside WSL this folder is `/mnt/c/cfd/Venturi-Redesign`. The scripts copy each
case to `~/turitap_cfd/...` inside Linux to run (fast disk, no spaces) and copy
the results back into the repository folder.

### A6. Smoke test (about 5 minutes)

In Ubuntu:

```bash
cd /mnt/c/cfd/Venturi-Redesign/cfd/assembly3d
source /opt/openfoam13/etc/bashrc
ASM_ITERS=20 python3 setup_asm_case.py /mnt/c/cfd/Venturi-Redesign/cfd/asm34_smoke/Q30 /mnt/c/cfd/Venturi-Redesign/cfd/asm34_current/fluid 30 coarse
bash run_asm.sh /mnt/c/cfd/Venturi-Redesign/cfd/asm34_smoke/Q30
cat /mnt/c/cfd/Venturi-Redesign/cfd/asm34_smoke/Q30/STATUS
```

It should end with `done (iteration 20)` and a `mesh_summary.txt` saying
`Mesh OK`. (Twenty iterations is not a result; this only proves the tools
work.) Then delete the test: `rm -rf /mnt/c/cfd/Venturi-Redesign/cfd/asm34_smoke ~/turitap_cfd/asm34_smoke`.

---

## Part B: the overnight sweep

### B1. What runs

The baseline: 3/4" E-Sig (flat side down) + 3/4" threaded elbow + current
Venturi, at 6, 8, 10, 20, 40, 50 and 60 L/min (30 L/min is already done and is
skipped), then the mesh check at 30 L/min on a coarse and a fine mesh.

`cfd/assembly3d/run_sweep.sh` does every case one after another without
supervision: it meshes once and reuses that mesh for every flow, starts each
flow from the nearest finished one, runs checkMesh before every solve (stops
the case if it fails or non-orthogonality is above 65 degrees), solves 2000
iterations on 8 cores, post-processes, builds the 3-D viewer file and updates
`cfd/results/assembly_3d_stage5b.csv`. A case that already has `results.json`
is skipped, so running the same command again resumes after a crash or reboot.

Expected time: about 45-60 min per flow, so about 6 h for the flows, plus
about 30 min (coarse) and 2-3 h (fine) for the mesh check.

### B2. Start it

Open **Windows Terminal** (or PowerShell), so the run keeps going when nothing
else is open, and run (one line):

```powershell
wsl -d Ubuntu-22.04 -- bash -lc "cd /mnt/c/cfd/Venturi-Redesign/cfd/assembly3d && export PUSH_STATUS=1 LABEL='3/4in E-Sig + 3/4in elbow + current Venturi' && bash run_sweep.sh ../asm34_current '20 10 8 6 40 50 60' && bash run_sweep.sh ../asm34_current '30' coarse && bash run_sweep.sh ../asm34_current '30' fine"
```

The order 20, 10, 8, 6, 40, 50, 60 is deliberate: each flow starts from the
nearest finished one, so stepping outward converges faster.

Leave this window open. Closing it stops the sweep (re-running the same
command resumes it).

### B3. Watch it on the computer

Open a second terminal and run:

```powershell
wsl -d Ubuntu-22.04 -- bash /mnt/c/cfd/Venturi-Redesign/cfd/assembly3d/live_status.sh
```

It refreshes every 15 s: the current case, mesh check, iteration and
residuals, inlet pressure, flow in/out, and the sweep log.

### B4. Check it from a phone

With `PUSH_STATUS=1` the sweep pushes a small progress page every 10 minutes
and after every case to the branch `sweep-status` (it never touches `main`):

https://github.com/kamronsoltani/Venturi-Redesign/blob/sweep-status/STATUS.md

The first time it pushes, Windows may open a GitHub sign-in window on the lab
computer; sign in once with an account that can push to the repository.

It lists each flow with its status, iteration, residual, water level and port
pressure, plus the last lines of the sweep log. If the page stops updating for
more than ~15 minutes, the machine has slept, rebooted or lost the network.

### B5. In the morning: check and publish

1. Read `cfd/asm34_current/sweep.log` and `cfd/results/assembly_3d_stage5b.csv`.
   For each run check: `STATUS` says `done`; `inlet_level_drift_last500_cm`
   is small (under ~0.5 cm); `mass_imbalance_pct` under ~0.5 %;
   `inlet_static_negative_flag` (YES means the real part would draw air at the
   inlet, so the pipe-full answer is not physical; report it, do not hide it).
2. Commit and push (raw fields, meshes and logs are ignored by `.gitignore`;
   the case dictionaries, results, CSV and viewer files are kept):

   ```powershell
   cd C:\cfd\Venturi-Redesign
   git pull --rebase
   git add -A
   git commit -m "Baseline sweep: current Venturi behind 3/4in E-Sig + elbow, 6-60 L/min, plus mesh check"
   git push
   ```

3. A minute later the 3-D viewer shows every flow:
   https://kamronsoltani.github.io/Venturi-Redesign/3d/

Two computers can work at the same time as long as each one runs different
assembly folders (`cfd/asm34_current`, `cfd/asm1e_V1`, ...); always
`git pull --rebase` before pushing.

---

## Part C: next assemblies (after the baseline)

See "Stage 5b" in `PROJECT_NOTES.md`. Each new assembly needs its water volume
first: `cfd/assembly3d/build_fluid.py` currently builds only the 3/4" E-Sig +
3/4" elbow + current Venturi (Venturi bore from `geometry/bore_stations.csv`).
For the other Venturis behind the same parts, change the Venturi tag it reads
and the output folder; for the 1" parts, the E-Sig and elbow sections need the
1" STEP files in `geometry/upstream/`. Check every new volume with
`showfluid.py` (picture of the symmetry plane) before meshing, then run it with
`run_sweep.sh` exactly as above.

---

## Rules kept throughout the project

- Never print full solver logs; read the last ~50 lines or `postProcessing`.
- checkMesh before every solve; stop on failure or non-orthogonality above 65 degrees.
- Do not tune settings to match the 1-D model; report disagreements.
- One case folder per geometry and flow rate under `cfd/`.
- Summarise results in a CSV under `cfd/results/` and say plainly which runs converged.
- Ask before installing anything beyond Part A.

---

## Prompt to paste into the coding agent on the lab computer

```text
You are continuing a CFD project in C:\cfd\Venturi-Redesign (GitHub:
kamronsoltani/Venturi-Redesign). Read LAB_SETUP.md and PROJECT_NOTES.md fully
before doing anything, and follow the rules at the end of LAB_SETUP.md.

1. Check whether Part A of LAB_SETUP.md is done on this machine (WSL
   Ubuntu-22.04, OpenFOAM 13 at /opt/openfoam13, gmsh 4.15.2 for python3 in
   WSL, the repo cloned at C:\cfd\Venturi-Redesign, .wslconfig memory). For
   anything missing, tell me what is missing and ask before installing it,
   then walk me through it step by step. Run the smoke test in A6 and show me
   the result.
2. Start the overnight sweep from Part B2 in a separate terminal window that
   keeps running on its own (not inside your session), and open the live
   status view from B3 in a second terminal. Confirm the status page in B4
   updates (give me the link).
3. Do not change any physics, mesh or solver settings. Do not commit or push
   until I ask. If a case fails, say which, show the last 30 lines of its log,
   and stop; do not change settings to make it pass without asking me.
4. When the sweep finishes, do Part B5 step 1: give me a table per flow
   (water level above the inlet centreline, static pressure at the inlet with
   any negative-pressure flag, port wall pressure mean / E-Sig side / opposite
   side, losses per part, outlet KE factor, diffuser reversed fraction,
   residuals and convergence), the mesh-check spread at 30 L/min, and say
   plainly which runs converged. Then ask me before committing and pushing.
```
