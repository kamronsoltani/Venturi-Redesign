# Builds the review picture for one assembly whose fluid_mm.stl is not in git.
# Works on a temp copy; nothing in the repo is changed.
A=${1:-asm34_current}
SRC=/mnt/c/cfd/Venturi-Redesign/cfd/$A/fluid
T=/tmp/preview_$A/fluid
rm -rf /tmp/preview_$A; mkdir -p $T
cp $SRC/assembly_info.json $T/
python3 - "$SRC/fluid.stl" "$T/fluid_mm.stl" <<'EOF'
import sys
src, dst = sys.argv[1], sys.argv[2]
lines = open(src).read().splitlines()
vals = [abs(float(v)) for l in lines if l.strip().startswith("vertex") for v in l.split()[1:]]
s = 1000.0 if max(vals) < 5 else 1.0     # metres -> mm if needed
print("max coord", max(vals), "scale", s)
with open(dst, "w") as f:
    for l in lines:
        t = l.strip()
        if t.startswith("vertex"):
            x, y, z = (float(v) * s for v in t.split()[1:])
            f.write(f"      vertex {x:.6g} {y:.6g} {z:.6g}\n")
        else:
            f.write(l + "\n")
EOF
cd /mnt/c/cfd/Venturi-Redesign/cfd/assembly3d
nice -n 19 python3 review_assembly.py $T /mnt/c/cfd/lab_agent/previews/$A.png $A 2>&1 | tail -3
ls -l /mnt/c/cfd/lab_agent/previews/$A.png
