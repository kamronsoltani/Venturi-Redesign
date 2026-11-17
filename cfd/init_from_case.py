#!/usr/bin/env python3
"""
init_from_case.py -- replace the internalField of a case's 0/ fields with the
converged internalField of another case on the SAME mesh (same geometry and
mesh level). Boundary conditions, solver settings and schemes are untouched;
only the starting guess changes. U is scaled by the flow-rate ratio so the
initial field is consistent with the target inlet flow.

Used when a run diverges from the uniform start (e.g. B_3-4in Q8, FPE at
iteration 162 on 2026-11-13).

Usage: python3 init_from_case.py <target case> <source case> <q_target> <q_source>
"""
import os
import re
import sys

FIELDS = ["U", "p", "k", "omega", "nut"]


def internal_block(text):
    a = text.index("internalField")
    b = text.index("\nboundaryField")
    return a, b


def latest_time(case):
    times = [d for d in os.listdir(case) if re.fullmatch(r"[0-9.]+", d) and float(d) > 0]
    return max(times, key=float)


def main():
    tgt, src, qt, qs = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
    st = latest_time(src)
    scale = qt / qs
    for f in FIELDS:
        s_text = open(os.path.join(src, st, f)).read()
        a, b = internal_block(s_text)
        block = s_text[a:b]
        if f == "U":
            block = re.sub(r"\(([-0-9.eE+]+) ([-0-9.eE+]+) ([-0-9.eE+]+)\)",
                           lambda m: "(%.10g %.10g %.10g)" % tuple(float(x) * scale for x in m.groups()),
                           block)
        t_path = os.path.join(tgt, "0", f)
        t_text = open(t_path).read()
        ta, tb = internal_block(t_text)
        with open(t_path, "w", newline="\n") as fh:
            fh.write(t_text[:ta] + block + t_text[tb:])
    with open(os.path.join(tgt, "case_info.txt"), "a", newline="\n") as fh:
        fh.write(f"initialised_from {os.path.basename(os.path.dirname(os.path.abspath(src)))}/"
                 f"{os.path.basename(os.path.abspath(src))} time {st} (U scaled x{scale:.4g})\n")
    print(f"{tgt}: internal fields from {src}/{st}, U x{scale:.4g}")


if __name__ == "__main__":
    main()
