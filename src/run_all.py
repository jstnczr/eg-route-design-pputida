#!/usr/bin/env python3
"""
Runs the whole analysis in the order the manuscript presents it.

    python src/run_all.py

Each script writes its own CSV into results/ and prints what it found.
Scripts are independent: you can run any one on its own.
"""
import subprocess, sys, os, glob
HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(os.path.dirname(HERE), "results"), exist_ok=True)
scripts = sorted(glob.glob(os.path.join(HERE, "0[1-9]_*.py")))
fails = 0
for s in scripts:
    name = os.path.basename(s)
    print(f"\n{'#'*78}\n#  {name}\n{'#'*78}")
    fails += subprocess.call([sys.executable, s])
print(f"\n{'='*78}")
print(f"{len(scripts)} steps run, {fails} reported a problem.")
print("Outputs in results/.  Each step's own output says what it checked.")
sys.exit(1 if fails else 0)
