"""Shared reporting: every script states computed vs manuscript, and PASS/FAIL."""
import os, csv
RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
os.makedirs(RESULTS, exist_ok=True)
_checks = []

def banner(n, title):
    print(f"\n{'='*78}\n{n}. {title}\n{'='*78}")
    print(f"    {'quantity':<46s} {'computed':>10s} {'claimed':>10s}")

def check(label, computed, claimed, tol=5e-4, unit=""):
    ok = claimed is None or abs(computed - claimed) <= tol
    verdict = "       " if claimed is None else ("PASS" if ok else ">>> FAIL")
    shown = " " * 10 if claimed is None else f"{claimed:10.4f}"
    print(f"    {label:<46s} {computed:10.4f} {shown} {unit:<12s} {verdict}")
    _checks.append((label, round(float(computed), 6),
                    "" if claimed is None else claimed,
                    "" if claimed is None else ("PASS" if ok else "FAIL")))
    return ok

def write_csv(name, header, rows):
    with open(os.path.join(RESULTS, name), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(header); w.writerows(rows)
    print(f"    -> results/{name}")

def finish(tag):
    path = os.path.join(RESULTS, "00_checks.csv")
    new = not os.path.exists(path)
    with open(path, "a", newline="") as fh:
        w = csv.writer(fh)
        if new: w.writerow(["script","quantity","computed","manuscript_claim","verdict"])
        for row in _checks: w.writerow([tag] + list(row))
    fails = [c for c in _checks if c[3] == "FAIL"]
    print(f"\n    {len([c for c in _checks if c[3]=='PASS'])} passed, {len(fails)} failed")
    for label, computed, claimed, _ in fails:
        print(f"    FAIL  {label}: computed {computed}, manuscript says {claimed}")
    return 1 if fails else 0
