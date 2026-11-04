"""
Compare several Neware GCPL files at one cycle (sodium-ion half-cells, hard carbon vs Na).

Figures
-------
compare_discharge_profiles_cycle<N>.png   sodiation (discharge) curves of all samples on one graph
compare_slope_plateau_cycle<N>.png        sodiation capacity split into slope (> threshold) and
                                          plateau (< threshold) as a fraction of the total
plus compare_cycle<N>.xlsx with the numbers and the profile data (for OriginPro).

Cycle counting follows analyse_gcpl.py (CYCLE_START): with "desodiation" the initial sodiation is
not counted, and cycle N = desodiation N + the sodiation after it. The discharge of cycle N is that
sodiation.

Usage
-----
    python compare_samples.py                                   # all Neware .xlsx files in this folder
    python compare_samples.py file1.xlsx file2.xlsx --cycle 10
    python compare_samples.py --cycle 10 --threshold 0.1 --mass C=0.50   # correct one sample's mass

Needs analyse_gcpl.py in the same folder. Requirements: pandas, numpy, matplotlib, openpyxl
"""

import argparse
import contextlib
import io
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import analyse_gcpl as gcpl

# ============================ USER SETTINGS ==================================
CYCLE = 10                     # cycle to compare
PLATEAU_THRESHOLD_V = 0.1      # below = plateau, above = slope
CYCLE_START = "desodiation"    # same meaning as in analyse_gcpl.py
SAMPLE_ORDER = None            # e.g. ["A", "Q", "O", "F", "C"]; None = alphabetical
OUT_DIR = "comparison"
# =============================================================================

COLORS = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"]


def load_sample(path, cycle, start, mass_override=None):
    with contextlib.redirect_stdout(io.StringIO()) as log:
        mass, name = gcpl.mass_and_name(path)
    mass = mass_override or mass
    if not mass:
        raise ValueError(f"{path.name}: no mass found - give it with --mass {name or 'X'}=<mg>")
    mass_g = mass / 1000
    halves = gcpl.build_half_cycles(gcpl.load_detail(path), mass_g)
    cycles = gcpl.build_cycles(halves, start)
    c = next((cc for cc in cycles if cc["cycle"] == cycle), None)
    if c is None or c.get("sod") is None:
        raise ValueError(f"{path.name}: cycle {cycle} has no sodiation (the test has {len(cycles)} cycles)")
    h = c["sod"]
    if not h["complete"]:
        print(f"WARNING: {path.name}: the sodiation of cycle {cycle} did not reach the cut-off")
    q, v = h["Q"], h["V"]
    q_slope = gcpl.crossing_capacity(q, v, PLATEAU_THRESHOLD_V, "down")
    total = q[-1] / mass_g
    slope = q_slope / mass_g
    i_med = np.median([hh["I_mA"] for hh in halves])
    return {
        "name": name or path.stem, "file": path.name, "mass_mg": mass, "Q": q / mass_g, "V": v,
        "total": total, "slope": slope, "plateau": total - slope,
        "steps": ", ".join(str(s) for s in h["steps"]), "current_uA": i_med * 1000,
        "c_rate": i_med / (mass_g * gcpl.C_RATE_BASIS_MAH_G),
    }


def fig_profiles(samples, cycle, out):
    fig, ax = plt.subplots(figsize=(6.4, 4.5))
    for k, s in enumerate(samples):
        ax.plot(s["Q"], s["V"], color=COLORS[k % len(COLORS)], lw=1.6,
                label=f"Sample {s['name']} ({s['total']:.0f} mAh g$^{{-1}}$)")
    ax.axhline(PLATEAU_THRESHOLD_V, color="0.5", lw=0.8, ls=":")
    ax.text(ax.get_xlim()[1] * 0.99 if ax.get_xlim()[1] else 1, PLATEAU_THRESHOLD_V + 0.03,
            f"{PLATEAU_THRESHOLD_V:g} V", ha="right", va="bottom", fontsize=8, color="0.4")
    ax.set_xlim(left=0)
    ax.set_ylim(0, gcpl.DESODIATION_CUTOFF_V + 0.05)
    ax.set_xlabel(gcpl.CAP_LABEL)
    ax.set_ylabel(gcpl.POT_LABEL)
    ax.legend(fontsize=8.5, loc="upper right")
    ax.set_title(f"Sodiation (discharge) profiles, cycle {cycle}", fontsize=10)
    f = out / f"compare_discharge_profiles_cycle{cycle}.png"
    fig.savefig(f)
    plt.close(fig)
    return f


def fig_fraction(samples, cycle, out):
    names = [f"Sample {s['name']}" for s in samples]
    x = np.arange(len(samples))
    fs = np.array([100 * s["slope"] / s["total"] for s in samples])
    fp = 100 - fs
    fig, ax = plt.subplots(figsize=(1.25 * len(samples) + 2.2, 4.6))
    ax.bar(x, fs, 0.62, color="#9ecae1", edgecolor="0.3", lw=0.6, label=f"Slope (> {PLATEAU_THRESHOLD_V:g} V)")
    ax.bar(x, fp, 0.62, bottom=fs, color="#3182bd", edgecolor="0.3", lw=0.6,
           label=f"Plateau (< {PLATEAU_THRESHOLD_V:g} V)")
    for xi, s, a, b in zip(x, samples, fs, fp):
        ax.text(xi, a / 2, f"{a:.0f} %\n{s['slope']:.0f}", ha="center", va="center", fontsize=8)
        ax.text(xi, a + b / 2, f"{b:.0f} %\n{s['plateau']:.0f}", ha="center", va="center", fontsize=8,
                color="white")
        ax.text(xi, 101.5, f"{s['total']:.0f} mAh g$^{{-1}}$", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x, names)
    ax.set_ylim(0, 112)
    ax.set_yticks(range(0, 101, 20))
    ax.set_ylabel("Share of sodiation capacity (%)")
    ax.tick_params(top=False)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2)
    ax.set_title(f"Slope vs plateau capacity, cycle {cycle} (numbers in segments: mAh g$^{{-1}}$)", fontsize=10)
    f = out / f"compare_slope_plateau_cycle{cycle}.png"
    fig.savefig(f)
    plt.close(fig)
    return f


def main():
    global PLATEAU_THRESHOLD_V
    p = argparse.ArgumentParser(description="Compare GCPL files at one cycle")
    p.add_argument("files", nargs="*", help="Neware .xlsx files (default: all .xlsx in this folder)")
    p.add_argument("--cycle", type=int, default=CYCLE)
    p.add_argument("--threshold", type=float, default=PLATEAU_THRESHOLD_V)
    p.add_argument("--start", choices=["desodiation", "sodiation"], default=CYCLE_START)
    p.add_argument("--mass", nargs="*", default=[], help="override masses, e.g. C=0.50 A=0.85")
    p.add_argument("--out", default=OUT_DIR)
    args = p.parse_args()
    PLATEAU_THRESHOLD_V = args.threshold
    overrides = {k.strip().upper(): float(v) for k, v in (m.split("=") for m in args.mass)}

    files = [Path(f) for f in args.files] or sorted(Path(".").glob("*.xlsx"))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    samples = []
    for f in files:
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                _, nm = gcpl.mass_and_name(f)
            samples.append(load_sample(f, args.cycle, args.start, overrides.get((nm or "").upper())))
        except Exception as e:
            print(f"Skipped {f.name}: {e}")
    if not samples:
        raise SystemExit("No usable files.")
    order = SAMPLE_ORDER or sorted(s["name"] for s in samples)
    samples.sort(key=lambda s: order.index(s["name"]) if s["name"] in order else len(order))

    print(f"Cycle {args.cycle} ({'initial sodiation not counted' if args.start == 'desodiation' else 'standard counting'}), "
          f"threshold {PLATEAU_THRESHOLD_V:g} V\n")
    rows = []
    for s in samples:
        rows.append({
            "sample": s["name"], "file": s["file"], "mass (mg)": s["mass_mg"],
            "current (uA)": round(s["current_uA"], 2), "C-rate (1C = 250 mAh/g)": f"C/{1 / s['c_rate']:.1f}",
            "Neware step": s["steps"], "sodiation (mAh/g)": s["total"],
            f"slope > {PLATEAU_THRESHOLD_V:g} V (mAh/g)": s["slope"],
            f"plateau < {PLATEAU_THRESHOLD_V:g} V (mAh/g)": s["plateau"],
            "slope share (%)": 100 * s["slope"] / s["total"], "plateau share (%)": 100 * s["plateau"] / s["total"],
        })
        print(f"  Sample {s['name']:<3} {s['mass_mg']:.2f} mg  C/{1 / s['c_rate']:<5.1f} "
              f"sodiation {s['total']:6.1f} mAh/g   slope {s['slope']:6.1f} ({100 * s['slope'] / s['total']:.0f} %)   "
              f"plateau {s['plateau']:6.1f} ({100 * s['plateau'] / s['total']:.0f} %)")
        if abs(1 / s["c_rate"] - 10) > 1:
            print(f"    NOTE: the current ({s['current_uA']:.2f} uA) is not C/10 for {s['mass_mg']} mg - "
                  f"check the mass of sample {s['name']}")
    table = pd.DataFrame(rows)
    profiles = {}
    for s in samples:
        profiles[f"Sample {s['name']} Q (mAh/g)"] = pd.Series(s["Q"])
        profiles[f"Sample {s['name']} V (V)"] = pd.Series(s["V"])

    f1 = fig_profiles(samples, args.cycle, out)
    f2 = fig_fraction(samples, args.cycle, out)
    gcpl.write_excel(out / f"compare_cycle{args.cycle}.xlsx",
                     {"Summary": table, "Sodiation profiles": pd.DataFrame(profiles)})
    print(f"\nSaved: {f1.name}, {f2.name}, compare_cycle{args.cycle}.xlsx in {out.resolve()}")


if __name__ == "__main__":
    main()
