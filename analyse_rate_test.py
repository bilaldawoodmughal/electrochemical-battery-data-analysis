"""
Rate-capability analysis for Na-ion half-cells (hard carbon vs Na metal)
exported from a Neware/BTS cycler as .xlsx (sheets: Info, Cycle_*, Statis_*, Detail_*).

What it does
------------
1. Reads the Detail sheet (every data point) and splits it into half-cycles:
   every sodiation (discharge, voltage going down) and every desodiation
   (charge, voltage going up). The Neware "Cycle" sheet is not used.
2. Counts the cycles the way set in CYCLE_START:
     "desodiation": the first sodiation is the initial (formation) sodiation and
                    is not counted. Cycle 1 = first desodiation + the sodiation
                    after it, cycle 2 = next desodiation + next sodiation, ...
     "sodiation":   standard counting, cycle 1 = first sodiation + first desodiation.
3. Labels every half-cycle with its own C-rate (detected from the current).
4. Calculates capacities (mAh, mAh/g), coulombic efficiency, average potentials,
   voltage hysteresis and the slope/plateau capacity split.
   Coulombic efficiency is always: desodiation / the sodiation just before it
   (sodium taken out / sodium put in). For the first desodiation this is the ICE.
5. Summarises every C-rate block (mean +/- SD, retention, recovery).
6. Saves 4 figures (PNG, 300 dpi) and an Excel file with all numbers
   (per cycle, per C-rate, voltage profiles) ready for OriginPro.

Naming used in a half-cell
--------------------------
Discharge ("DChg") = sodiation of the hard carbon  (voltage goes down to ~5 mV)
Charge    ("Chg")  = desodiation of the hard carbon (voltage goes up to 2.5 V)
The desodiation capacity is the reversible capacity.

Usage
-----
    python analyse_rate_test.py                       # uses the settings below
    python analyse_rate_test.py myfile.xlsx --mass 1.09
    python analyse_rate_test.py myfile.xlsx --mass 1.09 --rates 0.05 0.1 0.2 0.5 1 2 0.05

Requirements: pandas, numpy, matplotlib, openpyxl
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

# ============================ USER SETTINGS ==================================
FILE = "rate test_Bilal_sample F_0.66mg-240040-2-1-0.xlsx"
ACTIVE_MASS_MG = None          # None = read from file name ("1_09mg", "1,09mg" or "0.66mg")
RATE_SEQUENCE = [0.05, 0.1, 0.2, 0.5, 1, 2, 0.05]   # C-rates in the order you ran them
C_RATE_BASIS_MAH_G = 250       # only used if RATE_SEQUENCE does not match the data
CYCLE_START = "desodiation"    # "desodiation": initial sodiation not counted, cycle = desod + next sod
                               # "sodiation":   standard, cycle = sod + desod
SODIATION_CUTOFF_V = 0.005     # lower cut-off (V vs Na+/Na)
DESODIATION_CUTOFF_V = 2.5     # upper cut-off (V vs Na+/Na)
PLATEAU_THRESHOLD_V = 0.1      # below this potential = plateau, above = slope
OUT_DIR = None                 # None = "<file name>_analysis" next to the file
SAMPLE_NAME = None             # None = taken from the file name (e.g. "Sample A, 1.09 mg")
# =============================================================================

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 10,
    "axes.linewidth": 1.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "legend.frameon": False,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})
CAP_LABEL = r"Specific capacity (mAh g$^{-1}$)"
POT_LABEL = r"Potential (V vs Na$^{+}$/Na)"
NAMES = {"sod": "sodiation", "desod": "desodiation"}


# ------------------------------- helpers -------------------------------------
def find_col(df, *keys):
    """Return the first column whose name contains all keys (case-insensitive)."""
    for c in df.columns:
        name = str(c).lower()
        if all(k.lower() in name for k in keys):
            return c
    raise KeyError(f"No column containing {keys} in {list(df.columns)}")


def current_to_mA(series, col_name):
    """Neware exports current in uA or mA depending on settings."""
    name = col_name.lower().replace("\u00b5", "u").replace("\u03bc", "u")
    factor = 1e-3 if "(ua)" in name else 1.0
    return series.astype(float) * factor


def mass_from_filename(path):
    m = re.search(r"(\d+)[_,.](\d+)\s*mg", Path(path).name, re.IGNORECASE)
    if m:
        return float(f"{m.group(1)}.{m.group(2)}")
    m = re.search(r"(\d+)\s*mg", Path(path).name, re.IGNORECASE)
    return float(m.group(1)) if m else None


def nice_sample_name(path, mass_mg):
    m = re.search(r"sample[_ ]?([A-Za-z0-9]+)", Path(path).name, re.IGNORECASE)
    name = f"Sample {m.group(1)}" if m else Path(path).stem.split("-")[0]
    return f"{name}, {mass_mg:g} mg"


def rate_label(rate):
    return f"{rate:g}C"


def crossing_capacity(q, v, threshold, direction):
    """Capacity at which the potential first crosses `threshold`.
    direction='down' for sodiation, 'up' for desodiation. Linear interpolation."""
    idx = np.where(v <= threshold)[0] if direction == "down" else np.where(v >= threshold)[0]
    if len(idx) == 0:
        return q[-1]            # never crossed: all capacity on one side
    i = idx[0]
    if i == 0:
        return 0.0
    v0, v1, q0, q1 = v[i - 1], v[i], q[i - 1], q[i]
    if v1 == v0:
        return q1
    return q0 + (threshold - v0) * (q1 - q0) / (v1 - v0)


def average_potential(q, v):
    """Capacity-weighted average potential (integral V dQ / Q)."""
    if q[-1] <= 0:
        return np.nan
    integrate = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return integrate(v, q) / q[-1]


# ------------------------------ loading --------------------------------------
def load_detail(path):
    xl = pd.ExcelFile(path)
    sheet = next(s for s in xl.sheet_names if s.lower().startswith("detail"))
    det = pd.read_excel(path, sheet_name=sheet)
    c_status = find_col(det, "status")
    c_step = find_col(det, "steps")
    c_volt = find_col(det, "voltage")
    c_cur = find_col(det, "current")
    c_cap = next(c for c in det.columns if str(c).lower().startswith("capacity"))
    c_rec = find_col(det, "record")
    out = pd.DataFrame({
        "record": det[c_rec],
        "status": det[c_status].astype(str),
        "step": det[c_step].astype(int),
        "V": det[c_volt].astype(float),
        "I_mA": current_to_mA(det[c_cur], str(c_cur)),
        "Q_mAh": det[c_cap].astype(float),
    }).sort_values("record").reset_index(drop=True)
    s = out["status"].str.lower()
    out["direction"] = np.where(s.str.contains("dchg"), "sod",
                                np.where(s.str.contains("chg"), "desod", "rest"))
    return out


def build_half_cycles(det, mass_g):
    """One entry per sodiation or desodiation, in time order. Consecutive steps with
    the same direction (e.g. CC + CV) are merged into one half-cycle."""
    halves = []
    for step, g in det.groupby("step", sort=True):
        d = g["direction"].iloc[0]
        if d == "rest":
            continue
        q = g["Q_mAh"].to_numpy(dtype=float)
        v = g["V"].to_numpy(dtype=float)
        i_med = float(np.median(np.abs(g["I_mA"].to_numpy())))
        if halves and halves[-1]["dir"] == d:          # merge CC + CV etc.
            h = halves[-1]
            h["Q"] = np.concatenate([h["Q"], q + h["Q"][-1]])
            h["V"] = np.concatenate([h["V"], v])
            h["steps"].append(step)
        else:
            halves.append({"dir": d, "Q": q, "V": v, "I_mA": i_med, "steps": [step]})

    for h in halves:                                    # per half-cycle numbers
        q, v = h["Q"], h["V"]
        h["mAh"] = q[-1]
        h["spec"] = q[-1] / mass_g
        h["avgV"] = average_potential(q, v)
        if h["dir"] == "sod":
            q_slope = crossing_capacity(q, v, PLATEAU_THRESHOLD_V, "down")
            h["slope"], h["plateau"] = q_slope / mass_g, (q[-1] - q_slope) / mass_g
            h["complete"] = bool(v.min() <= SODIATION_CUTOFF_V + 0.003)
        else:
            q_plat = crossing_capacity(q, v, PLATEAU_THRESHOLD_V, "up")
            h["plateau"], h["slope"] = q_plat / mass_g, (q[-1] - q_plat) / mass_g
            h["complete"] = bool(v.max() >= DESODIATION_CUTOFF_V - 0.01)

    first_sod = next((h for h in halves if h["dir"] == "sod"), None)
    for j, h in enumerate(halves):                      # CE = desodiation / sodiation before it
        h["is_first_sod"] = h is first_sod
        h["CE"], h["hyst"], h["is_ice"] = np.nan, np.nan, False
        if h["dir"] == "desod" and j > 0 and halves[j - 1]["dir"] == "sod":
            prev = halves[j - 1]
            if h["complete"] and prev["complete"]:
                h["CE"] = 100 * h["mAh"] / prev["mAh"]
                h["hyst"] = h["avgV"] - prev["avgV"]
            h["is_ice"] = prev is first_sod
    return halves


def assign_rates(halves, mass_g, rate_sequence, basis_mAh_g):
    """Split half-cycles into blocks of constant current and give each block a C-rate."""
    block = 0
    for k, h in enumerate(halves):
        if k > 0 and abs(h["I_mA"] - halves[k - 1]["I_mA"]) / halves[k - 1]["I_mA"] > 0.10:
            block += 1
        h["block"] = block + 1
    n_blocks = block + 1
    block_current = {b: float(np.median([h["I_mA"] for h in halves if h["block"] == b]))
                     for b in range(1, n_blocks + 1)}

    if rate_sequence and len(rate_sequence) == n_blocks:
        labels = {b: float(rate_sequence[b - 1]) for b in block_current}
        source = "rate sequence given in settings"
    else:
        one_c = mass_g * basis_mAh_g            # mA
        labels = {b: round(block_current[b] / one_c, 3) for b in block_current}
        source = f"current / (mass x {basis_mAh_g} mAh/g)"
        if rate_sequence:
            print(f"WARNING: {len(rate_sequence)} rates given but {n_blocks} current blocks found; "
                  f"labels calculated from the current instead.")
    for h in halves:
        h["rate"] = labels[h["block"]]
    implied = np.mean([block_current[b] / labels[b] / mass_g for b in block_current])
    print(f"C-rate labels from: {source}")
    print(f"Implied 1C basis: {implied:.0f} mAh/g  (1C = {implied * mass_g:.4f} mA)")
    return n_blocks


def build_cycles(halves, start):
    """Group half-cycles into numbered cycles.
    start='desodiation': sodiation(s) before the first desodiation are the initial
    sodiation and are not counted; each cycle = desodiation + the next sodiation.
    start='sodiation': each cycle = sodiation + the next desodiation."""
    first, second = ("desod", "sod") if start == "desodiation" else ("sod", "desod")
    i0 = next((i for i, h in enumerate(halves) if h["dir"] == first), len(halves))
    initial, seq = halves[:i0], halves[i0:]
    cycles, i = [], 0
    while i < len(seq):
        h = seq[i]
        if h["dir"] == first:
            nxt = seq[i + 1] if i + 1 < len(seq) and seq[i + 1]["dir"] == second else None
            cycles.append({first: h, second: nxt})
            i += 2 if nxt else 1
        else:
            cycles.append({first: None, second: h})
            i += 1
    for n, c in enumerate(cycles, start=1):
        c["cycle"] = n
        lead = c[first] or c[second]
        c["rate"], c["block"] = lead["rate"], lead["block"]     # cycle belongs to its first half
        for key in ("sod", "desod"):
            if c.get(key) is not None:
                c[key]["cycle"] = n
    for h in initial:
        h["cycle"] = 0
    return initial, cycles


# ------------------------------ analysis -------------------------------------
def per_cycle_table(cycles, start):
    order = ["desod", "sod"] if start == "desodiation" else ["sod", "desod"]
    rows = []
    for c in cycles:
        r = {"cycle": c["cycle"], "cycle C-rate": c["rate"]}
        for key in order:
            h, nm = c.get(key), NAMES[key]
            if h is None:
                continue
            r[f"{nm} C-rate"] = h["rate"]
            r[f"{nm} (mAh)"] = h["mAh"]
            r[f"{nm} (mAh/g)"] = h["spec"]
            r[f"{nm} slope >thr (mAh/g)"] = h["slope"]
            r[f"{nm} plateau <thr (mAh/g)"] = h["plateau"]
            r[f"avg. {nm} potential (V)"] = h["avgV"]
            r[f"{nm} complete"] = h["complete"]
            if key == "desod":
                r["CE (%) = desod / previous sod"] = h["CE"]
                r["hysteresis (V)"] = h["hyst"]
        r["Neware steps"] = ", ".join(str(s) for key in order if c.get(key) is not None
                                      for s in c[key]["steps"])
        rows.append(r)
    return pd.DataFrame(rows)


def per_rate_table(halves):
    blocks, first_mean = [], None
    for b in sorted({h["block"] for h in halves}):
        des = [h for h in halves if h["dir"] == "desod" and h["block"] == b and h["complete"]]
        sod = [h for h in halves if h["dir"] == "sod" and h["block"] == b and h["complete"]
               and not h["is_first_sod"]]
        if not des:
            continue
        ce = [h["CE"] for h in des if not h["is_ice"] and np.isfinite(h["CE"])]
        hy = [h["hyst"] for h in des if not h["is_ice"] and np.isfinite(h["hyst"])]
        dv = np.array([h["spec"] for h in des])
        if first_mean is None:
            first_mean = dv.mean()
        blocks.append({
            "block": b,
            "C-rate": des[0]["rate"],
            "cycles (desodiation)": f"{min(h['cycle'] for h in des)}-{max(h['cycle'] for h in des)}",
            "n desodiations": len(des),
            "desodiation mean (mAh/g)": dv.mean(),
            "desodiation std (mAh/g)": dv.std(ddof=1) if len(dv) > 1 else 0.0,
            "desodiation last (mAh/g)": dv[-1],
            "n sodiations": len(sod),
            "sodiation mean (mAh/g)": np.mean([h["spec"] for h in sod]) if sod else np.nan,
            "sod. slope mean (mAh/g)": np.mean([h["slope"] for h in sod]) if sod else np.nan,
            "sod. plateau mean (mAh/g)": np.mean([h["plateau"] for h in sod]) if sod else np.nan,
            "CE mean (%)": np.mean(ce) if ce else np.nan,
            "hysteresis mean (V)": np.mean(hy) if hy else np.nan,
            "retention vs block 1 (%)": 100 * dv.mean() / first_mean,
        })
    return pd.DataFrame(blocks)


def formation_info(halves, cycles, start):
    first_sod = next(h for h in halves if h["is_first_sod"])
    ice_des = next((h for h in halves if h["is_ice"]), None)
    return {"first_sod": first_sod, "first_desod": ice_des,
            "ice": ice_des["CE"] if ice_des is not None else np.nan,
            "cycle1": cycles[0] if cycles else None, "start": start}


# ------------------------------- figures -------------------------------------
def rate_colors(halves):
    rates = sorted({h["rate"] for h in halves})
    cmap = plt.get_cmap("viridis", max(len(rates), 2))
    return {r: cmap(i / max(len(rates) - 1, 1) * 0.9) for i, r in enumerate(rates)}


def block_names(per_rate):
    names, seen = [], set()
    for r in per_rate["C-rate"]:
        names.append(rate_label(r) + (" (rec.)" if r in seen else ""))
        seen.add(r)
    return names


def fig_rate_capability(cycles, colors, info, sample, out):
    counted = [c[k] for c in cycles for k in ("sod", "desod") if c.get(k) is not None and c[k]["complete"]]
    fig, (ax_ce, ax) = plt.subplots(2, 1, figsize=(7.0, 5.0), sharex=True,
                                    gridspec_kw={"height_ratios": [1, 3], "hspace": 0.06})

    # background bands and labels: a cycle belongs to the C-rate of its first half-cycle
    blocks = {}
    for c in cycles:
        blocks.setdefault(c["block"], []).append(c)
    for i, b in enumerate(sorted(blocks)):
        cs = blocks[b]
        x0, x1 = cs[0]["cycle"] - 0.5, cs[-1]["cycle"] + 0.5
        if i % 2 == 0:
            for a in (ax_ce, ax):
                a.axvspan(x0, x1, color="0.93", zorder=0)
        ax_ce.text((x0 + x1) / 2, 1.04, rate_label(cs[0]["rate"]), transform=ax_ce.get_xaxis_transform(),
                   ha="center", va="bottom", fontsize=9)

    # capacities; every marker has the colour of its own C-rate
    normal = [h["spec"] for h in counted if not h["is_first_sod"]]
    y_top = 1.3 * max(normal)
    for h in counted:
        col = colors[h["rate"]]
        y = min(h["spec"], y_top * 0.97)
        if h["dir"] == "desod":
            ax.plot(h["cycle"], y, "o", color=col, ms=6, zorder=3)
        else:
            ax.plot(h["cycle"], y, "o", mfc="white", mec=col, ms=6, zorder=3)
        if h["is_first_sod"] and h["spec"] > y_top * 0.97:
            ax.annotate(f"1st sodiation: {h['spec']:.0f} mAh/g (off scale)", xy=(h["cycle"], y),
                        xytext=(h["cycle"] + 1.2, y_top * 0.99), fontsize=8, va="top",
                        arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.set_ylim(0, y_top)
    ax.set_xlim(0.3, cycles[-1]["cycle"] + 0.7)
    ax.set_xlabel("Cycle number")
    ax.set_ylabel(CAP_LABEL)
    handles = [plt.Line2D([], [], ls="", marker="o", color="0.2", label="Desodiation (reversible)"),
               plt.Line2D([], [], ls="", marker="o", mfc="white", mec="0.2", label="Sodiation")]
    ax.legend(handles=handles, loc="lower left", fontsize=8)

    # coulombic efficiency (the ICE is written as text because it is far lower)
    ce = [(h["cycle"], h["CE"]) for h in counted if h["dir"] == "desod" and not h["is_ice"]
          and np.isfinite(h["CE"])]
    if ce:
        x, y = zip(*ce)
        ax_ce.plot(x, y, "D", color="0.3", ms=3.5, zorder=3)
        ax_ce.set_ylim(np.floor(min(y) / 5) * 5 - 2, 102)
    ax_ce.set_ylabel("CE (%)")
    ax_ce.text(0.988, 0.08, f"ICE = {info['ice']:.1f} %", transform=ax_ce.transAxes,
               fontsize=8, va="bottom", ha="right")
    ax_ce.set_title(sample, fontsize=10, pad=20)
    fig.savefig(out / "01_rate_capability.png")
    plt.close(fig)


def profile_pair(halves, block):
    """Last complete desodiation of a block and the sodiation just before it."""
    des = [h for h in halves if h["dir"] == "desod" and h["block"] == block and h["complete"]]
    if not des:
        return None, None
    d = des[-1]
    j = next(i for i, h in enumerate(halves) if h is d)
    s = halves[j - 1] if j > 0 and halves[j - 1]["dir"] == "sod" and halves[j - 1]["complete"] else None
    return s, d


def fig_profiles_by_rate(halves, per_rate, colors, mass_g, sample, out):
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    for name, (_, row) in zip(block_names(per_rate), per_rate.iterrows()):
        s, d = profile_pair(halves, row["block"])
        ls = "--" if "rec." in name else "-"
        col = colors[row["C-rate"]]
        ax.plot(d["Q"] / mass_g, d["V"], ls, color=col, lw=1.4, label=f"{name} - desod. cycle {d['cycle']}")
        if s is not None:
            ax.plot(s["Q"] / mass_g, s["V"], ls, color=col, lw=1.4)
    ax.axhline(PLATEAU_THRESHOLD_V, color="0.6", lw=0.7, ls=":")
    ax.set_xlim(left=0)
    ax.set_ylim(0, DESODIATION_CUTOFF_V + 0.05)
    ax.set_xlabel(CAP_LABEL)
    ax.set_ylabel(POT_LABEL)
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ax.set_title(f"{sample}: sodiation + desodiation at each C-rate", fontsize=10)
    fig.savefig(out / "02_voltage_profiles_by_rate.png")
    plt.close(fig)


def fig_first_cycles(cycles, info, mass_g, sample, out):
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    cmap = plt.get_cmap("plasma", 4)
    first_sod = info["first_sod"]
    if info["start"] == "desodiation":
        ax.plot(first_sod["Q"] / mass_g, first_sod["V"], "--", color=cmap(0), lw=1.4,
                label="Initial sodiation (not counted)")
        shown = cycles[:2]
    else:
        shown = cycles[:3]
    for k, c in enumerate(shown, start=1 if info["start"] == "desodiation" else 0):
        for key in ("sod", "desod"):
            h = c.get(key)
            if h is not None:
                ax.plot(h["Q"] / mass_g, h["V"], color=cmap(k), lw=1.4,
                        label=f"Cycle {c['cycle']}" if key == ("desod" if info["start"] == "desodiation" else "sod") else None)
    c1 = info["cycle1"]
    if info["start"] == "desodiation":
        sod1 = c1["sod"]["spec"] if c1.get("sod") is not None else np.nan
        txt = (f"Initial sodiation (not counted): {first_sod['spec']:.1f} mAh/g\n"
               f"Cycle 1 desodiation: {c1['desod']['spec']:.1f} mAh/g\n"
               f"Cycle 1 sodiation: {sod1:.1f} mAh/g\n"
               f"ICE: {info['ice']:.1f} %")
    else:
        txt = (f"Cycle 1 ({rate_label(c1['rate'])})\n"
               f"Sodiation: {c1['sod']['spec']:.1f} mAh/g\n"
               f"Desodiation: {c1['desod']['spec']:.1f} mAh/g\n"
               f"ICE: {info['ice']:.1f} %")
    ax.text(0.97, 0.60, txt, transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    ax.set_xlim(left=0)
    ax.set_ylim(0, DESODIATION_CUTOFF_V + 0.05)
    ax.set_xlabel(CAP_LABEL)
    ax.set_ylabel(POT_LABEL)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title(f"{sample}: first cycles", fontsize=10)
    fig.savefig(out / "03_first_cycles.png")
    plt.close(fig)


def fig_rate_summary(per_rate, colors, sample, out):
    names = block_names(per_rate)
    x = np.arange(len(per_rate))
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 3.9), gridspec_kw={"width_ratios": [1.5, 1]})

    a.bar(x, per_rate["sod. slope mean (mAh/g)"], 0.6, color="#9ecae1", edgecolor="0.3", lw=0.6,
          label=f"Sodiation slope (> {PLATEAU_THRESHOLD_V} V)")
    a.bar(x, per_rate["sod. plateau mean (mAh/g)"], 0.6, bottom=per_rate["sod. slope mean (mAh/g)"],
          color="#3182bd", edgecolor="0.3", lw=0.6, label=f"Sodiation plateau (< {PLATEAU_THRESHOLD_V} V)")
    a.errorbar(x, per_rate["desodiation mean (mAh/g)"], yerr=per_rate["desodiation std (mAh/g)"],
               fmt="D", color="black", ms=5, capsize=3, label="Desodiation (reversible), mean ± SD")
    top = per_rate[["sodiation mean (mAh/g)", "desodiation mean (mAh/g)"]].max(axis=1)
    for xi, t, ret in zip(x, top, per_rate["retention vs block 1 (%)"]):
        a.text(xi, t * 1.03 + 2, f"{ret:.0f} %", ha="center", va="bottom", fontsize=8)
    a.set_xticks(x, names, rotation=30, ha="right")
    a.set_ylabel(CAP_LABEL)
    a.set_ylim(0, top.max() * 1.35)
    a.legend(fontsize=7.5, loc="upper right")
    a.set_title("(a) Capacity per C-rate (label = retention vs first block)", fontsize=9)
    a.tick_params(top=False)

    seen = set()
    for _, row in per_rate.iterrows():
        rec = row["C-rate"] in seen
        seen.add(row["C-rate"])
        b.errorbar(row["C-rate"], row["desodiation mean (mAh/g)"], yerr=row["desodiation std (mAh/g)"],
                   fmt="s" if rec else "o", ms=6, capsize=3, color=colors[row["C-rate"]],
                   mfc="white" if rec else colors[row["C-rate"]], zorder=3)
    steps = per_rate.drop_duplicates("C-rate", keep="first")
    b.plot(steps["C-rate"], steps["desodiation mean (mAh/g)"], "-", color="0.6", lw=0.8, zorder=1)
    b.set_xscale("log")
    b.set_xticks(sorted(per_rate["C-rate"].unique()))
    b.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    b.set_xlabel("C-rate (C)")
    b.set_ylabel(CAP_LABEL)
    b.set_ylim(0, per_rate["desodiation mean (mAh/g)"].max() * 1.2)
    b.legend(handles=[plt.Line2D([], [], ls="", marker="o", color="0.3", label="Rate steps"),
                      plt.Line2D([], [], ls="", marker="s", mfc="white", mec="0.3", label="Recovery")],
             fontsize=8, loc="lower left")
    b.set_title("(b) Reversible capacity vs C-rate", fontsize=9)
    fig.suptitle(sample, fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "04_rate_summary.png")
    plt.close(fig)


# ------------------------------ Excel export ---------------------------------
def profiles_table(halves, per_rate, mass_g):
    """Sodiation + desodiation profile of the last complete desodiation of every block."""
    cols = {}
    for name, (_, row) in zip(block_names(per_rate), per_rate.iterrows()):
        s, d = profile_pair(halves, row["block"])
        if s is not None:
            cols[f"{name} sod (c{s['cycle']}) Q (mAh/g)"] = pd.Series(s["Q"] / mass_g)
            cols[f"{name} sod (c{s['cycle']}) V (V)"] = pd.Series(s["V"])
        cols[f"{name} desod (c{d['cycle']}) Q (mAh/g)"] = pd.Series(d["Q"] / mass_g)
        cols[f"{name} desod (c{d['cycle']}) V (V)"] = pd.Series(d["V"])
    return pd.DataFrame(cols)


def write_excel(path, per_cycle, per_rate, profiles, settings):
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        per_rate.to_excel(xw, sheet_name="Per C-rate", index=False)
        per_cycle.to_excel(xw, sheet_name="Per cycle", index=False)
        profiles.to_excel(xw, sheet_name="Profiles", index=False)
        pd.DataFrame(list(settings.items()), columns=["Setting", "Value"]).to_excel(
            xw, sheet_name="Settings", index=False)
        for ws in xw.book.worksheets:
            for col in ws.columns:
                width = max(len(str(col[0].value)), 10)
                ws.column_dimensions[col[0].column_letter].width = min(width + 2, 40)
            for row in ws.iter_rows(min_row=2):
                for cell in row:
                    if isinstance(cell.value, float):
                        cell.number_format = "0.000"


# --------------------------------- main --------------------------------------
def run(path, mass_mg=None, rates=None, start=None, out=None, sample=None, verbose=True):
    """Full analysis of one file. Returns a dict with the tables and the saved figure paths.
    Used by main() below and by app.py."""
    log = print if verbose else (lambda *a, **k: None)
    path = Path(path)
    rates = RATE_SEQUENCE if rates is None else rates
    start = start or CYCLE_START
    mass_mg = mass_mg or mass_from_filename(path)
    if not mass_mg:
        raise ValueError("Active mass not found in the file name - give it with --mass (in mg).")
    mass_g = mass_mg / 1000
    sample = sample or SAMPLE_NAME or nice_sample_name(path, mass_mg)
    out = Path(out) if out else path.with_name(path.stem + "_analysis")
    out.mkdir(parents=True, exist_ok=True)

    log(f"File: {path.name}\nActive mass: {mass_mg} mg\nCycle counting starts with: {start}")
    halves = build_half_cycles(load_detail(path), mass_g)
    assign_rates(halves, mass_g, rates, C_RATE_BASIS_MAH_G)
    initial, cycles = build_cycles(halves, start)
    info = formation_info(halves, cycles, start)
    per_cycle = per_cycle_table(cycles, start)
    per_rate = per_rate_table(halves)
    colors = rate_colors(halves)

    incomplete = [f"{NAMES[h['dir']]} (cycle {h['cycle']})" for h in halves if not h["complete"]]
    if incomplete:
        log(f"Not used (cut-off not reached): {', '.join(incomplete)}")

    fig_rate_capability(cycles, colors, info, sample, out)
    fig_profiles_by_rate(halves, per_rate, colors, mass_g, sample, out)
    fig_first_cycles(cycles, info, mass_g, sample, out)
    fig_rate_summary(per_rate, colors, sample, out)

    c1 = info["cycle1"]
    settings = {
        "Source file": path.name, "Active mass (mg)": mass_mg,
        "Cycle counting": ("initial sodiation not counted; cycle = desodiation + next sodiation"
                           if start == "desodiation" else "cycle = sodiation + desodiation"),
        "C-rate sequence": ", ".join(rate_label(r) for r in per_rate["C-rate"]),
        "Initial sodiation (mAh/g)": round(info["first_sod"]["spec"], 2),
        "First desodiation (mAh/g)": round(info["first_desod"]["spec"], 2) if info["first_desod"] else None,
        "ICE (%)": round(info["ice"], 2),
        "Sodiation cut-off (V)": SODIATION_CUTOFF_V, "Desodiation cut-off (V)": DESODIATION_CUTOFF_V,
        "Slope/plateau threshold (V)": PLATEAU_THRESHOLD_V,
        "CE definition": "desodiation / the sodiation just before it",
        "Means per C-rate": "all complete half-cycles at that C-rate; the initial sodiation and the ICE are "
                            "not included in the sodiation and CE means",
        "Note": "Discharge = sodiation, charge = desodiation (half-cell vs Na metal). "
                "Cycles rebuilt from the Detail sheet; the Neware Cycle sheet is not used.",
    }
    excel = out / "rate_test_results.xlsx"
    write_excel(excel, per_cycle, per_rate, profiles_table(halves, per_rate, mass_g), settings)

    log(f"\nInitial sodiation: {info['first_sod']['spec']:.1f} mAh/g"
        f"{' (not counted as a cycle)' if start == 'desodiation' else ''}")
    if c1 is not None and start == "desodiation":
        sod1 = c1["sod"]["spec"] if c1.get("sod") is not None else np.nan
        log(f"Cycle 1: desodiation {c1['desod']['spec']:.1f} mAh/g, sodiation {sod1:.1f} mAh/g")
    log(f"ICE: {info['ice']:.1f} %")
    show = per_rate[["C-rate", "cycles (desodiation)", "desodiation mean (mAh/g)", "desodiation std (mAh/g)",
                     "CE mean (%)", "retention vs block 1 (%)"]]
    log("\n" + show.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    log(f"\nFigures and Excel file saved in: {out.resolve()}")
    return {
        "sample": sample, "mass_mg": mass_mg, "out": out, "excel": excel, "ice": info["ice"],
        "summary": show, "per_rate": per_rate, "per_cycle": per_cycle, "notes": incomplete,
        "figures": [out / f for f in ("01_rate_capability.png", "02_voltage_profiles_by_rate.png",
                                      "03_first_cycles.png", "04_rate_summary.png")],
    }


def main():
    p = argparse.ArgumentParser(description="Rate test analysis for Na-ion half-cells (Neware .xlsx)")
    p.add_argument("file", nargs="?", default=FILE)
    p.add_argument("--mass", type=float, default=ACTIVE_MASS_MG, help="active mass in mg")
    p.add_argument("--rates", type=float, nargs="*", default=RATE_SEQUENCE,
                   help="C-rates in the order they were run, e.g. 0.05 0.1 0.2 0.5 1 2 0.05")
    p.add_argument("--start", choices=["desodiation", "sodiation"], default=CYCLE_START,
                   help="how cycles are counted (see CYCLE_START)")
    p.add_argument("--out", default=OUT_DIR)
    args = p.parse_args()
    try:
        run(args.file, args.mass, args.rates, args.start, args.out)
    except ValueError as e:
        raise SystemExit(str(e))


if __name__ == "__main__":
    main()
