"""
GCPL (galvanostatic cycling) analysis for Na-ion half-cells (hard carbon vs Na metal)
exported from a Neware/BTS cycler as .xlsx (sheets: Info, Cycle_*, Statis_*, Detail_*).

Figures
-------
01_cycling_performance.png  capacity (sodiation/desodiation) and CE vs cycle number
02_CE_stabilisation.png     CE vs cycle with thresholds + irreversible capacity per cycle
03_voltage_profiles.png     potential vs capacity for selected cycles
04_dQdV.png                 dQ/dV of the cycles in DQDV_CYCLES on one graph (sodiation by default)
plus gcpl_results.xlsx with every number (per cycle, CE summary, profiles, dQ/dV) for OriginPro.

Cycle counting (CYCLE_START)
----------------------------
"desodiation": the first sodiation is the initial (formation) sodiation and is not counted.
               Cycle 1 = first desodiation + the sodiation after it, and so on.
"sodiation":   standard counting, cycle 1 = first sodiation + first desodiation.
Coulombic efficiency is always: desodiation / the sodiation just before it
(sodium taken out / sodium put in). For the first desodiation this is the ICE.

Naming used in a half-cell
--------------------------
Discharge ("DChg") = sodiation (voltage goes down), charge ("Chg") = desodiation (voltage goes up).

Usage
-----
    python analyse_gcpl.py                               # uses the settings below
    python analyse_gcpl.py "Number 1-1-4-2818576166.xlsx" --mass 0.66

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

# ============================ USER SETTINGS ==================================
FILE = "Number 1-1-4-2818576166.xlsx"
ACTIVE_MASS_MG = None          # None = read from the file name or the Info sheet (e.g. "Sample_F_0,66mg")
SAMPLE_NAME = None             # None = read from the file name or the Info sheet
CYCLE_START = "desodiation"    # "desodiation" (initial sodiation not counted) or "sodiation" (standard)
SELECTED_CYCLES = [1, 2, 3, 5, 10, 20, "last"]   # cycles shown in the profile and dQ/dV plots
CE_THRESHOLDS = [95, 99]       # CE levels (%) drawn in the CE plot and reported
N_LAST_FOR_STABLE_CE = 5       # "stabilised CE" = mean of the last N complete cycles
C_RATE_BASIS_MAH_G = 250       # only used to report the C-rate (1C = mass x this)
SODIATION_CUTOFF_V = 0.005     # lower cut-off (V vs Na+/Na)
DESODIATION_CUTOFF_V = 2.5     # upper cut-off (V vs Na+/Na)
PLATEAU_THRESHOLD_V = 0.1      # below this potential = plateau, above = slope
DQDV_CYCLES = [2, 15]          # cycles compared in the dQ/dV plot (same graph)
DQDV_BRANCH = "sodiation"      # "sodiation" (= discharge), "desodiation" (= charge) or "both"
DQDV_STEP_V = 0.002            # potential grid for dQ/dV (V)
DQDV_SMOOTH_POINTS = 7         # moving-average window for dQ/dV (odd number, 1 = no smoothing)
DQDV_INSET_RANGE_V = (0.12, 1.3)  # potential window of the magnified inset (slope region)
OUT_DIR = None                 # None = "<file name>_GCPL_analysis" next to the file
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
DQDV_LABEL = r"dQ/dV (mAh g$^{-1}$ V$^{-1}$)"
NAMES = {"sod": "sodiation", "desod": "desodiation"}


# ------------------------------- helpers -------------------------------------
def find_col(df, *keys):
    for c in df.columns:
        name = str(c).lower()
        if all(k.lower() in name for k in keys):
            return c
    raise KeyError(f"No column containing {keys} in {list(df.columns)}")


def current_to_mA(series, col_name):
    name = col_name.lower().replace("\u00b5", "u").replace("\u03bc", "u")
    return series.astype(float) * (1e-3 if "(ua)" in name else 1.0)


MASS_RE = re.compile(r"(\d+)[_,.](\d+)\s*mg", re.IGNORECASE)
SAMPLE_RE = re.compile(r"sample[_ ]?([A-Za-z0-9]+)", re.IGNORECASE)


def mass_and_name(path):
    """Look for '<x>,<yy>mg' and 'Sample_<name>' in the file name, then in the Info sheet."""
    texts = [Path(path).name]
    try:
        info = pd.read_excel(path, sheet_name=0, header=None)
        texts += [str(v) for v in info.to_numpy().ravel() if isinstance(v, str)]
    except Exception:
        pass
    mass, name, mass_text = None, None, None
    for t in texts:
        m = MASS_RE.search(t)
        if m and mass is None:
            mass, mass_text = float(f"{m.group(1)}.{m.group(2)}"), t
    names = {m.group(1) for t in texts for m in [SAMPLE_RE.search(t)] if m}
    if mass_text and SAMPLE_RE.search(mass_text):
        name = SAMPLE_RE.search(mass_text).group(1)     # the name written next to the mass
    elif names:
        name = sorted(names)[0]
    if len(names) > 1:
        print(f"WARNING: This file mentions more than one sample ({', '.join(sorted(names))}). "
              f"Using {name}, which is written next to the mass in '{mass_text}'. Check that this is the right cell.")
    return mass, name


def crossing_capacity(q, v, threshold, direction):
    idx = np.where(v <= threshold)[0] if direction == "down" else np.where(v >= threshold)[0]
    if len(idx) == 0:
        return q[-1]
    i = idx[0]
    if i == 0:
        return 0.0
    v0, v1, q0, q1 = v[i - 1], v[i], q[i - 1], q[i]
    return q1 if v1 == v0 else q0 + (threshold - v0) * (q1 - q0) / (v1 - v0)


def average_potential(q, v):
    if q[-1] <= 0:
        return np.nan
    integrate = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return integrate(v, q) / q[-1]


def moving_average(y, n):
    if n <= 1 or len(y) < n:
        return y
    pad = n // 2
    yp = np.pad(y, pad, mode="edge")
    return np.convolve(yp, np.ones(n) / n, mode="valid")


def dqdv(half, mass_g, step, smooth):
    """dQ/dV on an even potential grid. Sodiation comes out negative, desodiation positive."""
    q, v = half["Q"], half["V"]
    v_mono = np.minimum.accumulate(v) if half["dir"] == "sod" else np.maximum.accumulate(v)
    order = np.argsort(v_mono, kind="stable")
    x, y = v_mono[order], q[order]
    xu, inv = np.unique(x, return_inverse=True)            # average Q at identical potentials
    yu = np.bincount(inv, weights=y) / np.bincount(inv)
    if len(xu) < 5:
        return np.array([]), np.array([])
    grid = np.arange(np.ceil(xu[0] / step) * step, xu[-1], step)
    if len(grid) < 5:
        return np.array([]), np.array([])
    qg = np.interp(grid, xu, yu)
    d = np.gradient(qg, grid) / mass_g
    return grid, moving_average(d, smooth)


# ------------------------------ loading --------------------------------------
def load_detail(path):
    xl = pd.ExcelFile(path)
    sheet = next(s for s in xl.sheet_names if s.lower().startswith("detail"))
    det = pd.read_excel(path, sheet_name=sheet)
    c_cur = find_col(det, "current")
    out = pd.DataFrame({
        "record": det[find_col(det, "record")],
        "status": det[find_col(det, "status")].astype(str),
        "step": det[find_col(det, "steps")].astype(int),
        "V": det[find_col(det, "voltage")].astype(float),
        "I_mA": current_to_mA(det[c_cur], str(c_cur)),
        "Q_mAh": det[next(c for c in det.columns if str(c).lower().startswith("capacity"))].astype(float),
    }).sort_values("record").reset_index(drop=True)
    s = out["status"].str.lower()
    out["direction"] = np.where(s.str.contains("dchg"), "sod",
                                np.where(s.str.contains("chg"), "desod", "rest"))
    return out


def build_half_cycles(det, mass_g):
    halves = []
    for step, g in det.groupby("step", sort=True):
        d = g["direction"].iloc[0]
        if d == "rest":
            continue
        q, v = g["Q_mAh"].to_numpy(dtype=float), g["V"].to_numpy(dtype=float)
        i_med = float(np.median(np.abs(g["I_mA"].to_numpy())))
        if halves and halves[-1]["dir"] == d:                 # merge CC + CV etc.
            h = halves[-1]
            h["Q"] = np.concatenate([h["Q"], q + h["Q"][-1]])
            h["V"] = np.concatenate([h["V"], v])
            h["steps"].append(step)
        else:
            halves.append({"dir": d, "Q": q, "V": v, "I_mA": i_med, "steps": [step]})

    for h in halves:
        q, v = h["Q"], h["V"]
        h["mAh"], h["spec"], h["avgV"] = q[-1], q[-1] / mass_g, average_potential(q, v)
        if h["dir"] == "sod":
            qs = crossing_capacity(q, v, PLATEAU_THRESHOLD_V, "down")
            h["slope"], h["plateau"] = qs / mass_g, (q[-1] - qs) / mass_g
            h["complete"] = bool(v.min() <= SODIATION_CUTOFF_V + 0.003)
        else:
            qp = crossing_capacity(q, v, PLATEAU_THRESHOLD_V, "up")
            h["plateau"], h["slope"] = qp / mass_g, (q[-1] - qp) / mass_g
            h["complete"] = bool(v.max() >= DESODIATION_CUTOFF_V - 0.01)

    first_sod = next((h for h in halves if h["dir"] == "sod"), None)
    for j, h in enumerate(halves):
        h["is_first_sod"] = h is first_sod
        h["CE"] = h["irrev"] = h["hyst"] = np.nan
        h["is_ice"] = False
        if h["dir"] == "desod" and j > 0 and halves[j - 1]["dir"] == "sod":
            prev = halves[j - 1]
            if h["complete"] and prev["complete"]:
                h["CE"] = 100 * h["mAh"] / prev["mAh"]
                h["irrev"] = prev["spec"] - h["spec"]
                h["hyst"] = h["avgV"] - prev["avgV"]
            h["is_ice"] = prev is first_sod
    return halves


def build_cycles(halves, start):
    first, second = ("desod", "sod") if start == "desodiation" else ("sod", "desod")
    i0 = next((i for i, h in enumerate(halves) if h["dir"] == first), len(halves))
    seq = halves[i0:]
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
        for key in ("sod", "desod"):
            if c.get(key) is not None:
                c[key]["cycle"] = n
    for h in halves[:i0]:
        h["cycle"] = 0
    return cycles


# ------------------------------ analysis -------------------------------------
def per_cycle_table(cycles, start, ref_desod):
    order = ["desod", "sod"] if start == "desodiation" else ["sod", "desod"]
    rows, cum = [], 0.0
    for c in cycles:
        r = {"cycle": c["cycle"]}
        for key in order:
            h, nm = c.get(key), NAMES[key]
            if h is None:
                continue
            r[f"{nm} (mAh)"] = h["mAh"]
            r[f"{nm} (mAh/g)"] = h["spec"]
            r[f"{nm} slope >thr (mAh/g)"] = h["slope"]
            r[f"{nm} plateau <thr (mAh/g)"] = h["plateau"]
            r[f"avg. {nm} potential (V)"] = h["avgV"]
            r[f"{nm} complete"] = h["complete"]
            if key == "desod":
                r["CE (%) = desod / previous sod"] = h["CE"]
                r["irreversible capacity (mAh/g)"] = h["irrev"]
                if np.isfinite(h["irrev"]):
                    cum += h["irrev"]
                r["cumulative irreversible (mAh/g)"] = cum
                r["hysteresis (V)"] = h["hyst"]
                r["desodiation retention vs cycle 1 (%)"] = (100 * h["spec"] / ref_desod
                                                             if h["complete"] else np.nan)
        r["Neware steps"] = ", ".join(str(s) for key in order if c.get(key) is not None
                                      for s in c[key]["steps"])
        rows.append(r)
    return pd.DataFrame(rows)


def ce_summary(halves, cycles):
    des = [h for h in halves if h["dir"] == "desod" and np.isfinite(h["CE"])]
    ce = {h["cycle"]: h["CE"] for h in des}
    ice = next((h["CE"] for h in des if h["is_ice"]), np.nan)
    counted = [n for n in sorted(ce) if n >= 1 and not any(h["is_ice"] and h["cycle"] == n for h in des)]
    last_n = counted[-N_LAST_FOR_STABLE_CE:]
    rows = [("ICE (%)", ice)]
    for n in (2, 3, 5, 10, 20):
        if n in ce:
            rows.append((f"CE cycle {n} (%)", ce[n]))
    if counted:
        rows.append((f"CE last complete cycle ({counted[-1]}) (%)", ce[counted[-1]]))
        vals = np.array([ce[n] for n in last_n])
        rows.append((f"Stabilised CE: mean of cycles {last_n[0]}-{last_n[-1]} (%)", vals.mean()))
        rows.append((f"Stabilised CE: SD of cycles {last_n[0]}-{last_n[-1]} (%)",
                     vals.std(ddof=1) if len(vals) > 1 else 0.0))
    for t in CE_THRESHOLDS:
        reached = next((n for n in counted if ce[n] >= t), None)
        stable = next((n for k, n in enumerate(counted) if all(ce[m] >= t for m in counted[k:])), None)
        rows.append((f"First cycle with CE >= {t} %", reached if reached else "not reached"))
        rows.append((f"CE >= {t} % from cycle on (and stays)", stable if stable else "not reached"))
    des_ok = [c["desod"] for c in cycles if c.get("desod") is not None and c["desod"]["complete"]]
    if des_ok:
        rows.append(("Desodiation cycle 1 (mAh/g)", des_ok[0]["spec"]))
        rows.append((f"Desodiation last complete cycle ({des_ok[-1]['cycle']}) (mAh/g)", des_ok[-1]["spec"]))
        rows.append(("Desodiation retention last vs cycle 1 (%)", 100 * des_ok[-1]["spec"] / des_ok[0]["spec"]))
    return pd.DataFrame(rows, columns=["Quantity", "Value"]), ce, counted


def resolve_selected(cycles):
    complete = [c["cycle"] for c in cycles
                if all(c.get(k) is not None and c[k]["complete"] for k in ("sod", "desod"))]
    if not complete:
        return []
    chosen = []
    for s in SELECTED_CYCLES:
        n = complete[-1] if s == "last" else int(s)
        if n in complete and n not in chosen:
            chosen.append(n)
    return chosen


# ------------------------------- figures -------------------------------------
def fig_cycling(cycles, info, sample, out):
    counted = [c[k] for c in cycles for k in ("sod", "desod") if c.get(k) is not None and c[k]["complete"]]
    fig, (ax_ce, ax) = plt.subplots(2, 1, figsize=(7.0, 5.0), sharex=True,
                                    gridspec_kw={"height_ratios": [1, 3], "hspace": 0.06})
    normal = [h["spec"] for h in counted if not h["is_first_sod"]]
    y_top = 1.3 * max(normal)
    for h in counted:
        y = min(h["spec"], y_top * 0.97)
        if h["dir"] == "desod":
            ax.plot(h["cycle"], y, "o", color="#2c7fb8", ms=5.5, zorder=3)
        else:
            ax.plot(h["cycle"], y, "o", mfc="white", mec="#d95f0e", ms=5.5, zorder=3)
        if h["is_first_sod"] and h["spec"] > y_top * 0.97:
            ax.annotate(f"1st sodiation: {h['spec']:.0f} mAh/g (off scale)", xy=(h["cycle"], y),
                        xytext=(h["cycle"] + 1.5, y_top * 0.99), fontsize=8, va="top",
                        arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.set_ylim(0, y_top)
    ax.set_xlim(0.3, cycles[-1]["cycle"] + 0.7)
    ax.set_xlabel("Cycle number")
    ax.set_ylabel(CAP_LABEL)
    ax.legend(handles=[plt.Line2D([], [], ls="", marker="o", color="#2c7fb8", label="Desodiation (reversible)"),
                       plt.Line2D([], [], ls="", marker="o", mfc="white", mec="#d95f0e", label="Sodiation")],
              loc="lower left", fontsize=8)
    ce = [(h["cycle"], h["CE"]) for h in counted if h["dir"] == "desod" and not h["is_ice"] and np.isfinite(h["CE"])]
    if ce:
        x, y = zip(*ce)
        ax_ce.plot(x, y, "D-", color="0.3", ms=3.5, lw=0.6, zorder=3)
        ax_ce.set_ylim(np.floor(min(y) / 5) * 5 - 2, 101)
    ax_ce.set_ylabel("CE (%)")
    ax_ce.text(0.988, 0.08, f"ICE = {info['ice']:.1f} %", transform=ax_ce.transAxes, fontsize=8,
               va="bottom", ha="right")
    ax_ce.set_title(sample, fontsize=10)
    fig.savefig(out / "01_cycling_performance.png")
    plt.close(fig)


def fig_ce_stabilisation(halves, ce, counted, info, sample, out):
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.5, 3.9))
    if counted:
        x = np.array(counted)
        y = np.array([ce[n] for n in counted])
        a.plot(x, y, "D-", color="0.25", ms=4, lw=0.8, zorder=3, label="CE")
        last = y[-N_LAST_FOR_STABLE_CE:]
        a.axhspan(last.mean() - last.std(ddof=1) if len(last) > 1 else last.mean(),
                  last.mean() + last.std(ddof=1) if len(last) > 1 else last.mean(),
                  color="#9ecae1", alpha=0.5, zorder=1,
                  label=f"Stabilised CE {last.mean():.1f} ± {last.std(ddof=1) if len(last) > 1 else 0:.1f} % "
                        f"(last {len(last)} cycles)")
        for t, ls in zip(CE_THRESHOLDS, ["--", ":"]):
            a.axhline(t, color="#d95f0e", lw=0.9, ls=ls, label=f"{t} %")
        a.set_ylim(min(np.floor(y.min() / 5) * 5 - 2, min(CE_THRESHOLDS) - 2), 101)
        a.set_xlim(0.3, x.max() + 0.7)
    a.set_xlabel("Cycle number")
    a.set_ylabel("Coulombic efficiency (%)")
    a.legend(fontsize=7.5, loc="lower right")
    a.set_title(f"(a) Coulombic efficiency (ICE = {info['ice']:.1f} %, not shown)", fontsize=9)

    des = [h for h in halves if h["dir"] == "desod" and np.isfinite(h["irrev"]) and h["cycle"] >= 1]
    if des:
        xi = np.array([h["cycle"] for h in des])
        yi = np.array([h["irrev"] for h in des])
        b.plot(xi, np.clip(yi, 1e-3, None), "o-", color="#d95f0e", ms=4, lw=0.8, label="per cycle")
        b.set_yscale("log")
        b.set_ylabel(r"Irreversible capacity per cycle (mAh g$^{-1}$)", color="#d95f0e")
        b.tick_params(axis="y", colors="#d95f0e", right=False)
        b2 = b.twinx()
        b2.plot(xi, np.cumsum(yi), "-", color="0.3", lw=1.4)
        b2.set_ylabel(r"Cumulative irreversible capacity (mAh g$^{-1}$)", color="0.3")
        b2.set_ylim(0, np.cumsum(yi)[-1] * 1.1)
        b.set_xlim(0.3, xi.max() + 0.7)
    b.set_xlabel("Cycle number")
    b.set_title("(b) Capacity lost per cycle (previous sodiation − desodiation)", fontsize=9)
    b.tick_params(top=False)
    fig.suptitle(sample, fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "02_CE_stabilisation.png")
    plt.close(fig)


def selected_curves(cycles, info, selected):
    curves = []
    if info["start"] == "desodiation" and info["first_sod"] is not None:
        curves.append(("Initial sodiation (not counted)", [info["first_sod"]], True))
    for n in selected:
        c = next(cc for cc in cycles if cc["cycle"] == n)
        curves.append((f"Cycle {n}", [c[k] for k in ("sod", "desod") if c.get(k) is not None], False))
    return curves


def fig_profiles(curves, mass_g, sample, out):
    """Charge/discharge curves of the selected cycles. If the initial sodiation is much longer than the
    later cycles, it gets its own panel so the counted cycles are not squeezed."""
    sel = [(lab, hs) for lab, hs, init in curves if not init]
    init = [(lab, hs) for lab, hs, ini in curves if ini]
    max_sel = max((h["Q"][-1] / mass_g for _, hs in sel for h in hs), default=1.0)
    max_init = max((h["Q"][-1] / mass_g for _, hs in init for h in hs), default=0.0)
    split = bool(init) and max_init > 1.5 * max_sel
    cmap = plt.get_cmap("viridis", max(len(sel), 2))
    col = {lab: cmap(k / max(len(sel) - 1, 1) * 0.92) for k, (lab, _) in enumerate(sel)}

    def draw(ax, items, dashed_init):
        for lab, hs, is_init in items:
            for j, h in enumerate(hs):
                ax.plot(h["Q"] / mass_g, h["V"], "--" if is_init else "-",
                        color="0.35" if is_init else col[lab], lw=1.3, label=lab if j == 0 else None)
        ax.axhline(PLATEAU_THRESHOLD_V, color="0.6", lw=0.7, ls=":")
        ax.set_xlim(left=0)
        ax.set_ylim(0, DESODIATION_CUTOFF_V + 0.05)
        ax.set_xlabel(CAP_LABEL)
        ax.set_ylabel(POT_LABEL)

    if split:
        fig, (a, b) = plt.subplots(1, 2, figsize=(11.0, 4.4), gridspec_kw={"width_ratios": [1, 1.25]})
        first = [(init[0][0], init[0][1], True)] + ([(sel[0][0], sel[0][1], False)] if sel else [])
        draw(a, first, True)
        a.legend(fontsize=8, loc="upper right")
        a.set_title(f"(a) Initial sodiation ({max_init:.0f} mAh/g) and cycle 1", fontsize=9)
        draw(b, [(lab, hs, False) for lab, hs in sel], False)
        b.set_xlim(0, max_sel * 1.08)
        b.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1.0))
        b.set_title("(b) Charge/discharge curves of the counted cycles", fontsize=9)
        fig.suptitle(f"{sample}: galvanostatic profiles", fontsize=10)
        fig.tight_layout()
    else:
        fig, ax = plt.subplots(figsize=(6.2, 4.4))
        draw(ax, curves, True)
        ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1.0))
        ax.set_title(f"{sample}: galvanostatic profiles", fontsize=10)
    fig.savefig(out / "03_voltage_profiles.png")
    plt.close(fig)


def fig_dqdv(cycles, mass_g, sample, out):
    """dQ/dV of the cycles in DQDV_CYCLES on one graph, for the branch in DQDV_BRANCH."""
    branches = {"sodiation": ["sod"], "desodiation": ["desod"], "both": ["sod", "desod"]}[DQDV_BRANCH]
    colors = ["#2c7fb8", "#d95f0e", "#1b9e77", "#7570b3", "#e7298a", "#66a61e"]
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ins = ax.inset_axes([0.42, 0.10, 0.54, 0.45])
    table, peaks, ins_vals = {}, [], []
    for k, n in enumerate(DQDV_CYCLES):
        c = next((cc for cc in cycles if cc["cycle"] == n), None)
        if c is None:
            print(f"WARNING: cycle {n} does not exist - skipped in dQ/dV")
            continue
        for key in branches:
            h = c.get(key)
            if h is None or not h["complete"]:
                print(f"WARNING: {NAMES[key]} of cycle {n} is missing or incomplete - skipped in dQ/dV")
                continue
            x, y = dqdv(h, mass_g, DQDV_STEP_V, DQDV_SMOOTH_POINTS)
            if not len(x):
                continue
            ls = "-" if key == "sod" else "--"
            lab = f"Cycle {n}" + ("" if len(branches) == 1 else f" {NAMES[key]}")
            ax.plot(x, y, ls, color=colors[k % len(colors)], lw=1.4, label=lab)
            m = (x >= DQDV_INSET_RANGE_V[0]) & (x <= DQDV_INSET_RANGE_V[1])
            ins.plot(x[m], y[m], ls, color=colors[k % len(colors)], lw=1.2)
            ins_vals.append(y[m])
            table[f"{lab} V (V)"] = pd.Series(x)
            table[f"{lab} dQ/dV (mAh/g/V)"] = pd.Series(y)
            low = x < 0.15                                  # main plateau peak
            if low.any():
                i = np.argmax(np.abs(y[low]))
                peaks.append((f"Cycle {n} {NAMES[key]}", float(x[low][i]), float(y[low][i])))

    ax.axhline(0, color="0.5", lw=0.6)
    ax.set_xlim(0, 1.5 if branches == ["sod"] else DESODIATION_CUTOFF_V)
    ax.set_xlabel(POT_LABEL)
    ax.set_ylabel(DQDV_LABEL)
    ax.legend(fontsize=9, loc="center left", bbox_to_anchor=(0.14, 0.45))
    ins.axhline(0, color="0.5", lw=0.5)
    ins.set_xlim(*DQDV_INSET_RANGE_V)
    if ins_vals:
        v = np.concatenate(ins_vals)
        lo, hi = np.nanmin(v), np.nanmax(v)
        pad = 0.1 * (hi - lo if hi > lo else 1)
        ins.set_ylim(min(lo - pad, 0), max(hi + pad, 0))
    ins.tick_params(labelsize=7, top=False, right=False)
    ins.set_title(f"Slope region {DQDV_INSET_RANGE_V[0]:g}-{DQDV_INSET_RANGE_V[1]:g} V (magnified)", fontsize=7.5)
    ax.indicate_inset_zoom(ins, edgecolor="0.6")
    branch_txt = {"sodiation": "sodiation (discharge)", "desodiation": "desodiation (charge)",
                  "both": "sodiation and desodiation"}[DQDV_BRANCH]
    ax.set_title(f"{sample}: dQ/dV, {branch_txt}", fontsize=10)
    fig.savefig(out / "04_dQdV.png")
    plt.close(fig)
    for name, vpk, ypk in peaks:
        print(f"dQ/dV main low-voltage peak, {name}: {vpk * 1000:.0f} mV ({ypk:.0f} mAh/g/V)")
    peak_df = pd.DataFrame(peaks, columns=["curve", "peak potential (V)", "peak dQ/dV (mAh/g/V)"])
    return pd.DataFrame(table), peak_df


def profiles_table(curves, mass_g):
    cols = {}
    for label, hs, _ in curves:
        for h in hs:
            tag = f"{label} {NAMES[h['dir']]}"
            cols[f"{tag} Q (mAh/g)"] = pd.Series(h["Q"] / mass_g)
            cols[f"{tag} V (V)"] = pd.Series(h["V"])
    return pd.DataFrame(cols)


def write_excel(path, sheets):
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for name, df in sheets.items():
            df.to_excel(xw, sheet_name=name, index=False)
        for ws in xw.book.worksheets:
            for col in ws.columns:
                width = max(len(str(col[0].value)), 10)
                ws.column_dimensions[col[0].column_letter].width = min(width + 2, 45)
            for row in ws.iter_rows(min_row=2):
                for cell in row:
                    if isinstance(cell.value, float):
                        cell.number_format = "0.000"


# --------------------------------- main --------------------------------------
def run(path, mass_mg=None, start=None, out=None, sample=None, verbose=True):
    """Full analysis of one file. Returns a dict with the tables and the saved figure paths.
    Used by main() below and by app.py."""
    log = print if verbose else (lambda *a, **k: None)
    path = Path(path)
    start = start or CYCLE_START
    auto_mass, auto_name = mass_and_name(path)
    mass_mg = mass_mg or auto_mass
    if not mass_mg:
        raise ValueError("Active mass not found in the file name or Info sheet - give it with --mass (in mg).")
    mass_g = mass_mg / 1000
    sample = sample or SAMPLE_NAME or (f"Sample {auto_name}, {mass_mg:g} mg" if auto_name
                                       else f"{path.stem}, {mass_mg:g} mg")
    out = Path(out) if out else path.with_name(path.stem + "_GCPL_analysis")
    out.mkdir(parents=True, exist_ok=True)

    log(f"File: {path.name}\nSample: {sample}\nCycle counting starts with: {start}")
    halves = build_half_cycles(load_detail(path), mass_g)
    cycles = build_cycles(halves, start)
    i_med = np.median([h["I_mA"] for h in halves])
    c_rate = i_med / (mass_g * C_RATE_BASIS_MAH_G)
    log(f"Current: {i_med * 1000:.2f} uA  ->  C/{1 / c_rate:.1f} "
        f"(1C = {C_RATE_BASIS_MAH_G} mAh/g). If this is not the C-rate you set, check the mass!")

    first_sod = next((h for h in halves if h["is_first_sod"]), None)
    ice_des = next((h for h in halves if h["is_ice"]), None)
    info = {"start": start, "first_sod": first_sod, "ice": ice_des["CE"] if ice_des else np.nan}

    des_ok = [c["desod"] for c in cycles if c.get("desod") is not None and c["desod"]["complete"]]
    per_cycle = per_cycle_table(cycles, start, des_ok[0]["spec"] if des_ok else np.nan)
    summary, ce, counted = ce_summary(halves, cycles)
    selected = resolve_selected(cycles)
    curves = selected_curves(cycles, info, selected)

    incomplete = [f"{NAMES[h['dir']]} (cycle {h['cycle']})" for h in halves if not h["complete"]]
    if incomplete:
        log(f"Not used (cut-off not reached): {', '.join(incomplete)}")

    fig_cycling(cycles, info, sample, out)
    fig_ce_stabilisation(halves, ce, counted, info, sample, out)
    fig_profiles(curves, mass_g, sample, out)
    dq_table, dq_peaks = fig_dqdv(cycles, mass_g, sample, out)

    settings = pd.DataFrame([
        ("Source file", path.name), ("Sample", sample), ("Active mass (mg)", mass_mg),
        ("Current (uA)", round(i_med * 1000, 3)),
        ("Cycle counting", "initial sodiation not counted; cycle = desodiation + next sodiation"
         if start == "desodiation" else "cycle = sodiation + desodiation"),
        ("Initial sodiation (mAh/g)", round(first_sod["spec"], 2) if first_sod else None),
        ("CE definition", "desodiation / the sodiation just before it"),
        ("Irreversible capacity", "previous sodiation - desodiation (mAh/g)"),
        ("Cut-offs (V)", f"{SODIATION_CUTOFF_V} - {DESODIATION_CUTOFF_V}"),
        ("Slope/plateau threshold (V)", PLATEAU_THRESHOLD_V),
        ("dQ/dV", f"cycles {', '.join(str(n) for n in DQDV_CYCLES)}, {DQDV_BRANCH}; grid {DQDV_STEP_V} V, "
                  f"moving average {DQDV_SMOOTH_POINTS} pts; sodiation negative, desodiation positive"),
        ("Selected cycles", ", ".join(str(n) for n in selected)),
    ], columns=["Setting", "Value"])
    excel = out / "gcpl_results.xlsx"
    write_excel(excel, {
        "CE summary": summary, "Per cycle": per_cycle, "Profiles": profiles_table(curves, mass_g),
        "dQdV": dq_table, "dQdV peaks": dq_peaks, "Settings": settings})

    log()
    for q, v in summary.itertuples(index=False):
        log(f"  {q:<48} {v:.2f}" if isinstance(v, (float, np.floating)) else f"  {q:<48} {v}")
    log(f"\nFigures and Excel file saved in: {out.resolve()}")
    return {
        "sample": sample, "mass_mg": mass_mg, "out": out, "excel": excel, "ice": info["ice"],
        "current_uA": i_med * 1000, "c_rate": c_rate, "summary": summary, "per_cycle": per_cycle,
        "dqdv_peaks": dq_peaks, "notes": incomplete,
        "figures": [out / f for f in ("01_cycling_performance.png", "02_CE_stabilisation.png",
                                      "03_voltage_profiles.png", "04_dQdV.png")],
    }


def main():
    p = argparse.ArgumentParser(description="GCPL analysis for Na-ion half-cells (Neware .xlsx)")
    p.add_argument("file", nargs="?", default=FILE)
    p.add_argument("--mass", type=float, default=ACTIVE_MASS_MG, help="active mass in mg")
    p.add_argument("--start", choices=["desodiation", "sodiation"], default=CYCLE_START)
    p.add_argument("--out", default=OUT_DIR)
    args = p.parse_args()
    try:
        run(args.file, args.mass, args.start, args.out)
    except ValueError as e:
        raise SystemExit(str(e))


if __name__ == "__main__":
    main()
