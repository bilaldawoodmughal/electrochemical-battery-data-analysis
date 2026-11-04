"""
GITT analysis for BioLogic EC-Lab files (.mpr or .mpt text export) - Na half-cells
================================================================================
Bilal Dawood, TU Darmstadt (FAME-AIS thesis)

What it does
  1. Finds every current pulse and the rest after it (sodiation = negative current,
     desodiation = positive current).
  2. Checks the data: rest drift at the end of every rest, IR drop and sudden jumps in it,
     pulses cut short by the voltage limit, the first pulse after an unrelaxed rest.
  3. Calculates the Na diffusion coefficient per pulse with the Weppner-Huggins equation
         D = 4/(pi*tau) * L^2 * (dEs/dEt)^2          valid for tau << L^2/D
     L is the diffusion length (see LENGTH_MODE below).
  4. Saves figures (PNG, 300 dpi) and an Excel table.

Run on its own
  pip install numpy pandas matplotlib openpyxl   (galvani is optional: galvani_biologic.py is used)
  python analyse_gitt.py                (uses the settings below)
  python analyse_gitt.py myfile.mpr     (other file, same settings)
The same code runs inside the Neware by Bilal app.
"""
import struct
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ============================== SETTINGS =====================================
FILE = "example_data/Sample_C1200_GITT_C-20_2.414mg.mpr"
SAMPLE_NAME = None             # None = electrode name stored in the .mpr, or the file name
MASS_MG = None                 # None = active mass stored in the .mpr settings (or in the file name)
DIAMETER_MM = 10.0             # electrode diameter (mm)
THICKNESS_MM = 0.75            # electrode thickness (mm)

# Diffusion length L used in the Weppner-Huggins equation:
#   "thickness" : L = electrode thickness (THICKNESS_MM)
#   "mass"      : L = m * V_M / (M * S) with S = disc area  (= m / (density * area),
#                 the thickness the carbon would have with no pores; Ledwoch et al. 2022 form)
#   "fibre"     : L = V/S of one fibre    = r/2  (needs SEM_DIAMETER_UM)
#   "particle"  : L = V/S of one particle = r/3  (spheres, needs SEM_DIAMETER_UM)
LENGTH_MODE = "mass"
SEM_DIAMETER_UM = None         # mean fibre or particle diameter from SEM (um)
V_M_CM3_MOL = 48.0             # molar volume of C6 (Ledwoch et al. 2022: 72 g/mol / 1.5 g/cm3)
M_G_MOL = 72.0                 # molar mass of C6

C_RATE_BASIS_MAH_G = 250       # 1C definition (only used to report the C-rate)
IR_TIME_S = 1.0                # the voltage change in the first second of a pulse is
                               # counted as the IR jump; dEt is measured from this point
MIN_PULSE_FRACTION = 0.95      # pulses shorter than this x nominal length (stopped by the
                               # voltage limit) are not used for D
SKIP_FIRST_PULSE = True        # first pulse comes after an unrelaxed rest -> not used
EXCLUDE_PULSES = {"sodiation": [], "desodiation": []}   # pulse numbers to leave out by hand
REST_DRIFT_LIMIT_MV_H = 3.6    # Ledwoch criterion: dE/dt <= 0.001 mV/s at the end of a rest
DRIFT_WINDOW_S = 600           # last part of each rest used to measure the drift
APPARENT_LIMIT = 0.3           # tau*D/L^2 above this = assumption tau << L^2/D not met
RELAX_AVG_S = 30.0             # relaxed potential = median of the last 30 s of each rest
STEP_MIN_MV = 10.0             # a sudden potential step larger than this (and STEP_RATIO x the normal change)
STEP_RATIO = 8.0               # inside a pulse or the rest after it = contact / counter-electrode disturbance;
                               # that pulse is not used
OUT_DIR = None                 # None = "<file name>_GITT_results" next to the file
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
COL = {"sodiation": "#2a78d6", "desodiation": "#eb6834"}
MRK = {"sodiation": "o", "desodiation": "s"}
POT_LABEL = r"Potential (V vs Na$^{+}$/Na)"
CAP_LABEL = r"Specific capacity (mAh g$^{-1}$)"
D_LABEL = r"$D_{Na^+}$ (cm$^{2}$ s$^{-1}$)"
DIRECTIONS = ("sodiation", "desodiation")


def _biologic():
    """BioLogic file reader: the copy bundled with the app, or an installed galvani."""
    try:
        import galvani_biologic as BioLogic
    except ImportError:
        from galvani import BioLogic
    return BioLogic


# ----------------------------------------------------------------------------- reading
def read_cell_info(path):
    """Active mass (mg) and electrode name stored in the .mpr settings, else from the file name."""
    import re
    path = Path(path)
    mass = name = None
    source = None
    if path.suffix.lower() == ".mpr":
        try:
            BioLogic = _biologic()
            with open(path, "rb") as f:
                f.read(len(BioLogic.MPR_MAGIC))
                for mod in BioLogic.read_VMP_modules(f):
                    if mod["shortname"] == b"VMP Set   ":
                        data = mod["data"]
                        m = struct.unpack_from("<f", data, 0x107)[0]
                        if 0.001 < m < 10000:
                            mass, source = round(float(m), 4), "the EC-Lab settings in the file"
                        n = data[0x11E]
                        txt = data[0x11F:0x11F + n].decode("latin-1", "ignore").strip()
                        if txt and txt.isprintable():
                            name = txt
                    break
        except Exception:
            pass
    if mass is None:
        mm = re.search(r"(\d+)[_,.](\d+)\s*mg", path.name, re.I)
        if mm:
            mass, source = float(f"{mm.group(1)}.{mm.group(2)}"), "the file name"
    if name is None:
        mm = re.search(r"sample[_ ]?([A-Za-z0-9]+)", path.stem, re.I)
        name = mm.group(1) if mm else None
    return mass, name, source


def load(path):
    BioLogic = _biologic()
    path = Path(path)
    if path.suffix.lower() == ".mpr":
        df = pd.DataFrame(BioLogic.MPRfile(str(path)).data)
    else:
        df = pd.DataFrame(BioLogic.MPTfile(str(path), encoding="latin-1")[0])
    cols = set(df.columns)
    t = df["time/s"].to_numpy(float)
    E = df["Ewe/V"].to_numpy(float)
    if "control/V/mA" in cols:          # GCPL: the set current, clean zero during rests
        I = df["control/V/mA"].to_numpy(float) * 1000
    else:                               # measured current: treat noise around zero as rest
        I = df["I/mA"].to_numpy(float) * 1000
        I[np.abs(I) < 0.02 * np.nanmax(np.abs(I))] = 0.0
    if "dq/mA.h" in cols:
        dq = df["dq/mA.h"].to_numpy(float) * 1000
    else:
        dq = np.r_[0, np.diff(t)] * I / 3600
    ns = df["Ns"].to_numpy(int) if "Ns" in cols else np.zeros(len(t), int)
    return pd.DataFrame({"t": t, "E": E, "I_uA": np.round(I, 3), "dq_uAh": dq, "Ns": ns})


def segments(d):
    """Split into constant-current blocks (pulse) and zero-current blocks (rest)."""
    sign = np.sign(d["I_uA"].to_numpy())
    change = np.r_[False, (np.diff(sign) != 0) | (np.diff(d["Ns"].to_numpy()) != 0)]
    d = d.assign(seg=np.cumsum(change))
    out = []
    for s, g in d.groupby("seg"):
        if len(g) < 2:
            continue
        out.append(dict(seg=s, kind="rest" if sign[g.index[0]] == 0 else "pulse",
                        I_uA=g["I_uA"].median(), t0=g["t"].iloc[0], t1=g["t"].iloc[-1],
                        idx0=g.index[0], idx1=g.index[-1], Q_uAh=g["dq_uAh"].sum()))
    s = pd.DataFrame(out)
    # merge neighbouring blocks of the same kind (a 1-point split at a sequence change)
    keep = []
    for _, row in s.iterrows():
        if keep and keep[-1]["kind"] == row["kind"] and np.sign(keep[-1]["I_uA"]) == np.sign(row["I_uA"]):
            keep[-1].update(t1=row["t1"], idx1=row["idx1"], Q_uAh=keep[-1]["Q_uAh"] + row["Q_uAh"])
        else:
            keep.append(row.to_dict())
    return d, pd.DataFrame(keep)


def count_pulses(segs):
    """Number of pulses that have a rest before and after them, per direction."""
    n = {"sodiation": 0, "desodiation": 0}
    kinds = segs["kind"].tolist()
    for k in range(1, len(segs) - 1):
        if kinds[k] == "pulse" and kinds[k - 1] == "rest" and kinds[k + 1] == "rest":
            n["sodiation" if segs.iloc[k]["I_uA"] < 0 else "desodiation"] += 1
    return n


# ----------------------------------------------------------------------------- analysis
def disc_area_cm2():
    return np.pi * (DIAMETER_MM / 20) ** 2


def lengths_cm():
    """All diffusion lengths that can be calculated from the given electrode data (cm)."""
    L = {"thickness": THICKNESS_MM / 10 if THICKNESS_MM else None,
         "mass": (MASS_MG / 1000) * V_M_CM3_MOL / (M_G_MOL * disc_area_cm2()) if DIAMETER_MM else None,
         "fibre": SEM_DIAMETER_UM / 4 * 1e-4 if SEM_DIAMETER_UM else None,
         "particle": SEM_DIAMETER_UM / 6 * 1e-4 if SEM_DIAMETER_UM else None}
    return {k: v for k, v in L.items() if v}


def diffusion_length_cm():
    L = lengths_cm()
    if LENGTH_MODE not in ("thickness", "mass", "fibre", "particle"):
        raise ValueError("LENGTH_MODE must be 'thickness', 'mass', 'fibre' or 'particle'")
    if LENGTH_MODE not in L:
        need = {"thickness": "the electrode thickness", "mass": "the electrode diameter",
                "fibre": "the fibre diameter from SEM", "particle": "the particle diameter from SEM"}
        raise ValueError(f"Enter {need[LENGTH_MODE]} to use this diffusion length.")
    v = L[LENGTH_MODE]
    text = {"thickness": f"L = electrode thickness = {v * 1e4:.0f} µm",
            "mass": f"L = m V$_M$/(M S) = {v * 1e4:.1f} µm",
            "fibre": f"L = r/2 of a fibre = {v * 1e4:.3g} µm",
            "particle": f"L = r/3 of a particle = {v * 1e4:.3g} µm"}[LENGTH_MODE]
    return v, text


def step_mV(t, e, t_from, win):
    """Largest potential change within `win` seconds after `t_from`, and how many times larger it is than
    the typical change in the same block. A clean pulse or rest changes smoothly; a contact problem shows
    up as a sudden step."""
    m = t > t_from
    t, e = t[m], e[m]
    if len(t) < 10:
        return 0.0, 0.0
    idx = np.clip(np.searchsorted(t, t + win), 0, len(t) - 1)
    ok = idx > np.arange(len(t))
    ch = np.abs(e[idx] - e)[ok] * 1000
    if len(ch) == 0:
        return 0.0, 0.0
    return float(ch.max()), float(ch.max() / max(np.median(ch), 0.3))


def disturbed(t, e, t_from, win):
    big, ratio = step_mV(t, e, t_from, win)
    return big > STEP_MIN_MV and ratio > STEP_RATIO


def relaxed_E(d, r):
    """Relaxed potential at the end of a rest: median of the last RELAX_AVG_S seconds (less noisy than the
    single last point)."""
    g = d.loc[r["idx0"]:r["idx1"]]
    tail = g.loc[g["t"] >= r["t1"] - RELAX_AVG_S, "E"]
    return float(tail.median()) if len(tail) else float(g["E"].iloc[-1])


def rest_drift_mV_h(d, r):
    m = (d["t"] >= r["t1"] - DRIFT_WINDOW_S) & (d.index >= r["idx0"]) & (d.index <= r["idx1"])
    g = d[m]
    if len(g) < 5:
        return np.nan
    return np.polyfit(g["t"], g["E"], 1)[0] * 1000 * 3600


def analyse(d, segs, L_cm):
    mass_g = MASS_MG / 1000
    rows = []
    counter = {"sodiation": 0, "desodiation": 0}
    cum = {"sodiation": 0.0, "desodiation": 0.0}
    for k in range(len(segs)):
        p = segs.iloc[k]
        if p["kind"] != "pulse":
            continue
        direction = "sodiation" if p["I_uA"] < 0 else "desodiation"
        before = segs.iloc[k - 1] if k > 0 and segs.iloc[k - 1]["kind"] == "rest" else None
        after = segs.iloc[k + 1] if k + 1 < len(segs) and segs.iloc[k + 1]["kind"] == "rest" else None
        if before is None or after is None:
            continue
        counter[direction] += 1
        g = d.loc[p["idx0"]:p["idx1"]]
        tt = g["t"].to_numpy() - p["t0"]
        ee = g["E"].to_numpy()
        tau = tt[-1]
        q = abs(p["Q_uAh"]) / 1000 / mass_g  # mAh/g in this pulse
        cap_mid = cum[direction] + q / 2
        cum[direction] += q
        Es0 = relaxed_E(d, before)
        Es1 = relaxed_E(d, after)
        E_ir = np.interp(IR_TIME_S, tt, ee) if tau > IR_TIME_S else ee[-1]
        E_end = ee[-1]
        dEs = Es1 - Es0
        dEt = E_end - E_ir
        ratio = dEs / dEt if dEt != 0 else np.nan
        sel = tt >= IR_TIME_S
        r2 = np.nan
        if sel.sum() > 10:
            x = np.sqrt(tt[sel]); y = ee[sel]
            cf = np.polyfit(x, y, 1)
            res = y - np.polyval(cf, x)
            r2 = 1 - (res ** 2).sum() / ((y - y.mean()) ** 2).sum()
        rows.append(dict(
            direction=direction, pulse=counter[direction], t_start_h=p["t0"] / 3600,
            I_uA=p["I_uA"], tau_s=tau, Q_pulse_mAh_g=q, capacity_mid_mAh_g=cap_mid,
            capacity_end_mAh_g=cum[direction], E_relaxed_before_V=Es0, E_relaxed_after_V=Es1,
            E_mid_V=(Es0 + Es1) / 2, E_pulse_end_V=E_end, IR_drop_mV=(E_ir - Es0) * 1000,
            dEs_mV=dEs * 1000, dEt_mV=dEt * 1000, dEs_over_dEt=ratio, sqrt_t_R2=r2,
            rest_h=(after["t1"] - after["t0"]) / 3600, rest_end_drift_mV_h=rest_drift_mV_h(d, after),
            disturbed=disturbed(tt, ee, 60, 3) or disturbed(
                d.loc[after["idx0"]:after["idx1"], "t"].to_numpy() - after["t0"],
                d.loc[after["idx0"]:after["idx1"], "E"].to_numpy(), 600, 30),
        ))
    res = pd.DataFrame(rows)
    if res.empty:
        raise ValueError("No current pulses with a rest before and after them were found.")
    nominal = res.groupby("direction")["tau_s"].transform("max")
    res["nominal_tau_s"] = nominal
    # D/L^2 (1/s): multiply by any L^2 (cm2) to get D for another diffusion length
    res["D_over_L2_1_s"] = 4 / (np.pi * res["tau_s"]) * res["dEs_over_dEt"] ** 2
    res["D_cm2_s"] = res["D_over_L2_1_s"] * L_cm ** 2
    for key, L in lengths_cm().items():
        res[f"D_{key}_cm2_s"] = res["D_over_L2_1_s"] * L ** 2
    # tau*D/L^2 = 4/pi*(dEs/dEt)^2 : must be << 1 for the equation to hold (independent of L)
    res["tauD_over_L2"] = 4 / np.pi * res["dEs_over_dEt"] ** 2

    reasons = []
    for i, r in res.iterrows():
        why = []
        if SKIP_FIRST_PULSE and r["pulse"] == 1:
            why.append("first pulse (rest before not relaxed)")
        if r["pulse"] in EXCLUDE_PULSES.get(r["direction"], []):
            why.append("left out by hand")
        if r["disturbed"] and r["tau_s"] >= MIN_PULSE_FRACTION * nominal[i]:
            why.append("sudden potential step during the pulse or rest")
        prev = res[(res["direction"] == r["direction"]) & (res["pulse"] == r["pulse"] - 1)]
        if len(prev) and bool(prev["disturbed"].iloc[0]) and not r["disturbed"]:
            why.append("starts from the rest of a disturbed pulse")
        if r["tau_s"] < MIN_PULSE_FRACTION * nominal[i]:
            why.append("pulse stopped early by voltage limit")
        expected = -1 if r["direction"] == "sodiation" else 1
        if np.sign(r["dEs_mV"]) != expected or np.sign(r["dEt_mV"]) != expected:
            why.append("dEs or dEt has the wrong sign")
        reasons.append("; ".join(why))
    res["excluded_because"] = reasons
    res["used"] = res["excluded_because"] == ""
    notes = []
    for _, r in res.iterrows():
        n = []
        if r["tauD_over_L2"] > APPARENT_LIMIT:
            n.append("dEs not << dEt (apparent D)")
        if r["tau_s"] < 0.98 * r["nominal_tau_s"]:
            n.append("cut short by voltage limit")
        if abs(r["rest_end_drift_mV_h"]) > REST_DRIFT_LIMIT_MV_H:
            n.append("rest not fully relaxed")
        notes.append("; ".join(n))
    res["note"] = notes
    return res


def find_problems(res):
    """Plain-language list of problems in the measurement."""
    out = []
    for direction in DIRECTIONS:
        r = res[res["direction"] == direction]
        if r.empty:
            continue
        u = r[r["used"]]
        D = direction.capitalize()
        # 1. loop kept running after the voltage limit
        cut = r[r["tau_s"] < MIN_PULSE_FRACTION * r["nominal_tau_s"]]
        if len(cut) > 1:
            end_h = r["t_start_h"].iloc[-1] + r["tau_s"].iloc[-1] / 3600 + r["rest_h"].iloc[-1]
            out.append(f"{D}: {len(cut)} pulses stopped early at the voltage limit, and the loop kept going for "
                       f"{end_h - cut['t_start_h'].min():.1f} h after the limit was first reached. These pulses "
                       f"are not used. In EC-Lab set EL = pass and nc to the number of pulses you need.")
        if u.empty:
            out.append(f"{D}: no pulse could be used for D.")
            continue
        # 2. rest too short
        drift = u["rest_end_drift_mV_h"].abs()
        n_bad = int((drift > REST_DRIFT_LIMIT_MV_H).sum())
        if n_bad > 0.25 * len(u):
            out.append(f"{D}: {n_bad} of {len(u)} rests were still drifting more than {REST_DRIFT_LIMIT_MV_H:g} mV/h "
                       f"when the next pulse started (up to {drift.max():.0f} mV/h). The relaxed potentials there are not at "
                       f"equilibrium, so D is less certain. A longer rest or a dE/dt stop condition helps.")
        # 3. sudden jump in the IR drop
        ir = u.sort_values("pulse")[["pulse", "IR_drop_mV"]].assign(IR=lambda x: x["IR_drop_mV"].abs())
        diff = ir["IR"].diff().abs()
        thr = max(8.0, 5 * np.nanmedian(diff)) if diff.notna().sum() > 3 else np.inf
        jumps = ir[diff > thr]
        if len(jumps):
            first = ir.index.get_loc(jumps.index[0])
            last = ir.index.get_loc(jumps.index[-1])
            p0, p1 = ir.iloc[max(first - 1, 0)], ir.iloc[last]
            out.append(f"{D}: the IR drop jumped from {p0['IR']:.0f} to {p1['IR']:.0f} mV between pulses "
                       f"{int(p0['pulse'])} and {int(p1['pulse'])} while the relaxed potential stayed smooth. This "
                       f"usually comes from the Na counter electrode or a contact, not from your electrode; "
                       f"D may show a step there.")
        # 3b. pulses with a sudden potential step
        full = r[r["tau_s"] >= MIN_PULSE_FRACTION * r["nominal_tau_s"]]
        steps = full.loc[full["disturbed"], "pulse"].astype(int).tolist()
        if steps:
            out.append(f"{D}: pulse{'s' if len(steps) > 1 else ''} {', '.join(map(str, steps))} show a sudden "
                       f"potential step during the pulse or the rest after it (a contact or counter-electrode "
                       f"disturbance). {'They are' if len(steps) > 1 else 'It is'} not used.")
        # 4. plateau not reached (sodiation only)
        if direction == "sodiation" and (r["tau_s"] < 0.98 * r["nominal_tau_s"]).any():
            last_full = r[r["tau_s"] >= MIN_PULSE_FRACTION * r["nominal_tau_s"]].iloc[-1]
            lowest_relaxed = r["E_relaxed_after_V"].min()
            lower_limit = r["E_pulse_end_V"].min()
            if lowest_relaxed - lower_limit > 0.04:
                over = (last_full["E_relaxed_after_V"] - last_full["E_pulse_end_V"]) * 1000
                out.append(f"Sodiation: the relaxed potential never went below {lowest_relaxed * 1000:.0f} mV. Under "
                           f"current the cell reached the lower limit about {over:.0f} mV earlier, so the "
                           f"low-voltage plateau was not measured.")
        # 5. Weppner-Huggins assumption weak
        n_app = int((u["tauD_over_L2"] > APPARENT_LIMIT).sum())
        if n_app:
            out.append(f"{D}: {n_app} of {len(u)} points have dEs not much smaller than dEt, so the "
                       f"Weppner-Huggins assumption is weak there (open symbols, apparent D).")
    return out


def build_summary(res):
    mass_g = MASS_MG / 1000
    rows = {}
    for direction in DIRECTIONS:
        r = res[res["direction"] == direction]
        if r.empty:
            continue
        u = r[r["used"]]
        I = r["I_uA"].abs().median()
        rows[direction.capitalize()] = {
            "Pulses": len(r),
            "Pulses used for D": len(u),
            "Pulses cut short by the voltage limit": int((r["tau_s"] < 0.98 * r["nominal_tau_s"]).sum()),
            "Current (µA)": f"{I:.3f}",
            "C-rate": f"C/{C_RATE_BASIS_MAH_G * mass_g * 1000 / I:.0f} (1C = {C_RATE_BASIS_MAH_G} mAh/g)",
            "Pulse length (s)": f"{r['nominal_tau_s'].iloc[0]:.0f}",
            "Rest length (h)": f"{r['rest_h'].median():.2f}",
            "Capacity (mAh/g)": f"{r['capacity_end_mAh_g'].iloc[-1]:.1f}",
            "Capacity in used pulses (mAh/g)": f"{u['capacity_end_mAh_g'].max():.1f}" if len(u) else "-",
            "Relaxed potential range (V)": f"{r['E_relaxed_after_V'].min():.3f} - {r['E_relaxed_after_V'].max():.3f}",
            "Median IR drop, used pulses (mV)": f"{u['IR_drop_mV'].abs().median():.1f}" if len(u) else "-",
            f"Rests drifting > {REST_DRIFT_LIMIT_MV_H:g} mV/h (used)":
                f"{int((u['rest_end_drift_mV_h'].abs() > REST_DRIFT_LIMIT_MV_H).sum())} of {len(u)}",
            "D min (cm²/s)": f"{u['D_cm2_s'].min():.2e}" if len(u) else "-",
            "D median (cm²/s)": f"{u['D_cm2_s'].median():.2e}" if len(u) else "-",
            "D max (cm²/s)": f"{u['D_cm2_s'].max():.2e}" if len(u) else "-",
        }
    s = pd.DataFrame(rows)
    if {"Sodiation", "Desodiation"} <= set(s.columns):
        sod = res.loc[res["direction"] == "sodiation", "capacity_end_mAh_g"].iloc[-1]
        des = res.loc[res["direction"] == "desodiation", "capacity_end_mAh_g"].iloc[-1]
        s.loc["Coulombic efficiency of this cycle (%)"] = ["", f"{100 * des / sod:.1f}"]
    return s.reset_index().rename(columns={"index": "Quantity"})


# ----------------------------------------------------------------------------- figures
def fig_overview(d, res, out, sample):
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.2))
    a.plot(d["t"] / 3600, d["E"], color="0.25", lw=0.6)
    top = float(np.nanmax(d["E"])) + 0.35
    for direction in DIRECTIONS:
        r = res[res["direction"] == direction]
        if len(r):
            a.axvspan(r["t_start_h"].min(), r["t_start_h"].max() + r["rest_h"].iloc[-1] + r["tau_s"].iloc[-1] / 3600,
                      color=COL[direction], alpha=0.10, lw=0)
            a.text((r["t_start_h"].min() + r["t_start_h"].max()) / 2, top - 0.23,
                   direction.capitalize() + f" ({len(r)} pulses)", ha="center", va="bottom", fontsize=9, color="0.2")
    a.set_xlabel("Time (h)"); a.set_ylabel(POT_LABEL); a.set_ylim(-0.05, top)
    a.set_title("(a) Full GITT run", loc="left", fontsize=10)

    for direction in DIRECTIONS:
        r = res[res["direction"] == direction]
        if not len(r):
            continue
        b.plot(r["capacity_end_mAh_g"], r["E_relaxed_after_V"], MRK[direction] + "-", color=COL[direction],
               ms=3.5, lw=1.0, label=f"{direction.capitalize()}: relaxed (quasi-OCV)")
        b.plot(r["capacity_end_mAh_g"], r["E_pulse_end_V"], "-", color=COL[direction], lw=0.8, alpha=0.55,
               label=f"{direction.capitalize()}: end of pulse")
    b.set_xlabel(CAP_LABEL); b.set_ylabel(POT_LABEL); b.set_ylim(-0.05, top - 0.25)
    b.legend(fontsize=8, loc="upper center")
    b.set_title("(b) Quasi-OCV and pulse-end potential", loc="left", fontsize=10)
    fig.suptitle(f"{sample}: GITT overview", fontsize=11)
    fig.tight_layout()
    path = out / "02_GITT_overview.png"
    fig.savefig(path); plt.close(fig)
    return path


def fig_diffusion(res, out, sample, L_text):
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True)
    handles = []
    any_open = False
    for direction in DIRECTIONS:
        r = res[(res["direction"] == direction) & res["used"]]
        if not len(r):
            continue
        ok = r["tauD_over_L2"] <= APPARENT_LIMIT
        any_open |= bool((~ok).any())
        for ax, x in ((a, "E_mid_V"), (b, "capacity_mid_mAh_g")):
            ax.plot(r[x], r["D_cm2_s"], "-", color=COL[direction], lw=0.8, alpha=0.6, zorder=1)
            ax.plot(r.loc[ok, x], r.loc[ok, "D_cm2_s"], MRK[direction], color=COL[direction], ms=5, zorder=3)
            ax.plot(r.loc[~ok, x], r.loc[~ok, "D_cm2_s"], MRK[direction], mfc="white", mec=COL[direction],
                    ms=5, zorder=3)
        handles.append(Line2D([], [], color=COL[direction], marker=MRK[direction], ms=5, lw=0.8,
                              label=direction.capitalize()))
    if any_open:
        handles.append(Line2D([], [], color="0.4", marker="o", mfc="white", ls="", ms=5,
                              label=r"open: $\Delta E_s$ not $\ll \Delta E_\tau$ (apparent)"))
    for ax in (a, b):
        ax.set_yscale("log")
    a.set_xlabel(POT_LABEL); a.set_ylabel(D_LABEL)
    b.set_xlabel(CAP_LABEL)
    a.set_title("(a) vs potential", loc="left", fontsize=10)
    b.set_title("(b) vs specific capacity", loc="left", fontsize=10)
    if handles:
        a.legend(handles=handles, fontsize=8, loc="best")
    fig.suptitle(f"{sample}: Na$^+$ diffusion coefficient from GITT   ({L_text})", fontsize=10.5)
    fig.tight_layout()
    path = out / "01_diffusion_coefficient.png"
    fig.savefig(path); plt.close(fig)
    return path


def fig_checks(res, out, sample):
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.2))
    items = [
        ("rest_end_drift_mV_h", "Drift at end of rest (mV h$^{-1}$)", "(a) Is the rest long enough?"),
        ("IR_drop_mV", "IR drop in first second (mV)", "(b) IR drop"),
        ("tau_s", "Pulse length (s)", "(c) Pulses cut short by the voltage limit"),
        ("tauD_over_L2", r"$\tau D/L^2 = (4/\pi)(\Delta E_s/\Delta E_\tau)^2$", r"(d) Validity: must be $\ll$ 1"),
    ]
    absolute = ("rest_end_drift_mV_h", "IR_drop_mV")
    for ax, (col, ylab, title) in zip(axs.flat, items):
        for direction in DIRECTIONS:
            r = res[res["direction"] == direction]
            y = r[col].abs() if col in absolute else r[col]
            ax.plot(r["pulse"], y, MRK[direction], color=COL[direction], ms=3.5, label=direction.capitalize())
            bad = r[~r["used"]]
            yb = bad[col].abs() if col in absolute else bad[col]
            ax.plot(bad["pulse"], yb, "x", color="0.15", ms=4.5, mew=0.9)
        ax.set_xlabel("Pulse number"); ax.set_ylabel(ylab); ax.set_title(title, loc="left", fontsize=10)
    axs[0, 0].axhline(REST_DRIFT_LIMIT_MV_H, color="0.3", ls="--", lw=0.8)
    axs[0, 0].text(0.98, REST_DRIFT_LIMIT_MV_H * 0.75, f"limit {REST_DRIFT_LIMIT_MV_H:g} mV/h", fontsize=8,
                   color="0.3", ha="right", va="top", transform=axs[0, 0].get_yaxis_transform())
    axs[1, 1].axhline(APPARENT_LIMIT, color="0.3", ls="--", lw=0.8)
    for ax in (axs[0, 0], axs[1, 0], axs[1, 1]):
        ax.set_yscale("log")
    h = [Line2D([], [], color=COL[x], marker=MRK[x], ls="", ms=4, label=x.capitalize()) for x in DIRECTIONS]
    h.append(Line2D([], [], color="0.15", marker="x", ls="", ms=5, label="not used for D"))
    axs[0, 1].legend(handles=h, fontsize=8)
    fig.suptitle(f"{sample}: GITT data checks", fontsize=11)
    fig.tight_layout()
    path = out / "03_GITT_data_checks.png"
    fig.savefig(path); plt.close(fig)
    return path


def write_excel(res, out, L_text, L_cm, summary, problems, src_name):
    base = ["direction", "pulse", "used", "excluded_because", "note", "t_start_h", "I_uA", "tau_s",
            "Q_pulse_mAh_g", "capacity_mid_mAh_g", "capacity_end_mAh_g", "E_relaxed_before_V",
            "E_relaxed_after_V", "E_mid_V", "E_pulse_end_V", "IR_drop_mV", "dEs_mV", "dEt_mV",
            "dEs_over_dEt", "tauD_over_L2", "D_cm2_s", "D_over_L2_1_s"]
    extra = [c for c in res.columns if c.startswith("D_") and c not in base]
    cols = base + extra + ["sqrt_t_R2", "rest_h", "rest_end_drift_mV_h"]
    settings = pd.DataFrame({
        "setting": ["file", "sample", "mass (mg)", "diameter (mm)", "thickness (mm)", "SEM diameter (um)",
                    "length mode", "diffusion length L (cm)", "L description", "IR time (s)",
                    "rest drift limit (mV/h)", "pulses left out by hand", "equation"],
        "value": [src_name, SAMPLE_NAME, MASS_MG, DIAMETER_MM, THICKNESS_MM, SEM_DIAMETER_UM, LENGTH_MODE, L_cm,
                  L_text.replace("$_M$", "_M"), IR_TIME_S, REST_DRIFT_LIMIT_MV_H, str(EXCLUDE_PULSES),
                  "D = 4/(pi*tau) * L^2 * (dEs/dEt)^2  (Weppner & Huggins 1977)"]})
    path = out / "gitt_results.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        summary.to_excel(xw, sheet_name="Summary", index=False)
        pd.DataFrame({"Problems found": problems or ["None"]}).to_excel(xw, sheet_name="Data checks", index=False)
        res[cols].to_excel(xw, sheet_name="Per pulse", index=False)
        settings.to_excel(xw, sheet_name="Settings", index=False)
        for ws in xw.book.worksheets:
            for colcells in ws.columns:
                if str(colcells[0].value).startswith("D_"):
                    for c in colcells[1:]:
                        c.number_format = "0.00E+00"
                w = max(len(str(c.value)) if c.value is not None else 0 for c in colcells[:60])
                ws.column_dimensions[colcells[0].column_letter].width = min(max(10, w + 2), 90)
    return path


# ----------------------------------------------------------------------------- main
def run(path=FILE, data=None, out=None, verbose=True):
    """Analyse one GITT file. `data` = already loaded table (from load()), to skip reading again."""
    global MASS_MG, SAMPLE_NAME
    path = Path(path)
    out = Path(out or OUT_DIR or path.with_name(path.stem + "_GITT_results"))
    out.mkdir(parents=True, exist_ok=True)
    if MASS_MG is None or SAMPLE_NAME is None:
        mass, name, _ = read_cell_info(path)
        MASS_MG = MASS_MG or mass
        SAMPLE_NAME = SAMPLE_NAME or name or path.stem
    if not MASS_MG:
        raise ValueError("No active mass found in the file. Set MASS_MG.")
    d = load(path) if data is None else data
    d, segs = segments(d)
    L_cm, L_text = diffusion_length_cm()
    res = analyse(d, segs, L_cm)
    problems = find_problems(res)
    summary = build_summary(res)
    figures = [fig_diffusion(res, out, SAMPLE_NAME, L_text),
               fig_overview(d, res, out, SAMPLE_NAME),
               fig_checks(res, out, SAMPLE_NAME)]
    excel = write_excel(res, out, L_text, L_cm, summary, problems, path.name)
    if verbose:
        pd.set_option("display.width", 200)
        print(summary.to_string(index=False))
        print("\nData checks:")
        for p in problems or ["No problems found."]:
            print(" - " + p)
        print(f"\nSaved to {out}")
    return {"per_pulse": res, "summary": summary, "problems": problems, "figures": figures, "excel": excel,
            "out": out, "L_text": L_text.replace("$_M$", "_M"), "L_cm": L_cm, "lengths": lengths_cm()}


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else FILE)
