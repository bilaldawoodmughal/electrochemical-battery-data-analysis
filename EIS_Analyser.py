"""
EIS Analyser for EC-Lab files
=============================

Upload EC-Lab impedance files (.mpr, or .mpt text exports) and get

  * a Kramers-Kronig quality check (Lin-KK, Schoenleber et al., Electrochim. Acta 131 (2014) 20)
    that also suggests how far down in frequency each spectrum can be trusted,
  * an equivalent-circuit fit with the circuit of Jin et al., Electrochim. Acta 141 (2014) 302:
        Re - (Rf || CPE1) - (CPE2 || (Rct - Zw))
    (or a one-arc Randles circuit), with parameter uncertainties,
  * report-ready Nyquist and Bode figures (PNG, PDF, SVG),
  * the fit table (CSV and LaTeX) and the data with fitted curves (CSV, e.g. for Origin).

Run it on its own:        streamlit run EIS_Analyser.py
Inside Neware by Bilal:   page_eis.py calls main(embedded=True).
"""
import io
import logging
import re
import zipfile

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.lines import Line2D                 # noqa: E402
from matplotlib.patches import Rectangle            # noqa: E402

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

APP_TITLE = "EIS Analyser"

# =====================================================================================
# 1. Reading EC-Lab files
# =====================================================================================
# EC-Lab names files  <name>_<NN>_<TECHNIQUE>_<channel>.mpr , e.g. "Sample A_1,06mg_05_PEIS_03.mpr"
FNAME_RE = re.compile(r"^(?P<base>.+?)_(?P<seq>\d{1,3})_(?P<tech>[A-Za-z]+)(?:_C?(?P<ch>\d{1,3}))?$")


def _find_col(columns, *prefixes):
    for c in columns:
        name = str(c).strip().lower()
        if any(name.startswith(p.lower()) for p in prefixes):
            return c
    return None


def _read_mpr(data: bytes) -> pd.DataFrame:
    try:
        import galvani_biologic as BioLogic          # copy bundled with Neware by Bilal
    except ImportError:
        try:
            from galvani import BioLogic
        except ImportError as exc:
            raise RuntimeError("Reading .mpr files needs the 'galvani' package. "
                               "Install it with:  pip install galvani") from exc
    try:
        mpr = BioLogic.MPRfile(io.BytesIO(data))
    except Exception as exc:  # unknown column types in newer EC-Lab versions etc.
        raise RuntimeError(f"Could not read this .mpr file ({exc}). In EC-Lab, use "
                           "'Export as Text' and upload the .mpt file instead.") from exc
    return pd.DataFrame(mpr.data)


def _read_mpt(data: bytes) -> pd.DataFrame:
    """EC-Lab text export (.mpt / .txt), with or without the 'EC-Lab ASCII FILE' header.
    Works with decimal points and decimal commas."""
    lines = data.decode("latin-1").splitlines()
    skip = 0
    if lines and lines[0].startswith("EC-Lab ASCII FILE"):
        m = re.search(r"Nb header lines\s*:\s*(\d+)", "\n".join(lines[:5]))
        if m:
            skip = int(m.group(1)) - 1
    header = [h.strip() for h in lines[skip].split("\t")]
    rows = [ln.split("\t") for ln in lines[skip + 1:] if ln.strip()]
    df = pd.DataFrame(rows)
    df = df.iloc[:, :len(header)]
    df.columns = header[:df.shape[1]]
    df = df.loc[:, [c for c in df.columns if c]]
    return df.apply(lambda s: pd.to_numeric(s.astype(str).str.strip().str.replace(",", ".", regex=False),
                                            errors="coerce"))


def read_eclab(name: str, data: bytes):
    """Return (group, seq, technique, spectra). Each spectrum: dict(f, Z, E), f from high to low.
    'group' is the sample/run name taken from the EC-Lab file name; non-EIS files give no spectra."""
    stem = re.sub(r"\.(mpr|mpt|txt)$", "", name, flags=re.IGNORECASE)
    m = FNAME_RE.match(stem)
    group = m.group("base").strip() if m else stem
    seq = int(m.group("seq")) if m else 0
    tech = m.group("tech").upper() if m else ""
    df = _read_mpr(data) if name.lower().endswith(".mpr") else _read_mpt(data)
    fcol = _find_col(df.columns, "freq")
    recol = _find_col(df.columns, "re(z)")
    imcol = _find_col(df.columns, "-im(z)")
    if fcol is None or recol is None or imcol is None:
        return group, seq, tech, []
    ecol = _find_col(df.columns, "<ewe>", "ewe")
    ccol = _find_col(df.columns, "cycle number")
    f = df[fcol].to_numpy(float)
    Z = df[recol].to_numpy(float) - 1j * df[imcol].to_numpy(float)
    E = df[ecol].to_numpy(float) if ecol is not None else np.full(len(f), np.nan)
    ok = np.isfinite(f) & (f > 0) & np.isfinite(Z.real) & np.isfinite(Z.imag)
    # one file can hold several sweeps: split on the cycle number and where the frequency jumps back up
    cyc = df[ccol].to_numpy(float) if ccol is not None else np.zeros(len(f))
    brk = np.zeros(len(f), bool)
    brk[1:] = (np.diff(cyc) != 0) | (f[1:] > 2 * f[:-1])
    seg = np.cumsum(brk)
    spectra = []
    for s in np.unique(seg[ok]):
        idx = np.where(ok & (seg == s))[0]
        if len(idx) < 10:
            continue
        o = idx[np.argsort(-f[idx])]
        e = E[o]
        spectra.append(dict(f=f[o], Z=Z[o], E=float(np.nanmean(e)) if np.isfinite(e).any() else np.nan))
    return group, seq, tech, spectra


# =====================================================================================
# 2. Kramers-Kronig test (Lin-KK) and automatic fit window
# =====================================================================================
def lin_kk(f, Z, c=0.85, max_M=80):
    """Lin-KK (Schoenleber et al. 2014): fit R0 + jwL + 1/(jwC) + sum Rk/(1+jw tau_k) by weighted
    linear least squares; M chosen with the mu < c criterion. Returns residuals |dZ|/|Z| in %."""
    f = np.asarray(f, float); Z = np.asarray(Z, complex)
    w = 2 * np.pi * f
    wgt = 1 / np.abs(Z)
    for M in range(2, max_M + 1):
        taus = np.logspace(np.log10(1 / w.max()), np.log10(1 / w.min()), M)
        A = np.column_stack([np.ones_like(w, dtype=complex), 1j * w, 1 / (1j * w)]
                            + [1 / (1 + 1j * w * t) for t in taus])
        Aw = A * wgt[:, None]
        p, *_ = np.linalg.lstsq(np.vstack([Aw.real, Aw.imag]),
                                np.concatenate([(Z * wgt).real, (Z * wgt).imag]), rcond=None)
        Rk = p[3:]
        pos, neg = Rk[Rk >= 0].sum(), -Rk[Rk < 0].sum()
        mu = 1 - neg / pos if pos > 0 else 1
        if mu <= c:
            break
    return 100 * np.abs(Z - A @ p) / np.abs(Z)


def suggest_fmin(f, res, thr=1.0, gap_tol=2, min_bad=3):
    """Lowest frequency to fit. Starting at the low-frequency end, follow the run of points whose
    KK residual exceeds `thr` % (allowing up to `gap_tol` good points inside the run). If that run
    contains at least `min_bad` bad points, everything up to its top is cut (drift / non-stationarity).
    Isolated stray points at higher frequency are left in."""
    order = np.argsort(f)
    fs, bad = np.asarray(f)[order], np.asarray(res)[order] > thr
    nbad, good_run, top = 0, 0, None
    for i in range(len(fs)):
        if bad[i]:
            nbad += 1; good_run = 0; top = i
        else:
            good_run += 1
            if good_run > gap_tol:
                break
    if top is None or nbad < min_bad:
        return float(fs[0])
    return float(fs[min(top + 1, len(fs) - 1)])


# =====================================================================================
# 3. Equivalent-circuit model and fit
# =====================================================================================
CIRCUITS = {
    "jin": "Re–(Rf‖CPE1)–(CPE2‖(Rct–Zw))  [Jin et al. 2014, two arcs]",
    "randles": "Re–(CPE‖(Rct–Zw))  [Randles, one arc]",
}
PARAMS = {"jin": ["Re", "Rf", "Q1", "n1", "Rct", "Q2", "n2", "Aw", "alpha"],
          "randles": ["Re", "Rct", "Q2", "n2", "Aw", "alpha"]}
LOG_PARAMS = {"Re", "Rf", "Q1", "Rct", "Q2", "Aw"}          # fitted as log10
BOUNDS = {"Rf": (-2, 6), "Q1": (-10, -1), "n1": (0.4, 1.0), "Rct": (-1, 7), "Q2": (-10, 0),
          "n2": (0.3, 1.0), "Aw": (-2, 8), "alpha": (0.3, 1.0)}
START = {"Rf": (0, 3.5), "Q1": (-7, -3), "n1": (0.6, 0.95), "Rct": (0.5, 4.5), "Q2": (-6, -2),
         "n2": (0.5, 0.95), "Aw": (1, 5), "alpha": (0.4, 0.8)}


def circuit_Z(model, p, f):
    """Impedance of the circuit. Zw = Aw (jw)^-alpha (alpha = 0.5: ideal 45-degree Warburg)."""
    jw = 1j * 2 * np.pi * np.asarray(f, float)
    Zw = p["Aw"] / jw ** p["alpha"]
    Z = p["Re"] + 1 / (p["Q2"] * jw ** p["n2"] + 1 / (p["Rct"] + Zw))
    if model == "jin":
        Z = Z + p["Rf"] / (1 + p["Rf"] * p["Q1"] * jw ** p["n1"])
    return Z


def _tau(R, Q, n):
    return (R * Q) ** (1 / n)


def hf_intercept(Z):
    k = np.polyfit(-Z.imag[:3], Z.real[:3], 1)
    v = float(np.polyval(k, 0))
    return v if v > 0 else float(max(np.min(Z.real), 1e-3))


def fit_spectrum(f, Z, model="jin", ideal_w=False, n_starts=80, prev=None, consistency=0.15, seed=0):
    """Complex non-linear least squares, modulus weighting, many random starts.
    For 'jin' the arc with the shorter time constant is always Rf (high-frequency arc), as in the paper.
    If `prev` (parameters of the previous spectrum of the same sample) is given, the solution closest
    to it is used when its error is within `consistency` (fraction) of the best one; this keeps the
    same arc assignment through a series of near-identical spectra."""
    f = np.asarray(f, float); Z = np.asarray(Z, complex)
    names = [n for n in PARAMS[model] if not (ideal_w and n == "alpha")]
    fixed = {"alpha": 0.5} if ideal_w else {}
    Re0 = hf_intercept(Z)
    lo = np.array([np.log10(Re0) - 0.5 if n == "Re" else BOUNDS[n][0] for n in names])
    hi = np.array([np.log10(Re0) + 0.5 if n == "Re" else BOUNDS[n][1] for n in names])

    def unpack(x):
        p = dict(fixed)
        for n, v in zip(names, x):
            p[n] = 10 ** v if n in LOG_PARAMS else v
        return p

    def pack(p):
        return np.array([np.log10(p[n]) if n in LOG_PARAMS else p[n] for n in names])

    wgt = 1 / np.abs(Z)

    def resid(x):
        d = (circuit_Z(model, unpack(x), f) - Z) * wgt
        return np.concatenate([d.real, d.imag])

    def valid(p):
        return model != "jin" or _tau(p["Rf"], p["Q1"], p["n1"]) < _tau(p["Rct"], p["Q2"], p["n2"])

    rng = np.random.default_rng(seed)
    starts = [pack(prev)] if prev is not None else []
    for _ in range(n_starts):
        starts.append(np.array([np.log10(Re0) if n == "Re" else rng.uniform(*START[n]) for n in names]))
    sols = []
    with np.errstate(all="ignore"):
        for x0 in starts:
            x0 = np.clip(x0, lo + 1e-6, hi - 1e-6)
            try:
                r = least_squares(resid, x0, bounds=(lo, hi), method="trf", max_nfev=1500)
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                continue
            if np.isfinite(r.cost) and valid(unpack(r.x)):
                sols.append((r.cost, r.x))
    if not sols:
        raise RuntimeError("The fit did not converge for this spectrum.")
    sols.sort(key=lambda s: s[0])
    best_cost, chosen = sols[0]
    if prev is not None and consistency > 0:
        keys = ["Rf", "Rct"] if model == "jin" else ["Rct"]
        target = np.log10([prev[k] for k in keys])
        cands = [x for c, x in sols if c <= best_cost * (1 + consistency)]
        chosen = min(cands, key=lambda x: np.linalg.norm(np.log10([unpack(x)[k] for k in keys]) - target))
    with np.errstate(all="ignore"):
        r = least_squares(resid, chosen, bounds=(lo, hi), method="trf", max_nfev=5000, xtol=1e-12, ftol=1e-12)
    x = r.x
    p = unpack(x)
    r0 = resid(x)
    npts = len(f)
    chi2 = float(r0 @ r0 / (2 * npts - len(x)))
    # uncertainties from the Jacobian; parameters sitting on a bound are held fixed
    free = [i for i in range(len(x)) if (x[i] - lo[i]) > 1e-4 and (hi[i] - x[i]) > 1e-4]
    err = {n: np.nan for n in names}
    if free:
        h = 1e-6
        J = np.column_stack([(resid(np.where(np.arange(len(x)) == i, x + h, x)) - r0) / h for i in free])
        try:
            sd = np.sqrt(np.abs(np.diag(np.linalg.inv(J.T @ J) * chi2)))
            for i, s in zip(free, sd):
                n = names[i]
                err[n] = float(np.log(10) * p[n] * s) if n in LOG_PARAMS else float(s)
        except np.linalg.LinAlgError:
            pass
    res = 100 * np.hypot(r0[:npts], r0[npts:])
    out = dict(model=model, ideal_w=ideal_w, params=p, errors=err, chi2=chi2,
               rms=float(np.sqrt(np.mean(res ** 2))), max_res=float(res.max()), residuals=res,
               cost_ratio=float(r.cost / best_cost), n_local=len(sols))
    # derived quantities: characteristic frequency and effective (Brug) capacitance of each arc
    if model == "jin":
        t1 = _tau(p["Rf"], p["Q1"], p["n1"])
        out["f1"], out["C1"] = 1 / (2 * np.pi * t1), t1 / p["Rf"]
    t2 = _tau(p["Rct"], p["Q2"], p["n2"])
    out["f2"], out["C2"] = 1 / (2 * np.pi * t2), t2 / p["Rct"]
    # plain-language flags
    flags = []
    for i, n in enumerate(names):
        if i not in free:
            if n in ("n1", "n2") and x[i] > 0.999:
                flags.append(f"{n} = 1 (ideal capacitor)")
            elif n in ("Rf", "Rct"):
                flags.append(f"{n} at its fit limit: this arc is not resolved")
            else:
                flags.append(f"{n} at its fit limit")
    for n in ("Rf", "Rct"):
        if n in p and np.isfinite(err.get(n, np.nan)) and err[n] / p[n] > 0.5:
            flags.append(f"{n} poorly determined (±{100 * err[n] / p[n]:.0f} %)")
    out["flags"] = flags
    return out


# =====================================================================================
# 4. Figures
# =====================================================================================
ORIGIN_COLORS = ["#000000", "#e8000b", "#0000ff", "#c000c0", "#808000", "#000080", "#8000ff", "#800000",
                 "#008080", "#ff8000"]
ORIGIN_MARKERS = ["s", "o", "^", "v", "D", "<", ">", "p", "h", "*"]
SERIF = ["Times New Roman", "TeX Gyre Termes", "Liberation Serif", "Nimbus Roman", "STIXGeneral", "DejaVu Serif"]


def _style():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": SERIF, "mathtext.fontset": "stix", "font.size": 11,
        "axes.linewidth": 1.2, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "xtick.major.size": 5, "ytick.major.size": 5,
        "xtick.minor.size": 2.5, "ytick.minor.size": 2.5, "xtick.minor.visible": True,
        "ytick.minor.visible": True, "xtick.major.width": 1.0, "ytick.major.width": 1.0,
        "savefig.facecolor": "white", "figure.facecolor": "white", "axes.facecolor": "white"})


def nice_ceil(v):
    """Round up to 1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6 or 8 x 10^k."""
    if not np.isfinite(v) or v <= 0:
        return 1.0
    k = np.floor(np.log10(v))
    for m in (1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10):
        if m * 10 ** k >= v * (1 - 1e-9):
            return float(m * 10 ** k)
    return float(10 ** (k + 1))


def nice_ticks(lim):
    for step in (1, 2, 2.5, 5):
        s = step * 10 ** np.floor(np.log10(lim / 4))
        if lim / s <= 5.01:
            return np.arange(0, lim + s / 2, s)
    return np.linspace(0, lim, 5)


def draw_circuit(ax, model):
    ax.set_xlim(0, 10); ax.set_ylim(0, 3.2); ax.axis("off")
    k = dict(color="k", lw=1.2, solid_capstyle="round")

    def wire(x0, y0, x1, y1):
        ax.plot([x0, x1], [y0, y1], **k)

    def resistor(x0, x1, y, label, below=False):
        xs = np.r_[x0, np.linspace(x0 + 0.08, x1 - 0.08, 12), x1]
        ys = np.r_[y, y + 0.17 * np.tile([1, -1], 6), y]
        ax.plot(xs, ys, **k)
        ax.text((x0 + x1) / 2, y + (-0.33 if below else 0.3), label, ha="center",
                va="top" if below else "bottom", fontsize=10.5)

    def cpe(xc, y, label):
        for dx in (-0.1, 0.1):
            ax.plot([xc + dx - 0.12, xc + dx + 0.06, xc + dx - 0.12], [y + 0.22, y, y - 0.22], **k)
        wire(xc - 0.6, y, xc - 0.16, y); wire(xc + 0.16, y, xc + 0.6, y)
        ax.text(xc, y + 0.3, label, ha="center", va="bottom", fontsize=10.5)

    def warburg(x0, x1, y, label):
        ax.add_patch(Rectangle((x0, y - 0.2), x1 - x0, 0.4, fill=False, lw=1.2))
        ax.text((x0 + x1) / 2, y, "W", ha="center", va="center", fontsize=9)
        ax.text((x0 + x1) / 2, y - 0.33, label, ha="center", va="top", fontsize=10.5)

    def node(x, y):
        ax.plot(x, y, "o", color="k", ms=4.5)

    yT, yB = 2.2, 1.15
    wire(0.0, yT, 0.4, yT); resistor(0.4, 1.6, yT, r"$R_\mathrm{e}$"); wire(1.6, yT, 2.2, yT)
    if model == "jin":
        node(2.2, yT); wire(2.2, yT, 2.2, yB); wire(4.6, yT, 4.6, yB)
        wire(2.2, yT, 2.8, yT); cpe(3.4, yT, r"CPE$_1$"); wire(4.0, yT, 4.6, yT)
        wire(2.2, yB, 2.8, yB); resistor(2.8, 4.0, yB, r"$R_\mathrm{f}$", below=True); wire(4.0, yB, 4.6, yB)
        node(4.6, yT); wire(4.6, yT, 5.2, yT)
        x0, cpe_label = 5.2, r"CPE$_2$"
    else:
        x0, cpe_label = 2.2, "CPE"
    node(x0, yT); wire(x0, yT, x0, yB); wire(x0 + 4.0, yT, x0 + 4.0, yB)
    wire(x0, yT, x0 + 1.4, yT); cpe(x0 + 2.0, yT, cpe_label); wire(x0 + 2.6, yT, x0 + 4.0, yT)
    wire(x0, yB, x0 + 0.3, yB); resistor(x0 + 0.3, x0 + 1.5, yB, r"$R_\mathrm{ct}$", below=True)
    wire(x0 + 1.5, yB, x0 + 2.4, yB); warburg(x0 + 2.4, x0 + 3.4, yB, r"$Z_\mathrm{w}$")
    wire(x0 + 3.4, yB, x0 + 4.0, yB); node(x0 + 4.0, yT); wire(x0 + 4.0, yT, x0 + 4.6, yT)


def _arc_end(item):
    """Re + Rf + Rct of the fit (real-axis end of the arcs); falls back to the data."""
    fit = item.get("fit")
    if fit:
        p = fit["params"]
        return p["Re"] + p.get("Rf", 0) + p["Rct"]
    return float(np.percentile(item["Z"].real, 50))


def auto_limits(items):
    ends = [_arc_end(it) for it in items]
    lim = nice_ceil(2.0 * max(ends))
    inset = nice_ceil(1.25 * min(ends))
    return lim, (inset if lim / inset >= 3 else None)


def _plot_series(ax, items, ms, lw, clip=True):
    for it in items:
        Z, m = it["Z"], it["mask"]
        ax.plot(Z.real[m], -Z.imag[m], it["marker"], color=it["color"], ms=ms, mew=0.6, zorder=3, clip_on=clip)
        if (~m).any():
            ax.plot(Z.real[~m], -Z.imag[~m], it["marker"], mfc="white", mec=it["color"], ms=ms, mew=0.9,
                    zorder=3, clip_on=clip)
        if it.get("fit"):
            fd = np.logspace(np.log10(it["f"][m].max()), np.log10(it["f"][m].min()), 600)
            Zf = circuit_Z(it["fit"]["model"], it["fit"]["params"], fd)
            ax.plot(Zf.real, -Zf.imag, "-", color=it["color"], lw=lw, zorder=4, clip_on=clip)


def nyquist_figure(items, title="", model="jin", lim=None, inset_lim=None, show_circuit=True,
                   show_inset=True, legend_pos="upper left", inset_pos="lower right", size=(6.2, 6.0)):
    """Journal-style Nyquist plot (after Jin et al. 2014, Fig. 8). items: list of dicts with
    label, f, Z, mask (fitted points), fit (result of fit_spectrum or None), color, marker."""
    _style()
    a_lim, a_inset = auto_limits(items)
    lim = lim or a_lim
    inset_lim = inset_lim or a_inset
    fig = plt.figure(figsize=size, dpi=100)
    ax = fig.add_axes([0.14, 0.10, 0.82, 0.85])
    _plot_series(ax, items, ms=4.2, lw=1.3)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal", adjustable="box")
    t = nice_ticks(lim); ax.set_xticks(t); ax.set_yticks(t)
    ax.set_xlabel(r"$Z'$ ($\Omega$)", fontsize=14, fontweight="bold", labelpad=4)
    ax.set_ylabel(r"$-Z''$ ($\Omega$)", fontsize=14, fontweight="bold", labelpad=4)
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontsize(12)
    if show_circuit:
        draw_circuit(ax.inset_axes([0.02, 0.765, 0.50 if model == "jin" else 0.36, 0.215]), model)
    handles = [Line2D([], [], color=it["color"], marker=it["marker"], ms=5, lw=1.3, label=it["label"]) for it in items]
    if any((~it["mask"]).any() for it in items):
        handles.append(Line2D([], [], color="0.35", marker="o", mfc="white", mec="0.35", ls="none", ms=5,
                              label="Not fitted"))
    if legend_pos == "upper left":
        legend_kw = dict(loc="upper left", bbox_to_anchor=(0.015, 0.745 if show_circuit else 0.985))
    else:
        legend_kw = dict(loc=legend_pos)
    leg = ax.legend(handles=handles, frameon=False, fontsize=10.5, title=title or None,
                    title_fontproperties={"weight": "bold", "size": 12}, handlelength=2.2,
                    labelspacing=0.35, borderaxespad=0.4, **legend_kw)
    leg._legend_box.align = "left"
    if show_inset and inset_lim:
        box = {"lower right": [0.67, 0.10, 0.29, 0.29], "upper right": [0.67, 0.63, 0.29, 0.29],
               "center right": [0.67, 0.36, 0.29, 0.29]}[inset_pos]
        iax = ax.inset_axes(box)
        iax.set_zorder(5)
        _plot_series(iax, items, ms=3.3, lw=1.1)
        iax.set_xlim(0, inset_lim); iax.set_ylim(0, inset_lim); iax.set_aspect("equal", adjustable="box")
        ti = nice_ticks(inset_lim)[::2] if len(nice_ticks(inset_lim)) > 4 else nice_ticks(inset_lim)
        iax.set_xticks(ti); iax.set_yticks(ti)
        iax.tick_params(labelsize=9, length=3.5); iax.tick_params(which="minor", length=2)
        iax.set_xlabel(r"$Z'$ ($\Omega$)", fontsize=10, fontweight="bold", labelpad=1)
        iax.set_ylabel(r"$-Z''$ ($\Omega$)", fontsize=10, fontweight="bold", labelpad=1)
    return fig


def bode_figure(items, title=""):
    _style()
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.2, 6.0), dpi=100, sharex=True)
    for it in items:
        f, Z, m = it["f"], it["Z"], it["mask"]
        for ax, y in ((a1, np.abs(Z)), (a2, -np.degrees(np.angle(Z)))):
            ax.plot(f[m], y[m], it["marker"], color=it["color"], ms=4, mew=0.6)
            if (~m).any():
                ax.plot(f[~m], y[~m], it["marker"], mfc="white", mec=it["color"], ms=4, mew=0.9)
        if it.get("fit"):
            fd = np.logspace(np.log10(f[m].max()), np.log10(f[m].min()), 400)
            Zf = circuit_Z(it["fit"]["model"], it["fit"]["params"], fd)
            a1.plot(fd, np.abs(Zf), "-", color=it["color"], lw=1.2)
            a2.plot(fd, -np.degrees(np.angle(Zf)), "-", color=it["color"], lw=1.2)
    a1.set_xscale("log"); a1.set_yscale("log")
    a1.set_ylabel(r"$|Z|$ ($\Omega$)", fontsize=13, fontweight="bold")
    a2.set_ylabel(r"$-$Phase ($^\circ$)", fontsize=13, fontweight="bold")
    a2.set_xlabel("Frequency (Hz)", fontsize=13, fontweight="bold")
    a2.set_ylim(bottom=min(0, a2.get_ylim()[0]))
    a1.legend(handles=[Line2D([], [], color=it["color"], marker=it["marker"], ms=5, lw=1.2, label=it["label"])
                       for it in items], frameon=False, fontsize=10, title=title or None,
              title_fontproperties={"weight": "bold", "size": 11})
    fig.tight_layout()
    return fig


def kk_figure(items, title=""):
    _style()
    fig, ax = plt.subplots(figsize=(6.2, 3.6), dpi=100)
    for it in items:
        ax.plot(it["f"], it["kk"], "-" + it["marker"], color=it["color"], ms=3.5, lw=1.0, label=it["label"])
        ax.axvline(it["fmin"], color=it["color"], lw=0.9, ls=":")
    ax.axhline(1.0, color="0.5", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("Frequency (Hz)", fontsize=12, fontweight="bold")
    ax.set_ylabel("KK residual (%)", fontsize=12, fontweight="bold")
    ax.set_ylim(0, max(3.0, 1.1 * max(float(np.max(it["kk"])) for it in items)))
    ax.legend(frameon=False, fontsize=9, title=title or None)
    ax.text(0.99, 0.97, "dotted lines: lowest fitted frequency", transform=ax.transAxes, ha="right", va="top",
            fontsize=8.5, color="0.35")
    fig.tight_layout()
    return fig


def fig_bytes(fig, fmt, dpi=300):
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=dpi, bbox_inches="tight", pad_inches=0.05)
    return buf.getvalue()


# =====================================================================================
# 5. Tables
# =====================================================================================
def _fmt(v, e=None):
    if v is None or not np.isfinite(v):
        return ""
    a = abs(v)
    d = 0 if a >= 100 else 1 if a >= 10 else 2 if a >= 0.1 else 3
    s = f"{v:.{d}f}"
    if e is not None and np.isfinite(e):
        s += f" ± {e:.{d}f}"
    return s


def fit_table(rows):
    """rows: list of (sample, label, E, fmin, fit) -> DataFrame (numbers) for CSV."""
    recs = []
    for sample, label, E, fmin, fit in rows:
        p, e = fit["params"], fit["errors"]
        rec = {"Sample": sample, "Spectrum": label, "OCV (V)": E, "Fitted down to (Hz)": fmin,
               "Re (Ohm)": p["Re"], "Re err": e.get("Re")}
        if fit["model"] == "jin":
            rec.update({"Rf (Ohm)": p["Rf"], "Rf err": e.get("Rf"), "CPE1 Q (F s^(n-1))": p["Q1"],
                        "CPE1 n": p["n1"], "f1 (Hz)": fit["f1"], "C1 (uF)": 1e6 * fit["C1"]})
        rec.update({"Rct (Ohm)": p["Rct"], "Rct err": e.get("Rct"), "CPE2 Q (F s^(n-1))": p["Q2"],
                    "CPE2 n": p["n2"], "f2 (Hz)": fit["f2"], "C2 (uF)": 1e6 * fit["C2"],
                    "Aw (Ohm s^-alpha)": p["Aw"], "alpha": p["alpha"],
                    "Rf+Rct (Ohm)": p.get("Rf", 0) + p["Rct"],
                    "RMS residual (%)": fit["rms"], "Max residual (%)": fit["max_res"],
                    "Notes": "; ".join(fit["flags"])})
        recs.append(rec)
    return pd.DataFrame(recs)


def display_table(rows):
    recs = []
    for sample, label, E, fmin, fit in rows:
        p, e = fit["params"], fit["errors"]
        rec = {"Sample": sample, "Spectrum": label, "OCV (V)": _fmt(E), "Re (Ω)": _fmt(p["Re"])}
        if fit["model"] == "jin":
            rec["Rf (Ω)"] = _fmt(p["Rf"], e.get("Rf"))
        rec["Rct (Ω)"] = _fmt(p["Rct"], e.get("Rct"))
        rec["α"] = f"{p['alpha']:.2f}"
        rec["Fit to (Hz)"] = f"{fmin:.3g}"
        rec["RMS (%)"] = f"{fit['rms']:.1f}"
        rec["Notes"] = "; ".join(fit["flags"])
        recs.append(rec)
    return pd.DataFrame(recs)


def latex_table(rows, model):
    jin = model == "jin"
    cols = "ll" + "c" * (5 if jin else 4)
    head = r"Sample & Spectrum & $R_\mathrm{e}$ ($\Omega$) & " + (r"$R_\mathrm{f}$ ($\Omega$) & " if jin else "") \
        + r"$R_\mathrm{ct}$ ($\Omega$) & $\alpha$ \\"
    lines = [r"\begin{table}[htbp]", r"\centering",
             r"\caption{Fitted EIS parameters (equivalent circuit in Fig.~X).}", r"\label{tab:eis}",
             rf"\begin{{tabular}}{{{cols}}}", r"\toprule", head, r"\midrule"]
    last = None
    for sample, label, E, fmin, fit in rows:
        p, e = fit["params"], fit["errors"]
        cells = [sample if sample != last else "", label, _fmt(p["Re"])]
        if jin:
            cells.append("$" + _fmt(p["Rf"], e.get("Rf")).replace("±", r"\pm") + "$")
        cells.append("$" + _fmt(p["Rct"], e.get("Rct")).replace("±", r"\pm") + "$")
        cells.append(f"{p['alpha']:.2f}")
        lines.append(" & ".join(cells) + r" \\")
        last = sample
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(lines) + "\n% needs \\usepackage{booktabs}\n"


def data_csv(items, sample):
    frames = []
    for it in items:
        f, Z = it["f"], it["Z"]
        df = pd.DataFrame({"Sample": sample, "Spectrum": it["label"], "Frequency (Hz)": f,
                           "Z' (Ohm)": Z.real, "-Z'' (Ohm)": -Z.imag, "|Z| (Ohm)": np.abs(Z),
                           "-Phase (deg)": -np.degrees(np.angle(Z)), "KK residual (%)": it["kk"],
                           "Fitted": it["mask"]})
        if it.get("fit"):
            Zf = circuit_Z(it["fit"]["model"], it["fit"]["params"], f)
            df["Z' fit (Ohm)"] = np.where(it["mask"], Zf.real, np.nan)
            df["-Z'' fit (Ohm)"] = np.where(it["mask"], -Zf.imag, np.nan)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


# =====================================================================================
# 6. Streamlit app
# =====================================================================================
def main(embedded=False):
    """embedded=True: running as a page of another app (that app sets the page and draws the title)."""
    import streamlit as st

    if not embedded:
        try:
            st.set_page_config(page_title=APP_TITLE, page_icon="〰️", layout="wide")
        except Exception:  # already set by a parent multipage app
            pass
    side_header = st.subheader if embedded else st.header

    @st.cache_data(show_spinner=False)
    def cached_read(name, data):
        return read_eclab(name, data)

    @st.cache_data(show_spinner=False)
    def cached_kk(f, Z):
        return lin_kk(f, Z)

    @st.cache_data(show_spinner=False)
    def cached_fit(f, Z, model, ideal_w, n_starts, prev, consistency):
        return fit_spectrum(f, Z, model=model, ideal_w=ideal_w, n_starts=n_starts,
                            prev=dict(prev) if prev else None, consistency=consistency)

    if not embedded:
        st.title(APP_TITLE)
        st.caption("Drop in EC-Lab impedance files (.mpr or exported .mpt). The app checks the data with a "
                   "Kramers–Kronig test, fits the equivalent circuit of Jin et al. (Electrochim. Acta 141, 2014) "
                   "and draws report-ready Nyquist plots.")

    # ---------------- sidebar: settings ----------------
    with st.sidebar:
        side_header("Fit settings")
        model = st.radio("Equivalent circuit", list(CIRCUITS),
                         format_func=lambda k: {"jin": "Two arcs + Warburg (Jin et al. 2014)",
                                                "randles": "One arc + Warburg (Randles)"}[k])
        st.caption(CIRCUITS[model].split("  [")[0])
        w_choice = st.radio("Warburg element Zw", ["Generalised: Aw·(jω)^−α, α free",
                                                   "Ideal 45° line: α = 0.5"])
        ideal_w = w_choice.startswith("Ideal")
        st.subheader("Fit window")
        win_mode = st.radio("Lowest frequency used in the fit",
                            ["Automatic (Kramers–Kronig)", "Same value for all spectra", "Use all points"])
        kk_thr = 1.0
        fixed_fmin = 0.1
        if win_mode.startswith("Automatic"):
            kk_thr = st.slider("KK residual limit (%)", 0.5, 3.0, 1.0, 0.1,
                               help="Low-frequency points are dropped where the KK residuals stay above "
                                    "this limit (drift during the measurement).")
        elif win_mode.startswith("Same"):
            fixed_fmin = st.number_input("Lowest frequency (Hz)", min_value=1e-4, value=0.15, format="%.4f")
        with st.expander("Advanced"):
            consistent = st.checkbox("Keep the same arc assignment within a sample", value=True,
                                     help="When two arc splits fit almost equally well, use the one closest "
                                          "to the previous spectrum of the same sample (within 15 % of the best fit).")
            n_starts = st.slider("Random starts per fit", 20, 200, 80, 10,
                                 help="More starts = more reliable, but slower.")
        side_header("Figure")
        show_circuit = st.checkbox("Draw the circuit", value=True)
        show_inset = st.checkbox("High-frequency zoom inset", value=True)
        inset_pos = st.selectbox("Inset position", ["lower right", "upper right", "center right"])
        legend_pos = st.selectbox("Legend position", ["upper left", "upper right", "lower right", "center right"])
        add_ocv = st.checkbox("Add OCV to legend labels", value=False)
        auto_axes = st.checkbox("Automatic axis range", value=True)
        man_lim = man_inset = None
        if not auto_axes:
            man_lim = st.number_input("Main plot: max Z′ and −Z″ (Ω)", min_value=1.0, value=4000.0, step=500.0)
            man_inset = st.number_input("Inset: max Z′ and −Z″ (Ω)", min_value=1.0, value=250.0, step=50.0)

    # ---------------- upload ----------------
    files = st.file_uploader("EC-Lab files", type=["mpr", "mpt", "txt"], accept_multiple_files=True,
                             help="You can drop all files of a run; only impedance (PEIS/GEIS) files are used.")
    if not files:
        st.info("Upload one or more EC-Lab impedance files to start. Files from the same run "
                "(same name before _NN_PEIS) are grouped as one sample.")
        return

    spectra, skipped, errors = [], [], []
    for uf in files:
        try:
            group, seq, tech, specs = cached_read(uf.name, uf.getvalue())
        except Exception as exc:
            errors.append(f"**{uf.name}**: {exc}")
            continue
        if not specs:
            skipped.append(uf.name)
            continue
        for k, s in enumerate(specs):
            spectra.append(dict(id=f"{uf.name}#{k}", file=uf.name, group=group, seq=seq, sub=k, **s))
    for e in errors:
        st.error(e)
    if skipped:
        with st.expander(f"{len(skipped)} file(s) skipped because they contain no impedance data (OCV, GCPL, …)"):
            st.write(", ".join(sorted(skipped)))
    if not spectra:
        st.warning("None of the files contains impedance data.")
        return
    spectra.sort(key=lambda s: (s["group"].lower(), s["seq"], s["sub"]))
    count = {}
    for s in spectra:
        count[s["group"]] = count.get(s["group"], 0) + 1
        s["default_label"] = f"EIS {count[s['group']]}"

    st.subheader("Spectra")
    st.caption("Edit the labels (e.g. 'Fresh', 'After 3rd desodiation') and the sample names; untick spectra "
               "you don't want. Spectra with the same sample name are fitted and plotted together.")
    table = pd.DataFrame({
        "Use": True,
        "Sample": [s["group"] for s in spectra],
        "Label": [s["default_label"] for s in spectra],
        "File": [s["file"] + (f" (sweep {s['sub'] + 1})" if s["sub"] else "") for s in spectra],
        "OCV (V)": [round(s["E"], 3) if np.isfinite(s["E"]) else None for s in spectra],
        "Points": [len(s["f"]) for s in spectra],
        "Frequency range": [f"{s['f'].max():.3g} Hz – {s['f'].min():.3g} Hz" for s in spectra],
    })
    edited = st.data_editor(table, hide_index=True,
                            disabled=["File", "OCV (V)", "Points", "Frequency range"],
                            key="spectra_table_" + "|".join(s["id"] for s in spectra))
    for s, (_, row) in zip(spectra, edited.iterrows()):
        s["use"] = bool(row["Use"]); s["sample"] = str(row["Sample"]).strip() or s["group"]
        s["label"] = str(row["Label"]).strip() or s["default_label"]

    run = st.button("Run analysis", type="primary")
    if run:
        st.session_state["eis_run"] = True
    if not st.session_state.get("eis_run"):
        return

    chosen = [s for s in spectra if s["use"]]
    if not chosen:
        st.warning("Tick at least one spectrum.")
        return

    # ---------------- analysis ----------------
    progress = st.progress(0.0, text="Fitting…")
    samples = {}
    prev_by_sample = {}
    for i, s in enumerate(chosen):
        s["kk"] = cached_kk(s["f"], s["Z"])
        if win_mode.startswith("Automatic"):
            s["fmin"] = suggest_fmin(s["f"], s["kk"], thr=kk_thr)
        elif win_mode.startswith("Same"):
            s["fmin"] = float(fixed_fmin)
        else:
            s["fmin"] = float(s["f"].min())
        s["mask"] = s["f"] >= s["fmin"] * (1 - 1e-9)
        prev = prev_by_sample.get(s["sample"]) if consistent else None
        progress.progress(i / len(chosen), text=f"Fitting {s['sample']} – {s['label']} ({i + 1}/{len(chosen)})")
        if s["mask"].sum() < 12:
            s["fit"] = None
            s["fit_error"] = "fewer than 12 points left in the fit window"
        else:
            try:
                s["fit"] = cached_fit(s["f"][s["mask"]], s["Z"][s["mask"]], model, ideal_w, n_starts,
                                      tuple(sorted(prev.items())) if prev else None,
                                      0.15 if consistent else 0.0)
                prev_by_sample[s["sample"]] = s["fit"]["params"]
            except RuntimeError as exc:
                s["fit"] = None
                s["fit_error"] = str(exc)
        samples.setdefault(s["sample"], []).append(s)
    progress.empty()

    for sample, items in samples.items():
        for k, it in enumerate(items):
            it["color"] = ORIGIN_COLORS[k % len(ORIGIN_COLORS)]
            it["marker"] = ORIGIN_MARKERS[k % len(ORIGIN_MARKERS)]
            it["legend"] = it["label"] + (f" ({it['E']:.2f} V)" if add_ocv and np.isfinite(it["E"]) else "")
    rows = [(s["sample"], s["label"], s["E"], s["fmin"], s["fit"]) for s in chosen if s.get("fit")]
    trimmed = []
    for s in chosen:
        if not s.get("fit"):
            st.error(f"{s['sample']} – {s['label']}: fit failed ({s.get('fit_error', 'unknown error')}).")
        f_cut = s["f"][s["mask"]].min() if s["mask"].any() else np.nan
        if win_mode.startswith("Automatic") and np.isfinite(f_cut) and f_cut > 5 * s["f"].min():
            trimmed.append(f"- {s['sample']} – {s['label']}: fitted down to {f_cut:.3g} Hz")
    if trimmed:
        st.warning("Low-frequency points failed the Kramers–Kronig test and were left out of the fit "
                   "(the cell was probably still relaxing; a longer rest before the EIS would help):\n"
                   + "\n".join(trimmed))

    downloads = {}
    tab_fig, tab_table, tab_kk, tab_bode, tab_dl = st.tabs(
        ["Nyquist plot", "Fit results", "Quality check (KK)", "Bode plot", "Downloads"])

    with tab_fig:
        cols = st.columns(min(2, len(samples)))
        for j, (sample, items) in enumerate(samples.items()):
            plot_items = [dict(it, label=it["legend"]) for it in items]
            fig = nyquist_figure(plot_items, title=sample, model=model, lim=man_lim, inset_lim=man_inset,
                                 show_circuit=show_circuit, show_inset=show_inset, legend_pos=legend_pos,
                                 inset_pos=inset_pos)
            png = fig_bytes(fig, "png")
            for fmt in ("png", "pdf", "svg"):
                downloads[f"{_safe(sample)}_Nyquist.{fmt}"] = png if fmt == "png" else fig_bytes(fig, fmt)
            plt.close(fig)
            with cols[j % len(cols)]:
                _show_image(st, png, caption=sample)
                st.download_button("PNG (300 dpi)", png, f"{_safe(sample)}_Nyquist.png", "image/png",
                                   key=f"dl_png_{j}")
        st.caption("Symbols: data. Lines: fit. Open symbols: low-frequency points not used in the fit. "
                   "If a label or the inset covers data, change the positions or the axis range in the sidebar.")

    with tab_table:
        if rows:
            st.dataframe(display_table(rows), hide_index=True)
            st.caption("± values are 1σ fit uncertainties. f = characteristic frequency and C = effective "
                       "capacitance of each arc are in the CSV.")
            csv = fit_table(rows).to_csv(index=False).encode()
            tex = latex_table(rows, model)
            downloads["EIS_fit_results.csv"] = csv
            downloads["EIS_fit_table.tex"] = tex.encode()
            c1, c2 = st.columns(2)
            c1.download_button("Fit table (CSV)", csv, "EIS_fit_results.csv", "text/csv")
            c2.download_button("Fit table (LaTeX)", tex, "EIS_fit_table.tex", "text/plain")
            with st.expander("LaTeX code"):
                st.code(tex, language="latex")
        else:
            st.info("No successful fits.")

    with tab_kk:
        st.markdown("Lin-KK residuals for every spectrum. Residuals below ~1 % mean the data are consistent; "
                    "a rise at low frequency means the cell changed during the measurement. The dotted line "
                    "marks the lowest frequency used in the fit.")
        for j, (sample, items) in enumerate(samples.items()):
            fig = kk_figure([dict(it, label=it["label"]) for it in items], title=sample)
            png = fig_bytes(fig, "png")
            downloads[f"{_safe(sample)}_KK_check.png"] = png
            plt.close(fig)
            _show_image(st, png)

    with tab_bode:
        for j, (sample, items) in enumerate(samples.items()):
            fig = bode_figure([dict(it, label=it["legend"]) for it in items], title=sample)
            png = fig_bytes(fig, "png")
            downloads[f"{_safe(sample)}_Bode.png"] = png
            downloads[f"{_safe(sample)}_Bode.pdf"] = fig_bytes(fig, "pdf")
            plt.close(fig)
            _show_image(st, png)

    with tab_dl:
        for sample, items in samples.items():
            downloads[f"{_safe(sample)}_data_and_fit.csv"] = data_csv(items, sample).to_csv(index=False).encode()
        zbuf = io.BytesIO()
        with zipfile.ZipFile(zbuf, "w", zipfile.ZIP_DEFLATED) as z:
            for name, data in downloads.items():
                z.writestr(name, data)
        st.download_button("Download everything (ZIP)", zbuf.getvalue(), "EIS_results.zip", "application/zip",
                           type="primary")
        st.markdown("The ZIP holds the Nyquist plots (PNG 300 dpi, PDF, SVG), Bode plots, KK checks, the fit "
                    "table (CSV + LaTeX) and, per sample, the data with the fitted curves (CSV, for Origin).")
        with st.expander(f"{len(downloads)} files in the ZIP"):
            st.markdown("\n".join(f"- {name}" for name in sorted(downloads)))

    with st.expander("How the analysis works"):
        st.markdown(
            "- **Kramers–Kronig check (Lin-KK, Schönleber et al. 2014).** Each spectrum is fitted with a chain of "
            "RC elements that is always KK-consistent; the residuals show which points are trustworthy. Points at "
            "the low-frequency end are dropped automatically where the residuals stay above the limit.\n"
            "- **Circuit (Jin et al. 2014).** Re: electrolyte resistance. Rf‖CPE1: high-frequency arc "
            "(surface film / SEI). CPE2‖(Rct–Zw): middle-frequency arc (charge transfer) and the low-frequency "
            "diffusion line. Zw = Aw·(jω)^−α; α = 0.5 is the ideal 45° Warburg line.\n"
            "- **Fit.** Complex non-linear least squares with modulus weighting and many random starts. The arc "
            "with the shorter time constant is always called Rf. Uncertainties are 1σ from the Jacobian; "
            "parameters sitting at a limit (e.g. n = 1) are held fixed and noted.\n"
            "- **Two-electrode half cells.** The sodium counter electrode also contributes to the arcs, so compare "
            "samples measured the same way rather than reading the absolute values too literally.")


def _show_image(st, data, caption=None):
    try:
        st.image(data, caption=caption, width="stretch")          # Streamlit >= 1.46
    except Exception:
        st.image(data, caption=caption, use_container_width=True)  # older Streamlit


def _safe(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name)).strip("_") or "sample"


if __name__ == "__main__":
    main()
