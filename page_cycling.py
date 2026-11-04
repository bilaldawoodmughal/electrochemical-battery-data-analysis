"""
Neware by Bilal - page "Cycling and GITT":
rate tests and GCPL from Neware (BTS .xlsx exports or raw .nda files),
GITT from BioLogic EC-Lab (.mpr files or .mpt text exports).
Started from app.py (the menu); not meant to be run on its own.
"""

import contextlib
import hashlib
import io
import os
import re
import tempfile
import warnings
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

import analyse_gcpl as gcpl
import analyse_gitt as gitt
import analyse_rate_test as rate
import nda_to_excel
import ui

BIOLOGIC = (".mpr", ".mpt")
LENGTH_LABELS = {"thickness": "Electrode thickness",
                 "mass": "Dense-carbon thickness",
                 "fibre": "Fibre radius (SEM)",
                 "particle": "Particle radius (SEM)"}
SHORT_LABELS = {"thickness": "Electrode thickness", "mass": "Dense carbon",
                "fibre": "Fibre, r/2", "particle": "Particle, r/3"}

warnings.filterwarnings("ignore", message="Workbook contains no default style")
WORK_DIR = Path(tempfile.gettempdir()) / "neware_analyser"
ui.apply_style()


# ------------------------------- helpers -------------------------------------
def parse_numbers(text, cast=float):
    """'0.05, 0.1 0.2' -> [0.05, 0.1, 0.2]. Returns [] for an empty field."""
    parts = [p for p in re.split(r"[,\s;]+", text.strip()) if p]
    return [cast(p) for p in parts]


def parse_cycles(text):
    out = []
    for p in [p for p in re.split(r"[,\s;]+", text.strip()) if p]:
        out.append("last" if p.lower() == "last" else int(p))
    return out


@st.cache_data(show_spinner=False)
def inspect_file(path_str, _digest):
    """Test type, current levels, mass and sample name found in the file."""
    path = Path(path_str)
    xl = pd.ExcelFile(path)
    statis = next((s for s in xl.sheet_names if s.lower().startswith("statis")), None)
    detail = next((s for s in xl.sheet_names if s.lower().startswith("detail")), None)
    if detail is None:
        raise ValueError("No 'Detail' sheet found. Export the file from Neware with the detail (record) data.")
    levels = []
    if statis:
        df = pd.read_excel(path, sheet_name=statis)
        status = df[gcpl.find_col(df, "status")].astype(str).str.lower()
        cur_col = next((c for c in df.columns if "termination current" in str(c).lower()),
                       gcpl.find_col(df, "current"))
        cur = np.abs(gcpl.current_to_mA(df[cur_col], str(cur_col)))[~status.str.contains("rest")]
        for c in np.sort(cur[cur > 0].to_numpy()):
            if not levels or c > levels[-1] * 1.10:
                levels.append(c)
    kind = "Rate test" if len(levels) >= 3 else "GCPL"
    mass, name = gcpl.mass_and_name(path)
    return kind, [lv * 1000 for lv in levels], mass, name


@st.cache_data(show_spinner=False, max_entries=6)
def load_gitt(path_str, _digest):
    """Read a BioLogic file once; the table is reused for every analysis of that file."""
    return gitt.load(path_str)


@st.cache_data(show_spinner=False)
def inspect_gitt(path_str, _digest):
    d = load_gitt(path_str, _digest)
    _, segs = gitt.segments(d)
    n = gitt.count_pulses(segs)
    pulses = segs[segs["kind"] == "pulse"]
    mass, name, source = gitt.read_cell_info(path_str)
    return {"n": n, "mass": mass, "name": name, "source": source,
            "current_uA": float(pulses["I_uA"].abs().median()) if len(pulses) else np.nan,
            "pulse_s": float((pulses["t1"] - pulses["t0"]).max()) if len(pulses) else np.nan,
            "hours": float((d["t"].iloc[-1] - d["t"].iloc[0]) / 3600)}


SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def sci(v):
    """2.4e-07 -> '2.4×10⁻⁷'"""
    if v is None or not np.isfinite(v) or v <= 0:
        return "-"
    e = int(np.floor(np.log10(v)))
    return f"{v / 10 ** e:.1f}×10{str(e).translate(SUP)}"


def run_gitt(path, d, mass_mg, sample, geom, gopts, basis, out):
    """Apply the GITT inputs to the analysis module and run one file."""
    gitt.MASS_MG = mass_mg
    gitt.SAMPLE_NAME = sample
    gitt.DIAMETER_MM = geom["diameter"] or None
    gitt.THICKNESS_MM = geom["thickness"] or None
    gitt.SEM_DIAMETER_UM = geom["sem"] or None
    gitt.LENGTH_MODE = gopts["mode"]
    gitt.IR_TIME_S = gopts["ir_time"]
    gitt.REST_DRIFT_LIMIT_MV_H = gopts["drift"]
    gitt.SKIP_FIRST_PULSE = gopts["skip_first"]
    gitt.EXCLUDE_PULSES = {"sodiation": gopts["excl_sod"], "desodiation": gopts["excl_desod"]}
    gitt.M_G_MOL = 72.0
    gitt.V_M_CM3_MOL = 72.0 / gopts["density"]
    gitt.C_RATE_BASIS_MAH_G = basis
    return gitt.run(path, data=d, out=out, verbose=False)


@st.cache_data(show_spinner=False)
def is_impedance_file(path_str, _digest):
    """True if a BioLogic file holds impedance spectra (PEIS/GEIS)."""
    try:
        import EIS_Analyser as eis
        path = Path(path_str)
        return bool(eis.read_eclab(path.name, path.read_bytes())[3])
    except Exception:
        return False


def eis_hint():
    st.info("This is an impedance (EIS) file. It is analysed on the EIS page.")
    st.page_link("page_eis.py", label="Open the EIS page", icon=":material/ssid_chart:")


def zip_folder(folder):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(Path(folder).glob("*")):
            z.write(f, arcname=f.name)
    return buf.getvalue()


def run_analysis(path, kind, mass_mg, sample, opts, out):
    """Apply the sidebar settings to the analysis modules and run one file."""
    for mod in (rate, gcpl):
        mod.SODIATION_CUTOFF_V = opts["sod_cutoff"]
        mod.DESODIATION_CUTOFF_V = opts["desod_cutoff"]
        mod.PLATEAU_THRESHOLD_V = opts["plateau"]
        mod.C_RATE_BASIS_MAH_G = opts["basis"]
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        if kind == "Rate test":
            res = rate.run(path, mass_mg, opts["rates"], opts["start"], out, sample)
        else:
            gcpl.DQDV_CYCLES = opts["dqdv_cycles"]
            gcpl.DQDV_BRANCH = opts["dqdv_branch"]
            gcpl.SELECTED_CYCLES = opts["profile_cycles"]
            gcpl.CE_THRESHOLDS = opts["ce_thresholds"]
            res = gcpl.run(path, mass_mg, opts["start"], out, sample)
    res["log"] = log.getvalue()
    res["kind"] = kind
    return res


# ------------------------------- sidebar -------------------------------------
with st.sidebar:
    st.subheader("Analysis settings")
    with st.expander("Cycle counting and limits", expanded=True):
        start = st.radio(
            "Cycle counting",
            ["desodiation", "sodiation"],
            format_func=lambda s: ("First sodiation not counted (desodiation + next sodiation)"
                                   if s == "desodiation" else "Standard (sodiation + desodiation)"))
        c1, c2 = st.columns(2)
        sod_cutoff = c1.number_input("Lower cut-off (V)", value=0.005, step=0.001, format="%.3f")
        desod_cutoff = c2.number_input("Upper cut-off (V)", value=2.5, step=0.1, format="%.2f")
        plateau = st.number_input("Slope / plateau boundary (V)", value=0.10, step=0.01, format="%.2f")
        basis = st.number_input("1C defined as (mAh/g)", value=250, step=10)
    with st.expander("Rate test"):
        rates_txt = st.text_input("C-rates in the order you ran them", "0.05, 0.1, 0.2, 0.5, 1, 2, 0.05",
                                  help="Leave empty to calculate the C-rates from the current and the 1C definition.")
    with st.expander("GCPL"):
        dqdv_txt = st.text_input("dQ/dV cycles, compared on one graph", "2, 15")
        dqdv_branch = st.selectbox("dQ/dV curve", ["sodiation", "desodiation", "both"],
                                   format_func=lambda b: {"sodiation": "Sodiation (discharge)",
                                                          "desodiation": "Desodiation (charge)",
                                                          "both": "Both"}[b])
        profiles_txt = st.text_input("Cycles in the voltage-profile plot", "1, 2, 3, 5, 10, 20, last")
        ce_txt = st.text_input("CE reference lines (%)", "95, 99")
    with st.expander("GITT (BioLogic)"):
        length_mode = st.selectbox(
            "Diffusion length L", list(LENGTH_LABELS), format_func=LENGTH_LABELS.get,
            help="Electrode thickness: L = thickness. Dense-carbon thickness: L = mass / (density × disc area), "
                 "the form used by Ledwoch et al. 2022. Fibre radius: L = r/2. Particle radius: L = r/3 (spheres). "
                 "D scales with L², so another L moves the whole curve up or down; the shape stays the same. "
                 "Use the same choice for every sample.")
        density = st.number_input("Carbon density for the dense-carbon thickness (g/cm³)", value=1.5, step=0.1,
                                  format="%.2f", help="Ledwoch et al. 2022: V_M = 48 cm³/mol for C6 = 72 g/mol / 1.5 g/cm³.")
        g1, g2 = st.columns(2)
        ir_time = g1.number_input("IR time (s)", value=1.0, min_value=0.0, step=0.5, format="%.1f",
                                  help="The potential change in this first part of each pulse counts as the IR drop.")
        drift_limit = g2.number_input("Rest drift limit (mV/h)", value=3.6, min_value=0.0, step=0.5, format="%.1f",
                                      help="Ledwoch et al. 2022: rest until dE/dt ≤ 0.001 mV/s = 3.6 mV/h.")
        skip_first = st.checkbox("Leave out the first pulse of each direction", value=True,
                                 help="The rest before it starts from an unrelaxed potential.")
        excl_sod_txt = st.text_input("Leave out sodiation pulses", "", placeholder="e.g. 22")
        excl_desod_txt = st.text_input("Leave out desodiation pulses", "", placeholder="e.g. 2, 3")

try:
    opts = {
        "start": start, "sod_cutoff": sod_cutoff, "desod_cutoff": desod_cutoff, "plateau": plateau,
        "basis": basis, "rates": parse_numbers(rates_txt), "dqdv_cycles": parse_cycles(dqdv_txt),
        "dqdv_branch": dqdv_branch, "profile_cycles": parse_cycles(profiles_txt),
        "ce_thresholds": parse_numbers(ce_txt),
    }
    gopts = {
        "mode": length_mode, "density": density, "ir_time": ir_time, "drift": drift_limit,
        "skip_first": skip_first, "excl_sod": parse_numbers(excl_sod_txt, int),
        "excl_desod": parse_numbers(excl_desod_txt, int),
    }
    settings_ok = True
except ValueError:
    st.sidebar.error("Use numbers separated by commas in the text fields, and 'last' for the profile cycles.")
    settings_ok = False
    gopts = {}


# --------------------------------- GITT --------------------------------------
def gitt_tab(up, path, folder, digest):
    """Inputs, analysis and results for one BioLogic GITT file."""
    try:
        with st.spinner("Reading the BioLogic file..."):
            d = load_gitt(str(path), digest)
            info = inspect_gitt(str(path), digest)
    except Exception as e:
        if is_impedance_file(str(path), digest):
            eis_hint()
        else:
            st.error(f"This file could not be read as a BioLogic EC-Lab file: {e}. "
                     "If it is a text export, try the original .mpr file instead.")
        return
    n_sod, n_des = info["n"]["sodiation"], info["n"]["desodiation"]
    if n_sod + n_des < 5:
        if is_impedance_file(str(path), digest):
            eis_hint()
            return
        st.error(f"Only {n_sod + n_des} current pulses with a rest before and after them were found, "
                 "so this does not look like a GITT measurement.")
        return

    with st.container(border=True):
        a, b, c, go = st.columns([1, 1, 1.6, 0.8], vertical_alignment="bottom")
        a.selectbox("Test type", ["GITT (BioLogic)"], key=f"kind_{digest}", disabled=True)
        mass = b.number_input("Active mass (mg)", value=float(info["mass"] or 0.0), min_value=0.0, step=0.01,
                              format="%.3f", key=f"mass_{digest}")
        default_name = f"Sample {info['name']}" if info["name"] else Path(up.name).stem
        name = c.text_input("Sample name on the figures", default_name, key=f"name_{digest}")
        e, f, g, h = st.columns([1, 1, 1.6, 0.8], vertical_alignment="bottom")
        diameter = e.number_input("Electrode diameter (mm)", value=10.0, min_value=0.0, step=0.5, format="%.2f",
                                  key=f"diam_{digest}")
        thickness = f.number_input("Electrode thickness (mm)", value=0.0, min_value=0.0, step=0.01, format="%.3f",
                                   key=f"thick_{digest}", help="Needed only for L = electrode thickness.")
        sem = g.number_input("Fibre or particle diameter from SEM (µm)", value=0.0, min_value=0.0, step=0.1,
                             format="%.3f", key=f"sem_{digest}",
                             help="Needed only for L = fibre or particle radius.")
        geom = {"diameter": diameter, "thickness": thickness, "sem": sem}
        missing = None
        if settings_ok:
            needs = {"thickness": (thickness, "the electrode thickness"), "mass": (diameter, "the electrode diameter"),
                     "fibre": (sem, "the fibre diameter from SEM"), "particle": (sem, "the particle diameter from SEM")}
            value, what = needs[gopts["mode"]]
            missing = None if value > 0 else what
        clicked = go.button("Analyse", key=f"run_{digest}", type="primary",
                            disabled=not settings_ok or mass <= 0 or missing is not None)
        src = f"Mass read from {info['source']}." if info["mass"] else "No mass in the file, enter it above."
        st.caption(f"Detected as GITT: {n_sod} sodiation and {n_des} desodiation pulses of "
                   f"{info['pulse_s']:.0f} s at {info['current_uA']:.2f} µA, {info['hours']:.0f} h in total. {src}"
                   + (f" Diffusion length set to {LENGTH_LABELS[gopts['mode']][0].lower() + LENGTH_LABELS[gopts['mode']][1:]} in the sidebar: "
                      f"enter {missing}, or pick another diffusion length." if missing else ""))

    if mass <= 0 or not settings_ok:
        return
    run_key = f"{digest}|GITT|{mass}|{name}|{sorted(geom.items())}|{sorted(gopts.items())}|{opts['basis']}"
    if clicked:
        with st.spinner("Finding the pulses, calculating D and making the figures..."):
            try:
                out = folder / f"{Path(up.name).stem}_GITT_results"
                res = run_gitt(path, d, mass, f"{name}, {mass:g} mg", geom, gopts, opts["basis"], out)
                st.session_state["results"][digest] = (run_key, res)
            except Exception as e:
                st.error(f"The analysis stopped: {e}")

    stored = st.session_state["results"].get(digest)
    if not stored:
        return
    key_used, res = stored
    if key_used != run_key:
        st.info("Settings changed since the last analysis. Click Analyse to update the figures.")

    if res["problems"]:
        nice = [p.replace("dEs", "ΔEs").replace("dEt", "ΔEτ") for p in res["problems"]]
        st.warning("**Data checks**\n\n" + "\n".join(f"- {p}" for p in nice))
    else:
        st.success("Data checks passed: rests relaxed, no pulses cut short, no jumps in the IR drop.")

    pp = res["per_pulse"]
    used = pp[pp["used"]]
    m = st.columns(4)
    cur = pp["I_uA"].abs().median()
    m[0].metric("Pulse current", f"{cur:.2f} µA (C/{opts['basis'] * mass / 1000 / cur * 1000:.0f})")
    caps = [pp.loc[pp["direction"] == x, "capacity_end_mAh_g"] for x in gitt.DIRECTIONS]
    m[1].metric("Capacity, sodiation / desodiation",
                " / ".join(f"{c.iloc[-1]:.1f}" if len(c) else "-" for c in caps) + " mAh/g")
    for col, direction in zip(m[2:], gitt.DIRECTIONS):
        u = used[used["direction"] == direction]
        if len(u):
            col.metric(f"Median D, {direction} (cm²/s)", sci(u["D_cm2_s"].median()),
                       help=f"Range {sci(u['D_cm2_s'].min())} to {sci(u['D_cm2_s'].max())} cm²/s "
                            f"from {len(u)} pulses.")
        else:
            col.metric(f"Median D, {direction}", "not measured")

    figs = res["figures"]
    with st.container(border=True):
        st.image(str(figs[0]))
        note = (f"{res['L_text']}. D scales with L², so another diffusion length moves the whole curve up or "
                "down without changing its shape. Open symbols: ΔEs is not much smaller than ΔEτ, so D is apparent.")
        if gopts["mode"] == "thickness":
            note += (" With the electrode thickness, D usually comes out far above literature values for hard "
                     "carbon (about 10⁻¹² to 10⁻⁹ cm²/s), because Na⁺ only diffuses in the solid inside each fibre "
                     "or particle.")
        st.caption(note)
    for fig in figs[1:]:
        with st.container(border=True):
            st.image(str(fig))

    left, right = st.columns([1.25, 1])
    with left:
        st.subheader("Summary")
        st.dataframe(res["summary"].astype(str), hide_index=True, height=34 * (len(res["summary"]) + 1) + 4)
    with right:
        st.subheader("Median D for each diffusion length")
        rows = []
        for key, L in res["lengths"].items():
            row = {"Diffusion length": SHORT_LABELS[key], "L (µm)": f"{L * 1e4:.3g}"}
            for direction in gitt.DIRECTIONS:
                u = used[used["direction"] == direction]
                row[f"D {direction} (cm²/s)"] = sci(u[f"D_{key}_cm2_s"].median()) if len(u) else "-"
            rows.append(row)
        st.dataframe(pd.DataFrame(rows), hide_index=True)
        st.caption("Lengths that need a value you have not entered are not shown.")

    with st.expander("Every pulse"):
        show = pp[["direction", "pulse", "used", "excluded_because", "note", "tau_s", "capacity_end_mAh_g",
                   "E_mid_V", "IR_drop_mV", "dEs_mV", "dEt_mV", "D_cm2_s", "rest_end_drift_mV_h"]]
        st.dataframe(show, hide_index=True, column_config={
            "D_cm2_s": st.column_config.NumberColumn("D (cm²/s)", format="%.2e"),
            "tau_s": st.column_config.NumberColumn("pulse (s)", format="%.0f"),
            "capacity_end_mAh_g": st.column_config.NumberColumn("capacity (mAh/g)", format="%.1f"),
            "E_mid_V": st.column_config.NumberColumn("potential (V)", format="%.3f"),
            "IR_drop_mV": st.column_config.NumberColumn("IR drop (mV)", format="%.1f"),
            "dEs_mV": st.column_config.NumberColumn("ΔEs (mV)", format="%.2f"),
            "dEt_mV": st.column_config.NumberColumn("ΔEτ (mV)", format="%.2f"),
            "rest_end_drift_mV_h": st.column_config.NumberColumn("rest drift (mV/h)", format="%.1f"),
        })

    d1, d2, _ = st.columns([1, 1, 2])
    d1.download_button("Download figures and Excel (.zip)", zip_folder(res["out"]),
                       file_name=f"{res['out'].name}.zip", mime="application/zip", key=f"zip_{digest}")
    d2.download_button("Download Excel", Path(res["excel"]).read_bytes(),
                       file_name=Path(res["excel"]).name, key=f"xlsx_{digest}",
                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")


# -------------------------------- main ---------------------------------------
ui.footer()
ui.header("Cycling data analysis",
          "Rate tests and galvanostatic cycling from Neware (BTS exports or raw .nda files), and GITT from "
          "BioLogic EC-Lab (.mpr). Drop one or more files to get the figures and an Excel summary.")

uploads = st.file_uploader("Neware files (.xlsx export or raw .nda) or BioLogic GITT files (.mpr, .mpt)",
                           type=["xlsx", "nda", "mpr", "mpt"], accept_multiple_files=True)
if "results" not in st.session_state:
    st.session_state["results"] = {}

if not uploads:
    with st.container(border=True):
        st.html('<div class="nb-empty"><h3>No files loaded</h3>'
                '<p>Drop a Neware .xlsx export (with the Detail sheet) or the raw .nda file in the box above; '
                '.nda files are converted to the Neware Excel layout first. '
                'For GITT, drop the BioLogic .mpr file to get the diffusion coefficient and a check of the data. '
                'The test type, active mass and sample name are read from the file, and you can change them '
                'before analysing.</p></div>')
    st.stop()

tabs = st.tabs([u.name for u in uploads])
for up, tab in zip(uploads, tabs):
    with tab:
        data = up.getvalue()
        digest = hashlib.md5(data).hexdigest()[:10]
        folder = WORK_DIR / digest
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / up.name
        if not path.exists():
            path.write_bytes(data)
        if path.suffix.lower() in BIOLOGIC:               # BioLogic GITT file
            gitt_tab(up, path, folder, digest)
            continue
        converted = None
        if path.suffix.lower() == ".nda":                 # raw Neware file: convert to the Excel layout
            xlsx = path.with_suffix(".xlsx")
            try:
                if not xlsx.exists():
                    with st.spinner("Converting the .nda file to the Neware Excel layout..."):
                        nda_to_excel.convert(path, xlsx)
                converted, path = xlsx, xlsx
            except Exception as e:
                st.error(f"The .nda file could not be converted: {e}")
                continue

        try:
            kind_auto, levels_uA, mass_auto, name_auto = inspect_file(str(path), digest)
        except Exception as e:
            st.error(f"This file could not be read as a Neware export: {e}")
            continue

        with st.container(border=True):
            a, b, c, d = st.columns([1, 1, 1.6, 0.8], vertical_alignment="bottom")
            kind = a.selectbox("Test type", ["Rate test", "GCPL"], index=0 if kind_auto == "Rate test" else 1,
                               key=f"kind_{digest}")
            mass = b.number_input("Active mass (mg)", value=float(mass_auto or 0.0), min_value=0.0, step=0.01,
                                  format="%.3f", key=f"mass_{digest}")
            default_name = f"Sample {name_auto}" if name_auto else Path(up.name).stem
            name = c.text_input("Sample name on the figures", default_name, key=f"name_{digest}")
            clicked = d.button("Analyse", key=f"run_{digest}", type="primary",
                               disabled=not settings_ok or mass <= 0)
            levels_txt = ", ".join(f"{lv:.2f}" for lv in levels_uA) if levels_uA else "none found"
            st.caption(f"Detected as {kind_auto}, current levels {levels_txt} µA. "
                       f"{'Mass read from the file.' if mass_auto else 'No mass in the file, enter it above.'}"
                       + (" Converted from .nda." if converted else ""))
            if converted:
                st.download_button("Download the converted Neware Excel file", converted.read_bytes(),
                                   file_name=converted.name, key=f"conv_{digest}",
                                   mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

        if mass <= 0:
            continue
        run_key = f"{digest}|{kind}|{mass}|{name}|{sorted(opts.items()) if settings_ok else ''}"
        if clicked:
            with st.spinner("Reading the Detail sheet and making the figures..."):
                try:
                    out = folder / f"{Path(up.name).stem}_{'rate' if kind == 'Rate test' else 'GCPL'}_results"
                    res = run_analysis(path, kind, mass, f"{name}, {mass:g} mg", opts, out)
                    st.session_state["results"][digest] = (run_key, res)
                except Exception as e:
                    st.error(f"The analysis stopped: {e}. Check the test type, the mass and the cut-offs.")

        stored = st.session_state["results"].get(digest)
        if not stored:
            continue
        key_used, res = stored
        if key_used != run_key:
            st.info("Settings changed since the last analysis. Click Analyse to update the figures.")

        for line in res["log"].splitlines():
            if line.startswith("WARNING"):
                msg = line.replace("WARNING: ", "")
                st.warning(msg[:1].upper() + msg[1:])

        m = st.columns(4)
        m[0].metric("Initial coulombic efficiency", f"{res['ice']:.1f} %")
        if res["kind"] == "Rate test":
            pr = res["per_rate"]
            m[1].metric(f"Reversible capacity at {pr['C-rate'].iloc[0]:g}C",
                        f"{pr['desodiation mean (mAh/g)'].iloc[0]:.1f} mAh/g")
            m[2].metric(f"Reversible capacity at {pr['C-rate'].max():g}C",
                        f"{pr.loc[pr['C-rate'].idxmax(), 'desodiation mean (mAh/g)']:.1f} mAh/g")
            if len(pr) > 1 and pr["C-rate"].iloc[-1] == pr["C-rate"].iloc[0]:
                m[3].metric("Recovery in the last block", f"{pr['retention vs block 1 (%)'].iloc[-1]:.0f} %")
            else:
                m[3].metric("Recovery block", "not finished yet")
            table = res["summary"]
        else:
            summ = dict(zip(res["summary"]["Quantity"], res["summary"]["Value"]))
            stab = next((v for k, v in summ.items() if k.startswith("Stabilised CE: mean")), np.nan)
            last = next((v for k, v in summ.items() if k.startswith("Desodiation last")), np.nan)
            m[1].metric("Stabilised coulombic efficiency", f"{stab:.1f} %")
            m[2].metric("Last reversible capacity", f"{last:.1f} mAh/g")
            m[3].metric("Current", f"{res['current_uA']:.2f} µA (C/{1 / res['c_rate']:.1f})")
            table = res["summary"].astype({"Value": str})
        if res["notes"]:
            st.caption("Not used because the cut-off was not reached: " + ", ".join(res["notes"]) + ".")

        figs = res["figures"]
        for row in range(0, len(figs), 2):                 # 2 x 2 grid, rows aligned
            cols = st.columns(2)
            for col, fig in zip(cols, figs[row:row + 2]):
                with col.container(border=True):
                    st.image(str(fig))

        st.subheader("Summary")
        st.dataframe(table, hide_index=True)
        d1, d2, _ = st.columns([1, 1, 2])
        d1.download_button("Download figures and Excel (.zip)", zip_folder(res["out"]),
                           file_name=f"{res['out'].name}.zip", mime="application/zip", key=f"zip_{digest}")
        d2.download_button("Download Excel", Path(res["excel"]).read_bytes(),
                           file_name=Path(res["excel"]).name, key=f"xlsx_{digest}",
                           mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        with st.expander("Analysis log"):
            st.code(res["log"] or "(empty)")
