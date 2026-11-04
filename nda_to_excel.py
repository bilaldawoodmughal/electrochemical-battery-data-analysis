"""
Convert Neware .nda files (BTS 7.x / 8.x, nda version 29) into the same Excel layout that
Neware BTS exports: sheets Info, Cycle_<dev>_<unit>_<chnl>, Statis_..., Detail_...

The result opens in Excel and works with analyse_rate_test.py, analyse_gcpl.py,
compare_samples.py and the "Neware by Bilal" app.

Usage
-----
    python nda_to_excel.py "rate test_Bilal_sample F_0.66mg-240040-2-1-0.nda"
    python nda_to_excel.py *.nda                      # several files
    python nda_to_excel.py "D:\\Neware data files"      # every .nda in a folder

Each .xlsx is written next to its .nda (an existing .xlsx with the same name is not overwritten
unless you add --overwrite).

Notes
-----
* Reads the binary records directly (record layout as documented by the open-source NewareNDA
  project). Damaged records, e.g. a half-written last record of a test that is still running,
  are skipped.
* nda version 130 (BTS 9.x) is handed to the NewareNDA package if it is installed
  (pip install NewareNDA); the "Jump" column is then left empty.

Requirements: pandas, numpy, openpyxl
"""

import argparse
import mmap
import re
import struct
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

# current-range code -> multiplier to get mA (from the NewareNDA project)
RANGE_MULTIPLIER = {
    -100000000: 10, -200000: 1e-2, -100000: 1e-2, -60000: 1e-2, -50000: 1e-2, -40000: 1e-2,
    -30000: 1e-2, -20000: 1e-2, -12000: 1e-2, -10000: 1e-2, -6000: 1e-2, -5000: 1e-2,
    -3000: 1e-2, -2000: 1e-2, -1000: 1e-2, -500: 1e-3, -100: 1e-3, -50: 1e-4, -25: 1e-4,
    -20: 1e-4, -10: 1e-4, -5: 1e-5, -2: 1e-5, -1: 1e-5, 0: 0, 1: 1e-4, 2: 1e-4, 5: 1e-4,
    10: 1e-3, 20: 1e-3, 25: 1e-3, 50: 1e-3, 100: 1e-2, 200: 1e-2, 250: 1e-2, 500: 1e-2,
    1000: 1e-1, 6000: 1e-1, 10000: 1e-1, 12000: 1e-1, 20000: 1e-1, 30000: 1e-1, 40000: 1e-1,
    50000: 1e-1, 60000: 1e-1, 100000: 1e-1, 200000: 1e-1,
}
STATUS = {1: "CC Chg", 2: "CC DChg", 3: "CV Chg", 4: "Rest", 5: "Cycle", 7: "CCCV Chg",
          8: "CP DChg", 9: "CP Chg", 10: "CR DChg", 13: "Pause", 16: "Pulse", 17: "SIM",
          19: "CV DChg", 20: "CCCV DChg", 21: "Control", 22: "OCV", 26: "CPCV DChg", 27: "CPCV Chg"}


# ------------------------------- reading -------------------------------------
def _header_info(mm):
    """Text fields from the file header (process name, remark, P/N, operator)."""
    head = mm[:min(len(mm), 8192)]
    texts = [(m.start(), m.group().decode("latin-1").strip())
             for m in re.finditer(rb"[\x20-\x7e]{4,}", head)]
    info = {"process": "", "remark": "", "pn": "", "start_text": ""}
    for pos, t in texts:
        if t.lower().endswith(".xml") and not info["process"]:
            info["process"] = t
        if re.fullmatch(r"\d{2}\.\d{2}\.\d{4} \d{2}:\d{2}:\d{2}", t) and not info["start_text"]:
            info["start_text"] = t
        elif re.fullmatch(r"\d{2}\.\d{2}\.\d{4}", t) and not info["pn"]:
            info["pn"] = t
    # the remark is the free text right before the process name
    if info["process"]:
        before = [t for pos, t in texts if pos < next(p for p, x in texts if x == info["process"])]
        cand = [t for t in before if not re.fullmatch(r"[\d.: ]+", t) and not t.startswith("NEWARE")]
        if cand:
            info["remark"] = cand[-1]
    return info


def _read_v29(mm):
    record_len, ident = 86, b"\x00\x00\x00\x00\x55\x00"
    size = mm.size()
    start = mm.find(ident)
    while start != -1 and start + 4 + record_len < size and (
            mm[start + 4 + record_len] != 0x55 or mm[start + 4 + 12] == 0):
        start = mm.find(ident, start + 4)
    if start == -1:
        raise ValueError("no data records found in the file")
    rows, skipped, pos = [], 0, start + 4
    while pos + record_len <= size:
        b = mm[pos:pos + record_len]
        pos += record_len
        if b[0:2] != b"\x55\x00" or b[82:86] != b"\x00\x00\x00\x00":
            continue
        index, cycle, step_index, status = struct.unpack("<IIHB", b[2:13])
        jump = b[13]
        [step_raw] = struct.unpack("<I", b[10:14])
        t_ms, volt, cur = struct.unpack("<Qii", b[14:30])
        cc, dc, ce, de = struct.unpack("<qqqq", b[38:70])
        y, mo, d, h, mi, s = struct.unpack("<HBBBBB", b[70:77])
        [rng] = struct.unpack("<i", b[78:82])
        if index == 0 or status == 0:
            continue
        mult = RANGE_MULTIPLIER.get(rng)
        try:
            ts = datetime(y, mo, d, h, mi, s)
        except ValueError:
            ts = None
        if mult is None or ts is None or status not in STATUS:
            skipped += 1
            continue
        rows.append((index, cycle + 1, step_raw, step_index, jump, status, t_ms, volt / 10000,
                     cur * mult, (cc if status in (1, 3, 7, 9, 27) else dc) * mult / 3600,
                     cc * mult / 3600, dc * mult / 3600, ce * mult / 3600, de * mult / 3600, ts))
    cols = ["index", "cycle", "step_raw", "step_index", "jump", "status", "t_ms", "V", "I_mA",
            "Q_mAh", "Qc_mAh", "Qd_mAh", "Ec_mWh", "Ed_mWh", "time"]
    return pd.DataFrame(rows, columns=cols), skipped


def _read_with_newarenda(path):
    import NewareNDA                                     # BTS 9 files
    df = NewareNDA.read(str(path))
    st_inv = {v.replace("_", " "): k for k, v in STATUS.items()}
    status_txt = df["Status"].astype(str).str.replace("_", " ")
    status = status_txt.map(st_inv).fillna(0).astype(int)
    q = np.where(status_txt.str.contains("DChg"), df["Discharge_Capacity(mAh)"], df["Charge_Capacity(mAh)"])
    out = pd.DataFrame({
        "index": df["Index"], "cycle": df["Cycle"], "step_raw": df["Step"], "step_index": df["Step_Index"],
        "jump": np.nan, "status": status, "t_ms": (df["Time"] * 1000).round(), "V": df["Voltage"],
        "I_mA": df["Current(mA)"], "Q_mAh": q, "Qc_mAh": df["Charge_Capacity(mAh)"],
        "Qd_mAh": df["Discharge_Capacity(mAh)"], "Ec_mWh": df["Charge_Energy(mWh)"],
        "Ed_mWh": df["Discharge_Energy(mWh)"], "time": pd.to_datetime(df["Timestamp"]),
    })
    return out, 0


def read_nda(path):
    path = Path(path)
    with open(path, "rb") as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        if mm[:6] != b"NEWARE":
            raise ValueError(f"{path.name} is not a Neware .nda file")
        version = mm[14]
        info = _header_info(mm)
        if version == 29:
            df, skipped = _read_v29(mm)
        elif version == 130:
            try:
                df, skipped = _read_with_newarenda(path)
            except ImportError:
                raise ValueError("this .nda is from BTS 9 (version 130): install the reader with "
                                 "'pip install NewareNDA' and run again")
        else:
            raise ValueError(f"nda version {version} is not supported")
    df = df.drop_duplicates("index").sort_values("index").reset_index(drop=True)
    df["steps"] = (df["step_raw"] != df["step_raw"].shift()).cumsum()     # sequential step number
    info.update({"version": version, "skipped": skipped})
    return df, info


# ------------------------------- building ------------------------------------
def fmt_duration(ms):
    ms = int(round(ms))
    h, rem = divmod(ms, 3_600_000)
    m, rem = divmod(rem, 60_000)
    s, ms = divmod(rem, 1000)
    return f"{h}:{m:02d}:{s:02d}.{ms:03d}"


def channel_ids(path):
    m = re.search(r"-(\d+)-(\d+)-(\d+)-\d+$", Path(path).stem)
    if not m:
        return 0, 0, 1
    return int(m.group(1)) % 1000, int(m.group(2)), int(m.group(3))


def build_sheets(df, info, path):
    dev, unit, chnl = channel_ids(path)
    status_txt = df["status"].map(STATUS)
    detail = pd.DataFrame({
        "Record number": df["index"],
        "status": status_txt,
        "Jump": df["jump"],
        "Cycle": df["cycle"],
        "Steps": df["steps"],
        "Current(μA)": (df["I_mA"] * 1000).round(2),
        "Voltage(V)": df["V"].round(4),
        "Capacity(mAh)": df["Q_mAh"].round(6),
        "Energy(Wh)": np.where(status_txt.str.contains("DChg"), df["Ed_mWh"], df["Ec_mWh"]) / 1000,
        "Relative Time(h:min:s.ms)": df["t_ms"].map(fmt_duration),
        "Real Time(h:min:s.ms)": df["time"].dt.strftime("%Y-%m-%d %H:%M:%S"),
    })

    g = df.groupby("steps", sort=True)
    first, last = g.first(), g.last()
    st = last["status"].map(STATUS)
    is_d = st.str.contains("DChg")
    is_c = st.str.contains("Chg") & ~is_d
    qc = np.where(is_c, last["Qc_mAh"], 0.0)
    qd = np.where(is_d, last["Qd_mAh"], 0.0)
    ec = np.where(is_c, last["Ec_mWh"], 0.0) / 1000
    ed = np.where(is_d, last["Ed_mWh"], 0.0) / 1000
    statis = pd.DataFrame({
        "Chnl": chnl,
        "Cycle": first["cycle"].values,
        "Steps": first.index.values,
        "Original step": first["step_index"].values,
        "status": st.values,
        "Start Volt(V)": first["V"].round(4).values,
        "End Volt(V)": last["V"].round(4).values,
        "Starting current(μA)": (first["I_mA"] * 1000).round(2).values,
        "Termination current(μA)": (last["I_mA"] * 1000).round(2).values,
        "Capacity(mAh)": np.round(qc + qd, 6),
        "Continuous Time(h:min:s.ms)": last["t_ms"].map(fmt_duration).values,
        "Relative Time(h:min:s.ms)": last["t_ms"].map(fmt_duration).values,
        "Real Time": last["time"].dt.strftime("%Y-%m-%d %H:%M:%S").values,
        "Net discharge capacity(mAh)": np.round(qd - qc, 6),
        "Chg Cap(mAh)": np.round(qc, 6),
        "DChg Cap(mAh)": np.round(qd, 6),
        "Net discharge energy(Wh)": np.round(ed - ec, 6),
        "Chg Eng(Wh)": np.round(ec, 6),
        "DChg Eng(Wh)": np.round(ed, 6),
    })

    cyc = statis.groupby("Cycle", sort=True)[["Chg Cap(mAh)", "DChg Cap(mAh)"]].sum().reset_index()
    first_d = cyc["DChg Cap(mAh)"].iloc[0] if len(cyc) else np.nan
    cycle = pd.DataFrame({
        "Chnl": chnl,
        "Cycle number": cyc["Cycle"],
        "Chg Cap(mAh)": cyc["Chg Cap(mAh)"].round(6),
        "DChg Cap(mAh)": cyc["DChg Cap(mAh)"].round(6),
        "Discharge capacity decay rate(%)": (100 * cyc["DChg Cap(mAh)"] / first_d).round(6) if first_d else np.nan,
    })

    xlsx_name = Path(path).with_suffix(".xlsx").name
    info_rows = [
        ["Test information"],
        [xlsx_name],
        [datetime.now().strftime("%Y-%m-%d %H:%M:%S")],
        [f"Converted from {Path(path).name} with nda_to_excel.py"],
        ["Circular statistics", "Default Step"],
        ["Process name:", info.get("process", "")],
        ["device", "Unit", "Chnl", "P/N", "Step File", "Starting time", "End Time", "Sorting files", "Class", "Remark"],
        [dev, unit, chnl, info.get("pn", ""), "-",
         df["time"].iloc[0].strftime("%Y-%m-%d %H:%M:%S"), df["time"].iloc[-1].strftime("%Y-%m-%d %H:%M:%S"),
         None, 0, info.get("remark", "")],
    ]
    width = max(len(r) for r in info_rows)
    info_df = pd.DataFrame([r + [None] * (width - len(r)) for r in info_rows])
    suffix = f"{dev}_{unit}_{chnl}"
    return {"Info": info_df, f"Cycle_{suffix}": cycle, f"Statis_{suffix}": statis, f"Detail_{suffix}": detail}


def convert(nda_path, out_path=None, overwrite=False):
    nda_path = Path(nda_path)
    out = Path(out_path) if out_path else nda_path.with_suffix(".xlsx")
    if out.exists() and not overwrite:
        raise FileExistsError(f"{out.name} already exists (use --overwrite to replace it)")
    df, info = read_nda(nda_path)
    sheets = build_sheets(df, info, nda_path)
    with pd.ExcelWriter(out, engine="openpyxl") as xw:
        for name, frame in sheets.items():
            frame.to_excel(xw, sheet_name=name, index=False, header=(name != "Info"))
        for ws in xw.book.worksheets:
            for col in ws.columns:
                width = max(len(str(col[0].value or "")), 10)
                ws.column_dimensions[col[0].column_letter].width = min(width + 2, 30)
    return out, {"records": len(df), "steps": int(df["steps"].max()), "cycles": int(df["cycle"].max()),
                 "skipped": info["skipped"], "remark": info.get("remark", ""),
                 "end": df["time"].iloc[-1]}


def main():
    p = argparse.ArgumentParser(description="Convert Neware .nda files to the Neware Excel export layout")
    p.add_argument("paths", nargs="+", help=".nda files or folders")
    p.add_argument("--overwrite", action="store_true", help="replace existing .xlsx files")
    args = p.parse_args()
    files = []
    for x in map(Path, args.paths):
        files += sorted(x.glob("*.nda")) if x.is_dir() else [x]
    for f in files:
        try:
            out, s = convert(f, overwrite=args.overwrite)
            print(f"{f.name} -> {out.name}: {s['records']} records, {s['steps']} steps, {s['cycles']} cycles, "
                  f"last record {s['end']:%Y-%m-%d %H:%M:%S}"
                  + (f", {s['skipped']} damaged record(s) skipped" if s["skipped"] else ""))
        except Exception as e:
            print(f"{f.name}: not converted - {e}")


if __name__ == "__main__":
    main()
