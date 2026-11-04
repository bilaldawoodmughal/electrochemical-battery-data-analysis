"""
Neware by Bilal - analysis of sodium-ion half-cells.

Pages (menu at the top of the sidebar):
    Cycling and GITT : rate tests and GCPL from Neware (BTS .xlsx exports or raw .nda files),
                       GITT from BioLogic EC-Lab (.mpr files or .mpt text exports)
    EIS              : impedance from BioLogic EC-Lab (PEIS .mpr / .mpt): Kramers-Kronig check,
                       equivalent-circuit fit (Jin et al. 2014), Nyquist and Bode plots

Start it with:
    python -m streamlit run app.py
(or double-click "Start Neware by Bilal.bat", or press Run in VS Code)

Folder layout (keep it together):
    app.py              this menu
    page_cycling.py     page "Cycling and GITT"
    page_eis.py         page "EIS" (runs EIS_Analyser.py)
    ui.py               shared look (style, header, footer)
    analyse_rate_test.py, analyse_gcpl.py, analyse_gitt.py, nda_to_excel.py, EIS_Analyser.py
    galvani_biologic.py reader for BioLogic files (from the galvani project, GPL-3.0-or-later)
    assets/   logo and icon
    static/   Barlow font files (served locally)
    licenses/ licences of the fonts and of the BioLogic reader
    .streamlit/config.toml   dark theme
Requirements: numpy, pandas, matplotlib, openpyxl, scipy, streamlit
"""
import os
import sys
from pathlib import Path

import streamlit as st

APP_DIR = Path(__file__).resolve().parent
ASSETS = APP_DIR / "assets"

# Started with the Run button or "python app.py"? Then start it properly as a Streamlit app.
if __name__ == "__main__":
    from streamlit import runtime
    if not runtime.exists():
        from streamlit.web import cli as stcli
        os.chdir(APP_DIR)                      # so the theme in .streamlit/config.toml is found
        sys.argv = ["streamlit", "run", str(APP_DIR / "app.py")]
        sys.exit(stcli.main())

st.set_page_config(page_title="Neware by Bilal", page_icon=str(ASSETS / "icon.png"), layout="wide")
st.logo(str(ASSETS / "logo.svg"), size="large", icon_image=str(ASSETS / "icon.svg"))

page = st.navigation([
    st.Page(str(APP_DIR / "page_cycling.py"), title="Cycling and GITT", icon=":material/battery_charging_full:",
            default=True),
    st.Page(str(APP_DIR / "page_eis.py"), title="EIS", icon=":material/ssid_chart:", url_path="eis"),
])
page.run()
