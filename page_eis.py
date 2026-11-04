"""
Neware by Bilal - page "EIS": impedance analysis of BioLogic EC-Lab files (PEIS .mpr / .mpt).
The analysis itself is in EIS_Analyser.py (it also runs on its own: streamlit run EIS_Analyser.py).
Started from app.py (the menu); not meant to be run on its own.
"""
import sys

import matplotlib.pyplot as plt
import streamlit as st

import ui

ui.apply_style()
ui.header("Impedance analysis",
          "EIS from BioLogic EC-Lab (.mpr or .mpt). Each spectrum gets a Kramers–Kronig check, an "
          "equivalent-circuit fit (Jin et al., Electrochim. Acta 141, 2014) and Nyquist and Bode plots ready "
          "for the report.")

try:
    import EIS_Analyser as eis
except ImportError as exc:
    st.error(f"The EIS page needs a package that is not installed: {exc.name}. Close the app, run this in a "
             f"terminal and start the app again:")
    st.code(f'"{sys.executable}" -m pip install {exc.name}', language="bat")
    st.stop()

with plt.rc_context():          # the EIS figure style must not change the figures on the other page
    eis.main(embedded=True)
ui.footer()
