"""Shared look of Neware by Bilal: page style, header with the cycling trace, sidebar footer."""
import streamlit as st

VERSION = "1.2"
MUTED, SURFACE, EDGE, AMBER = "#8B97A8", "#1A2230", "#2B3646", "#F2B233"


def apply_style():
    """CSS shared by all pages (call once at the top of every page)."""
    st.markdown(f"""
    <style>
    .block-container {{ padding-top: 2.4rem; max-width: 1400px; }}
    .nb-head h1 {{ font-family: 'Barlow Semi Condensed', Barlow, sans-serif; font-weight: 600;
                  font-size: 2.2rem; line-height: 1.1; margin: 0; padding: 0; }}
    .nb-head p {{ color: {MUTED}; margin: .45rem 0 0 0; max-width: 64ch; font-size: 1rem; }}
    .nb-trace {{ display: block; width: 100%; height: 26px; margin: 1.1rem 0 1.4rem 0; }}
    [data-testid="stMetric"] {{ background: {SURFACE}; border: 1px solid {EDGE}; border-radius: .4rem;
                               padding: .7rem 1rem .8rem 1rem; }}
    [data-testid="stMetricLabel"] p {{ color: {MUTED}; font-size: .85rem; }}
    [data-testid="stMetricValue"] {{ font-family: 'Barlow Semi Condensed', Barlow, sans-serif; font-weight: 600;
                                    font-size: 1.85rem; font-variant-numeric: tabular-nums; }}
    [data-testid="stFileUploaderDropzone"] {{ border: 1px dashed #3A475A; background: #161E2A; }}
    [data-testid="stSidebarHeader"] img, [data-testid="stLogo"] {{ height: 2.4rem !important; max-width: 100%; }}
    [data-testid="stBaseButton-primary"], [data-testid="stBaseButton-primary"] p {{ color: #121821 !important;
                                                                                  font-weight: 600; }}
    [data-testid="stBaseButton-primary"]:disabled, [data-testid="stBaseButton-primary"]:disabled p {{
        color: #6B778A !important; }}
    .nb-foot {{ color: {MUTED}; font-size: .8rem; line-height: 1.45; margin-top: 1.6rem; }}
    .nb-empty h3 {{ font-family: 'Barlow Semi Condensed', Barlow, sans-serif; font-weight: 600; margin: 0 0 .4rem 0; }}
    .nb-empty p {{ color: {MUTED}; margin: 0; max-width: 70ch; }}
    </style>
    """, unsafe_allow_html=True)


# decorative divider: a run of galvanostatic cycles, fading in from the left
TRACE_SVG = '<img class="nb-trace" alt="" src="data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAxMjAwIDI2IiBwcmVzZXJ2ZUFzcGVjdFJhdGlvPSJub25lIj48ZGVmcz48bGluZWFyR3JhZGllbnQgaWQ9ImYiIHgxPSIwIiB4Mj0iMSIgeTE9IjAiIHkyPSIwIj48c3RvcCBvZmZzZXQ9IjAiIHN0b3AtY29sb3I9IiNGMkIyMzMiIHN0b3Atb3BhY2l0eT0iMC4wOCIvPjxzdG9wIG9mZnNldD0iMC44IiBzdG9wLWNvbG9yPSIjRjJCMjMzIiBzdG9wLW9wYWNpdHk9IjAuNTUiLz48c3RvcCBvZmZzZXQ9IjEiIHN0b3AtY29sb3I9IiNGMkIyMzMiIHN0b3Atb3BhY2l0eT0iMSIvPjwvbGluZWFyR3JhZGllbnQ+PC9kZWZzPjxwYXRoIGQ9Ik0wLjAgMjMgQzguNiAyMiAxNS40IDE0IDI1LjcgMyBDMzAuOSAxMiA0Mi45IDIxIDY0LjMgMjIuNSBMODUuNyAyMyBMODUuNyAyMyBDOTQuMyAyMiAxMDEuMSAxNCAxMTEuNCAzIEMxMTYuNiAxMiAxMjguNiAyMSAxNTAuMCAyMi41IEwxNzEuNCAyMyBMMTcxLjQgMjMgQzE4MC4wIDIyIDE4Ni45IDE0IDE5Ny4xIDMgQzIwMi4zIDEyIDIxNC4zIDIxIDIzNS43IDIyLjUgTDI1Ny4xIDIzIEwyNTcuMSAyMyBDMjY1LjcgMjIgMjcyLjYgMTQgMjgyLjkgMyBDMjg4LjAgMTIgMzAwLjAgMjEgMzIxLjQgMjIuNSBMMzQyLjkgMjMgTDM0Mi45IDIzIEMzNTEuNCAyMiAzNTguMyAxNCAzNjguNiAzIEMzNzMuNyAxMiAzODUuNyAyMSA0MDcuMSAyMi41IEw0MjguNiAyMyBMNDI4LjYgMjMgQzQzNy4xIDIyIDQ0NC4wIDE0IDQ1NC4zIDMgQzQ1OS40IDEyIDQ3MS40IDIxIDQ5Mi45IDIyLjUgTDUxNC4zIDIzIEw1MTQuMyAyMyBDNTIyLjkgMjIgNTI5LjcgMTQgNTQwLjAgMyBDNTQ1LjEgMTIgNTU3LjEgMjEgNTc4LjYgMjIuNSBMNjAwLjAgMjMgTDYwMC4wIDIzIEM2MDguNiAyMiA2MTUuNCAxNCA2MjUuNyAzIEM2MzAuOSAxMiA2NDIuOSAyMSA2NjQuMyAyMi41IEw2ODUuNyAyMyBMNjg1LjcgMjMgQzY5NC4zIDIyIDcwMS4xIDE0IDcxMS40IDMgQzcxNi42IDEyIDcyOC42IDIxIDc1MC4wIDIyLjUgTDc3MS40IDIzIEw3NzEuNCAyMyBDNzgwLjAgMjIgNzg2LjkgMTQgNzk3LjEgMyBDODAyLjMgMTIgODE0LjMgMjEgODM1LjcgMjIuNSBMODU3LjEgMjMgTDg1Ny4xIDIzIEM4NjUuNyAyMiA4NzIuNiAxNCA4ODIuOSAzIEM4ODguMCAxMiA5MDAuMCAyMSA5MjEuNCAyMi41IEw5NDIuOSAyMyBMOTQyLjkgMjMgQzk1MS40IDIyIDk1OC4zIDE0IDk2OC42IDMgQzk3My43IDEyIDk4NS43IDIxIDEwMDcuMSAyMi41IEwxMDI4LjYgMjMgTDEwMjguNiAyMyBDMTAzNy4xIDIyIDEwNDQuMCAxNCAxMDU0LjMgMyBDMTA1OS40IDEyIDEwNzEuNCAyMSAxMDkyLjkgMjIuNSBMMTExNC4zIDIzIEwxMTE0LjMgMjMgQzExMjIuOSAyMiAxMTI5LjcgMTQgMTE0MC4wIDMgQzExNDUuMSAxMiAxMTU3LjEgMjEgMTE3OC42IDIyLjUgTDEyMDAuMCAyMyIgZmlsbD0ibm9uZSIgc3Ryb2tlPSJ1cmwoI2YpIiBzdHJva2Utd2lkdGg9IjEuNiIgdmVjdG9yLWVmZmVjdD0ibm9uLXNjYWxpbmctc3Ryb2tlIiBzdHJva2UtbGluZWpvaW49InJvdW5kIi8+PC9zdmc+">'


def header(title, text):
    st.html(f'<div class="nb-head"><h1>{title}</h1><p>{text}</p></div>' + TRACE_SVG)


def footer():
    with st.sidebar:
        st.markdown(f'<div class="nb-foot">Version {VERSION}. Made by Bilal Dawood for sodium-ion half-cell data.'
                    f'<br>Files are analysed on this computer and are not uploaded anywhere.</div>',
                    unsafe_allow_html=True)
