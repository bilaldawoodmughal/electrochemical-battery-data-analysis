@echo off
REM Starts the analyser in your browser. Keep this window open while you use the app.
REM If you use a virtual environment, activate it before running this file, or replace
REM "python" below with the full path to the python.exe inside it.
cd /d "%~dp0"
python -m streamlit run app.py
pause
