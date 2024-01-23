@echo off
title Molecule Odor Analyzer Initialization

echo Starting pgAdmin 4 data writing...
set /p password=Enter your postgres password:
setlocal enableextensions enabledelayedexpansion
set "search=PLACEHOLDER"
set "replace=!password!"
set "file=Molecule_Odor_Analyzer\Baza\Baza_account.py"
set "tempfile=Molecule_Odor_Analyzer\Baza\Baza_account_temp.py"
(for /f "delims=" %%i in ('type "%file%"') do (
    set "line=%%i"
    set "line=!line:%search%=%replace%!"
    echo(!line!
))>"%tempfile%"
move /y "%tempfile%" "%file%" >nul

echo Installing required libraries...
pip install -r requirements.txt
echo Starting initialization...
python -m Molecule_Odor_Analyzer.Baza.Baza_inicijalizacija
echo Starting database fill process...
python -m Molecule_Odor_Analyzer.Baza.Baza_punjenje
echo Script successfully executed
pause