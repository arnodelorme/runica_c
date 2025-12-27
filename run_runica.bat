@echo off
REM run_runica.bat - Run runica_simple on EEG data
REM Usage: run_runica.bat [dataset_basename]
REM
REM Example: run_runica.bat .\data\eeglab_data

setlocal enabledelayedexpansion

REM Default parameters
set DATASET=%~1
if "%DATASET%"=="" set DATASET=.\data\eeglab_data

REM Derived paths
set DATAFILE=%DATASET%.fdt
set WEIGHTSFILE=%DATASET%.wts_win
set SPHEREFILE=%DATASET%.sph_win

echo ========================================
echo Running runica_simple
echo ========================================
echo Dataset: %DATASET%
echo.

REM Check if data file exists
if not exist "%DATAFILE%" (
    echo Error: Data file %DATAFILE% not found
    exit /b 1
)

REM Check if executable exists
if not exist "runica.exe" (
    echo Error: runica.exe not found. Please compile the project first.
    exit /b 1
)

REM Run runica
echo Running runica.exe...
runica.exe "%DATAFILE%" "%WEIGHTSFILE%" "%SPHEREFILE%"

if errorlevel 1 (
    echo Error: runica.exe failed
    exit /b 1
)

echo.
echo ========================================
echo All done!
echo ========================================
echo Output files:
echo   Weights: %WEIGHTSFILE%
echo   Sphere:  %SPHEREFILE%
echo.

endlocal
