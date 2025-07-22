@echo off
echo ===============================================
echo    Exotic Options Derivative Pricer - C++
echo ===============================================
echo.

echo Compiling with optimizations...
g++ -std=c++17 -O3 -Wall -Wextra -march=native -o exotic_options.exe main.cpp

if %ERRORLEVEL% EQU 0 (
    echo.
    echo Compilation successful!
    echo.
    echo Starting Exotic Options Pricer...
    echo.
    exotic_options.exe
) else (
    echo.
    echo Compilation failed! Please check for errors.
    echo Make sure you have a C++17 compatible compiler installed.
    echo.
    pause
)

echo.
echo Program finished.
pause