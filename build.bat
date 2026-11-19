@echo off
rem Builds the solver and the unit tests with CMake (Visual Studio 2019/2022/2026 or MinGW g++).
rem Usage: build.bat [msvc|mingw]
rem The newest Visual Studio with the C++ tools is found with vswhere; its bundled CMake and Ninja are used.
rem (No PATH expansion inside parenthesised blocks: "Program Files (x86)" in PATH would break the parser.)
setlocal
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
set "VSPATH="
set "VSVER="
if not exist "%VSWHERE%" goto novswhere
for /f "usebackq delims=" %%i in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VSPATH=%%i"
for /f "usebackq tokens=1 delims=." %%v in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationVersion`) do set "VSVER=%%v"
:novswhere
if not defined VSPATH goto nocmake
set "VSCMAKE=%VSPATH%\Common7\IDE\CommonExtensions\Microsoft\CMake"
set "PATH=%VSCMAKE%\CMake\bin;%VSCMAKE%\Ninja;%PATH%"
:nocmake

if /I "%1"=="mingw" goto mingw

set "GEN="
if "%VSVER%"=="18" set "GEN=Visual Studio 18 2026"
if "%VSVER%"=="17" set "GEN=Visual Studio 17 2022"
if "%VSVER%"=="16" set "GEN=Visual Studio 16 2019"
if not defined GEN (
	echo No supported Visual Studio with C++ tools found. Install Visual Studio or use "build.bat mingw".
	exit /b 1
)
echo Using %GEN%
rem a build directory configured with another Visual Studio version cannot be reused
if not exist build-msvc\CMakeCache.txt goto configure
findstr /C:"CMAKE_GENERATOR:INTERNAL=%GEN%" build-msvc\CMakeCache.txt >nul && goto configure
echo build-msvc was configured with another generator, reconfiguring
rmdir /S /Q build-msvc
:configure
cmake -S . -B build-msvc -G "%GEN%" -A x64 || exit /b 1
cmake --build build-msvc --config Release || exit /b 1
copy /Y build-msvc\Release\cLBMnet.exe cLBMnet.exe >nul
goto done

:mingw
set "PATH=C:\msys64\mingw64\bin;%PATH%"
cmake -S . -B build-mingw -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ || exit /b 1
cmake --build build-mingw || exit /b 1
copy /Y build-mingw\cLBMnet.exe cLBMnet.exe >nul

:done
echo Build finished: cLBMnet.exe
