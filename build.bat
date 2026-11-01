@echo off
rem Builds the solver and the unit tests with CMake (Visual Studio 2022 or MinGW g++).
rem Usage: build.bat [msvc|mingw]
setlocal
set VSCMAKE=C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake
set PATH=%VSCMAKE%\CMake\bin;%VSCMAKE%\Ninja;%PATH%
if /I "%1"=="mingw" (
	set PATH=C:\msys64\mingw64\bin;%PATH%
	cmake -S . -B build-mingw -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ || exit /b 1
	cmake --build build-mingw || exit /b 1
	copy /Y build-mingw\cLBMnet.exe cLBMnet.exe >nul
) else (
	cmake -S . -B build-msvc -G "Visual Studio 17 2022" -A x64 || exit /b 1
	cmake --build build-msvc --config Release || exit /b 1
	copy /Y build-msvc\Release\cLBMnet.exe cLBMnet.exe >nul
)
echo Build finished: cLBMnet.exe
