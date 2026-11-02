@echo off
rem Builds LBMStudio.exe (.NET Framework 4.8) with the C# compiler of Visual Studio 2022.
rem With the .NET SDK installed you can use "dotnet build LBMStudio\LBMStudio.csproj" instead.
setlocal
cd /d "%~dp0"
set CSC=C:\Program Files\Microsoft Visual Studio\2022\Community\MSBuild\Current\Bin\Roslyn\csc.exe
if not exist "%CSC%" set CSC=C:\Program Files\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\Roslyn\csc.exe
if not exist "%CSC%" set CSC=C:\Program Files\Microsoft Visual Studio\2022\Enterprise\MSBuild\Current\Bin\Roslyn\csc.exe
set FW=%WINDIR%\Microsoft.NET\Framework64\v4.0.30319
set WPF=%FW%\WPF
if not exist bin mkdir bin
"%CSC%" /nologo /target:winexe /optimize+ /langversion:7.3 /out:bin\LBMStudio.exe /utf8output ^
 /r:"%FW%\System.dll" /r:"%FW%\System.Core.dll" /r:"%FW%\System.Drawing.dll" /r:"%FW%\System.Windows.Forms.dll" ^
 /r:"%FW%\System.Windows.Forms.DataVisualization.dll" /r:"%FW%\System.Xaml.dll" ^
 /r:"%WPF%\WindowsBase.dll" /r:"%WPF%\PresentationCore.dll" /r:"%WPF%\PresentationFramework.dll" /r:"%WPF%\WindowsFormsIntegration.dll" ^
 LBMStudio\*.cs || exit /b 1
echo Build finished: gui\bin\LBMStudio.exe
