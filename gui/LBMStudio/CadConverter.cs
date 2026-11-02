using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;

namespace LBMStudio
{
    /// <summary>
    /// Conversion of STEP files to STL with an external CAD kernel:
    /// gmsh (gmsh.exe or the Python module "gmsh") or FreeCAD (FreeCADCmd.exe).
    /// </summary>
    public static class CadConverter
    {
        /// <summary>
        /// Mesh size from the curvature for gmsh (Mesh.MeshSizeFromCurvature, -clcurv). Without it gmsh divides
        /// circles into 12-15 segments; 32 gives 64 segments on a full circle with gmsh 4.15 (chord error 0.12 % of the radius).
        /// </summary>
        const int GmshCurvatureElements = 32;
        /// <summary>FreeCAD angular deflection [rad], 0.1 gives about 64 segments on a full circle.</summary>
        const double FreeCadAngularDeflection = 0.1;

        public static bool IsStep(string path)
        {
            string ext = Path.GetExtension(path).ToLowerInvariant();
            return ext == ".step" || ext == ".stp";
        }

        static string FindOnPath(string exe)
        {
            string path = Environment.GetEnvironmentVariable("PATH") ?? "";
            foreach (var dir in path.Split(Path.PathSeparator))
            {
                try
                {
                    string f = Path.Combine(dir.Trim(), exe);
                    if (File.Exists(f))
                    {
                        return f;
                    }
                }
                catch (ArgumentException)
                {
                }
            }
            return null;
        }

        static string FindFreeCad(AppSettings s)
        {
            if (!string.IsNullOrEmpty(s.FreeCadExe) && File.Exists(s.FreeCadExe))
            {
                return s.FreeCadExe;
            }
            string onPath = FindOnPath("FreeCADCmd.exe");
            if (onPath != null)
            {
                return onPath;
            }
            foreach (var pf in new[] { Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles), Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData) + "\\Programs" })
            {
                try
                {
                    foreach (var d in Directory.GetDirectories(pf, "FreeCAD*"))
                    {
                        string exe = Path.Combine(d, "bin", "FreeCADCmd.exe");
                        if (File.Exists(exe))
                        {
                            return exe;
                        }
                    }
                }
                catch (Exception)
                {
                }
            }
            return null;
        }

        static string FindGmsh(AppSettings s)
        {
            if (!string.IsNullOrEmpty(s.GmshExe) && File.Exists(s.GmshExe))
            {
                return s.GmshExe;
            }
            return FindOnPath("gmsh.exe");
        }

        static string pythonWithGmsh;
        static bool pythonSearched;

        /// <summary>
        /// Full path of a Python interpreter with the module "gmsh", or null. Process.Start("python") resolves
        /// only python.exe on PATH and skips .bat shims (pyenv-win), so it can pick a different interpreter
        /// than the terminal does. The candidates are therefore started through cmd.exe and the py launcher
        /// and the interpreter reports its own path.
        /// </summary>
        static string FindPythonWithGmsh()
        {
            if (pythonSearched)
            {
                return pythonWithGmsh;
            }
            const string probe = "-c \"import gmsh, sys; print(sys.executable)\"";
            var candidates = new[]
            {
                new[] { "cmd.exe", "/c python " + probe },
                new[] { "py", "-3 " + probe },
                new[] { "python", probe }
            };
            foreach (var c in candidates)
            {
                try
                {
                    string output;
                    if (Run(c[0], c[1], 60000, out output) != 0)
                    {
                        continue;
                    }
                    string exe = output.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries)
                        .Select(l => l.Trim()).LastOrDefault(File.Exists);
                    if (exe != null)
                    {
                        pythonWithGmsh = exe;
                        break;
                    }
                }
                catch (Exception)
                {
                }
            }
            pythonSearched = true;
            return pythonWithGmsh;
        }

        /// <summary>Text describing the available converters.</summary>
        public static string Describe(AppSettings s)
        {
            var list = new List<string>();
            if (FindGmsh(s) != null)
            {
                list.Add("gmsh");
            }
            else if (FindPythonWithGmsh() != null)
            {
                list.Add("gmsh (Python)");
            }
            if (FindFreeCad(s) != null)
            {
                list.Add("FreeCAD");
            }
            return list.Count == 0 ? "not found (neither gmsh nor FreeCAD)" : string.Join(", ", list);
        }

        /// <summary>Converts STEP to STL, the result keeps the units of the STEP file (usually mm).</summary>
        /// <param name="maxElementSize">maximal triangle size in the units of the file (0 - default)</param>
        public static void StepToStl(AppSettings s, string step, string stl, double maxElementSize, Action<string> log)
        {
            if (File.Exists(stl))
            {
                File.Delete(stl);
            }
            string output;
            string gmsh = FindGmsh(s);
            if (gmsh != null)
            {
                string args = "\"" + step + "\" -2 -clcurv " + GmshCurvatureElements + " -format stl -o \"" + stl + "\"";
                if (maxElementSize > 0)
                {
                    args += " -clmax " + Fmt.D(maxElementSize);
                }
                log("gmsh " + args);
                Run(gmsh, args, 600000, out output);
                log(output);
                if (File.Exists(stl))
                {
                    return;
                }
            }
            string python = FindPythonWithGmsh();
            if (python != null)
            {
                string script = Path.Combine(Path.GetTempPath(), "lbmstudio_step2stl.py");
                File.WriteAllText(script,
                    "import sys, gmsh\n" +
                    "gmsh.initialize()\n" +
                    "gmsh.option.setNumber('General.Terminal', 1)\n" +
                    "gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', " + GmshCurvatureElements + ")\n" +
                    (maxElementSize > 0 ? "gmsh.option.setNumber('Mesh.MeshSizeMax', " + Fmt.D(maxElementSize) + ")\n" : "") +
                    "gmsh.open(sys.argv[1])\n" +
                    "gmsh.model.mesh.generate(2)\n" +
                    "gmsh.write(sys.argv[2])\n" +
                    "gmsh.finalize()\n");
                log(python + " (gmsh) " + step);
                Run(python,"\"" + script + "\" \"" + step + "\" \"" + stl + "\"", 600000, out output);
                log(output);
                if (File.Exists(stl))
                {
                    return;
                }
            }
            string freecad = FindFreeCad(s);
            if (freecad != null)
            {
                string script = Path.Combine(Path.GetTempPath(), "lbmstudio_step2stl_fc.py");
                double deflection = maxElementSize > 0 ? maxElementSize / 10 : 0.05;
                File.WriteAllText(script,
                    "import Part, MeshPart\n" +
                    "shape = Part.read(r'" + step + "')\n" +
                    "mesh = MeshPart.meshFromShape(Shape=shape, LinearDeflection=" + Fmt.D(deflection) + ", AngularDeflection=" + Fmt.D(FreeCadAngularDeflection) + ")\n" +
                    "mesh.write(r'" + stl + "')\n");
                log("FreeCADCmd " + script);
                Run(freecad, "\"" + script + "\"", 600000, out output);
                log(output);
                if (File.Exists(stl))
                {
                    return;
                }
            }
            throw new InvalidOperationException(
                "STEP → STL conversion failed. It needs gmsh (gmsh.exe in PATH or in Settings, or 'pip install gmsh') " +
                "or FreeCAD (FreeCADCmd.exe). Alternatively export the geometry from CAD directly to STL.");
        }

        static int Run(string exe, string args, int timeoutMs, out string output)
        {
            var psi = new ProcessStartInfo(exe, args)
            {
                UseShellExecute = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                CreateNoWindow = true
            };
            var sb = new StringBuilder();
            using (var p = new Process { StartInfo = psi })
            {
                p.OutputDataReceived += (o, e) => { if (e.Data != null) lock (sb) sb.AppendLine(e.Data); };
                p.ErrorDataReceived += (o, e) => { if (e.Data != null) lock (sb) sb.AppendLine(e.Data); };
                p.Start();
                p.BeginOutputReadLine();
                p.BeginErrorReadLine();
                if (!p.WaitForExit(timeoutMs))
                {
                    p.Kill();
                    output = sb + "\n(timeout)";
                    return -1;
                }
                p.WaitForExit();
                output = sb.ToString();
                return p.ExitCode;
            }
        }
    }
}
