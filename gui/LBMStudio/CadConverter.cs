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

        static bool PythonHasGmsh()
        {
            try
            {
                int code = Run("python", "-c \"import gmsh\"", 60000, out _);
                return code == 0;
            }
            catch (Exception)
            {
                return false;
            }
        }

        /// <summary>Text describing the available converters.</summary>
        public static string Describe(AppSettings s)
        {
            var list = new List<string>();
            if (FindGmsh(s) != null)
            {
                list.Add("gmsh");
            }
            if (FindFreeCad(s) != null)
            {
                list.Add("FreeCAD");
            }
            return list.Count == 0 ? "nenalezen (gmsh ani FreeCAD)" : string.Join(", ", list);
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
                string args = "\"" + step + "\" -2 -format stl -o \"" + stl + "\"";
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
            if (PythonHasGmsh())
            {
                string script = Path.Combine(Path.GetTempPath(), "lbmstudio_step2stl.py");
                File.WriteAllText(script,
                    "import sys, gmsh\n" +
                    "gmsh.initialize()\n" +
                    "gmsh.option.setNumber('General.Terminal', 1)\n" +
                    (maxElementSize > 0 ? "gmsh.option.setNumber('Mesh.MeshSizeMax', " + Fmt.D(maxElementSize) + ")\n" : "") +
                    "gmsh.open(sys.argv[1])\n" +
                    "gmsh.model.mesh.generate(2)\n" +
                    "gmsh.write(sys.argv[2])\n" +
                    "gmsh.finalize()\n");
                log("python gmsh " + step);
                Run("python", "\"" + script + "\" \"" + step + "\" \"" + stl + "\"", 600000, out output);
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
                    "mesh = MeshPart.meshFromShape(Shape=shape, LinearDeflection=" + Fmt.D(deflection) + ", AngularDeflection=0.3)\n" +
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
                "Převod STEP → STL se nepodařil. Je potřeba gmsh (gmsh.exe v PATH nebo v Nastavení, případně 'pip install gmsh') " +
                "nebo FreeCAD (FreeCADCmd.exe). Alternativně exportujte geometrii z CAD přímo do STL.");
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
