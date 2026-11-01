using System;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Text.RegularExpressions;

namespace LBMStudio
{
    public class ProgressInfo
    {
        public int Iteration;
        public double Time, UMax, Mass = double.NaN, Mlups;
    }

    /// <summary>Runs cLBMnet.exe for one simulation and reports its output.</summary>
    public class SolverRunner
    {
        Process process;
        string simDir;
        static readonly Regex IterRegex = new Regex(@"Iteration:\s*(\d+),\s*t\s*=\s*([^,]+),\s*uMax\s*=\s*([^,]+)(?:,\s*mass\s*=\s*([^,]+))?.*MLUPS\s*=\s*(\S+)", RegexOptions.Compiled);

        public event Action<string> Line;
        public event Action<ProgressInfo> Progress;
        public event Action<int> Exited;

        public bool Running { get { return process != null && !process.HasExited; } }
        public string SimulationDir { get { return simDir; } }

        public void Start(string solverExe, string geomDir, string simDir)
        {
            if (Running)
            {
                throw new InvalidOperationException("Výpočet už běží.");
            }
            if (string.IsNullOrEmpty(solverExe) || !File.Exists(solverExe))
            {
                throw new FileNotFoundException("Řešič cLBMnet.exe nebyl nalezen, nastavte cestu v Nastavení.");
            }
            this.simDir = simDir;
            // the solver appends the trailing separator itself
            var psi = new ProcessStartInfo(solverExe, "--paths \"" + geomDir.TrimEnd('\\', '/') + "\" \"" + simDir.TrimEnd('\\', '/') + "\"")
            {
                WorkingDirectory = simDir,
                UseShellExecute = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                CreateNoWindow = true
            };
            process = new Process { StartInfo = psi, EnableRaisingEvents = true };
            process.OutputDataReceived += (o, e) => { if (e.Data != null) Handle(e.Data); };
            process.ErrorDataReceived += (o, e) => { if (e.Data != null) Handle(e.Data); };
            process.Exited += (o, e) =>
            {
                var p = (Process)o;
                p.WaitForExit();
                Exited?.Invoke(p.ExitCode);
            };
            process.Start();
            process.BeginOutputReadLine();
            process.BeginErrorReadLine();
        }

        void Handle(string line)
        {
            Line?.Invoke(line);
            var m = IterRegex.Match(line);
            if (m.Success)
            {
                var info = new ProgressInfo
                {
                    Iteration = int.Parse(m.Groups[1].Value, CultureInfo.InvariantCulture),
                    Time = ParseD(m.Groups[2].Value),
                    UMax = ParseD(m.Groups[3].Value),
                    Mlups = ParseD(m.Groups[5].Value)
                };
                if (m.Groups[4].Success)
                {
                    info.Mass = ParseD(m.Groups[4].Value);
                }
                Progress?.Invoke(info);
            }
        }

        static double ParseD(string s)
        {
            double v;
            return double.TryParse(s.Trim(), NumberStyles.Float, CultureInfo.InvariantCulture, out v) ? v : double.NaN;
        }

        /// <summary>Asks the solver to stop at the next save (it writes the results).</summary>
        public void RequestStop()
        {
            if (Running)
            {
                File.WriteAllText(Path.Combine(simDir, "STOP"), "stop");
            }
        }

        public void Kill()
        {
            if (Running)
            {
                process.Kill();
            }
        }

        /// <summary>History from result/residuum.txt (iteration, uMax[, mass]).</summary>
        public static System.Collections.Generic.List<double[]> ReadResiduum(string simDir)
        {
            var list = new System.Collections.Generic.List<double[]>();
            string path = Path.Combine(simDir, "result", "residuum.txt");
            if (!File.Exists(path))
            {
                return list;
            }
            using (var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
            using (var r = new StreamReader(fs))
            {
                string line;
                while ((line = r.ReadLine()) != null)
                {
                    var t = line.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    var v = new double[t.Length];
                    bool ok = t.Length >= 2;
                    for (int i = 0; i < t.Length && ok; i++)
                    {
                        ok = double.TryParse(t[i], NumberStyles.Float, CultureInfo.InvariantCulture, out v[i]);
                    }
                    if (ok)
                    {
                        list.Add(v);
                    }
                }
            }
            return list;
        }
    }
}
