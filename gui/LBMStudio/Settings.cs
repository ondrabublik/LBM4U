using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace LBMStudio
{
    /// <summary>Ordered "key = value" file, the same format as parameters.txt of the solver.</summary>
    public class KeyValueFile
    {
        public readonly List<KeyValuePair<string, string>> Items = new List<KeyValuePair<string, string>>();

        public static KeyValueFile Load(string path)
        {
            var kv = new KeyValueFile();
            if (!File.Exists(path))
            {
                return kv;
            }
            foreach (var raw in File.ReadAllLines(path))
            {
                string line = raw;
                int c = line.IndexOfAny(new[] { '#', '%' });
                if (c >= 0)
                {
                    line = line.Substring(0, c);
                }
                int eq = line.IndexOf('=');
                if (eq <= 0)
                {
                    continue;
                }
                kv.Set(line.Substring(0, eq).Trim(), line.Substring(eq + 1).Trim());
            }
            return kv;
        }

        public void Save(string path)
        {
            var sb = new StringBuilder();
            foreach (var item in Items)
            {
                sb.Append(item.Key).Append(" = ").Append(item.Value).Append("\r\n");
            }
            File.WriteAllText(path, sb.ToString());
        }

        public bool Has(string key)
        {
            return Items.Any(i => i.Key == key);
        }

        public string Get(string key, string def = "")
        {
            foreach (var item in Items)
            {
                if (item.Key == key)
                {
                    return item.Value;
                }
            }
            return def;
        }

        public double GetDouble(string key, double def = 0)
        {
            double v;
            return double.TryParse(Get(key, ""), NumberStyles.Float, CultureInfo.InvariantCulture, out v) ? v : def;
        }

        public int GetInt(string key, int def = 0)
        {
            double v;
            return double.TryParse(Get(key, ""), NumberStyles.Float, CultureInfo.InvariantCulture, out v) ? (int)Math.Round(v) : def;
        }

        public void Set(string key, string value)
        {
            for (int i = 0; i < Items.Count; i++)
            {
                if (Items[i].Key == key)
                {
                    Items[i] = new KeyValuePair<string, string>(key, value);
                    return;
                }
            }
            Items.Add(new KeyValuePair<string, string>(key, value));
        }

        public void Set(string key, double value)
        {
            Set(key, Fmt.D(value));
        }

        public void Set(string key, int value)
        {
            Set(key, value.ToString(CultureInfo.InvariantCulture));
        }
    }

    public static class Fmt
    {
        public static string D(double v)
        {
            return v.ToString("R", CultureInfo.InvariantCulture);
        }

        public static string G(double v, int digits = 4)
        {
            return v.ToString("G" + digits, CultureInfo.InvariantCulture);
        }

        public static bool TryParse(string s, out double v)
        {
            s = (s ?? "").Trim().Replace(',', '.');
            return double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out v);
        }

        public static double Parse(string s, string what)
        {
            double v;
            if (!TryParse(s, out v))
            {
                throw new FormatException("Neplatná hodnota pole \"" + what + "\": " + s);
            }
            return v;
        }
    }

    /// <summary>Application settings stored in %APPDATA%\LBMStudio\settings.ini.</summary>
    public class AppSettings
    {
        public string Workspace = "";
        public string SolverExe = "";
        public string GmshExe = "";
        public string FreeCadExe = "";
        public string ParaViewExe = "";

        static string FilePath
        {
            get
            {
                string dir = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.ApplicationData), "LBMStudio");
                Directory.CreateDirectory(dir);
                return Path.Combine(dir, "settings.ini");
            }
        }

        public static AppSettings Load()
        {
            var kv = KeyValueFile.Load(FilePath);
            var s = new AppSettings
            {
                Workspace = kv.Get("workspace"),
                SolverExe = kv.Get("solver"),
                GmshExe = kv.Get("gmsh"),
                FreeCadExe = kv.Get("freecad"),
                ParaViewExe = kv.Get("paraview")
            };
            if (string.IsNullOrEmpty(s.SolverExe) || !File.Exists(s.SolverExe))
            {
                s.SolverExe = FindSolver();
            }
            if (string.IsNullOrEmpty(s.Workspace))
            {
                string root = RepositoryRoot();
                s.Workspace = Path.Combine(root ?? Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments), "simulations");
            }
            if (string.IsNullOrEmpty(s.ParaViewExe))
            {
                s.ParaViewExe = FindParaView();
            }
            return s;
        }

        public void Save()
        {
            var kv = new KeyValueFile();
            kv.Set("workspace", Workspace);
            kv.Set("solver", SolverExe);
            kv.Set("gmsh", GmshExe);
            kv.Set("freecad", FreeCadExe);
            kv.Set("paraview", ParaViewExe);
            kv.Save(FilePath);
        }

        /// <summary>Directory with CMakeLists.txt of the solver above the executable.</summary>
        public static string RepositoryRoot()
        {
            var dir = new DirectoryInfo(AppDomain.CurrentDomain.BaseDirectory);
            while (dir != null)
            {
                if (File.Exists(Path.Combine(dir.FullName, "CMakeLists.txt")) && File.Exists(Path.Combine(dir.FullName, "Mesh.cpp")))
                {
                    return dir.FullName;
                }
                dir = dir.Parent;
            }
            return null;
        }

        static string FindSolver()
        {
            var candidates = new List<string> { Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "cLBMnet.exe") };
            string root = RepositoryRoot();
            if (root != null)
            {
                candidates.Add(Path.Combine(root, "cLBMnet.exe"));
                candidates.Add(Path.Combine(root, "build-msvc", "Release", "cLBMnet.exe"));
                candidates.Add(Path.Combine(root, "build-mingw", "cLBMnet.exe"));
            }
            return candidates.FirstOrDefault(File.Exists) ?? "";
        }

        static string FindParaView()
        {
            string pf = Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles);
            try
            {
                foreach (var d in Directory.GetDirectories(pf, "ParaView*").OrderByDescending(x => x))
                {
                    string exe = Path.Combine(d, "bin", "paraview.exe");
                    if (File.Exists(exe))
                    {
                        return exe;
                    }
                }
            }
            catch (Exception)
            {
            }
            return "";
        }
    }
}
