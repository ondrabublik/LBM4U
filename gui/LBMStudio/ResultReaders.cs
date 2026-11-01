using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace LBMStudio
{
    /// <summary>Fields written by the solver to result/fields_&lt;iter&gt;.bin, assembled on the global grid.</summary>
    public class FieldsData
    {
        public int Iteration, ProblemId, Ni, Nj, Nk, GX, GY, GZ;
        public double Time, Ch;
        public string[] Names;
        /// <summary>[field][i + GX*(j + GY*k)], NaN outside of the blocks</summary>
        public float[][] Data;

        public int FieldIndex(string name)
        {
            return Array.IndexOf(Names, name);
        }

        public long Index(int i, int j, int k)
        {
            return i + (long)GX * (j + (long)GY * k);
        }

        public static FieldsData Read(string path, int gx, int gy, int gz)
        {
            using (var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
            using (var r = new BinaryReader(fs))
            {
                var magic = Encoding.ASCII.GetString(r.ReadBytes(4));
                if (magic != "LBMF")
                {
                    throw new InvalidDataException("Neplatný soubor výsledků " + path);
                }
                int version = r.ReadInt32();
                var d = new FieldsData
                {
                    Iteration = r.ReadInt32(),
                    Time = r.ReadDouble(),
                    ProblemId = r.ReadInt32()
                };
                int nBlock = r.ReadInt32();
                d.Ni = r.ReadInt32();
                d.Nj = r.ReadInt32();
                d.Nk = r.ReadInt32();
                d.Ch = r.ReadDouble();
                int nFields = r.ReadInt32();
                d.Names = new string[nFields];
                for (int f = 0; f < nFields; f++)
                {
                    d.Names[f] = Encoding.ASCII.GetString(r.ReadBytes(16)).TrimEnd('\0');
                }
                d.GX = gx;
                d.GY = gy;
                d.GZ = gz;
                long n = (long)gx * gy * gz;
                d.Data = new float[nFields][];
                for (int f = 0; f < nFields; f++)
                {
                    d.Data[f] = new float[n];
                    for (long i = 0; i < n; i++)
                    {
                        d.Data[f][i] = float.NaN;
                    }
                }
                int cells = d.Ni * d.Nj * d.Nk;
                var buf = new byte[cells * 4];
                var vals = new float[cells];
                for (int b = 0; b < nBlock; b++)
                {
                    int bi = r.ReadInt32(), bj = r.ReadInt32(), bk = r.ReadInt32();
                    for (int f = 0; f < nFields; f++)
                    {
                        int read = r.Read(buf, 0, buf.Length);
                        if (read != buf.Length)
                        {
                            throw new EndOfStreamException("Neúplný soubor " + path);
                        }
                        Buffer.BlockCopy(buf, 0, vals, 0, buf.Length);
                        for (int k = 0; k < d.Nk; k++)
                        {
                            for (int j = 0; j < d.Nj; j++)
                            {
                                for (int i = 0; i < d.Ni; i++)
                                {
                                    int gi = bi * d.Ni + i, gj = bj * d.Nj + j, gk = bk * d.Nk + k;
                                    if (gi < gx && gj < gy && gk < gz)
                                    {
                                        d.Data[f][d.Index(gi, gj, gk)] = vals[i + d.Ni * (j + d.Nj * k)];
                                    }
                                }
                            }
                        }
                    }
                }
                return d;
            }
        }

        /// <summary>Result files of a simulation sorted by iteration.</summary>
        public static List<KeyValuePair<int, string>> List(string simDir)
        {
            var list = new List<KeyValuePair<int, string>>();
            string dir = Path.Combine(simDir, "result");
            if (!Directory.Exists(dir))
            {
                return list;
            }
            var rx = new Regex(@"fields_(\d+)\.bin$");
            foreach (var f in Directory.GetFiles(dir, "fields_*.bin"))
            {
                var m = rx.Match(f);
                if (m.Success)
                {
                    list.Add(new KeyValuePair<int, string>(int.Parse(m.Groups[1].Value), f));
                }
            }
            return list.OrderBy(x => x.Key).ToList();
        }
    }

    /// <summary>Triangulated free surface from result/outputFS_*.vtk (legacy VTK, binary or ASCII).</summary>
    public static class SurfaceVtk
    {
        /// <summary>Returns 9 floats per triangle in lattice units.</summary>
        public static float[] Read(string path)
        {
            byte[] all = File.ReadAllBytes(path);
            int pos = 0;
            Func<string> readLine = () =>
            {
                int start = pos;
                while (pos < all.Length && all[pos] != '\n')
                {
                    pos++;
                }
                string s = Encoding.ASCII.GetString(all, start, pos - start).Trim();
                pos++;
                return s;
            };
            readLine();
            readLine();
            bool binary = readLine().ToUpperInvariant().StartsWith("BINARY");
            string line;
            do
            {
                line = readLine();
            } while (pos < all.Length && !line.StartsWith("POINTS"));
            if (!line.StartsWith("POINTS"))
            {
                return new float[0];
            }
            int n = int.Parse(line.Split(' ')[1], CultureInfo.InvariantCulture);
            var pts = new float[3 * n];
            if (binary)
            {
                for (int i = 0; i < 3 * n; i++)
                {
                    byte[] b = { all[pos + 4 * i + 3], all[pos + 4 * i + 2], all[pos + 4 * i + 1], all[pos + 4 * i] };
                    pts[i] = BitConverter.ToSingle(b, 0);
                }
            }
            else
            {
                int read = 0;
                while (read < 3 * n && pos < all.Length)
                {
                    foreach (var t in readLine().Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries))
                    {
                        if (read < 3 * n)
                        {
                            pts[read++] = float.Parse(t, CultureInfo.InvariantCulture);
                        }
                    }
                }
            }
            return pts;
        }

        public static string FileForIteration(string simDir, int iteration, double time)
        {
            // the solver names the files by time in microseconds
            long us = (long)Math.Round(time * 1e6);
            string path = Path.Combine(simDir, "result", "outputFS_" + us + "us.vtk");
            if (File.Exists(path))
            {
                return path;
            }
            // tolerate rounding differences
            string dir = Path.Combine(simDir, "result");
            if (!Directory.Exists(dir))
            {
                return null;
            }
            var rx = new Regex(@"outputFS_(\d+)us\.vtk$");
            string best = null;
            long bestDiff = long.MaxValue;
            foreach (var f in Directory.GetFiles(dir, "outputFS_*us.vtk"))
            {
                var m = rx.Match(f);
                if (m.Success)
                {
                    long diff = Math.Abs(long.Parse(m.Groups[1].Value) - us);
                    if (diff < bestDiff)
                    {
                        bestDiff = diff;
                        best = f;
                    }
                }
            }
            return bestDiff <= 2 ? best : null;
        }
    }
}
