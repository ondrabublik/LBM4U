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
        /// <summary>display cells per cell of the level 0 along an axis (2^display level of a refined mesh)</summary>
        public int Scale = 1;
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

        /// <summary>Finest level at which the results of the mesh are assembled (at most 64 million display cells).</summary>
        public static int DisplayLevel(GeometryDef g)
        {
            int d = 0;
            long n = (long)g.GX * g.GY * g.GZ;
            while (d + 1 < g.MeshLevels && n * (8L << (3 * d)) <= 64000000)
            {
                d++;
            }
            return d;
        }

        /// <param name="gx">dimensions of the level 0 grid</param>
        /// <param name="displayLevel">level of the assembled grid, coarser blocks are replicated, finer ones averaged</param>
        public static FieldsData Read(string path, int gx, int gy, int gz, int displayLevel = 0)
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
                int D = version >= 2 ? displayLevel : 0;
                d.Scale = 1 << D;
                gx <<= D;
                gy <<= D;
                gz <<= D;
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
                int[] count = null;       // finer blocks than the display level: values are averaged
                float[][] maxVal = null;
                int itype = Array.IndexOf(d.Names, "type"), iflag = Array.IndexOf(d.Names, "flag");
                for (int b = 0; b < nBlock; b++)
                {
                    int bi = r.ReadInt32(), bj = r.ReadInt32(), bk = r.ReadInt32();
                    int level = version >= 2 ? r.ReadInt32() : 0;
                    int up = Math.Max(0, D - level), down = Math.Max(0, level - D);
                    int rep = 1 << up;
                    if (down > 0 && count == null)
                    {
                        count = new int[n];
                        maxVal = new float[nFields][];
                    }
                    for (int f = 0; f < nFields; f++)
                    {
                        int read = r.Read(buf, 0, buf.Length);
                        if (read != buf.Length)
                        {
                            throw new EndOfStreamException("Neúplný soubor " + path);
                        }
                        Buffer.BlockCopy(buf, 0, vals, 0, buf.Length);
                        if (down > 0 && maxVal[f] == null)
                        {
                            maxVal[f] = new float[n];
                        }
                        for (int k = 0; k < d.Nk; k++)
                        {
                            for (int j = 0; j < d.Nj; j++)
                            {
                                for (int i = 0; i < d.Ni; i++)
                                {
                                    float v = vals[i + d.Ni * (j + d.Nj * k)];
                                    long gi = (long)bi * d.Ni + i, gj = (long)bj * d.Nj + j, gk = (long)bk * d.Nk + k;
                                    if (down > 0)
                                    {
                                        gi >>= down;
                                        gj >>= down;
                                        gk >>= down;
                                        if (gi >= gx || gj >= gy || gk >= gz)
                                        {
                                            continue;
                                        }
                                        long id = d.Index((int)gi, (int)gj, (int)gk);
                                        float old = d.Data[f][id];
                                        if (float.IsNaN(old))
                                        {
                                            d.Data[f][id] = v;
                                            maxVal[f][id] = v;
                                        }
                                        else if (f == itype || f == iflag)
                                        {
                                            d.Data[f][id] = Math.Min(old, v);
                                            maxVal[f][id] = Math.Max(maxVal[f][id], v);
                                        }
                                        else
                                        {
                                            d.Data[f][id] = old + v;
                                        }
                                        if (f == 0)
                                        {
                                            count[id]++;
                                        }
                                        continue;
                                    }
                                    gi <<= up;
                                    gj <<= up;
                                    gk <<= up;
                                    for (int c = 0; c < rep * rep * rep; c++)
                                    {
                                        long a = gi + c % rep, bb = gj + c / rep % rep, cc = gk + c / (rep * rep);
                                        if (a < gx && bb < gy && cc < gz)
                                        {
                                            d.Data[f][d.Index((int)a, (int)bb, (int)cc)] = v;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if (count != null)
                {
                    for (long id = 0; id < n; id++)
                    {
                        if (count[id] <= 1)
                        {
                            continue;
                        }
                        for (int f = 0; f < nFields; f++)
                        {
                            if (f == iflag)
                            {
                                // fluid and gas or interface inside: interface
                                if (d.Data[f][id] != maxVal[f][id])
                                {
                                    d.Data[f][id] = 3;
                                }
                            }
                            else if (f != itype)
                            {
                                d.Data[f][id] /= count[id];
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
