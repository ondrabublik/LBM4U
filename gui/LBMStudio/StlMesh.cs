using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;

namespace LBMStudio
{
    public struct Vec3
    {
        public double X, Y, Z;

        public Vec3(double x, double y, double z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        public double this[int i]
        {
            get { return i == 0 ? X : (i == 1 ? Y : Z); }
        }

        public static Vec3 operator +(Vec3 a, Vec3 b) { return new Vec3(a.X + b.X, a.Y + b.Y, a.Z + b.Z); }
        public static Vec3 operator -(Vec3 a, Vec3 b) { return new Vec3(a.X - b.X, a.Y - b.Y, a.Z - b.Z); }
        public static Vec3 operator *(Vec3 a, double s) { return new Vec3(a.X * s, a.Y * s, a.Z * s); }
        public double Length { get { return Math.Sqrt(X * X + Y * Y + Z * Z); } }

        public static Vec3 Cross(Vec3 a, Vec3 b)
        {
            return new Vec3(a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z, a.X * b.Y - a.Y * b.X);
        }
    }

    /// <summary>Triangle surface mesh (STL), coordinates in metres after scaling.</summary>
    public class StlMesh
    {
        /// <summary>9 values per triangle.</summary>
        public double[] V = new double[0];

        public int TriangleCount { get { return V.Length / 9; } }

        public Vec3 Min, Max;

        public static StlMesh Read(string path)
        {
            byte[] data = File.ReadAllBytes(path);
            var mesh = new StlMesh();
            bool binary = true;
            if (data.Length >= 84)
            {
                uint n = BitConverter.ToUInt32(data, 80);
                if (84 + 50L * n != data.Length)
                {
                    binary = false;
                }
            }
            else
            {
                binary = false;
            }
            if (binary)
            {
                uint n = BitConverter.ToUInt32(data, 80);
                mesh.V = new double[9 * n];
                for (int t = 0; t < n; t++)
                {
                    int off = 84 + 50 * t + 12;
                    for (int k = 0; k < 9; k++)
                    {
                        mesh.V[9 * t + k] = BitConverter.ToSingle(data, off + 4 * k);
                    }
                }
            }
            else
            {
                var list = new List<double>();
                string text = Encoding.ASCII.GetString(data);
                foreach (var raw in text.Split('\n'))
                {
                    string line = raw.Trim();
                    if (line.StartsWith("vertex", StringComparison.OrdinalIgnoreCase))
                    {
                        var p = line.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                        for (int k = 1; k <= 3; k++)
                        {
                            list.Add(double.Parse(p[k], NumberStyles.Float, CultureInfo.InvariantCulture));
                        }
                    }
                }
                if (list.Count == 0 || list.Count % 9 != 0)
                {
                    throw new InvalidDataException("File " + path + " is not a valid STL.");
                }
                mesh.V = list.ToArray();
            }
            if (mesh.TriangleCount == 0)
            {
                throw new InvalidDataException("The STL contains no triangles.");
            }
            mesh.UpdateBounds();
            return mesh;
        }

        public void WriteBinary(string path)
        {
            using (var w = new BinaryWriter(File.Create(path)))
            {
                w.Write(new byte[80]);
                w.Write((uint)TriangleCount);
                for (int t = 0; t < TriangleCount; t++)
                {
                    var a = Vertex(t, 0);
                    var n = Vec3.Cross(Vertex(t, 1) - a, Vertex(t, 2) - a);
                    double l = n.Length;
                    if (l > 0)
                    {
                        n = n * (1 / l);
                    }
                    w.Write((float)n.X);
                    w.Write((float)n.Y);
                    w.Write((float)n.Z);
                    for (int k = 0; k < 9; k++)
                    {
                        w.Write((float)V[9 * t + k]);
                    }
                    w.Write((ushort)0);
                }
            }
        }

        public Vec3 Vertex(int t, int v)
        {
            return new Vec3(V[9 * t + 3 * v], V[9 * t + 3 * v + 1], V[9 * t + 3 * v + 2]);
        }

        public void UpdateBounds()
        {
            double[] mn = { double.MaxValue, double.MaxValue, double.MaxValue };
            double[] mx = { double.MinValue, double.MinValue, double.MinValue };
            for (int i = 0; i < V.Length; i++)
            {
                int d = i % 3;
                mn[d] = Math.Min(mn[d], V[i]);
                mx[d] = Math.Max(mx[d], V[i]);
            }
            Min = new Vec3(mn[0], mn[1], mn[2]);
            Max = new Vec3(mx[0], mx[1], mx[2]);
        }

        public void Scale(double s)
        {
            for (int i = 0; i < V.Length; i++)
            {
                V[i] *= s;
            }
            UpdateBounds();
        }

        /// <summary>Volume by the divergence theorem (positive for outward normals).</summary>
        public double SignedVolume()
        {
            double vol = 0;
            for (int t = 0; t < TriangleCount; t++)
            {
                var a = Vertex(t, 0);
                var b = Vertex(t, 1);
                var c = Vertex(t, 2);
                vol += a.X * (b.Y * c.Z - b.Z * c.Y) - a.Y * (b.X * c.Z - b.Z * c.X) + a.Z * (b.X * c.Y - b.Y * c.X);
            }
            return vol / 6;
        }
    }
}
