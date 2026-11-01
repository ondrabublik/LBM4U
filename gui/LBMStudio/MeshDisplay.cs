using System;
using System.Collections.Generic;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;
using W = System.Windows;

namespace LBMStudio
{
    /// <summary>
    /// 3D display of the solver mesh: faces between fluid (inlet, outlet) cells and walls, every face textured
    /// with one cell so that the grid is visible, and the edges of the blocks.
    /// </summary>
    public static class MeshDisplay
    {
        /// <summary>Limit of displayed faces (memory and speed of the WPF viewport).</summary>
        public const long MaxFaces = 1500000;

        static readonly Color[] TypeColors = { Color.FromRgb(120, 160, 210), Colors.Transparent, Color.FromRgb(220, 60, 50), Color.FromRgb(60, 180, 75) };

        /// <summary>
        /// Boundary faces of the fluid region. With a cut (cutAxis 0-2, -1 = no cut) only the cells with the index
        /// along cutAxis up to cutPos (from cutPos with keepUpper) are shown, the cut plane then shows the cells
        /// inside. Returns null when the number of faces exceeds MaxFaces.
        /// </summary>
        public static Model3D Cells(GeometryDef g, byte[] types, int cutAxis, int cutPos, bool keepUpper, out long faceCount)
        {
            int GX = g.GX, GY = g.GY, GZ = g.GZ;
            double dx = g.Dx;
            Vec3 o = g.Origin;
            Func<int, int, int, bool> active = (i, j, k) =>
            {
                if (i < 0 || j < 0 || k < 0 || i >= GX || j >= GY || k >= GZ)
                {
                    return false;
                }
                if (cutAxis >= 0)
                {
                    int c = cutAxis == 0 ? i : (cutAxis == 1 ? j : k);
                    if (keepUpper ? c < cutPos : c > cutPos)
                    {
                        return false;
                    }
                }
                return types[i + (long)GX * (j + (long)GY * k)] != CaseGenerator.Wall;
            };
            int[,] dirs = { { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 } };

            faceCount = 0;
            for (int k = 0; k < GZ; k++)
                for (int j = 0; j < GY; j++)
                    for (int i = 0; i < GX; i++)
                    {
                        if (!active(i, j, k))
                        {
                            continue;
                        }
                        for (int d = 0; d < 6; d++)
                        {
                            if (!active(i + dirs[d, 0], j + dirs[d, 1], k + dirs[d, 2]))
                            {
                                faceCount++;
                            }
                        }
                    }
            if (faceCount > MaxFaces)
            {
                return null;
            }

            // collections are filled before they are attached to the mesh, otherwise every Add raises change events
            var positions = new Point3DCollection[TypeColors.Length];
            var uvs = new PointCollection[TypeColors.Length];
            var indices = new Int32Collection[TypeColors.Length];
            var corners = new Point3D[4];
            for (int k = 0; k < GZ; k++)
                for (int j = 0; j < GY; j++)
                    for (int i = 0; i < GX; i++)
                    {
                        if (!active(i, j, k))
                        {
                            continue;
                        }
                        byte type = types[i + (long)GX * (j + (long)GY * k)];
                        for (int d = 0; d < 6; d++)
                        {
                            if (active(i + dirs[d, 0], j + dirs[d, 1], k + dirs[d, 2]))
                            {
                                continue;
                            }
                            int axis = d / 2;
                            int a1 = axis == 0 ? 1 : 0, a2 = axis == 2 ? 1 : 2;
                            double[] lo = { o.X + i * dx, o.Y + j * dx, o.Z + k * dx };
                            var p = new double[3];
                            p[axis] = lo[axis] + (d % 2 == 1 ? dx : 0);
                            for (int q = 0; q < 4; q++)
                            {
                                p[a1] = lo[a1] + (q == 1 || q == 2 ? dx : 0);
                                p[a2] = lo[a2] + (q >= 2 ? dx : 0);
                                corners[q] = new Point3D(p[0], p[1], p[2]);
                            }
                            if (positions[type] == null)
                            {
                                positions[type] = new Point3DCollection();
                                uvs[type] = new PointCollection();
                                indices[type] = new Int32Collection();
                            }
                            var pos = positions[type];
                            var uv = uvs[type];
                            var idx = indices[type];
                            int b = pos.Count;
                            foreach (var pt in corners)
                            {
                                pos.Add(pt);
                            }
                            uv.Add(new W.Point(0, 0));
                            uv.Add(new W.Point(1, 0));
                            uv.Add(new W.Point(1, 1));
                            uv.Add(new W.Point(0, 1));
                            idx.Add(b);
                            idx.Add(b + 1);
                            idx.Add(b + 2);
                            idx.Add(b);
                            idx.Add(b + 2);
                            idx.Add(b + 3);
                        }
                    }

            var group = new Model3DGroup();
            for (int t = 0; t < positions.Length; t++)
            {
                if (positions[t] == null)
                {
                    continue;
                }
                var mesh = new MeshGeometry3D { Positions = positions[t], TextureCoordinates = uvs[t], TriangleIndices = indices[t] };
                mesh.Freeze();
                var mat = new DiffuseMaterial(CellBrush(TypeColors[t]));
                mat.Freeze();
                group.Children.Add(new GeometryModel3D(mesh, mat) { BackMaterial = mat });
            }
            group.Freeze();
            return group;
        }

        /// <summary>Texture of one cell: the colour with a darker border.</summary>
        static Brush CellBrush(Color color)
        {
            const int n = 32;
            var px = new int[n * n];
            Color edge = Color.FromRgb((byte)(color.R * 0.45), (byte)(color.G * 0.45), (byte)(color.B * 0.45));
            for (int y = 0; y < n; y++)
                for (int x = 0; x < n; x++)
                {
                    Color c = x == 0 || y == 0 || x == n - 1 || y == n - 1 ? edge : color;
                    px[y * n + x] = (255 << 24) | (c.R << 16) | (c.G << 8) | c.B;
                }
            var bmp = BitmapSource.Create(n, n, 96, 96, PixelFormats.Bgra32, null, px, n * 4);
            bmp.Freeze();
            var brush = new ImageBrush(bmp) { ViewportUnits = BrushMappingMode.Absolute, Viewport = new W.Rect(0, 0, 1, 1), TileMode = TileMode.None };
            brush.Freeze();
            return brush;
        }

        /// <summary>Edges of the blocks that contain other cells than walls.</summary>
        public static Model3D BlockEdges(GeometryDef g, byte[] types, Color color)
        {
            int bs = g.BlockSize, GX = g.GX, GY = g.GY;
            double h = bs * g.Dx;
            Vec3 o = g.Origin;
            var edges = new HashSet<long>(); // node index * 3 + axis, shared edges are drawn once
            Func<int, int, int, long> node = (i, j, k) => i + (long)(g.NBi + 1) * (j + (long)(g.NBj + 1) * k);
            for (int bk = 0; bk < g.NBk; bk++)
                for (int bj = 0; bj < g.NBj; bj++)
                    for (int bi = 0; bi < g.NBi; bi++)
                    {
                        bool used = false;
                        for (int k = 0; k < bs && !used; k++)
                            for (int j = 0; j < bs && !used; j++)
                                for (int i = 0; i < bs && !used; i++)
                                {
                                    used = types[bi * bs + i + (long)GX * (bj * bs + j + (long)GY * (bk * bs + k))] != CaseGenerator.Wall;
                                }
                        if (!used)
                        {
                            continue;
                        }
                        for (int e = 0; e < 4; e++)
                        {
                            int p = e & 1, q = e >> 1;
                            edges.Add(node(bi, bj + p, bk + q) * 3 + 0);
                            edges.Add(node(bi + p, bj, bk + q) * 3 + 1);
                            edges.Add(node(bi + p, bj + q, bk) * 3 + 2);
                        }
                    }
            var xyz = new List<double>();
            double r = h * 0.012;
            foreach (long key in edges)
            {
                int axis = (int)(key % 3);
                long n = key / 3;
                int i = (int)(n % (g.NBi + 1)), j = (int)(n / (g.NBi + 1) % (g.NBj + 1)), k = (int)(n / ((long)(g.NBi + 1) * (g.NBj + 1)));
                var a = new Vec3(o.X + i * h, o.Y + j * h, o.Z + k * h);
                var b = a + new Vec3(axis == 0 ? h : 0, axis == 1 ? h : 0, axis == 2 ? h : 0);
                Viewer3D.AddCylinder(xyz, a, b, r, r, 4);
            }
            return xyz.Count == 0 ? null : Viewer3D.Triangles(xyz, color, 1);
        }
    }
}
