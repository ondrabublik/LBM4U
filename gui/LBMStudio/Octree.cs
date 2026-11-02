using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace LBMStudio
{
    /// <summary>Block of the octree mesh: BlockSize^3 cells of the size Dx / 2^Level, Types[i + bs*(j + bs*k)].</summary>
    public class OctBlock
    {
        public int Bi, Bj, Bk, Level;
        public byte[] Types;
    }

    /// <summary>
    /// Static octree of blocks (mesh/octree.bin). A block of the level L with the coordinates (bi, bj, bk) covers the
    /// cells (bi*bs .. bi*bs+bs-1, ...) of the level L, the cell size is Dx / 2^L; neighbouring blocks differ by one level at most.
    /// </summary>
    public class OctreeMesh
    {
        public int BlockSize, NBi, NBj, NBk, MaxLevel;
        public List<OctBlock> Blocks = new List<OctBlock>();
        Dictionary<long, OctBlock> map = new Dictionary<long, OctBlock>();

        public static long Key(int level, int bi, int bj, int bk)
        {
            return ((long)level << 60) | ((long)bi << 40) | ((long)bj << 20) | (long)bk;
        }

        public void BuildIndex()
        {
            map = new Dictionary<long, OctBlock>();
            MaxLevel = 0;
            foreach (var b in Blocks)
            {
                map[Key(b.Level, b.Bi, b.Bj, b.Bk)] = b;
                MaxLevel = Math.Max(MaxLevel, b.Level);
            }
        }

        public OctBlock Find(int level, int bi, int bj, int bk)
        {
            if (bi < 0 || bj < 0 || bk < 0)
            {
                return null;
            }
            OctBlock b;
            return map.TryGetValue(Key(level, bi, bj, bk), out b) ? b : null;
        }

        /// <summary>Identification of the block in the file names of the solver.</summary>
        public string Tag(OctBlock b)
        {
            if (b.Level == 0)
            {
                return (NBi * NBj * b.Bk + NBi * b.Bj + b.Bi).ToString();
            }
            return "L" + b.Level + "_" + b.Bi + "_" + b.Bj + "_" + b.Bk;
        }

        public long CellCount { get { return (long)Blocks.Count * BlockSize * BlockSize * BlockSize; } }

        /// <summary>Block holding the cell (gi, gj, gk) of the level and the cell index in it, at the same, a coarser or a finer level.</summary>
        public OctBlock Locate(int level, int gi, int gj, int gk, out int index)
        {
            index = 0;
            if (gi < 0 || gj < 0 || gk < 0)
            {
                return null;
            }
            int bs = BlockSize;
            for (int l = level; l >= 0; l--)
            {
                int s = level - l;
                int ci = gi >> s, cj = gj >> s, ck = gk >> s;
                var b = Find(l, ci / bs, cj / bs, ck / bs);
                if (b != null)
                {
                    index = ci % bs + bs * (cj % bs + bs * (ck % bs));
                    return b;
                }
            }
            for (int l = level + 1; l <= MaxLevel; l++)
            {
                int s = l - level;
                long ci = (long)gi << s, cj = (long)gj << s, ck = (long)gk << s;
                var b = Find(l, (int)(ci / bs), (int)(cj / bs), (int)(ck / bs));
                if (b != null)
                {
                    index = (int)(ci % bs + bs * (cj % bs + bs * (ck % bs)));
                    return b;
                }
            }
            return null;
        }

        /// <summary>Cell type at the cell of the level (wall outside of the blocks).</summary>
        public byte TypeAt(int level, int gi, int gj, int gk)
        {
            int index;
            var b = Locate(level, gi, gj, gk, out index);
            return b == null ? CaseGenerator.Wall : b.Types[index];
        }

        public void Save(string path)
        {
            using (var w = new BinaryWriter(File.Create(path)))
            {
                w.Write(new[] { (byte)'L', (byte)'B', (byte)'M', (byte)'O' });
                w.Write(1);
                w.Write(BlockSize);
                w.Write(NBi);
                w.Write(NBj);
                w.Write(NBk);
                w.Write(Blocks.Count);
                foreach (var b in Blocks)
                {
                    w.Write(b.Bi);
                    w.Write(b.Bj);
                    w.Write(b.Bk);
                    w.Write(b.Level);
                    w.Write(b.Types);
                }
            }
        }

        public static OctreeMesh Load(string path)
        {
            using (var r = new BinaryReader(File.OpenRead(path)))
            {
                if (new string(r.ReadChars(4)) != "LBMO" || r.ReadInt32() != 1)
                {
                    throw new InvalidDataException("Neplatný soubor " + path);
                }
                var o = new OctreeMesh { BlockSize = r.ReadInt32(), NBi = r.ReadInt32(), NBj = r.ReadInt32(), NBk = r.ReadInt32() };
                int n = r.ReadInt32();
                int cells = o.BlockSize * o.BlockSize * o.BlockSize;
                for (int i = 0; i < n; i++)
                {
                    var b = new OctBlock { Bi = r.ReadInt32(), Bj = r.ReadInt32(), Bk = r.ReadInt32(), Level = r.ReadInt32() };
                    b.Types = r.ReadBytes(cells);
                    o.Blocks.Add(b);
                }
                o.BuildIndex();
                return o;
            }
        }
    }

    /// <summary>Uniform grid of triangles for ray queries.</summary>
    class TriangleGrid
    {
        readonly StlMesh mesh;
        readonly Vec3 lo;
        readonly double h;
        readonly int nx, ny, nz;
        readonly int[][] cells;
        public readonly Vec3[] Normal;   // unit normals (zero for degenerate triangles)

        public TriangleGrid(StlMesh mesh, Vec3 min, Vec3 max, double cellSize)
        {
            this.mesh = mesh;
            Vec3 pad = new Vec3(cellSize, cellSize, cellSize) * 0.5;
            lo = min - pad;
            Vec3 size = max + pad - lo;
            h = cellSize;
            nx = Math.Max(1, (int)Math.Ceiling(size.X / h));
            ny = Math.Max(1, (int)Math.Ceiling(size.Y / h));
            nz = Math.Max(1, (int)Math.Ceiling(size.Z / h));
            var lists = new List<int>[(long)nx * ny * nz];
            int nt = mesh.TriangleCount;
            Normal = new Vec3[nt];
            for (int t = 0; t < nt; t++)
            {
                Vec3 a = mesh.Vertex(t, 0), b = mesh.Vertex(t, 1), c = mesh.Vertex(t, 2);
                Vec3 n = Vec3.Cross(b - a, c - a);
                double len = n.Length;
                Normal[t] = len > 0 ? n * (1 / len) : new Vec3(0, 0, 0);
                int i0 = Clamp((int)Math.Floor((Math.Min(a.X, Math.Min(b.X, c.X)) - lo.X) / h), nx);
                int i1 = Clamp((int)Math.Floor((Math.Max(a.X, Math.Max(b.X, c.X)) - lo.X) / h), nx);
                int j0 = Clamp((int)Math.Floor((Math.Min(a.Y, Math.Min(b.Y, c.Y)) - lo.Y) / h), ny);
                int j1 = Clamp((int)Math.Floor((Math.Max(a.Y, Math.Max(b.Y, c.Y)) - lo.Y) / h), ny);
                int k0 = Clamp((int)Math.Floor((Math.Min(a.Z, Math.Min(b.Z, c.Z)) - lo.Z) / h), nz);
                int k1 = Clamp((int)Math.Floor((Math.Max(a.Z, Math.Max(b.Z, c.Z)) - lo.Z) / h), nz);
                for (int k = k0; k <= k1; k++)
                    for (int j = j0; j <= j1; j++)
                        for (int i = i0; i <= i1; i++)
                        {
                            long id = i + (long)nx * (j + (long)ny * k);
                            if (lists[id] == null)
                            {
                                lists[id] = new List<int>();
                            }
                            lists[id].Add(t);
                        }
            }
            cells = new int[lists.Length][];
            for (long i = 0; i < lists.Length; i++)
            {
                cells[i] = lists[i] == null ? null : lists[i].ToArray();
            }
        }

        static int Clamp(int v, int n)
        {
            return v < 0 ? 0 : (v >= n ? n - 1 : v);
        }

        /// <summary>Distance to the nearest triangle hit by the ray (unit direction d) beyond tMin, infinity if none.</summary>
        public double Cast(Vec3 p, Vec3 d, double tMin, double tMax, int skip)
        {
            double best = double.PositiveInfinity;
            // 3D DDA through the cells
            double fx = (p.X - lo.X) / h, fy = (p.Y - lo.Y) / h, fz = (p.Z - lo.Z) / h;
            int i = (int)Math.Floor(fx), j = (int)Math.Floor(fy), k = (int)Math.Floor(fz);
            int si = d.X > 0 ? 1 : -1, sj = d.Y > 0 ? 1 : -1, sk = d.Z > 0 ? 1 : -1;
            double tdx = d.X != 0 ? h / Math.Abs(d.X) : double.PositiveInfinity;
            double tdy = d.Y != 0 ? h / Math.Abs(d.Y) : double.PositiveInfinity;
            double tdz = d.Z != 0 ? h / Math.Abs(d.Z) : double.PositiveInfinity;
            double tx = d.X != 0 ? ((si > 0 ? i + 1 - fx : fx - i) * h) / Math.Abs(d.X) : double.PositiveInfinity;
            double ty = d.Y != 0 ? ((sj > 0 ? j + 1 - fy : fy - j) * h) / Math.Abs(d.Y) : double.PositiveInfinity;
            double tz = d.Z != 0 ? ((sk > 0 ? k + 1 - fz : fz - k) * h) / Math.Abs(d.Z) : double.PositiveInfinity;
            double tCell = 0;
            while (tCell <= Math.Min(best, tMax))
            {
                if (i >= 0 && j >= 0 && k >= 0 && i < nx && j < ny && k < nz)
                {
                    var list = cells[i + (long)nx * (j + (long)ny * k)];
                    if (list != null)
                    {
                        foreach (int t in list)
                        {
                            if (t == skip)
                            {
                                continue;
                            }
                            double tt = Intersect(t, p, d);
                            if (tt > tMin && tt < best)
                            {
                                best = tt;
                            }
                        }
                    }
                }
                else if ((si > 0 ? i >= nx : i < 0) || (sj > 0 ? j >= ny : j < 0) || (sk > 0 ? k >= nz : k < 0))
                {
                    break; // left the grid
                }
                if (tx <= ty && tx <= tz)
                {
                    tCell = tx;
                    tx += tdx;
                    i += si;
                }
                else if (ty <= tz)
                {
                    tCell = ty;
                    ty += tdy;
                    j += sj;
                }
                else
                {
                    tCell = tz;
                    tz += tdz;
                    k += sk;
                }
            }
            return best <= tMax ? best : double.PositiveInfinity;
        }

        /// <summary>Moeller-Trumbore, ray parameter of the hit or -1.</summary>
        double Intersect(int t, Vec3 p, Vec3 d)
        {
            Vec3 a = mesh.Vertex(t, 0), b = mesh.Vertex(t, 1), c = mesh.Vertex(t, 2);
            Vec3 e1 = b - a, e2 = c - a;
            Vec3 pv = Vec3.Cross(d, e2);
            double det = Dot(e1, pv);
            if (Math.Abs(det) < 1e-30)
            {
                return -1;
            }
            double inv = 1 / det;
            Vec3 tv = p - a;
            double u = Dot(tv, pv) * inv;
            if (u < 0 || u > 1)
            {
                return -1;
            }
            Vec3 qv = Vec3.Cross(tv, e1);
            double v = Dot(d, qv) * inv;
            if (v < 0 || u + v > 1)
            {
                return -1;
            }
            return Dot(e2, qv) * inv;
        }

        public static double Dot(Vec3 a, Vec3 b)
        {
            return a.X * b.X + a.Y * b.Y + a.Z * b.Z;
        }
    }

    /// <summary>Settings of the refinement driven by the geometry.</summary>
    public class RefinementSpec
    {
        public int MaxLevel = 2;
        /// <summary>minimum number of cells across a gap of the fluid (0 - off)</summary>
        public double GapCells = 10;
        /// <summary>minimum number of cells per radius of curvature of the surface (0 - off)</summary>
        public double CurvatureCells = 0;
        /// <summary>dihedral angles above this value are sharp edges, not curvature [deg]</summary>
        public double SharpAngle = 60;
        /// <summary>minimum level of the blocks containing the surface</summary>
        public int WallLevel = 0;
    }

    /// <summary>Builds the octree of blocks from the geometry criteria and voxelizes its blocks.</summary>
    public static class OctreeBuilder
    {
        /// <param name="baseTypes">types of the level 0 grid nbi*bs x nbj*bs x nbk*bs (cells outside of the voxel grid are walls)</param>
        /// <param name="vg">voxel grid of the level 0 (origin, dimensions, outer layer)</param>
        public static OctreeMesh Build(GeometryDef def, StlMesh mesh, VoxelGrid vg, byte[] baseTypes, int nbi, int nbj, int nbk, RefinementSpec spec, Action<string> log)
        {
            int bs = def.BlockSize;
            int L = Math.Max(0, Math.Min(spec.MaxLevel, 6));
            double dx = def.Dx;
            Vec3 origin = vg.Origin;
            double finestBlock = bs * dx / (1 << L);

            // required levels at the blocks of the finest level
            var req = new Dictionary<long, int>();
            Action<Vec3, int> mark = (p, level) =>
            {
                if (level <= 0)
                {
                    return;
                }
                int ci = (int)Math.Floor((p.X - origin.X) / finestBlock), cj = (int)Math.Floor((p.Y - origin.Y) / finestBlock), ck = (int)Math.Floor((p.Z - origin.Z) / finestBlock);
                if (ci < 0 || cj < 0 || ck < 0 || ci >= (nbi << L) || cj >= (nbj << L) || ck >= (nbk << L))
                {
                    return;
                }
                long key = OctreeMesh.Key(L, ci, cj, ck);
                lock (req)
                {
                    int old;
                    if (!req.TryGetValue(key, out old) || old < level)
                    {
                        req[key] = level;
                    }
                }
            };
            Func<double, int> levelFor = dxReq =>
            {
                if (dxReq <= 0)
                {
                    return L;
                }
                int l = (int)Math.Ceiling(Math.Log(dx / dxReq, 2) - 1e-9);
                return Math.Max(0, Math.Min(L, l));
            };

            Vec3 domMax = origin + new Vec3(nbi * bs, nbj * bs, nbk * bs) * dx;
            var grid = new TriangleGrid(mesh, origin, domMax, Math.Max(finestBlock, 2 * dx));
            int nt = mesh.TriangleCount;
            // outward normals for a positive signed volume, the fluid lies against them inside the body
            double sign = mesh.SignedVolume() >= 0 ? 1 : -1;
            double toFluid = def.Region == FluidRegion.InsideBody ? -sign : sign;
            double sample = finestBlock / 2;
            double eps = 1e-4 * dx / (1 << L);
            long rays = 0;

            if (spec.GapCells > 0 || spec.WallLevel > 0)
            {
                Parallel.For(0, nt, t =>
                {
                    Vec3 n = grid.Normal[t];
                    if (n.Length == 0)
                    {
                        return;
                    }
                    Vec3 a = mesh.Vertex(t, 0), b = mesh.Vertex(t, 1), c = mesh.Vertex(t, 2);
                    double longest = Math.Max((b - a).Length, Math.Max((c - a).Length, (c - b).Length));
                    int ns = Math.Max(1, (int)Math.Ceiling(longest / sample));
                    Vec3 d = n * toFluid;
                    double maxGap = spec.GapCells > 0 ? spec.GapCells * dx : 0; // wider gaps need no refinement
                    // barycentric grid of points inside the triangle, spacing about sample
                    for (int iu = 0; iu < ns; iu++)
                    {
                        for (int iv = 0; iu + iv < ns; iv++)
                        {
                            double u = (iu + 1.0 / 3) / ns, v = (iv + 1.0 / 3) / ns;
                            Vec3 p = a + (b - a) * u + (c - a) * v;
                            mark(p, spec.WallLevel);
                            if (maxGap <= 0)
                            {
                                continue;
                            }
                            System.Threading.Interlocked.Increment(ref rays);
                            double gap = grid.Cast(p + d * eps, d, eps, maxGap, t);
                            if (double.IsInfinity(gap))
                            {
                                continue;
                            }
                            int level = levelFor(gap / spec.GapCells);
                            if (level == 0)
                            {
                                continue;
                            }
                            // the whole gap between the walls gets the level
                            int steps = Math.Max(1, (int)Math.Ceiling(gap / sample));
                            for (int s = 0; s <= steps; s++)
                            {
                                mark(p + d * (gap * s / steps), level);
                            }
                        }
                    }
                });
                if (spec.GapCells > 0)
                {
                    log("Zahušťování: " + rays + " paprsků přes mezery");
                }
            }

            if (spec.CurvatureCells > 0)
            {
                MarkCurvature(mesh, grid, spec, dx, sample, levelFor, mark, log);
            }

            // required levels of the coarser blocks
            var need = new Dictionary<long, int>(req);
            var current = req;
            for (int l = L - 1; l >= 0; l--)
            {
                var next = new Dictionary<long, int>();
                foreach (var kv in current)
                {
                    long key = kv.Key;
                    int ci = (int)((key >> 40) & 0xFFFFF), cj = (int)((key >> 20) & 0xFFFFF), ck = (int)(key & 0xFFFFF);
                    long pk = OctreeMesh.Key(l, ci >> 1, cj >> 1, ck >> 1);
                    int old;
                    if (!next.TryGetValue(pk, out old) || old < kv.Value)
                    {
                        next[pk] = kv.Value;
                    }
                }
                foreach (var kv in next)
                {
                    need[kv.Key] = kv.Value;
                }
                current = next;
            }

            // leaves: refine the blocks of the level 0 (with fluid or with a required level) recursively
            var leaves = new HashSet<long>();
            Func<long, int> needOf = key => { int v; return need.TryGetValue(key, out v) ? v : 0; };
            Action<int, int, int, int> add = null;
            add = (l, i, j, k) =>
            {
                if (l < L && needOf(OctreeMesh.Key(l, i, j, k)) > l)
                {
                    for (int o = 0; o < 8; o++)
                    {
                        add(l + 1, 2 * i + (o & 1), 2 * j + ((o >> 1) & 1), 2 * k + (o >> 2));
                    }
                }
                else
                {
                    leaves.Add(OctreeMesh.Key(l, i, j, k));
                }
            };
            int GX = nbi * bs, GY = nbj * bs;
            for (int bk = 0; bk < nbk; bk++)
                for (int bj = 0; bj < nbj; bj++)
                    for (int bi = 0; bi < nbi; bi++)
                    {
                        bool keep = needOf(OctreeMesh.Key(0, bi, bj, bk)) > 0;
                        for (int k = 0; k < bs && !keep; k++)
                            for (int j = 0; j < bs && !keep; j++)
                                for (int i = 0; i < bs && !keep; i++)
                                {
                                    keep = baseTypes[bi * bs + i + (long)GX * (bj * bs + j + (long)GY * (bk * bs + k))] != CaseGenerator.Wall;
                                }
                        if (keep)
                        {
                            add(0, bi, bj, bk);
                        }
                    }

            Balance(leaves, nbi, nbj, nbk);

            // voxelization of the blocks at their levels
            var columns = new ColumnIndex(mesh, origin, domMax, finestBlock);
            var list = leaves.ToList();
            var blocks = new OctBlock[list.Count];
            Parallel.For(0, list.Count, n =>
            {
                long key = list[n];
                int l = (int)(key >> 60), bi = (int)((key >> 40) & 0xFFFFF), bj = (int)((key >> 20) & 0xFFFFF), bk = (int)(key & 0xFFFFF);
                blocks[n] = VoxelizeBlock(def, vg, columns, l, bi, bj, bk);
            });
            var tree = new OctreeMesh { BlockSize = bs };
            foreach (var b in blocks)
            {
                if (b.Types.Any(t => t != CaseGenerator.Wall))
                {
                    tree.Blocks.Add(b);
                }
            }
            if (tree.Blocks.Count == 0)
            {
                throw new InvalidOperationException("Síť neobsahuje žádnou buňku tekutiny. Zkontrolujte jednotky, velikost buňky a uzavřenost geometrie.");
            }
            // block grid of the level 0 cropped to the used blocks (the solver derives it from the largest indexes)
            tree.NBi = tree.Blocks.Max(b => (b.Bi >> b.Level) + 1);
            tree.NBj = tree.Blocks.Max(b => (b.Bj >> b.Level) + 1);
            tree.NBk = tree.Blocks.Max(b => (b.Bk >> b.Level) + 1);
            tree.Blocks.Sort((x, y) => x.Level != y.Level ? x.Level.CompareTo(y.Level) : (x.Bk != y.Bk ? x.Bk.CompareTo(y.Bk) : (x.Bj != y.Bj ? x.Bj.CompareTo(y.Bj) : x.Bi.CompareTo(y.Bi))));
            tree.BuildIndex();
            return tree;
        }

        /// <summary>Radius of curvature from the dihedral angle of neighbouring triangles, R = d / (2 sin(theta / 2)).</summary>
        static void MarkCurvature(StlMesh mesh, TriangleGrid grid, RefinementSpec spec, double dx, double sample, Func<double, int> levelFor, Action<Vec3, int> mark, Action<string> log)
        {
            int nt = mesh.TriangleCount;
            double q = Math.Max(1e-12, dx * 1e-6);
            Func<Vec3, long> vkey = v => ((long)Math.Round(v.X / q) * 73856093L) ^ ((long)Math.Round(v.Y / q) * 19349663L) ^ ((long)Math.Round(v.Z / q) * 83492791L);
            var edges = new Dictionary<Tuple<long, long>, int>();
            long marked = 0;
            double sharp = spec.SharpAngle * Math.PI / 180;
            for (int t = 0; t < nt; t++)
            {
                if (grid.Normal[t].Length == 0)
                {
                    continue;
                }
                for (int e = 0; e < 3; e++)
                {
                    long ka = vkey(mesh.Vertex(t, e)), kb = vkey(mesh.Vertex(t, (e + 1) % 3));
                    var key = ka < kb ? Tuple.Create(ka, kb) : Tuple.Create(kb, ka);
                    int other;
                    if (!edges.TryGetValue(key, out other))
                    {
                        edges[key] = t;
                        continue;
                    }
                    edges.Remove(key);
                    Vec3 n1 = grid.Normal[t], n2 = grid.Normal[other];
                    double cos = Math.Max(-1, Math.Min(1, TriangleGrid.Dot(n1, n2)));
                    double theta = Math.Acos(cos);
                    if (theta < 1e-3 || theta > sharp)
                    {
                        continue;
                    }
                    Vec3 c1 = Centroid(mesh, t), c2 = Centroid(mesh, other);
                    double r = (c1 - c2).Length / (2 * Math.Sin(theta / 2));
                    int level = levelFor(r / spec.CurvatureCells);
                    if (level == 0)
                    {
                        continue;
                    }
                    marked++;
                    // the neighbourhood of the curved surface up to the radius gets the level
                    Vec3 m = (c1 + c2) * 0.5;
                    double ext = Math.Min(r, 4 * sample);
                    int ns = (int)Math.Ceiling(ext / sample);
                    for (int a = -ns; a <= ns; a++)
                        for (int b = -ns; b <= ns; b++)
                            for (int c = -ns; c <= ns; c++)
                            {
                                mark(m + new Vec3(a, b, c) * sample, level);
                            }
                }
            }
            log("Zahušťování: " + marked + " zakřivených hran");
        }

        static Vec3 Centroid(StlMesh mesh, int t)
        {
            return (mesh.Vertex(t, 0) + mesh.Vertex(t, 1) + mesh.Vertex(t, 2)) * (1.0 / 3);
        }

        /// <summary>Splits the leaves until the neighbouring leaves (26-neighbourhood) differ by one level at most.</summary>
        static void Balance(HashSet<long> leaves, int nbi, int nbj, int nbk)
        {
            Func<int, int, int, int, long> cover = (l, i, j, k) =>
            {
                for (int lc = l; lc >= 0; lc--)
                {
                    int s = l - lc;
                    long key = OctreeMesh.Key(lc, i >> s, j >> s, k >> s);
                    if (leaves.Contains(key))
                    {
                        return key;
                    }
                }
                return -1;
            };
            bool changed = true;
            while (changed)
            {
                changed = false;
                foreach (long key in leaves.ToList())
                {
                    if (!leaves.Contains(key))
                    {
                        continue;
                    }
                    int l = (int)(key >> 60), bi = (int)((key >> 40) & 0xFFFFF), bj = (int)((key >> 20) & 0xFFFFF), bk = (int)(key & 0xFFFFF);
                    if (l < 2)
                    {
                        continue;
                    }
                    for (int dk = -1; dk <= 1; dk++)
                        for (int dj = -1; dj <= 1; dj++)
                            for (int di = -1; di <= 1; di++)
                            {
                                int i = bi + di, j = bj + dj, k = bk + dk;
                                if (i < 0 || j < 0 || k < 0 || i >= (nbi << l) || j >= (nbj << l) || k >= (nbk << l))
                                {
                                    continue;
                                }
                                long c = cover(l, i, j, k);
                                if (c < 0)
                                {
                                    continue;
                                }
                                int lc = (int)(c >> 60);
                                if (lc >= l - 1)
                                {
                                    continue;
                                }
                                // too coarse neighbour: split it
                                int ci = (int)((c >> 40) & 0xFFFFF), cj = (int)((c >> 20) & 0xFFFFF), ck = (int)(c & 0xFFFFF);
                                leaves.Remove(c);
                                for (int o = 0; o < 8; o++)
                                {
                                    leaves.Add(OctreeMesh.Key(lc + 1, 2 * ci + (o & 1), 2 * cj + ((o >> 1) & 1), 2 * ck + (o >> 2)));
                                }
                                changed = true;
                            }
                }
            }
        }

        /// <summary>Cell types of a block at its level (inside/outside by z rays, region and outer walls as for the level 0).</summary>
        static OctBlock VoxelizeBlock(GeometryDef def, VoxelGrid vg, ColumnIndex columns, int l, int bi, int bj, int bk)
        {
            int bs = def.BlockSize;
            double h = def.Dx / (1 << l);
            var b = new OctBlock { Bi = bi, Bj = bj, Bk = bk, Level = l, Types = new byte[bs * bs * bs] };
            double x0 = vg.Origin.X + bi * bs * h, y0 = vg.Origin.Y + bj * bs * h, z0 = vg.Origin.Z + bk * bs * h;
            double sx = 0.5 + 1.3e-4, sy = 0.5 + 2.7e-4;
            var zs = new List<double>();
            for (int j = 0; j < bs; j++)
            {
                for (int i = 0; i < bs; i++)
                {
                    double x = x0 + (i + sx) * h, y = y0 + (j + sy) * h;
                    columns.Hits(x, y, zs);
                    for (int k = 0; k < bs; k++)
                    {
                        double z = z0 + (k + 0.5) * h;
                        int below = 0;
                        foreach (double zz in zs)
                        {
                            if (zz < z)
                            {
                                below++;
                            }
                        }
                        bool inside = (below & 1) == 1;
                        // cell of the level 0 holding this cell
                        int s = l;
                        int i0 = (bi * bs + i) >> s, j0 = (bj * bs + j) >> s, k0 = (bk * bs + k) >> s;
                        bool f = false;
                        if (i0 < vg.Nx && j0 < vg.Ny && k0 < vg.Nz)
                        {
                            if (def.Region == FluidRegion.InsideBody)
                            {
                                f = inside;
                            }
                            else
                            {
                                bool outer = i0 == 0 || j0 == 0 || k0 == 0 || i0 == vg.Nx - 1 || j0 == vg.Ny - 1 || k0 == vg.Nz - 1;
                                f = !inside && !outer;
                            }
                        }
                        b.Types[i + bs * (j + bs * k)] = f ? CaseGenerator.Fluid : CaseGenerator.Wall;
                    }
                }
            }
            return b;
        }

        /// <summary>Marks inlet/outlet cells of the octree like CaseGenerator.MarkPort marks them on the level 0 grid.</summary>
        public static long MarkPort(GeometryDef def, PortSpec port, OctreeMesh tree, byte value)
        {
            int axis = PortSpec.Axis(port.Face);
            int sign = PortSpec.InwardSign(port.Face);
            int a1 = axis == 0 ? 1 : 0, a2 = axis == 2 ? 1 : 2;
            int[] n0 = { def.Nx, def.Ny, def.Nz };
            int bs = tree.BlockSize;
            int maxDepth0 = Math.Max(4, (int)Math.Ceiling(def.Margin / def.Dx) + 4);
            long marked = 0;
            var toMark = new List<KeyValuePair<OctBlock, int>>();
            foreach (var b in tree.Blocks)
            {
                int l = b.Level;
                double h = def.Dx / (1 << l);
                int nAxis = n0[axis] << l, maxDepth = maxDepth0 << l;
                int[] bc = { b.Bi, b.Bj, b.Bk };
                for (int k = 0; k < bs; k++)
                    for (int j = 0; j < bs; j++)
                        for (int i = 0; i < bs; i++)
                        {
                            int index = i + bs * (j + bs * k);
                            if (b.Types[index] != CaseGenerator.Wall)
                            {
                                continue;
                            }
                            int[] g = { bc[0] * bs + i, bc[1] * bs + j, bc[2] * bs + k };
                            if (g[a1] >= (n0[a1] << l) || g[a2] >= (n0[a2] << l))
                            {
                                continue;
                            }
                            int depth = sign > 0 ? g[axis] : nAxis - 1 - g[axis];
                            if (depth < 0 || depth + 1 >= maxDepth || depth + 1 >= nAxis)
                            {
                                continue;
                            }
                            if (port.Radius > 0)
                            {
                                double c1 = def.Origin[a1] + (g[a1] + 0.5) * h - port.U, c2 = def.Origin[a2] + (g[a2] + 0.5) * h - port.V;
                                if (c1 * c1 + c2 * c2 > port.Radius * port.Radius)
                                {
                                    continue;
                                }
                            }
                            // the next cell inwards is fluid and no fluid lies between the face and this cell
                            int[] c = (int[])g.Clone();
                            c[axis] += sign;
                            if (tree.TypeAt(l, c[0], c[1], c[2]) != CaseGenerator.Fluid)
                            {
                                continue;
                            }
                            bool clear = true;
                            for (int d = depth - 1; d >= 0 && clear; d--)
                            {
                                c[axis] = sign > 0 ? d : nAxis - 1 - d;
                                clear = tree.TypeAt(l, c[0], c[1], c[2]) != CaseGenerator.Fluid;
                            }
                            if (clear)
                            {
                                toMark.Add(new KeyValuePair<OctBlock, int>(b, index));
                            }
                        }
            }
            foreach (var kv in toMark)
            {
                if (kv.Key.Types[kv.Value] == CaseGenerator.Wall)
                {
                    kv.Key.Types[kv.Value] = value;
                    marked++;
                }
            }
            return marked;
        }
    }

    /// <summary>Triangles bucketed by xy columns, intersections of vertical lines with the surface.</summary>
    class ColumnIndex
    {
        readonly StlMesh mesh;
        readonly Vec3 lo;
        readonly double h;
        readonly int nx, ny;
        readonly int[][] cols;

        public ColumnIndex(StlMesh mesh, Vec3 min, Vec3 max, double cellSize)
        {
            this.mesh = mesh;
            lo = min;
            h = cellSize;
            nx = Math.Max(1, (int)Math.Ceiling((max.X - min.X) / h) + 1);
            ny = Math.Max(1, (int)Math.Ceiling((max.Y - min.Y) / h) + 1);
            var lists = new List<int>[nx * ny];
            for (int t = 0; t < mesh.TriangleCount; t++)
            {
                Vec3 a = mesh.Vertex(t, 0), b = mesh.Vertex(t, 1), c = mesh.Vertex(t, 2);
                int i0 = Clamp((int)Math.Floor((Math.Min(a.X, Math.Min(b.X, c.X)) - lo.X) / h), nx);
                int i1 = Clamp((int)Math.Floor((Math.Max(a.X, Math.Max(b.X, c.X)) - lo.X) / h), nx);
                int j0 = Clamp((int)Math.Floor((Math.Min(a.Y, Math.Min(b.Y, c.Y)) - lo.Y) / h), ny);
                int j1 = Clamp((int)Math.Floor((Math.Max(a.Y, Math.Max(b.Y, c.Y)) - lo.Y) / h), ny);
                for (int j = j0; j <= j1; j++)
                    for (int i = i0; i <= i1; i++)
                    {
                        if (lists[i + nx * j] == null)
                        {
                            lists[i + nx * j] = new List<int>();
                        }
                        lists[i + nx * j].Add(t);
                    }
            }
            cols = lists.Select(x => x == null ? null : x.ToArray()).ToArray();
        }

        static int Clamp(int v, int n)
        {
            return v < 0 ? 0 : (v >= n ? n - 1 : v);
        }

        /// <summary>z coordinates where the vertical line through (x, y) crosses the surface, sorted.</summary>
        public void Hits(double x, double y, List<double> zs)
        {
            zs.Clear();
            int i = (int)Math.Floor((x - lo.X) / h), j = (int)Math.Floor((y - lo.Y) / h);
            if (i < 0 || j < 0 || i >= nx || j >= ny || cols[i + nx * j] == null)
            {
                return;
            }
            foreach (int t in cols[i + nx * j])
            {
                Vec3 a = mesh.Vertex(t, 0), b = mesh.Vertex(t, 1), c = mesh.Vertex(t, 2);
                double d = (b.X - a.X) * (c.Y - a.Y) - (c.X - a.X) * (b.Y - a.Y);
                if (Math.Abs(d) < 1e-30)
                {
                    continue;
                }
                double l1 = ((x - a.X) * (c.Y - a.Y) - (c.X - a.X) * (y - a.Y)) / d;
                double l2 = ((b.X - a.X) * (y - a.Y) - (x - a.X) * (b.Y - a.Y)) / d;
                double l0 = 1 - l1 - l2;
                if (l0 < 0 || l1 < 0 || l2 < 0)
                {
                    continue;
                }
                zs.Add(l0 * a.Z + l1 * b.Z + l2 * c.Z);
            }
            zs.Sort();
        }
    }
}
