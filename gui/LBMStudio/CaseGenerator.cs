using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;

namespace LBMStudio
{
    public enum Face { XMin, XMax, YMin, YMax, ZMin, ZMax }

    public enum FluidRegion { InsideBody, OutsideBody }

    /// <summary>Inlet or outlet: circular (or whole) part of a face of the bounding box.</summary>
    public class PortSpec
    {
        public bool Enabled;
        public Face Face = Face.ZMax;
        /// <summary>centre in the face plane [m] (y,z for X faces, x,z for Y faces, x,y for Z faces)</summary>
        public double U, V;
        /// <summary>radius [m], 0 - the whole face</summary>
        public double Radius;
        /// <summary>inlet velocity [m/s] or outlet pressure [Pa]</summary>
        public double Value;

        public static int Axis(Face f) { return (int)f / 2; }
        public static int InwardSign(Face f) { return ((int)f % 2 == 0) ? 1 : -1; }

        public static Vec3 InwardNormal(Face f)
        {
            double s = InwardSign(f);
            switch (Axis(f))
            {
                case 0: return new Vec3(s, 0, 0);
                case 1: return new Vec3(0, s, 0);
                default: return new Vec3(0, 0, s);
            }
        }

        public void Save(KeyValueFile kv, string prefix)
        {
            kv.Set(prefix + "Enabled", Enabled ? 1 : 0);
            kv.Set(prefix + "Face", Face.ToString());
            kv.Set(prefix + "U", U);
            kv.Set(prefix + "V", V);
            kv.Set(prefix + "Radius", Radius);
            kv.Set(prefix + "Value", Value);
        }

        public static PortSpec Load(KeyValueFile kv, string prefix)
        {
            var p = new PortSpec
            {
                Enabled = kv.GetInt(prefix + "Enabled") != 0,
                U = kv.GetDouble(prefix + "U"),
                V = kv.GetDouble(prefix + "V"),
                Radius = kv.GetDouble(prefix + "Radius"),
                Value = kv.GetDouble(prefix + "Value")
            };
            Face f;
            if (Enum.TryParse(kv.Get(prefix + "Face", "ZMax"), out f))
            {
                p.Face = f;
            }
            return p;
        }
    }

    /// <summary>Geometry and mesh definition, stored in simulations/&lt;name&gt;/geometry.ini.</summary>
    public class GeometryDef
    {
        public string Name = "";
        public string SourceFile = "";
        public double UnitScale = 0.001;
        public double Dx = 0.002;
        public int BlockSize = 16;
        public FluidRegion Region = FluidRegion.InsideBody;
        public double Margin = 0.0;
        public PortSpec Inlet = new PortSpec { Enabled = true, Face = Face.ZMax, Value = 0.5 };
        public PortSpec Outlet = new PortSpec { Enabled = false, Face = Face.ZMin };
        public Vec3 Gravity = new Vec3(0, 0, -9.81);

        // results of the mesh generation
        public bool MeshReady;
        public int Nx, Ny, Nz, NBi, NBj, NBk, BlockCount;
        public Vec3 Origin;
        public long FluidCells, InletCells, OutletCells;

        public string Dir(string workspace) { return Path.Combine(workspace, Name); }
        public string StlPath(string workspace) { return Path.Combine(Dir(workspace), "geometry.stl"); }

        public void Save(string workspace)
        {
            Directory.CreateDirectory(Dir(workspace));
            var kv = new KeyValueFile();
            kv.Set("name", Name);
            kv.Set("sourceFile", SourceFile);
            kv.Set("unitScale", UnitScale);
            kv.Set("dx", Dx);
            kv.Set("blockSize", BlockSize);
            kv.Set("region", Region.ToString());
            kv.Set("margin", Margin);
            Inlet.Save(kv, "inlet");
            Outlet.Save(kv, "outlet");
            kv.Set("gx", Gravity.X);
            kv.Set("gy", Gravity.Y);
            kv.Set("gz", Gravity.Z);
            kv.Set("meshReady", MeshReady ? 1 : 0);
            kv.Set("nx", Nx);
            kv.Set("ny", Ny);
            kv.Set("nz", Nz);
            kv.Set("nbi", NBi);
            kv.Set("nbj", NBj);
            kv.Set("nbk", NBk);
            kv.Set("blockCount", BlockCount);
            kv.Set("originX", Origin.X);
            kv.Set("originY", Origin.Y);
            kv.Set("originZ", Origin.Z);
            kv.Set("fluidCells", FluidCells.ToString(CultureInfo.InvariantCulture));
            kv.Set("inletCells", InletCells.ToString(CultureInfo.InvariantCulture));
            kv.Set("outletCells", OutletCells.ToString(CultureInfo.InvariantCulture));
            kv.Save(Path.Combine(Dir(workspace), "geometry.ini"));
        }

        public static GeometryDef Load(string dir)
        {
            var kv = KeyValueFile.Load(Path.Combine(dir, "geometry.ini"));
            var g = new GeometryDef
            {
                Name = Path.GetFileName(dir),
                SourceFile = kv.Get("sourceFile"),
                UnitScale = kv.GetDouble("unitScale", 0.001),
                Dx = kv.GetDouble("dx", 0.002),
                BlockSize = kv.GetInt("blockSize", 16),
                Margin = kv.GetDouble("margin"),
                Inlet = PortSpec.Load(kv, "inlet"),
                Outlet = PortSpec.Load(kv, "outlet"),
                Gravity = new Vec3(kv.GetDouble("gx"), kv.GetDouble("gy"), kv.GetDouble("gz", -9.81)),
                MeshReady = kv.GetInt("meshReady") != 0,
                Nx = kv.GetInt("nx"),
                Ny = kv.GetInt("ny"),
                Nz = kv.GetInt("nz"),
                NBi = kv.GetInt("nbi"),
                NBj = kv.GetInt("nbj"),
                NBk = kv.GetInt("nbk"),
                BlockCount = kv.GetInt("blockCount"),
                Origin = new Vec3(kv.GetDouble("originX"), kv.GetDouble("originY"), kv.GetDouble("originZ")),
                FluidCells = (long)kv.GetDouble("fluidCells"),
                InletCells = (long)kv.GetDouble("inletCells"),
                OutletCells = (long)kv.GetDouble("outletCells")
            };
            FluidRegion r;
            if (Enum.TryParse(kv.Get("region", "InsideBody"), out r))
            {
                g.Region = r;
            }
            return g;
        }

        /// <summary>Cell types of the whole (block aligned) grid, written by the mesh generator.</summary>
        public byte[] LoadCellTypes(string workspace)
        {
            string path = Path.Combine(Dir(workspace), "mesh", "celltypes.bin");
            return File.Exists(path) ? File.ReadAllBytes(path) : null;
        }

        public int GX { get { return NBi * BlockSize; } }
        public int GY { get { return NBj * BlockSize; } }
        public int GZ { get { return NBk * BlockSize; } }
    }

    /// <summary>Creates the solver mesh: blocks.txt, boundary files and the cell type map.</summary>
    public static class CaseGenerator
    {
        public const byte Fluid = 0, Wall = 1, Inlet = 2, Outlet = 3;

        public static void Generate(GeometryDef def, StlMesh mesh, string workspace, Action<string> log)
        {
            if (def.Dx <= 0)
            {
                throw new ArgumentException("Velikost buňky musí být kladná.");
            }
            int marginCells = def.Region == FluidRegion.OutsideBody ? (int)Math.Ceiling(def.Margin / def.Dx) : 0;
            int padding = 1 + marginCells;
            log("Voxelizace (" + mesh.TriangleCount + " trojúhelníků, dx = " + Fmt.G(def.Dx) + " m) ...");
            VoxelGrid vg = Voxelizer.Voxelize(mesh, def.Dx, padding);
            int bs = def.BlockSize;
            int nbi = (vg.Nx + bs - 1) / bs, nbj = (vg.Ny + bs - 1) / bs, nbk = (vg.Nz + bs - 1) / bs;
            int GX = nbi * bs, GY = nbj * bs, GZ = nbk * bs;
            var types = new byte[(long)GX * GY * GZ];
            Func<int, int, int, long> idx = (i, j, k) => i + (long)GX * (j + (long)GY * k);
            long fluid = 0;
            for (int k = 0; k < GZ; k++)
            {
                for (int j = 0; j < GY; j++)
                {
                    for (int i = 0; i < GX; i++)
                    {
                        bool f = false;
                        if (i < vg.Nx && j < vg.Ny && k < vg.Nz)
                        {
                            bool inside = vg.Inside[vg.Index(i, j, k)] != 0;
                            if (def.Region == FluidRegion.InsideBody)
                            {
                                f = inside;
                            }
                            else
                            {
                                bool outer = i == 0 || j == 0 || k == 0 || i == vg.Nx - 1 || j == vg.Ny - 1 || k == vg.Nz - 1;
                                f = !inside && !outer;
                            }
                        }
                        types[idx(i, j, k)] = f ? Fluid : Wall;
                        if (f)
                        {
                            fluid++;
                        }
                    }
                }
            }
            if (fluid == 0)
            {
                throw new InvalidOperationException("Síť neobsahuje žádnou buňku tekutiny. Zkontrolujte jednotky, velikost buňky a uzavřenost geometrie.");
            }
            log("Buněk tekutiny: " + fluid);

            def.Nx = vg.Nx;
            def.Ny = vg.Ny;
            def.Nz = vg.Nz;
            def.Origin = vg.Origin;
            def.InletCells = def.Inlet.Enabled ? MarkPort(def, def.Inlet, types, GX, GY, GZ, Inlet) : 0;
            def.OutletCells = def.Outlet.Enabled ? MarkPort(def, def.Outlet, types, GX, GY, GZ, Outlet) : 0;
            if (def.Inlet.Enabled)
            {
                log("Buněk vtoku: " + def.InletCells);
                if (def.InletCells == 0)
                {
                    log("Varování: vtok neleží na povrchu tekutiny – zkontrolujte stěnu, střed a poloměr.");
                }
            }
            if (def.Outlet.Enabled)
            {
                log("Buněk výtoku: " + def.OutletCells);
                if (def.OutletCells == 0)
                {
                    log("Varování: výtok neleží na povrchu tekutiny.");
                }
            }
            def.FluidCells = fluid;

            string dir = def.Dir(workspace);
            string meshDir = Path.Combine(dir, "mesh");
            string bdir = Path.Combine(meshDir, "boundary");
            if (Directory.Exists(bdir))
            {
                foreach (var f in Directory.GetFiles(bdir, "blockBoundary*.txt"))
                {
                    File.Delete(f);
                }
            }
            Directory.CreateDirectory(bdir);

            // blocks containing anything else than walls; the solver derives the block grid from the largest indexes
            var kept = new List<int[]>();
            int ebi = 0, ebj = 0, ebk = 0;
            for (int bk = 0; bk < nbk; bk++)
            {
                for (int bj = 0; bj < nbj; bj++)
                {
                    for (int bi = 0; bi < nbi; bi++)
                    {
                        bool keep = false;
                        for (int k = 0; k < bs && !keep; k++)
                            for (int j = 0; j < bs && !keep; j++)
                                for (int i = 0; i < bs && !keep; i++)
                                {
                                    keep = types[idx(bi * bs + i, bj * bs + j, bk * bs + k)] != Wall;
                                }
                        if (keep)
                        {
                            kept.Add(new[] { bi, bj, bk });
                            ebi = Math.Max(ebi, bi + 1);
                            ebj = Math.Max(ebj, bj + 1);
                            ebk = Math.Max(ebk, bk + 1);
                        }
                    }
                }
            }
            // crop the grid to the used blocks (the removed part contains only walls)
            if (ebi != nbi || ebj != nbj || ebk != nbk)
            {
                int ex = ebi * bs, ey = ebj * bs, ez = ebk * bs;
                var cropped = new byte[(long)ex * ey * ez];
                for (int k = 0; k < ez; k++)
                    for (int j = 0; j < ey; j++)
                        for (int i = 0; i < ex; i++)
                        {
                            cropped[i + (long)ex * (j + (long)ey * k)] = types[idx(i, j, k)];
                        }
                types = cropped;
                nbi = ebi;
                nbj = ebj;
                nbk = ebk;
                GX = ex;
                GY = ey;
                GZ = ez;
                idx = (i, j, k) => i + (long)GX * (j + (long)GY * k);
            }
            def.NBi = nbi;
            def.NBj = nbj;
            def.NBk = nbk;
            def.Nx = Math.Min(def.Nx, GX);
            def.Ny = Math.Min(def.Ny, GY);
            def.Nz = Math.Min(def.Nz, GZ);

            var blocks = new StringBuilder();
            int count = 0;
            var sb = new StringBuilder(bs * bs * bs * 3);
            foreach (var blk in kept)
            {
                int bi = blk[0], bj = blk[1], bk = blk[2];
                bool allFluid = true;
                for (int k = 0; k < bs && allFluid; k++)
                    for (int j = 0; j < bs && allFluid; j++)
                        for (int i = 0; i < bs && allFluid; i++)
                        {
                            allFluid = types[idx(bi * bs + i, bj * bs + j, bk * bs + k)] == Fluid;
                        }
                count++;
                blocks.Append(bi).Append(' ').Append(bj).Append(' ').Append(bk).Append("\r\n");
                if (allFluid)
                {
                    continue; // no boundary file: the whole block is fluid
                }
                int id = nbi * nbj * bk + nbi * bj + bi;
                sb.Clear();
                // order of the solver: i outer, then j, k fastest
                for (int i = 0; i < bs; i++)
                    for (int j = 0; j < bs; j++)
                        for (int k = 0; k < bs; k++)
                        {
                            sb.Append((char)('0' + types[idx(bi * bs + i, bj * bs + j, bk * bs + k)])).Append('\n');
                        }
                File.WriteAllText(Path.Combine(bdir, "blockBoundary" + id + ".txt"), sb.ToString());
            }
            File.WriteAllText(Path.Combine(meshDir, "blocks.txt"), blocks.ToString());
            File.WriteAllBytes(Path.Combine(meshDir, "celltypes.bin"), types);
            def.BlockCount = count;
            def.MeshReady = true;
            log("Bloků: " + count + " z " + (nbi * nbj * nbk) + " (" + bs + "³ buněk)");
            def.Save(workspace);
        }

        /// <summary>Marks the non-fluid cells in front of the fluid surface on the face as inlet/outlet.</summary>
        static long MarkPort(GeometryDef def, PortSpec port, byte[] types, int GX, int GY, int GZ, byte value)
        {
            int axis = PortSpec.Axis(port.Face);
            int sign = PortSpec.InwardSign(port.Face);
            int[] n = { def.Nx, def.Ny, def.Nz };
            int[] G = { GX, GY, GZ };
            int a1 = axis == 0 ? 1 : 0, a2 = axis == 2 ? 1 : 2;
            double o1 = def.Origin[a1], o2 = def.Origin[a2];
            long marked = 0;
            int maxDepth = Math.Max(4, (int)Math.Ceiling(def.Margin / def.Dx) + 4);
            for (int p = 0; p < n[a1]; p++)
            {
                for (int q = 0; q < n[a2]; q++)
                {
                    if (port.Radius > 0)
                    {
                        double c1 = o1 + (p + 0.5) * def.Dx - port.U, c2 = o2 + (q + 0.5) * def.Dx - port.V;
                        if (c1 * c1 + c2 * c2 > port.Radius * port.Radius)
                        {
                            continue;
                        }
                    }
                    int start = sign > 0 ? 0 : n[axis] - 1;
                    for (int d = 1; d < maxDepth && d < n[axis]; d++)
                    {
                        int[] c = new int[3];
                        c[axis] = start + sign * d;
                        c[a1] = p;
                        c[a2] = q;
                        long id = c[0] + (long)GX * (c[1] + (long)GY * c[2]);
                        if (types[id] == Fluid)
                        {
                            c[axis] -= sign;
                            long pid = c[0] + (long)GX * (c[1] + (long)GY * c[2]);
                            if (types[pid] == Wall)
                            {
                                types[pid] = value;
                                marked++;
                            }
                            break;
                        }
                    }
                }
            }
            return marked;
        }
    }
}
