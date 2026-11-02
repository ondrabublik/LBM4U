using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace LBMStudio
{
    public enum FillMode { Empty, Full, Level }

    /// <summary>One simulation of a geometry, stored in simulations/&lt;geometry&gt;/&lt;name&gt;/simulation.ini.</summary>
    public class SimulationDef
    {
        public string Name = "sim1";
        public int ProblemId = 1;          // 0 - single phase, 1 - free surface
        public double Density = 998;       // kg/m3 (water 20 °C)
        public double Viscosity = 1.004e-6; // kinematic viscosity m2/s
        public double SurfaceTension = 0.0728; // N/m
        public double ContactAngle = 0;    // deg, 0 - not used
        public int ViscosityModel = 0;     // 0 Newtonian, 1 power law, 2 Bingham, 3 Cross, 4 polynomial
        public string FluidParameters = "";
        public Vec3 Gravity = new Vec3(0, 0, -9.81);
        public double InletVelocity = 0.5;
        public double OutletPressure = 0;
        public int CollisionModel = 1;     // 0 BGK, 1 MRT, 2 TRT
        public double Smagorinsky = 0.15;
        public double LatticeVelocity = 0.1; // lattice velocity of the characteristic velocity (Mach number * sqrt(3))
        public double CharVelocity = 0;    // 0 - automatic
        public double EndTime = 1.0;
        public double SaveInterval = 0.05;
        public int Threads = Environment.ProcessorCount;
        public FillMode Fill = FillMode.Empty;
        public double FillLevel = 0;       // height of the initial liquid level above the lowest fluid cell [m]
        public double TauMax = 1.5;

        public string Dir(GeometryDef g, string workspace) { return Path.Combine(g.Dir(workspace), Name); }

        public void Save(string dir)
        {
            Directory.CreateDirectory(dir);
            var kv = new KeyValueFile();
            kv.Set("problemID", ProblemId);
            kv.Set("density", Density);
            kv.Set("viscosity", Viscosity);
            kv.Set("surfaceTension", SurfaceTension);
            kv.Set("contactAngle", ContactAngle);
            kv.Set("viscosityModel", ViscosityModel);
            kv.Set("fluidParameters", FluidParameters);
            kv.Set("gx", Gravity.X);
            kv.Set("gy", Gravity.Y);
            kv.Set("gz", Gravity.Z);
            kv.Set("inletVelocity", InletVelocity);
            kv.Set("outletPressure", OutletPressure);
            kv.Set("collisionModel", CollisionModel);
            kv.Set("smagorinsky", Smagorinsky);
            kv.Set("latticeVelocity", LatticeVelocity);
            kv.Set("charVelocity", CharVelocity);
            kv.Set("endTime", EndTime);
            kv.Set("saveInterval", SaveInterval);
            kv.Set("threads", Threads);
            kv.Set("fill", Fill.ToString());
            kv.Set("fillLevel", FillLevel);
            kv.Set("tauMax", TauMax);
            kv.Save(Path.Combine(dir, "simulation.ini"));
        }

        public static SimulationDef Load(string dir)
        {
            var kv = KeyValueFile.Load(Path.Combine(dir, "simulation.ini"));
            var s = new SimulationDef
            {
                Name = Path.GetFileName(dir),
                ProblemId = kv.GetInt("problemID", 1),
                Density = kv.GetDouble("density", 998),
                Viscosity = kv.GetDouble("viscosity", 1.004e-6),
                SurfaceTension = kv.GetDouble("surfaceTension"),
                ContactAngle = kv.GetDouble("contactAngle"),
                ViscosityModel = kv.GetInt("viscosityModel"),
                FluidParameters = kv.Get("fluidParameters"),
                Gravity = new Vec3(kv.GetDouble("gx"), kv.GetDouble("gy"), kv.GetDouble("gz", -9.81)),
                InletVelocity = kv.GetDouble("inletVelocity", 0.5),
                OutletPressure = kv.GetDouble("outletPressure"),
                CollisionModel = kv.GetInt("collisionModel", 1),
                Smagorinsky = kv.GetDouble("smagorinsky", 0.15),
                LatticeVelocity = kv.GetDouble("latticeVelocity", 0.1),
                CharVelocity = kv.GetDouble("charVelocity"),
                EndTime = kv.GetDouble("endTime", 1),
                SaveInterval = kv.GetDouble("saveInterval", 0.05),
                Threads = kv.GetInt("threads", Environment.ProcessorCount),
                FillLevel = kv.GetDouble("fillLevel"),
                TauMax = kv.GetDouble("tauMax", 1.5)
            };
            FillMode f;
            if (Enum.TryParse(kv.Get("fill", "Empty"), out f))
            {
                s.Fill = f;
            }
            return s;
        }

        public SimulationDef Clone()
        {
            return (SimulationDef)MemberwiseClone();
        }
    }

    /// <summary>Derived lattice quantities of a simulation.</summary>
    public class Derived
    {
        public double Uc, UMaxParam, Dt, Tau, TauEffMin, Mach, MemoryGB;
        public int MaxIter, SaveRate;
        public long Cells;
        public List<string> Warnings = new List<string>();

        public static double ExtentAlong(GeometryDef g, Vec3 dir)
        {
            double l = dir.Length;
            if (l == 0)
            {
                return 0;
            }
            return (Math.Abs(dir.X) * g.Nx + Math.Abs(dir.Y) * g.Ny + Math.Abs(dir.Z) * g.Nz) / l * g.Dx;
        }

        public static Derived Compute(GeometryDef g, SimulationDef s)
        {
            var d = new Derived();
            double gAbs = s.Gravity.Length;
            double fall = Math.Sqrt(2 * gAbs * ExtentAlong(g, s.Gravity));
            double uIn = (g.Inlet.Enabled ? Math.Abs(s.InletVelocity) : 0);
            d.Uc = s.CharVelocity > 0 ? s.CharVelocity : Math.Max(Math.Max(uIn, fall), 1e-3);
            double ul = s.LatticeVelocity > 0 ? s.LatticeVelocity : 0.1;
            d.UMaxParam = d.Uc * 0.1 / ul;          // solver maps uMax to lattice velocity 0.1
            d.Dt = g.Dx * ul / d.Uc;
            d.Tau = 0.5 + 3 * s.Viscosity * d.Dt / (g.Dx * g.Dx);
            d.Mach = ul * Math.Sqrt(3);
            d.MaxIter = (int)Math.Min(int.MaxValue - 1, Math.Ceiling(s.EndTime / d.Dt));
            d.SaveRate = Math.Max(1, (int)Math.Round(s.SaveInterval / d.Dt));
            d.Cells = g.MeshCells;
            double bytesPerCell = s.ProblemId == 0 ? 2 * 19 * 8 + 120 : 2 * 19 * 8 + 260;
            d.MemoryGB = d.Cells * bytesPerCell / 1e9;
            if (d.Tau < 0.501 && s.Smagorinsky <= 0)
            {
                d.Warnings.Add("τ = " + Fmt.G(d.Tau) + " is very close to 0.5 – the computation will be unstable. Enable LES (Smagorinsky ≈ 0.15), lower the lattice velocity or refine the mesh.");
            }
            else if (d.Tau < 0.51)
            {
                d.Warnings.Add("τ = " + Fmt.G(d.Tau) + " is small, the LES model keeps it stable.");
            }
            if (d.Tau > 2)
            {
                d.Warnings.Add("τ = " + Fmt.G(d.Tau) + " is large – accuracy drops, raise the lattice velocity or refine the mesh.");
            }
            if (ul > 0.15)
            {
                d.Warnings.Add("Lattice velocity above 0.15 – compressibility errors (Ma = " + Fmt.G(d.Mach, 3) + ").");
            }
            if (s.ProblemId == 0 && g.Inlet.Enabled && !g.Outlet.Enabled)
            {
                d.Warnings.Add("A single phase computation with an inlet needs an outlet, otherwise the pressure in the closed region grows.");
            }
            if (s.ProblemId == 1 && s.Fill == FillMode.Empty && !g.Inlet.Enabled)
            {
                d.Warnings.Add("Empty region without an inlet – nothing to compute.");
            }
            if (d.MaxIter > 5000000)
            {
                d.Warnings.Add("The number of iterations " + d.MaxIter + " is very high.");
            }
            return d;
        }
    }

    /// <summary>Writes parameters.txt and the initial condition of a simulation.</summary>
    public static class SimulationWriter
    {
        public static void Write(GeometryDef g, SimulationDef s, string workspace)
        {
            if (!g.MeshReady)
            {
                throw new InvalidOperationException("The geometry has no mesh.");
            }
            string dir = s.Dir(g, workspace);
            Directory.CreateDirectory(dir);
            s.Save(dir);
            Derived d = Derived.Compute(g, s);
            var p = new KeyValueFile();
            p.Set("problemID", s.ProblemId);
            p.Set("model", 319);
            p.Set("collisionModel", s.CollisionModel);
            p.Set("threads", Math.Max(1, s.Threads));
            p.Set("maxIter", d.MaxIter);
            p.Set("saveRate", d.SaveRate);
            p.Set("ni", g.BlockSize);
            p.Set("nj", g.BlockSize);
            p.Set("nk", g.BlockSize);
            p.Set("H", g.BlockSize * g.Dx);
            p.Set("uMax", d.UMaxParam);
            p.Set("density", s.Density);
            p.Set("vKin", s.Viscosity);
            p.Set("gx", s.Gravity.X);
            p.Set("gy", s.Gravity.Y);
            p.Set("gz", s.Gravity.Z);
            p.Set("smagorinsky", s.Smagorinsky);
            p.Set("tauMax", s.TauMax);
            // limits the relaxation time of the odd moments, required for stability when tau -> 0.5 (water)
            p.Set("tauMinusMax", 1.0);
            p.Set("xPeriodicity", 0);
            p.Set("yPeriodicity", 0);
            p.Set("zPeriodicity", 0);
            p.Set("movingLidTop", 0);
            if (g.Inlet.Enabled)
            {
                Vec3 u = PortSpec.InwardNormal(g.Inlet.Face) * s.InletVelocity;
                p.Set("uInlet", u.X);
                p.Set("vInlet", u.Y);
                p.Set("wInlet", u.Z);
            }
            if (g.Outlet.Enabled)
            {
                p.Set("pOutlet", s.OutletPressure);
            }
            if (s.ProblemId == 1)
            {
                p.Set("surfaceTension", s.SurfaceTension);
                p.Set("wetAngleEq", s.ContactAngle);
                p.Set("wetAngleAdv", s.ContactAngle);
                p.Set("wetAngleRec", s.ContactAngle);
                p.Set("initializeWithHydroPress", s.Fill == FillMode.Empty ? 0 : 1);
            }
            if (s.ViscosityModel > 0)
            {
                p.Set("viscosityModel", s.ViscosityModel);
                p.Set("fluidModelParameters", s.FluidParameters);
            }
            p.Set("vtkBinary", 1);
            p.Set("writeFields", 1);
            p.Save(Path.Combine(dir, "parameters.txt"));

            string initDir = Path.Combine(dir, "init");
            if (Directory.Exists(initDir))
            {
                foreach (var f in Directory.GetFiles(initDir, "blockInit*.txt"))
                {
                    File.Delete(f);
                }
            }
            if (s.ProblemId == 1)
            {
                OctreeMesh tree = g.LoadOctree(workspace);
                if (tree != null)
                {
                    WriteFreeSurfaceInit(g, s, tree, initDir);
                }
                else
                {
                    WriteFreeSurfaceInit(g, s, workspace, initDir);
                }
            }
        }

        /// <summary>Initial flags of a refined mesh, every block at its level (same rules as for the uniform mesh).</summary>
        static void WriteFreeSurfaceInit(GeometryDef g, SimulationDef s, OctreeMesh tree, string initDir)
        {
            Directory.CreateDirectory(initDir);
            int bs = tree.BlockSize;
            Vec3 up = s.Gravity.Length > 0 ? s.Gravity * (-1 / s.Gravity.Length) : new Vec3(0, 0, 1);
            // height of the cell centre above the domain origin measured against gravity
            Func<OctBlock, int, int, int, double> height = (b, i, j, k) =>
            {
                double h = g.Dx / (1 << b.Level);
                return ((b.Bi * bs + i + 0.5) * up.X + (b.Bj * bs + j + 0.5) * up.Y + (b.Bk * bs + k + 0.5) * up.Z) * h;
            };
            double hMin = double.MaxValue;
            foreach (var b in tree.Blocks)
                for (int k = 0; k < bs; k++)
                    for (int j = 0; j < bs; j++)
                        for (int i = 0; i < bs; i++)
                        {
                            if (b.Types[i + bs * (j + bs * k)] == CaseGenerator.Fluid)
                            {
                                hMin = Math.Min(hMin, height(b, i, j, k));
                            }
                        }
            int[,] nb = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
            var sb = new StringBuilder();
            foreach (var b in tree.Blocks)
            {
                sb.Clear();
                for (int i = 0; i < bs; i++)
                    for (int j = 0; j < bs; j++)
                        for (int k = 0; k < bs; k++)
                        {
                            bool liquid = false;
                            if (b.Types[i + bs * (j + bs * k)] == CaseGenerator.Fluid)
                            {
                                switch (s.Fill)
                                {
                                    case FillMode.Full:
                                        liquid = true;
                                        break;
                                    case FillMode.Level:
                                        liquid = height(b, i, j, k) - hMin < s.FillLevel;
                                        break;
                                }
                                // the first fluid layer in front of the inlet starts filled
                                for (int q = 0; q < 6 && !liquid; q++)
                                {
                                    liquid = tree.TypeAt(b.Level, b.Bi * bs + i + nb[q, 0], b.Bj * bs + j + nb[q, 1], b.Bk * bs + k + nb[q, 2]) == CaseGenerator.Inlet;
                                }
                            }
                            sb.Append(liquid ? "1 0 0 0 1\n" : "1 0 0 0 2\n");
                        }
                File.WriteAllText(Path.Combine(initDir, "blockInit" + tree.Tag(b) + ".txt"), sb.ToString());
            }
        }

        static void WriteFreeSurfaceInit(GeometryDef g, SimulationDef s, string workspace, string initDir)
        {
            Directory.CreateDirectory(initDir);
            byte[] types = g.LoadCellTypes(workspace);
            if (types == null)
            {
                throw new InvalidOperationException("mesh/celltypes.bin is missing, generate the mesh again.");
            }
            int GX = g.GX, GY = g.GY, GZ = g.GZ, bs = g.BlockSize;
            Func<int, int, int, long> idx = (i, j, k) => i + (long)GX * (j + (long)GY * k);
            // height above the lowest fluid cell measured against gravity
            Vec3 up = s.Gravity.Length > 0 ? s.Gravity * (-1 / s.Gravity.Length) : new Vec3(0, 0, 1);
            double hMin = double.MaxValue;
            for (int k = 0; k < GZ; k++)
                for (int j = 0; j < GY; j++)
                    for (int i = 0; i < GX; i++)
                    {
                        if (types[idx(i, j, k)] == CaseGenerator.Fluid)
                        {
                            hMin = Math.Min(hMin, (i * up.X + j * up.Y + k * up.Z) * g.Dx);
                        }
                    }
            var flags = new byte[types.Length];
            for (int k = 0; k < GZ; k++)
                for (int j = 0; j < GY; j++)
                    for (int i = 0; i < GX; i++)
                    {
                        long id = idx(i, j, k);
                        if (types[id] != CaseGenerator.Fluid)
                        {
                            flags[id] = 2;
                            continue;
                        }
                        bool liquid;
                        switch (s.Fill)
                        {
                            case FillMode.Full:
                                liquid = true;
                                break;
                            case FillMode.Level:
                                liquid = (i * up.X + j * up.Y + k * up.Z) * g.Dx - hMin < s.FillLevel;
                                break;
                            default:
                                liquid = false;
                                break;
                        }
                        // the first fluid layer in front of the inlet starts filled
                        if (!liquid)
                        {
                            int[,] nb = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
                            for (int q = 0; q < 6; q++)
                            {
                                int a = i + nb[q, 0], b = j + nb[q, 1], c = k + nb[q, 2];
                                if (a >= 0 && b >= 0 && c >= 0 && a < GX && b < GY && c < GZ && types[idx(a, b, c)] == CaseGenerator.Inlet)
                                {
                                    liquid = true;
                                }
                            }
                        }
                        flags[id] = (byte)(liquid ? 1 : 2);
                    }

            var blocks = File.ReadAllLines(Path.Combine(g.Dir(workspace), "mesh", "blocks.txt"));
            var sb = new StringBuilder();
            foreach (var line in blocks)
            {
                var t = line.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                if (t.Length < 3)
                {
                    continue;
                }
                int bi = int.Parse(t[0]), bj = int.Parse(t[1]), bk = int.Parse(t[2]);
                int id = g.NBi * g.NBj * bk + g.NBi * bj + bi;
                sb.Clear();
                for (int i = 0; i < bs; i++)
                    for (int j = 0; j < bs; j++)
                        for (int k = 0; k < bs; k++)
                        {
                            sb.Append(flags[idx(bi * bs + i, bj * bs + j, bk * bs + k)] == 1 ? "1 0 0 0 1\n" : "1 0 0 0 2\n");
                        }
                File.WriteAllText(Path.Combine(initDir, "blockInit" + id + ".txt"), sb.ToString());
            }
        }
    }
}
