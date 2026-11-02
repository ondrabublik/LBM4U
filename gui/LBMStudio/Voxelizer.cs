using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace LBMStudio
{
    /// <summary>Uniform voxel grid, cell (i, j, k) has its centre at Origin + (i + 0.5, j + 0.5, k + 0.5) * Dx.</summary>
    public class VoxelGrid
    {
        public int Nx, Ny, Nz;
        public double Dx;
        public Vec3 Origin;
        /// <summary>1 for cells inside the closed surface.</summary>
        public byte[] Inside;

        public int Index(int i, int j, int k)
        {
            return i + Nx * (j + Ny * k);
        }

        public Vec3 Center(int i, int j, int k)
        {
            return new Vec3(Origin.X + (i + 0.5) * Dx, Origin.Y + (j + 0.5) * Dx, Origin.Z + (k + 0.5) * Dx);
        }
    }

    /// <summary>Voxelization of a closed triangle surface by ray parity along z.</summary>
    public static class Voxelizer
    {
        /// <param name="padding">number of cells added around the bounding box</param>
        public static VoxelGrid Voxelize(StlMesh mesh, double dx, int padding)
        {
            var g = new VoxelGrid
            {
                Dx = dx,
                Origin = new Vec3(mesh.Min.X - padding * dx, mesh.Min.Y - padding * dx, mesh.Min.Z - padding * dx)
            };
            g.Nx = (int)Math.Ceiling((mesh.Max.X - mesh.Min.X) / dx - 1e-9) + 2 * padding;
            g.Ny = (int)Math.Ceiling((mesh.Max.Y - mesh.Min.Y) / dx - 1e-9) + 2 * padding;
            g.Nz = (int)Math.Ceiling((mesh.Max.Z - mesh.Min.Z) / dx - 1e-9) + 2 * padding;
            long n = (long)g.Nx * g.Ny * g.Nz;
            if (n > 400000000)
            {
                throw new InvalidOperationException("The mesh would have " + n + " cells, increase the cell size.");
            }
            g.Inside = new byte[n];

            // bucket the triangles by columns, rays are slightly shifted to avoid hitting edges exactly
            double sx = 0.5 + 1.3e-4, sy = 0.5 + 2.7e-4;
            var columns = new List<double>[g.Nx * g.Ny];
            int nt = mesh.TriangleCount;
            for (int t = 0; t < nt; t++)
            {
                Vec3 a = mesh.Vertex(t, 0), b = mesh.Vertex(t, 1), c = mesh.Vertex(t, 2);
                double minX = Math.Min(a.X, Math.Min(b.X, c.X)), maxX = Math.Max(a.X, Math.Max(b.X, c.X));
                double minY = Math.Min(a.Y, Math.Min(b.Y, c.Y)), maxY = Math.Max(a.Y, Math.Max(b.Y, c.Y));
                int i0 = Math.Max(0, (int)Math.Floor((minX - g.Origin.X) / dx - sx));
                int i1 = Math.Min(g.Nx - 1, (int)Math.Ceiling((maxX - g.Origin.X) / dx - sx));
                int j0 = Math.Max(0, (int)Math.Floor((minY - g.Origin.Y) / dx - sy));
                int j1 = Math.Min(g.Ny - 1, (int)Math.Ceiling((maxY - g.Origin.Y) / dx - sy));
                double d = (b.X - a.X) * (c.Y - a.Y) - (c.X - a.X) * (b.Y - a.Y);
                if (Math.Abs(d) < 1e-30)
                {
                    continue; // triangle parallel to the rays
                }
                for (int i = i0; i <= i1; i++)
                {
                    double x = g.Origin.X + (i + sx) * dx;
                    for (int j = j0; j <= j1; j++)
                    {
                        double y = g.Origin.Y + (j + sy) * dx;
                        // barycentric coordinates in the xy plane
                        double l1 = ((x - a.X) * (c.Y - a.Y) - (c.X - a.X) * (y - a.Y)) / d;
                        double l2 = ((b.X - a.X) * (y - a.Y) - (x - a.X) * (b.Y - a.Y)) / d;
                        double l0 = 1 - l1 - l2;
                        if (l0 < 0 || l1 < 0 || l2 < 0)
                        {
                            continue;
                        }
                        double z = l0 * a.Z + l1 * b.Z + l2 * c.Z;
                        int col = i + g.Nx * j;
                        if (columns[col] == null)
                        {
                            columns[col] = new List<double>();
                        }
                        columns[col].Add(z);
                    }
                }
            }

            Parallel.For(0, g.Nx * g.Ny, col =>
            {
                var zs = columns[col];
                if (zs == null || zs.Count < 2)
                {
                    return;
                }
                zs.Sort();
                int i = col % g.Nx, j = col / g.Nx;
                for (int s = 0; s + 1 < zs.Count; s += 2)
                {
                    int k0 = (int)Math.Ceiling((zs[s] - g.Origin.Z) / dx - 0.5);
                    int k1 = (int)Math.Floor((zs[s + 1] - g.Origin.Z) / dx - 0.5);
                    for (int k = Math.Max(k0, 0); k <= Math.Min(k1, g.Nz - 1); k++)
                    {
                        g.Inside[g.Index(i, j, k)] = 1;
                    }
                }
            });
            return g;
        }
    }
}
