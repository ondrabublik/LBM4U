using System;
using System.Collections.Generic;
using System.IO;
using System.Windows.Controls;
using System.Windows.Forms.Integration;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;
using W = System.Windows;

namespace LBMStudio
{
    /// <summary>Simple 3D view (WPF Viewport3D hosted in WinForms): left drag rotates, right drag pans, wheel zooms.</summary>
    public class Viewer3D : ElementHost
    {
        readonly Viewport3D viewport;
        readonly PerspectiveCamera camera;
        readonly Model3DGroup scene = new Model3DGroup();       // opaque models
        readonly Model3DGroup transparent = new Model3DGroup(); // drawn last, otherwise they hide what is behind them
        readonly Dictionary<string, Model3D> layers = new Dictionary<string, Model3D>();
        readonly Border border;
        Point3D target = new Point3D(0, 0, 0);
        double distance = 1, theta = -60, phi = 25;
        W.Point last;
        int dragButton; // 0 none, 1 rotate, 2 pan

        public Viewer3D()
        {
            viewport = new Viewport3D();
            camera = new PerspectiveCamera { FieldOfView = 40, UpDirection = new Vector3D(0, 0, 1) };
            viewport.Camera = camera;
            var lights = new Model3DGroup();
            lights.Children.Add(new AmbientLight(Color.FromRgb(90, 90, 90)));
            lights.Children.Add(new DirectionalLight(Color.FromRgb(200, 200, 200), new Vector3D(-1, -0.5, -1)));
            lights.Children.Add(new DirectionalLight(Color.FromRgb(110, 110, 120), new Vector3D(1, 0.8, 0.3)));
            var all = new Model3DGroup();
            all.Children.Add(lights);
            all.Children.Add(scene);
            all.Children.Add(transparent);
            viewport.Children.Add(new ModelVisual3D { Content = all });
            border = new Border
            {
                Background = new LinearGradientBrush(Color.FromRgb(250, 251, 253), Color.FromRgb(214, 221, 232), 90),
                Child = viewport
            };
            border.MouseDown += (s, e) =>
            {
                dragButton = e.ChangedButton == W.Input.MouseButton.Left ? 1 : 2;
                last = e.GetPosition(border);
                border.CaptureMouse();
            };
            border.MouseUp += (s, e) =>
            {
                dragButton = 0;
                border.ReleaseMouseCapture();
            };
            border.MouseMove += (s, e) =>
            {
                if (dragButton == 0)
                {
                    return;
                }
                var p = e.GetPosition(border);
                double dx = p.X - last.X, dy = p.Y - last.Y;
                last = p;
                if (dragButton == 1)
                {
                    theta -= dx * 0.4;
                    phi = Math.Max(-89, Math.Min(89, phi + dy * 0.4));
                }
                else
                {
                    Vector3D look = camera.LookDirection;
                    look.Normalize();
                    Vector3D right = Vector3D.CrossProduct(look, new Vector3D(0, 0, 1));
                    right.Normalize();
                    Vector3D up = Vector3D.CrossProduct(right, look);
                    double scale = distance * 0.0015;
                    target -= right * dx * scale;
                    target += up * dy * scale;
                }
                UpdateCamera();
            };
            border.MouseWheel += (s, e) =>
            {
                distance *= e.Delta > 0 ? 0.9 : 1.1;
                UpdateCamera();
            };
            Child = border;
            UpdateCamera();
        }

        void UpdateCamera()
        {
            double t = theta * Math.PI / 180, p = phi * Math.PI / 180;
            var dir = new Vector3D(Math.Cos(p) * Math.Cos(t), Math.Cos(p) * Math.Sin(t), Math.Sin(p));
            camera.Position = target + dir * distance;
            camera.LookDirection = -dir * distance;
            camera.NearPlaneDistance = distance * 0.01;
            camera.FarPlaneDistance = distance * 100;
        }

        public void FitView(Vec3 min, Vec3 max)
        {
            target = new Point3D((min.X + max.X) / 2, (min.Y + max.Y) / 2, (min.Z + max.Z) / 2);
            double size = (max - min).Length;
            distance = Math.Max(size, 1e-6) * 1.6;
            UpdateCamera();
        }

        public void SetLayer(string name, Model3D model, bool isTransparent = false)
        {
            Model3D old;
            if (layers.TryGetValue(name, out old))
            {
                scene.Children.Remove(old);
                transparent.Children.Remove(old);
                layers.Remove(name);
            }
            if (model != null)
            {
                layers[name] = model;
                (isTransparent ? transparent : scene).Children.Add(model);
            }
        }

        public void Clear()
        {
            scene.Children.Clear();
            transparent.Children.Clear();
            layers.Clear();
        }

        static Material MakeMaterial(Color color, double opacity)
        {
            var brush = new SolidColorBrush(color) { Opacity = opacity };
            brush.Freeze();
            var group = new MaterialGroup();
            group.Children.Add(new DiffuseMaterial(brush));
            if (opacity >= 1)
            {
                group.Children.Add(new SpecularMaterial(new SolidColorBrush(Color.FromArgb(60, 255, 255, 255)), 30));
            }
            return group;
        }

        /// <summary>Flat shaded triangles, 9 coordinates per triangle.</summary>
        public static Model3D Triangles(IList<double> xyz, Color color, double opacity)
        {
            var mesh = new MeshGeometry3D();
            var pos = new Point3DCollection(xyz.Count / 3);
            var idx = new Int32Collection(xyz.Count / 3);
            for (int i = 0; i + 2 < xyz.Count; i += 3)
            {
                pos.Add(new Point3D(xyz[i], xyz[i + 1], xyz[i + 2]));
                idx.Add(i / 3);
            }
            mesh.Positions = pos;
            mesh.TriangleIndices = idx;
            mesh.Freeze();
            var mat = MakeMaterial(color, opacity);
            var model = new GeometryModel3D(mesh, mat) { BackMaterial = mat };
            model.Freeze();
            return model;
        }

        /// <summary>Cylinder from a to b (used for arrows and edges).</summary>
        public static void AddCylinder(List<double> xyz, Vec3 a, Vec3 b, double r1, double r2, int seg = 16)
        {
            Vec3 axis = b - a;
            double len = axis.Length;
            if (len == 0)
            {
                return;
            }
            Vec3 n = axis * (1 / len);
            Vec3 helper = Math.Abs(n.Z) < 0.9 ? new Vec3(0, 0, 1) : new Vec3(1, 0, 0);
            Vec3 u = Vec3.Cross(n, helper);
            u = u * (1 / u.Length);
            Vec3 v = Vec3.Cross(n, u);
            for (int s = 0; s < seg; s++)
            {
                double a0 = 2 * Math.PI * s / seg, a1 = 2 * Math.PI * (s + 1) / seg;
                Vec3 d0 = u * Math.Cos(a0) + v * Math.Sin(a0), d1 = u * Math.Cos(a1) + v * Math.Sin(a1);
                Vec3 p0 = a + d0 * r1, p1 = a + d1 * r1, q0 = b + d0 * r2, q1 = b + d1 * r2;
                AddTri(xyz, p0, p1, q1);
                if (r2 > 0)
                {
                    AddTri(xyz, p0, q1, q0);
                }
                AddTri(xyz, a, p1, p0);
            }
        }

        public static void AddTri(List<double> xyz, Vec3 a, Vec3 b, Vec3 c)
        {
            xyz.Add(a.X); xyz.Add(a.Y); xyz.Add(a.Z);
            xyz.Add(b.X); xyz.Add(b.Y); xyz.Add(b.Z);
            xyz.Add(c.X); xyz.Add(c.Y); xyz.Add(c.Z);
        }

        public static Model3D Arrow(Vec3 start, Vec3 vector, double radius, Color color)
        {
            var xyz = new List<double>();
            Vec3 end = start + vector * 0.75;
            AddCylinder(xyz, start, end, radius, radius);
            AddCylinder(xyz, end, start + vector, radius * 2.5, 0);
            return Triangles(xyz, color, 1);
        }

        public static Model3D BoxEdges(Vec3 min, Vec3 max, double radius, Color color)
        {
            var xyz = new List<double>();
            var c = new Vec3[8];
            for (int i = 0; i < 8; i++)
            {
                c[i] = new Vec3((i & 1) == 0 ? min.X : max.X, (i & 2) == 0 ? min.Y : max.Y, (i & 4) == 0 ? min.Z : max.Z);
            }
            int[,] e = { { 0, 1 }, { 2, 3 }, { 4, 5 }, { 6, 7 }, { 0, 2 }, { 1, 3 }, { 4, 6 }, { 5, 7 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
            for (int i = 0; i < 12; i++)
            {
                AddCylinder(xyz, c[e[i, 0]], c[e[i, 1]], radius, radius, 6);
            }
            return Triangles(xyz, color, 1);
        }

        /// <summary>Disc or rectangle of a port on a face of the box.</summary>
        public static Model3D Port(PortSpec port, Vec3 min, Vec3 max, Color color)
        {
            int axis = PortSpec.Axis(port.Face);
            int a1 = axis == 0 ? 1 : 0, a2 = axis == 2 ? 1 : 2;
            double plane = PortSpec.InwardSign(port.Face) > 0 ? min[axis] : max[axis];
            double off = (max - min).Length * 0.002 * -PortSpec.InwardSign(port.Face);
            Func<double, double, Vec3> P = (p, q) =>
            {
                var a = new double[3];
                a[axis] = plane + off;
                a[a1] = p;
                a[a2] = q;
                return new Vec3(a[0], a[1], a[2]);
            };
            var xyz = new List<double>();
            if (port.Radius <= 0)
            {
                AddTri(xyz, P(min[a1], min[a2]), P(max[a1], min[a2]), P(max[a1], max[a2]));
                AddTri(xyz, P(min[a1], min[a2]), P(max[a1], max[a2]), P(min[a1], max[a2]));
            }
            else
            {
                int seg = 48;
                for (int s = 0; s < seg; s++)
                {
                    double t0 = 2 * Math.PI * s / seg, t1 = 2 * Math.PI * (s + 1) / seg;
                    AddTri(xyz, P(port.U, port.V), P(port.U + port.Radius * Math.Cos(t0), port.V + port.Radius * Math.Sin(t0)),
                        P(port.U + port.Radius * Math.Cos(t1), port.V + port.Radius * Math.Sin(t1)));
                }
            }
            return Triangles(xyz, color, 0.85);
        }

        /// <summary>Renders the current view to a PNG file.</summary>
        public void SaveImage(string path, int width, int height)
        {
            border.Measure(new W.Size(width, height));
            border.Arrange(new W.Rect(0, 0, width, height));
            border.UpdateLayout();
            var bmp = new RenderTargetBitmap(width, height, 96, 96, PixelFormats.Pbgra32);
            bmp.Render(border);
            var enc = new PngBitmapEncoder();
            enc.Frames.Add(BitmapFrame.Create(bmp));
            using (var fs = File.Create(path))
            {
                enc.Save(fs);
            }
        }
    }
}
