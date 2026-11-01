using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Windows.Forms;
using System.Windows.Media;
using Color = System.Drawing.Color;

namespace LBMStudio
{
    public partial class MainForm
    {
        ComboBox cbResSim, cbResField, cbResAxis;
        TrackBar tbResTime, tbResSlice;
        Label lblResTime, lblResSlice;
        CheckBox chkVectors, chkPlay, chkSurface3D;
        SliceView sliceView;
        Viewer3D resultView;
        Timer playTimer;
        List<KeyValuePair<int, string>> resultFiles = new List<KeyValuePair<int, string>>();
        FieldsData fields;
        string fieldsPath;
        byte[] cellTypes;

        static readonly string[] FieldLabels = { "Velikost rychlosti |u|", "Rychlost ux", "Rychlost uy", "Rychlost uz", "Přetlak p", "Objemový podíl VOF", "Stav buňky" };

        TabPage BuildResultsTab()
        {
            var page = new TabPage("3. Výsledky");
            var bar = new FlowLayoutPanel { Dock = DockStyle.Top, AutoSize = true, AutoSizeMode = AutoSizeMode.GrowAndShrink, Padding = new Padding(4), WrapContents = true };
            Func<string, Control, Control> labeled = (text, c) =>
            {
                var p = new FlowLayoutPanel { AutoSize = true, FlowDirection = FlowDirection.TopDown, Margin = new Padding(4, 0, 8, 0) };
                p.Controls.Add(new Label { Text = text, AutoSize = true, ForeColor = Color.DimGray });
                p.Controls.Add(c);
                return p;
            };
            cbResSim = Combo();
            cbResSim.Width = 150;
            cbResSim.SelectedIndexChanged += (s, e) => LoadResultList();
            cbResField = Combo(FieldLabels);
            cbResField.Width = 170;
            cbResField.SelectedIndexChanged += (s, e) => UpdateSlice();
            cbResAxis = Combo("řez kolmý na x", "řez kolmý na y", "řez kolmý na z");
            cbResAxis.SelectedIndex = 1;
            cbResAxis.Width = 120;
            cbResAxis.SelectedIndexChanged += (s, e) => { ResetSliceRange(); UpdateSlice(); };
            tbResTime = new TrackBar { Width = 260, Minimum = 0, Maximum = 0, TickStyle = TickStyle.None, Height = 30 };
            tbResTime.ValueChanged += (s, e) => LoadFields();
            lblResTime = new Label { AutoSize = true, Text = "—" };
            tbResSlice = new TrackBar { Width = 200, Minimum = 0, Maximum = 0, TickStyle = TickStyle.None, Height = 30 };
            tbResSlice.ValueChanged += (s, e) => UpdateSlice();
            lblResSlice = new Label { AutoSize = true, Text = "—" };
            chkVectors = new CheckBox { Text = "Vektory", AutoSize = true };
            chkVectors.CheckedChanged += (s, e) => UpdateSlice();
            chkSurface3D = new CheckBox { Text = "Hladina 3D", AutoSize = true, Checked = true };
            chkSurface3D.CheckedChanged += (s, e) => Update3DResult();
            chkPlay = new CheckBox { Text = "▶ Přehrát", Appearance = Appearance.Button, AutoSize = true };
            playTimer = new Timer { Interval = 400 };
            playTimer.Tick += (s, e) =>
            {
                if (tbResTime.Maximum == 0)
                {
                    return;
                }
                tbResTime.Value = tbResTime.Value >= tbResTime.Maximum ? 0 : tbResTime.Value + 1;
            };
            chkPlay.CheckedChanged += (s, e) => playTimer.Enabled = chkPlay.Checked;
            var btnRefresh = new Button { Text = "Obnovit", AutoSize = true };
            btnRefresh.Click += (s, e) => LoadResultList();
            var btnPv = new Button { Text = "Otevřít v ParaView", AutoSize = true };
            btnPv.Click += (s, e) => OpenParaView();
            var btnPng = new Button { Text = "Uložit obrázek…", AutoSize = true };
            btnPng.Click += (s, e) =>
            {
                using (var d = new SaveFileDialog { Filter = "PNG (*.png)|*.png", FileName = "rez.png" })
                {
                    if (d.ShowDialog(this) == DialogResult.OK)
                    {
                        sliceView.SaveImage(d.FileName);
                    }
                }
            };

            var timePanel = new FlowLayoutPanel { AutoSize = true, FlowDirection = FlowDirection.LeftToRight };
            timePanel.Controls.Add(tbResTime);
            timePanel.Controls.Add(lblResTime);
            var slicePanel = new FlowLayoutPanel { AutoSize = true, FlowDirection = FlowDirection.LeftToRight };
            slicePanel.Controls.Add(tbResSlice);
            slicePanel.Controls.Add(lblResSlice);
            bar.Controls.Add(labeled("Simulace", cbResSim));
            bar.Controls.Add(labeled("Veličina", cbResField));
            bar.Controls.Add(labeled("Řez", cbResAxis));
            bar.Controls.Add(labeled("Poloha řezu", slicePanel));
            bar.Controls.Add(labeled("Čas", timePanel));
            var opts = new FlowLayoutPanel { AutoSize = true, FlowDirection = FlowDirection.LeftToRight, Margin = new Padding(4, 14, 4, 0) };
            opts.Controls.AddRange(new Control[] { chkPlay, chkVectors, chkSurface3D, btnRefresh, btnPng, btnPv });
            bar.Controls.Add(opts);

            var split = new SplitContainer { Dock = DockStyle.Fill };
            page.Layout += (s, e) =>
            {
                if (!resultSplitInitialised && split.Width > 200)
                {
                    split.SplitterDistance = split.Width * 55 / 100;
                    resultSplitInitialised = true;
                }
            };
            sliceView = new SliceView { Dock = DockStyle.Fill };
            resultView = new Viewer3D { Dock = DockStyle.Fill };
            split.Panel1.Controls.Add(sliceView);
            split.Panel2.Controls.Add(resultView);
            page.Controls.Add(split);
            page.Controls.Add(bar);
            return page;
        }

        void RefreshResultsSimulations()
        {
            string current = cbResSim.SelectedItem as string;
            cbResSim.Items.Clear();
            if (geom == null)
            {
                LoadResultList();
                return;
            }
            foreach (var item in lbSims.Items)
            {
                cbResSim.Items.Add(item);
            }
            if (current != null && cbResSim.Items.Contains(current))
            {
                cbResSim.SelectedItem = current;
            }
            else if (lbSims.SelectedItem != null)
            {
                cbResSim.SelectedItem = lbSims.SelectedItem;
            }
            else if (cbResSim.Items.Count > 0)
            {
                cbResSim.SelectedIndex = 0;
            }
            cellTypes = geom.MeshReady ? geom.LoadCellTypes(settings.Workspace) : null;
            if (cbResSim.SelectedIndex < 0)
            {
                // no simulation selected, drop the fields of the previous geometry
                LoadResultList();
            }
        }

        string ResultSimDir
        {
            get
            {
                string name = cbResSim.SelectedItem as string;
                return geom == null || name == null ? null : Path.Combine(geom.Dir(settings.Workspace), name);
            }
        }

        void LoadResultList()
        {
            string dir = ResultSimDir;
            resultFiles = dir == null ? new List<KeyValuePair<int, string>>() : FieldsData.List(dir);
            fields = null;
            fieldsPath = null;
            tbResTime.Maximum = Math.Max(0, resultFiles.Count - 1);
            if (resultFiles.Count == 0)
            {
                lblResTime.Text = "žádné výsledky";
                sliceView.SetData(new float[1, 1] { { float.NaN } }, null, 1, 1, "", "", "Simulace zatím nemá výsledky", "");
                resultView.Clear();
                return;
            }
            ResetSliceRange();
            if (tbResTime.Value == tbResTime.Maximum)
            {
                LoadFields();
            }
            else
            {
                tbResTime.Value = tbResTime.Maximum;
            }
        }

        void ResetSliceRange()
        {
            if (geom == null || !geom.MeshReady)
            {
                return;
            }
            int axis = cbResAxis.SelectedIndex;
            int n = axis == 0 ? geom.Nx : (axis == 1 ? geom.Ny : geom.Nz);
            tbResSlice.Maximum = Math.Max(0, n - 1);
            tbResSlice.Value = n / 2;
        }

        void LoadFields()
        {
            if (resultFiles.Count == 0 || geom == null)
            {
                return;
            }
            var item = resultFiles[Math.Min(tbResTime.Value, resultFiles.Count - 1)];
            if (item.Value != fieldsPath)
            {
                try
                {
                    fields = FieldsData.Read(item.Value, geom.GX, geom.GY, geom.GZ);
                    fieldsPath = item.Value;
                }
                catch (Exception ex)
                {
                    SetStatus("Chyba čtení výsledků: " + ex.Message);
                    return;
                }
            }
            lblResTime.Text = "t = " + Fmt.G(fields.Time, 4) + " s (iterace " + fields.Iteration + ")";
            UpdateSlice();
            Update3DResult();
        }

        void UpdateSlice()
        {
            if (fields == null || geom == null || !geom.MeshReady)
            {
                return;
            }
            int axis = cbResAxis.SelectedIndex;
            int pos = tbResSlice.Value;
            int a1 = axis == 0 ? 1 : 0, a2 = axis == 2 ? 1 : 2;
            int[] n = { geom.Nx, geom.Ny, geom.Nz };
            string[] names = { "x", "y", "z" };
            lblResSlice.Text = names[axis] + " = " + Fmt.G(geom.Origin[axis] + (pos + 0.5) * geom.Dx, 4) + " m";
            int na = n[a1], nb = n[a2];
            var vals = new float[na, nb];
            var kind = new byte[na, nb];
            int field = cbResField.SelectedIndex;
            int iux = fields.FieldIndex("ux"), iuy = fields.FieldIndex("uy"), iuz = fields.FieldIndex("uz");
            int ip = fields.FieldIndex("p"), itype = fields.FieldIndex("type"), ivof = fields.FieldIndex("vof"), iflag = fields.FieldIndex("flag");
            int[] ucomp = { iux, iuy, iuz };
            float[,] vu = chkVectors.Checked ? new float[na, nb] : null;
            float[,] vv = chkVectors.Checked ? new float[na, nb] : null;
            string unit = "";
            for (int a = 0; a < na; a++)
            {
                for (int b = 0; b < nb; b++)
                {
                    int[] c = new int[3];
                    c[axis] = pos;
                    c[a1] = a;
                    c[a2] = b;
                    long id = fields.Index(c[0], c[1], c[2]);
                    float type = fields.Data[itype][id];
                    float flag = iflag >= 0 ? fields.Data[iflag][id] : 1;
                    if (float.IsNaN(type) || type != 0)
                    {
                        kind[a, b] = 1;
                        vals[a, b] = float.NaN;
                        continue;
                    }
                    if (flag == 2 && field != 5 && field != 6)
                    {
                        kind[a, b] = 2;
                    }
                    float ux = fields.Data[iux][id], uy = fields.Data[iuy][id], uz = fields.Data[iuz][id];
                    switch (field)
                    {
                        case 0:
                            vals[a, b] = (float)Math.Sqrt(ux * ux + uy * uy + uz * uz);
                            unit = "m/s";
                            break;
                        case 1:
                        case 2:
                        case 3:
                            vals[a, b] = fields.Data[ucomp[field - 1]][id];
                            unit = "m/s";
                            break;
                        case 4:
                            vals[a, b] = fields.Data[ip][id];
                            unit = "Pa";
                            break;
                        case 5:
                            vals[a, b] = ivof >= 0 ? (flag == 1 ? 1f : (flag == 2 ? 0f : Math.Max(0f, Math.Min(1f, fields.Data[ivof][id])))) : 1f;
                            unit = "–";
                            break;
                        case 6:
                            vals[a, b] = flag;
                            unit = "1 kap., 2 plyn, 3 rozhraní";
                            break;
                    }
                    if (vu != null)
                    {
                        vu[a, b] = fields.Data[ucomp[a1]][id];
                        vv[a, b] = fields.Data[ucomp[a2]][id];
                    }
                }
            }
            sliceView.FixedRange = field == 5;
            sliceView.RangeMin = 0;
            sliceView.RangeMax = 1;
            sliceView.SetData(vals, kind, na * geom.Dx, nb * geom.Dx, names[a1], names[a2],
                FieldLabels[field] + ", t = " + Fmt.G(fields.Time, 4) + " s", unit, vu, vv);
        }

        void Update3DResult()
        {
            resultView.Clear();
            if (geom == null || fields == null)
            {
                return;
            }
            Vec3 min = geom.Origin, max = geom.Origin + new Vec3(geom.Nx, geom.Ny, geom.Nz) * geom.Dx;
            string stl = geom.StlPath(settings.Workspace);
            if (File.Exists(stl))
            {
                try
                {
                    var body = StlMesh.Read(stl);
                    resultView.SetLayer("body", Viewer3D.Triangles(body.V, System.Windows.Media.Color.FromRgb(170, 175, 185), 0.25), true);
                }
                catch (Exception)
                {
                }
            }
            resultView.SetLayer("box", Viewer3D.BoxEdges(min, max, (max - min).Length * 0.0012, Colors.DimGray));
            // slice plane position
            int axis = cbResAxis.SelectedIndex;
            double plane = geom.Origin[axis] + (tbResSlice.Value + 0.5) * geom.Dx;
            var smin = new double[] { min.X, min.Y, min.Z };
            var smax = new double[] { max.X, max.Y, max.Z };
            smin[axis] = plane;
            smax[axis] = plane;
            var pp = new List<double>();
            int a1 = axis == 0 ? 1 : 0, a2 = axis == 2 ? 1 : 2;
            Func<double, double, Vec3> P = (u, v) => { var c = new double[3]; c[axis] = plane; c[a1] = u; c[a2] = v; return new Vec3(c[0], c[1], c[2]); };
            Viewer3D.AddTri(pp, P(smin[a1], smin[a2]), P(smax[a1], smin[a2]), P(smax[a1], smax[a2]));
            Viewer3D.AddTri(pp, P(smin[a1], smin[a2]), P(smax[a1], smax[a2]), P(smin[a1], smax[a2]));
            resultView.SetLayer("plane", Viewer3D.Triangles(pp, Colors.Orange, 0.25), true);
            if (chkSurface3D.Checked && fields.ProblemId != 0)
            {
                string fs = SurfaceVtk.FileForIteration(ResultSimDir, fields.Iteration, fields.Time);
                if (fs != null)
                {
                    try
                    {
                        float[] tri = SurfaceVtk.Read(fs);
                        var xyz = new double[tri.Length];
                        for (int i = 0; i < tri.Length; i++)
                        {
                            xyz[i] = geom.Origin[i % 3] + tri[i] * geom.Dx;
                        }
                        resultView.SetLayer("surface", Viewer3D.Triangles(xyz, System.Windows.Media.Color.FromRgb(40, 120, 220), 1.0));
                    }
                    catch (Exception ex)
                    {
                        SetStatus("Chyba čtení hladiny: " + ex.Message);
                    }
                }
            }
            if (!resultViewFitted)
            {
                resultView.FitView(min, max);
                resultViewFitted = true;
            }
        }

        bool resultViewFitted;
        bool resultSplitInitialised;

        void OpenParaView()
        {
            string dir = ResultSimDir;
            if (dir == null)
            {
                return;
            }
            string pvd = Path.Combine(dir, "result", "output.pvd");
            if (string.IsNullOrEmpty(settings.ParaViewExe) || !File.Exists(settings.ParaViewExe))
            {
                MessageBox.Show(this, "ParaView nebyl nalezen, nastavte cestu v Nastavení.", Text);
                return;
            }
            var args = new List<string>();
            if (File.Exists(pvd))
            {
                args.Add("\"" + pvd + "\"");
            }
            string pvdFs = Path.Combine(dir, "result", "outputFS.pvd");
            if (File.Exists(pvdFs))
            {
                args.Add("\"" + pvdFs + "\"");
            }
            Process.Start(settings.ParaViewExe, "--data=" + (args.Count > 0 ? args[0] : "\"" + dir + "\""));
        }

        /// <summary>Self test helper: shows the latest result and saves images.</summary>
        public void ShowResultsAndSave(string sliceFile, string view3dFile)
        {
            tabs.SelectedIndex = 2;
            RefreshResultsSimulations();
            Application.DoEvents();
            LoadResultList();
            Application.DoEvents();
            sliceView.SaveImage(sliceFile);
            resultView.SaveImage(view3dFile, 900, 700);
        }
    }
}
