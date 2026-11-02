using System;
using System.Drawing;
using System.IO;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Media;
using Color = System.Drawing.Color;

namespace LBMStudio
{
    public partial class MainForm
    {
        TextBox tbSource, tbMargin, tbDx;
        ComboBox cbUnits, cbBlock;
        CheckBox chkRefine;
        ComboBox cbMaxLevel, cbWallLevel;
        TextBox tbGapCells, tbCurvCells;
        RadioButton rbInside, rbOutside;
        Label lblGeomInfo, lblMeshEstimate;
        CheckBox chkInlet, chkOutlet;
        ComboBox cbInletFace, cbOutletFace;
        TextBox tbInletU, tbInletV, tbInletR, tbInletVel, tbOutletU, tbOutletV, tbOutletR, tbOutletP;
        TextBox tbGx, tbGy, tbGz;
        Button btnGenerate;
        TextBox geomLog;
        Viewer3D geomView;
        bool loadingGeometry;
        CheckBox chkShowModel, chkShowCells, chkShowBlocks, chkMeshCutFlip;
        ComboBox cbMeshCut;
        TrackBar tbMeshCut;
        Label lblMeshCut;
        byte[] meshTypes;               // cell types of the displayed mesh, loaded on demand
        OctreeMesh meshTree;            // blocks of a refined mesh, loaded on demand
        System.Windows.Media.Media3D.Model3D cellsModel, blocksModel; // cached, rebuilt when the mesh or the cut changes
        bool updatingMeshCut;

        static readonly string[] FaceNames = { "x min", "x max", "y min", "y max", "z min", "z max" };
        static readonly string[] UnitNames = { "mm", "m", "cm", "inches" };
        static readonly double[] UnitScales = { 0.001, 1, 0.01, 0.0254 };

        TabPage BuildGeometryTab()
        {
            var page = new TabPage("1. Geometry and mesh");

            // CAD
            var tCad = FieldTable();
            tbSource = new TextBox { ReadOnly = true };
            var btnBrowse = new Button { Text = "Load STEP / STL…", Height = 28 };
            btnBrowse.Click += (s, e) => BrowseGeometry();
            cbUnits = Combo(UnitNames);
            cbUnits.SelectedIndexChanged += (s, e) => { if (!loadingGeometry) ApplyUnits(); };
            tCad.Controls.Add(btnBrowse);
            tCad.SetColumnSpan(btnBrowse, 3);
            btnBrowse.Dock = DockStyle.Fill;
            AddRow(tCad, "File", tbSource);
            AddRow(tCad, "File units", cbUnits);
            lblGeomInfo = new Label { AutoSize = true, ForeColor = Color.DimGray, MaximumSize = new Size(390, 0) };
            tCad.Controls.Add(lblGeomInfo);
            tCad.SetColumnSpan(lblGeomInfo, 3);

            // fluid region
            var tReg = FieldTable();
            rbInside = new RadioButton { Text = "Fluid inside the body (cavity, mold)", AutoSize = true, Checked = true };
            rbOutside = new RadioButton { Text = "Flow around the body (box around it)", AutoSize = true };
            tReg.Controls.Add(rbInside);
            tReg.SetColumnSpan(rbInside, 3);
            tReg.Controls.Add(rbOutside);
            tReg.SetColumnSpan(rbOutside, 3);
            tbMargin = AddRow(tReg, "Margin around the body", Num("0.02"), "m");
            rbOutside.CheckedChanged += (s, e) => { tbMargin.Enabled = rbOutside.Checked; UpdateGeometryPreview(); };
            tbMargin.Enabled = false;

            // mesh
            var tMesh = FieldTable();
            tbDx = AddRow(tMesh, "Cell size dx", Num("0.002"), "m");
            cbBlock = AddRow(tMesh, "Block size", Combo("16", "8", "32"), "cells");
            chkRefine = new CheckBox { Text = "Refinement by geometry (octree)", AutoSize = true };
            tMesh.Controls.Add(chkRefine);
            tMesh.SetColumnSpan(chkRefine, 3);
            cbMaxLevel = AddRow(tMesh, "Max. refinement level", Combo("1", "2", "3", "4"), "dx / 2^L");
            tbGapCells = AddRow(tMesh, "Cells across a gap", Num("10"));
            tbCurvCells = AddRow(tMesh, "Cells per curvature radius", Num("0"), "0 = off");
            cbWallLevel = AddRow(tMesh, "Level at all walls", Combo("0", "1", "2", "3", "4"));
            cbMaxLevel.SelectedIndex = 1;
            cbWallLevel.SelectedIndex = 0;
            chkRefine.CheckedChanged += (s, e) => { UpdateRefineControls(); UpdateGeometryPreview(); };
            cbMaxLevel.SelectedIndexChanged += (s, e) => UpdateGeometryPreview();
            UpdateRefineControls();
            lblMeshEstimate = new Label { AutoSize = true, ForeColor = Color.DimGray, MaximumSize = new Size(390, 0) };
            tMesh.Controls.Add(lblMeshEstimate);
            tMesh.SetColumnSpan(lblMeshEstimate, 3);
            tbDx.TextChanged += (s, e) => UpdateGeometryPreview();
            tbMargin.TextChanged += (s, e) => UpdateGeometryPreview();

            // inlet
            var tIn = FieldTable();
            chkInlet = new CheckBox { Text = "Inlet (velocity)", AutoSize = true, Checked = true };
            tIn.Controls.Add(chkInlet);
            tIn.SetColumnSpan(chkInlet, 3);
            cbInletFace = AddRow(tIn, "Box face", Combo(FaceNames));
            tbInletU = AddRow(tIn, "Center – coordinate 1", Num("0"), "m");
            tbInletV = AddRow(tIn, "Center – coordinate 2", Num("0"), "m");
            tbInletR = AddRow(tIn, "Radius (0 = whole face)", Num("0"), "m");
            tbInletVel = AddRow(tIn, "Velocity (default)", Num("0.5"), "m/s");
            var btnCenterIn = new Button { Text = "Center on face", AutoSize = true };
            btnCenterIn.Click += (s, e) => CenterPort(cbInletFace, tbInletU, tbInletV);
            tIn.Controls.Add(btnCenterIn);
            tIn.SetColumnSpan(btnCenterIn, 3);
            cbInletFace.SelectedIndex = 5;

            // outlet
            var tOut = FieldTable();
            chkOutlet = new CheckBox { Text = "Outlet (pressure)", AutoSize = true };
            tOut.Controls.Add(chkOutlet);
            tOut.SetColumnSpan(chkOutlet, 3);
            cbOutletFace = AddRow(tOut, "Box face", Combo(FaceNames));
            tbOutletU = AddRow(tOut, "Center – coordinate 1", Num("0"), "m");
            tbOutletV = AddRow(tOut, "Center – coordinate 2", Num("0"), "m");
            tbOutletR = AddRow(tOut, "Radius (0 = whole face)", Num("0"), "m");
            tbOutletP = AddRow(tOut, "Gauge pressure (default)", Num("0"), "Pa");
            var btnCenterOut = new Button { Text = "Center on face", AutoSize = true };
            btnCenterOut.Click += (s, e) => CenterPort(cbOutletFace, tbOutletU, tbOutletV);
            tOut.Controls.Add(btnCenterOut);
            tOut.SetColumnSpan(btnCenterOut, 3);

            // gravity
            var tG = FieldTable();
            tbGx = AddRow(tG, "gx", Num("0"), "m/s²");
            tbGy = AddRow(tG, "gy", Num("0"), "m/s²");
            tbGz = AddRow(tG, "gz", Num("-9.81"), "m/s²");

            foreach (var c in new Control[] { chkInlet, chkOutlet })
            {
                ((CheckBox)c).CheckedChanged += (s, e) => UpdateGeometryPreview();
            }
            foreach (var c in new Control[] { cbInletFace, cbOutletFace })
            {
                ((ComboBox)c).SelectedIndexChanged += (s, e) => UpdateGeometryPreview();
            }
            foreach (var t in new[] { tbInletU, tbInletV, tbInletR, tbOutletU, tbOutletV, tbOutletR, tbGx, tbGy, tbGz })
            {
                t.TextChanged += (s, e) => UpdateGeometryPreview();
            }

            btnGenerate = new Button { Text = "Generate mesh", Dock = DockStyle.Top, Height = 38, Font = new Font(Font, FontStyle.Bold) };
            btnGenerate.Click += (s, e) => GenerateMesh();

            var column = Column(430,
                Group("CAD model", tCad),
                Group("Fluid region", tReg),
                Group("Mesh", tMesh),
                Group("Inlet", tIn),
                Group("Outlet", tOut),
                Group("Gravity (default for simulations)", tG),
                btnGenerate);

            geomView = new Viewer3D { Dock = DockStyle.Fill };
            geomLog = new TextBox { Dock = DockStyle.Bottom, Height = 120, Multiline = true, ScrollBars = ScrollBars.Vertical, ReadOnly = true, Font = new Font("Consolas", 8.5f) };
            var hint = new Label { Dock = DockStyle.Top, Height = 20, ForeColor = Color.DimGray, Text = "  Left button – rotate, right – pan, wheel – zoom. Red = inlet, green = outlet, blue arrow = gravity, blue cells = fluid." };

            // display of the solver mesh
            var meshBar = new FlowLayoutPanel { Dock = DockStyle.Top, AutoSize = true, AutoSizeMode = AutoSizeMode.GrowAndShrink, Padding = new Padding(4, 2, 4, 0), WrapContents = true };
            chkShowModel = new CheckBox { Text = "CAD model", AutoSize = true, Checked = true, Margin = new Padding(3, 6, 8, 0) };
            chkShowCells = new CheckBox { Text = "Mesh cells", AutoSize = true, Margin = new Padding(3, 6, 8, 0), Enabled = false };
            chkShowBlocks = new CheckBox { Text = "Blocks", AutoSize = true, Margin = new Padding(3, 6, 8, 0), Enabled = false };
            cbMeshCut = Combo("no cut", "cut normal to x", "cut normal to y", "cut normal to z");
            cbMeshCut.Width = 120;
            cbMeshCut.Enabled = false;
            tbMeshCut = new TrackBar { Width = 200, Minimum = 0, Maximum = 0, TickStyle = TickStyle.None, Height = 30, Enabled = false };
            lblMeshCut = new Label { AutoSize = true, ForeColor = Color.DimGray, Margin = new Padding(3, 6, 3, 0) };
            chkMeshCutFlip = new CheckBox { Text = "flip", AutoSize = true, Margin = new Padding(3, 6, 8, 0), Enabled = false };
            chkMeshCutFlip.CheckedChanged += (s, e) => { cellsModel = null; UpdateMeshLayers(); };
            chkShowModel.CheckedChanged += (s, e) => UpdateGeometryPreview();
            chkShowCells.CheckedChanged += (s, e) =>
            {
                if (chkShowCells.Checked)
                {
                    chkShowModel.Checked = false; // the model surface would cover the cells
                }
                UpdateMeshLayers();
            };
            chkShowBlocks.CheckedChanged += (s, e) => UpdateMeshLayers();
            cbMeshCut.SelectedIndexChanged += (s, e) => { ResetMeshCut(); cellsModel = null; UpdateMeshLayers(); };
            tbMeshCut.ValueChanged += (s, e) =>
            {
                if (!updatingMeshCut)
                {
                    cellsModel = null;
                    UpdateMeshLayers();
                }
            };
            meshBar.Controls.AddRange(new Control[] { new Label { Text = "Show:", AutoSize = true, Margin = new Padding(3, 6, 3, 0) }, chkShowModel, chkShowCells, chkShowBlocks, cbMeshCut, tbMeshCut, chkMeshCutFlip, lblMeshCut });

            var right = new Panel { Dock = DockStyle.Fill, Size = new Size(900, 800) };
            right.Controls.Add(geomView);
            right.Controls.Add(hint);
            right.Controls.Add(meshBar);
            right.Controls.Add(geomLog);
            page.Controls.Add(right);
            page.Controls.Add(column);
            return page;
        }

        void UpdateRefineControls()
        {
            foreach (Control c in new Control[] { cbMaxLevel, tbGapCells, tbCurvCells, cbWallLevel })
            {
                c.Enabled = chkRefine.Checked;
            }
        }

        void LogGeom(string text)
        {
            geomLog.AppendText(text + Environment.NewLine);
        }

        void LoadGeometryToUi()
        {
            loadingGeometry = true;
            tbSource.Text = geom.SourceFile;
            int ui = Array.FindIndex(UnitScales, s => Math.Abs(s - geom.UnitScale) < 1e-12);
            cbUnits.SelectedIndex = ui < 0 ? 0 : ui;
            rbInside.Checked = geom.Region == FluidRegion.InsideBody;
            rbOutside.Checked = geom.Region == FluidRegion.OutsideBody;
            tbMargin.Text = Fmt.D(geom.Margin);
            tbDx.Text = Fmt.D(geom.Dx);
            cbBlock.SelectedItem = geom.BlockSize.ToString();
            chkRefine.Checked = geom.Refine;
            cbMaxLevel.SelectedIndex = Math.Max(0, Math.Min(cbMaxLevel.Items.Count - 1, geom.MaxLevel - 1));
            tbGapCells.Text = Fmt.D(geom.GapCells);
            tbCurvCells.Text = Fmt.D(geom.CurvatureCells);
            cbWallLevel.SelectedIndex = Math.Max(0, Math.Min(cbWallLevel.Items.Count - 1, geom.WallLevel));
            UpdateRefineControls();
            chkInlet.Checked = geom.Inlet.Enabled;
            cbInletFace.SelectedIndex = (int)geom.Inlet.Face;
            tbInletU.Text = Fmt.D(geom.Inlet.U);
            tbInletV.Text = Fmt.D(geom.Inlet.V);
            tbInletR.Text = Fmt.D(geom.Inlet.Radius);
            tbInletVel.Text = Fmt.D(geom.Inlet.Value);
            chkOutlet.Checked = geom.Outlet.Enabled;
            cbOutletFace.SelectedIndex = (int)geom.Outlet.Face;
            tbOutletU.Text = Fmt.D(geom.Outlet.U);
            tbOutletV.Text = Fmt.D(geom.Outlet.V);
            tbOutletR.Text = Fmt.D(geom.Outlet.Radius);
            tbOutletP.Text = Fmt.D(geom.Outlet.Value);
            tbGx.Text = Fmt.D(geom.Gravity.X);
            tbGy.Text = Fmt.D(geom.Gravity.Y);
            tbGz.Text = Fmt.D(geom.Gravity.Z);
            loadingGeometry = false;
            geomLog.Clear();
            ResetMeshDisplay();
            ApplyUnits();
            if (mesh != null)
            {
                geomView.FitView(mesh.Min, mesh.Max);
            }
        }

        /// <summary>Reads the UI into geom, throws FormatException for invalid input.</summary>
        void ReadGeometryFromUi()
        {
            geom.UnitScale = UnitScales[Math.Max(0, cbUnits.SelectedIndex)];
            geom.Region = rbOutside.Checked ? FluidRegion.OutsideBody : FluidRegion.InsideBody;
            geom.Margin = D(tbMargin, "Margin");
            geom.Dx = D(tbDx, "Cell size");
            geom.BlockSize = int.Parse((string)cbBlock.SelectedItem);
            geom.Refine = chkRefine.Checked;
            geom.MaxLevel = cbMaxLevel.SelectedIndex + 1;
            geom.GapCells = D(tbGapCells, "Cells across a gap");
            geom.CurvatureCells = D(tbCurvCells, "Cells per curvature radius");
            geom.WallLevel = cbWallLevel.SelectedIndex;
            if (geom.Refine && geom.GapCells <= 0 && geom.CurvatureCells <= 0 && geom.WallLevel == 0)
            {
                throw new FormatException("Refinement has no criterion – set the cells across a gap, the cells per curvature radius or the level at walls.");
            }
            geom.Inlet = new PortSpec { Enabled = chkInlet.Checked, Face = (Face)cbInletFace.SelectedIndex, U = D(tbInletU, "Inlet – coordinate 1"), V = D(tbInletV, "Inlet – coordinate 2"), Radius = D(tbInletR, "Inlet – radius"), Value = D(tbInletVel, "Inlet – velocity") };
            geom.Outlet = new PortSpec { Enabled = chkOutlet.Checked, Face = (Face)cbOutletFace.SelectedIndex, U = D(tbOutletU, "Outlet – coordinate 1"), V = D(tbOutletV, "Outlet – coordinate 2"), Radius = D(tbOutletR, "Outlet – radius"), Value = D(tbOutletP, "Outlet – pressure") };
            geom.Gravity = new Vec3(D(tbGx, "gx"), D(tbGy, "gy"), D(tbGz, "gz"));
            if (geom.Dx <= 0)
            {
                throw new FormatException("Cell size must be positive.");
            }
        }

        void BrowseGeometry()
        {
            if (geom == null)
            {
                MessageBox.Show(this, "First create or select a geometry on the left.", Text);
                return;
            }
            using (var d = new OpenFileDialog { Filter = "CAD (*.step;*.stp;*.stl)|*.step;*.stp;*.stl|STEP (*.step;*.stp)|*.step;*.stp|STL (*.stl)|*.stl" })
            {
                if (d.ShowDialog(this) != DialogResult.OK)
                {
                    return;
                }
                ImportGeometryFile(d.FileName);
            }
        }

        /// <summary>Converts (STEP) or copies (STL) the file to source.stl of the geometry.</summary>
        public void ImportGeometryFile(string file)
        {
            string dst = Path.Combine(geom.Dir(settings.Workspace), "source.stl");
            Directory.CreateDirectory(geom.Dir(settings.Workspace));
            Cursor = Cursors.WaitCursor;
            try
            {
                if (CadConverter.IsStep(file))
                {
                    LogGeom("Converting STEP → STL …");
                    CadConverter.StepToStl(settings, file, dst, 0, LogGeom);
                    cbUnits.SelectedIndex = 0; // the converters write millimetres
                }
                else
                {
                    File.Copy(file, dst, true);
                }
                rawMesh = StlMesh.Read(dst);
                geom.SourceFile = file;
                tbSource.Text = file;
                LogGeom("Loaded " + rawMesh.TriangleCount + " triangles.");
                ApplyUnits();
                if (mesh != null)
                {
                    geomView.FitView(mesh.Min, mesh.Max);
                    // sensible defaults: about 60 cells along the longest side, inlet centred on the top face
                    Vec3 size = mesh.Max - mesh.Min;
                    double longest = Math.Max(size.X, Math.Max(size.Y, size.Z));
                    tbDx.Text = Fmt.G(longest / 60, 2);
                    CenterPort(cbInletFace, tbInletU, tbInletV);
                    CenterPort(cbOutletFace, tbOutletU, tbOutletV);
                    // a sprue-like inlet: 15 % of the smaller side of the face
                    int axis = cbInletFace.SelectedIndex / 2;
                    double side = Math.Min(size[axis == 0 ? 1 : 0], size[axis == 2 ? 1 : 2]);
                    tbInletR.Text = Fmt.G(0.15 * side, 2);
                }
                ReadGeometryFromUi();
                geom.MeshReady = false;
                geom.Save(settings.Workspace);
            }
            catch (Exception ex)
            {
                LogGeom("Error: " + ex.Message);
                MessageBox.Show(this, ex.Message, Text, MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
            finally
            {
                Cursor = Cursors.Default;
            }
        }

        void ApplyUnits()
        {
            if (rawMesh == null)
            {
                mesh = null;
                lblGeomInfo.Text = "No model – load a STEP or STL file.";
                UpdateGeometryPreview();
                return;
            }
            mesh = new StlMesh { V = (double[])rawMesh.V.Clone() };
            mesh.Scale(UnitScales[Math.Max(0, cbUnits.SelectedIndex)]);
            Vec3 s = mesh.Max - mesh.Min;
            lblGeomInfo.Text = string.Format("Size: {0} × {1} × {2} m, volume {3} m³, {4} triangles",
                Fmt.G(s.X), Fmt.G(s.Y), Fmt.G(s.Z), Fmt.G(Math.Abs(mesh.SignedVolume())), mesh.TriangleCount);
            UpdateGeometryPreview();
        }

        void CenterPort(ComboBox face, TextBox u, TextBox v)
        {
            if (mesh == null)
            {
                return;
            }
            int axis = face.SelectedIndex / 2;
            int a1 = axis == 0 ? 1 : 0, a2 = axis == 2 ? 1 : 2;
            u.Text = Fmt.G((mesh.Min[a1] + mesh.Max[a1]) / 2, 6);
            v.Text = Fmt.G((mesh.Min[a2] + mesh.Max[a2]) / 2, 6);
        }

        /// <summary>Domain box of the mesh (bounding box of the body plus the margin for flow around a body).</summary>
        void DomainBox(out Vec3 min, out Vec3 max)
        {
            double m = 0;
            if (rbOutside.Checked)
            {
                Fmt.TryParse(tbMargin.Text, out m);
            }
            min = mesh.Min - new Vec3(m, m, m);
            max = mesh.Max + new Vec3(m, m, m);
        }

        void UpdateGeometryPreview()
        {
            if (loadingGeometry || geomView == null)
            {
                return;
            }
            geomView.Clear();
            UpdateMeshLayers();
            if (mesh == null)
            {
                lblMeshEstimate.Text = "";
                return;
            }
            Vec3 min, max;
            DomainBox(out min, out max);
            double size = (max - min).Length;
            if (chkShowModel.Checked)
            {
                geomView.SetLayer("body", Viewer3D.Triangles(mesh.V, System.Windows.Media.Color.FromRgb(170, 175, 185), rbOutside.Checked ? 1.0 : 0.45), !rbOutside.Checked);
            }
            geomView.SetLayer("box", Viewer3D.BoxEdges(min, max, size * 0.0015, Colors.DimGray));
            double dx;
            if (Fmt.TryParse(tbDx.Text, out dx) && dx > 0)
            {
                Vec3 s = max - min;
                double cells = Math.Ceiling(s.X / dx + 2) * Math.Ceiling(s.Y / dx + 2) * Math.Ceiling(s.Z / dx + 2);
                lblMeshEstimate.Text = string.Format("Grid approx. {0} × {1} × {2} = {3:N0} cells",
                    Math.Ceiling(s.X / dx + 2), Math.Ceiling(s.Y / dx + 2), Math.Ceiling(s.Z / dx + 2), cells);
                if (chkRefine.Checked)
                {
                    int lmax = cbMaxLevel.SelectedIndex + 1;
                    lblMeshEstimate.Text += string.Format(" on level 0, smallest cell {0} m (level {1})", Fmt.G(dx / (1 << lmax)), lmax);
                }
            }
            try
            {
                var inlet = new PortSpec { Face = (Face)cbInletFace.SelectedIndex, U = D(tbInletU, ""), V = D(tbInletV, ""), Radius = D(tbInletR, "") };
                if (chkInlet.Checked)
                {
                    geomView.SetLayer("inlet", Viewer3D.Port(inlet, min, max, Colors.Red), true);
                }
                var outlet = new PortSpec { Face = (Face)cbOutletFace.SelectedIndex, U = D(tbOutletU, ""), V = D(tbOutletV, ""), Radius = D(tbOutletR, "") };
                if (chkOutlet.Checked)
                {
                    geomView.SetLayer("outlet", Viewer3D.Port(outlet, min, max, Colors.LimeGreen), true);
                }
                var g = new Vec3(D(tbGx, ""), D(tbGy, ""), D(tbGz, ""));
                if (g.Length > 0)
                {
                    Vec3 start = new Vec3(max.X + size * 0.08, (min.Y + max.Y) / 2, max.Z);
                    geomView.SetLayer("gravity", Viewer3D.Arrow(start, g * (size * 0.3 / g.Length), size * 0.006, Colors.RoyalBlue));
                }
            }
            catch (FormatException)
            {
            }
        }

        void GenerateMesh()
        {
            if (geom == null || mesh == null)
            {
                MessageBox.Show(this, "Load the geometry first.", Text);
                return;
            }
            try
            {
                ReadGeometryFromUi();
            }
            catch (FormatException ex)
            {
                MessageBox.Show(this, ex.Message, Text);
                return;
            }
            var g = geom;
            var m = mesh;
            string ws = settings.Workspace;
            btnGenerate.Enabled = false;
            Cursor = Cursors.WaitCursor;
            Action<string> log = t => BeginInvokeSafe(() => LogGeom(t));
            Task.Run(() =>
            {
                m.WriteBinary(g.StlPath(ws));
                CaseGenerator.Generate(g, m, ws, log);
            }).ContinueWith(t => BeginInvokeSafe(() =>
            {
                btnGenerate.Enabled = true;
                Cursor = Cursors.Default;
                if (t.Exception != null)
                {
                    string msg = t.Exception.InnerException.Message;
                    LogGeom("Error: " + msg);
                    MessageBox.Show(this, msg, Text, MessageBoxButtons.OK, MessageBoxIcon.Error);
                    return;
                }
                LogGeom("The mesh is ready. Continue on the Simulation tab." + (g.MeshLevels > 1 ? " Blocks are colored by level (0 gray, 1 blue, 2 green, 3 orange, 4 red)." : ""));
                SetStatus("Mesh: " + g.BlockCount + " blocks, " + g.FluidCells + " fluid cells");
                ResetMeshDisplay();
                UpdateMeshLayers();
                RefreshSimulationList(null);
            }));
        }

        /// <summary>Synchronous mesh generation (self test).</summary>
        public void GenerateMeshSync()
        {
            ReadGeometryFromUi();
            mesh.WriteBinary(geom.StlPath(settings.Workspace));
            CaseGenerator.Generate(geom, mesh, settings.Workspace, LogGeom);
            ResetMeshDisplay();
            UpdateMeshLayers();
        }

        // ---------------------------------------------------------------- mesh display

        /// <summary>Drops the cached mesh display, called when the geometry or its mesh changes.</summary>
        void ResetMeshDisplay()
        {
            meshTypes = null;
            meshTree = null;
            cellsModel = null;
            blocksModel = null;
            bool ready = geom != null && geom.MeshReady;
            chkShowCells.Enabled = ready;
            chkShowBlocks.Enabled = ready;
            cbMeshCut.Enabled = ready;
            ResetMeshCut();
        }

        void ResetMeshCut()
        {
            int axis = cbMeshCut.SelectedIndex - 1;
            bool on = axis >= 0 && geom != null && geom.MeshReady;
            tbMeshCut.Enabled = on;
            chkMeshCutFlip.Enabled = on;
            lblMeshCut.Text = "";
            if (!on)
            {
                return;
            }
            int n = axis == 0 ? geom.Nx : (axis == 1 ? geom.Ny : geom.Nz);
            updatingMeshCut = true;
            tbMeshCut.Maximum = Math.Max(0, n - 1);
            tbMeshCut.Value = n / 2;
            updatingMeshCut = false;
        }

        /// <summary>Shows or hides the cells and blocks of the mesh in the 3D view, builds them when needed.</summary>
        void UpdateMeshLayers()
        {
            if (geomView == null)
            {
                return;
            }
            bool ready = geom != null && geom.MeshReady;
            if (ready && (chkShowCells.Checked || chkShowBlocks.Checked) && meshTypes == null)
            {
                meshTypes = geom.LoadCellTypes(settings.Workspace);
                if (meshTypes != null && meshTypes.LongLength != (long)geom.GX * geom.GY * geom.GZ)
                {
                    LogGeom("File mesh/celltypes.bin does not match the mesh, generate the mesh again.");
                    meshTypes = null;
                }
                try
                {
                    meshTree = geom.LoadOctree(settings.Workspace);
                }
                catch (Exception ex)
                {
                    LogGeom("Error reading mesh/octree.bin: " + ex.Message);
                    meshTree = null;
                }
            }
            ready = ready && meshTypes != null;
            int axis = cbMeshCut.SelectedIndex - 1;
            // by default the half towards the initial camera position (+x, -y, +z) is removed
            bool keepUpper = (axis == 1) != chkMeshCutFlip.Checked;
            if (ready && axis >= 0)
            {
                string[] names = { "x", "y", "z" };
                int plane = keepUpper ? tbMeshCut.Value : tbMeshCut.Value + 1;
                lblMeshCut.Text = names[axis] + (keepUpper ? " ≥ " : " ≤ ") + Fmt.G(geom.Origin[axis] + plane * geom.Dx, 4) + " m";
            }
            if (ready && chkShowCells.Checked && cellsModel == null)
            {
                long faces;
                Cursor = Cursors.WaitCursor;
                try
                {
                    cellsModel = meshTree != null
                        ? MeshDisplay.OctreeCells(geom, meshTree, axis, tbMeshCut.Value, keepUpper, out faces)
                        : MeshDisplay.Cells(geom, meshTypes, axis, tbMeshCut.Value, keepUpper, out faces);
                }
                finally
                {
                    Cursor = Cursors.Default;
                }
                if (cellsModel == null)
                {
                    LogGeom(string.Format("The mesh surface has {0:N0} cell faces, at most {1:N0} can be shown. Use a cut or a coarser mesh.", faces, MeshDisplay.MaxFaces));
                }
            }
            if (ready && chkShowBlocks.Checked && blocksModel == null)
            {
                blocksModel = meshTree != null ? MeshDisplay.OctreeBlockEdges(geom, meshTree) : MeshDisplay.BlockEdges(geom, meshTypes, Colors.DarkSlateGray);
            }
            geomView.SetLayer("cells", ready && chkShowCells.Checked ? cellsModel : null);
            geomView.SetLayer("blocks", ready && chkShowBlocks.Checked ? blocksModel : null);
        }
    }
}
