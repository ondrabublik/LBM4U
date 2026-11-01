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
        System.Windows.Media.Media3D.Model3D cellsModel, blocksModel; // cached, rebuilt when the mesh or the cut changes
        bool updatingMeshCut;

        static readonly string[] FaceNames = { "x min", "x max", "y min", "y max", "z min", "z max" };
        static readonly string[] UnitNames = { "mm", "m", "cm", "palce" };
        static readonly double[] UnitScales = { 0.001, 1, 0.01, 0.0254 };

        TabPage BuildGeometryTab()
        {
            var page = new TabPage("1. Geometrie a síť");

            // CAD
            var tCad = FieldTable();
            tbSource = new TextBox { ReadOnly = true };
            var btnBrowse = new Button { Text = "Načíst STEP / STL…", Height = 28 };
            btnBrowse.Click += (s, e) => BrowseGeometry();
            cbUnits = Combo(UnitNames);
            cbUnits.SelectedIndexChanged += (s, e) => { if (!loadingGeometry) ApplyUnits(); };
            tCad.Controls.Add(btnBrowse);
            tCad.SetColumnSpan(btnBrowse, 3);
            btnBrowse.Dock = DockStyle.Fill;
            AddRow(tCad, "Soubor", tbSource);
            AddRow(tCad, "Jednotky souboru", cbUnits);
            lblGeomInfo = new Label { AutoSize = true, ForeColor = Color.DimGray, MaximumSize = new Size(390, 0) };
            tCad.Controls.Add(lblGeomInfo);
            tCad.SetColumnSpan(lblGeomInfo, 3);

            // fluid region
            var tReg = FieldTable();
            rbInside = new RadioButton { Text = "Tekutina uvnitř tělesa (dutina, forma)", AutoSize = true, Checked = true };
            rbOutside = new RadioButton { Text = "Těleso obtékané tekutinou (kvádr kolem)", AutoSize = true };
            tReg.Controls.Add(rbInside);
            tReg.SetColumnSpan(rbInside, 3);
            tReg.Controls.Add(rbOutside);
            tReg.SetColumnSpan(rbOutside, 3);
            tbMargin = AddRow(tReg, "Okraj kolem tělesa", Num("0.02"), "m");
            rbOutside.CheckedChanged += (s, e) => { tbMargin.Enabled = rbOutside.Checked; UpdateGeometryPreview(); };
            tbMargin.Enabled = false;

            // mesh
            var tMesh = FieldTable();
            tbDx = AddRow(tMesh, "Velikost buňky dx", Num("0.002"), "m");
            cbBlock = AddRow(tMesh, "Velikost bloku", Combo("16", "8", "32"), "buněk");
            lblMeshEstimate = new Label { AutoSize = true, ForeColor = Color.DimGray, MaximumSize = new Size(390, 0) };
            tMesh.Controls.Add(lblMeshEstimate);
            tMesh.SetColumnSpan(lblMeshEstimate, 3);
            tbDx.TextChanged += (s, e) => UpdateGeometryPreview();
            tbMargin.TextChanged += (s, e) => UpdateGeometryPreview();

            // inlet
            var tIn = FieldTable();
            chkInlet = new CheckBox { Text = "Vtok (rychlostní)", AutoSize = true, Checked = true };
            tIn.Controls.Add(chkInlet);
            tIn.SetColumnSpan(chkInlet, 3);
            cbInletFace = AddRow(tIn, "Stěna kvádru", Combo(FaceNames));
            tbInletU = AddRow(tIn, "Střed – souřadnice 1", Num("0"), "m");
            tbInletV = AddRow(tIn, "Střed – souřadnice 2", Num("0"), "m");
            tbInletR = AddRow(tIn, "Poloměr (0 = celá stěna)", Num("0"), "m");
            tbInletVel = AddRow(tIn, "Rychlost (výchozí)", Num("0.5"), "m/s");
            var btnCenterIn = new Button { Text = "Střed na stěnu", AutoSize = true };
            btnCenterIn.Click += (s, e) => CenterPort(cbInletFace, tbInletU, tbInletV);
            tIn.Controls.Add(btnCenterIn);
            tIn.SetColumnSpan(btnCenterIn, 3);
            cbInletFace.SelectedIndex = 5;

            // outlet
            var tOut = FieldTable();
            chkOutlet = new CheckBox { Text = "Výtok (tlakový)", AutoSize = true };
            tOut.Controls.Add(chkOutlet);
            tOut.SetColumnSpan(chkOutlet, 3);
            cbOutletFace = AddRow(tOut, "Stěna kvádru", Combo(FaceNames));
            tbOutletU = AddRow(tOut, "Střed – souřadnice 1", Num("0"), "m");
            tbOutletV = AddRow(tOut, "Střed – souřadnice 2", Num("0"), "m");
            tbOutletR = AddRow(tOut, "Poloměr (0 = celá stěna)", Num("0"), "m");
            tbOutletP = AddRow(tOut, "Přetlak (výchozí)", Num("0"), "Pa");
            var btnCenterOut = new Button { Text = "Střed na stěnu", AutoSize = true };
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

            btnGenerate = new Button { Text = "Vytvořit výpočetní síť", Dock = DockStyle.Top, Height = 38, Font = new Font(Font, FontStyle.Bold) };
            btnGenerate.Click += (s, e) => GenerateMesh();

            var column = Column(430,
                Group("CAD model", tCad),
                Group("Oblast tekutiny", tReg),
                Group("Síť", tMesh),
                Group("Vtok", tIn),
                Group("Výtok", tOut),
                Group("Gravitace (výchozí pro simulace)", tG),
                btnGenerate);

            geomView = new Viewer3D { Dock = DockStyle.Fill };
            geomLog = new TextBox { Dock = DockStyle.Bottom, Height = 120, Multiline = true, ScrollBars = ScrollBars.Vertical, ReadOnly = true, Font = new Font("Consolas", 8.5f) };
            var hint = new Label { Dock = DockStyle.Top, Height = 20, ForeColor = Color.DimGray, Text = "  Levé tlačítko – otáčení, pravé – posun, kolečko – zoom. Červená = vtok, zelená = výtok, modrá šipka = gravitace, modré buňky = tekutina." };

            // display of the solver mesh
            var meshBar = new FlowLayoutPanel { Dock = DockStyle.Top, AutoSize = true, AutoSizeMode = AutoSizeMode.GrowAndShrink, Padding = new Padding(4, 2, 4, 0), WrapContents = true };
            chkShowModel = new CheckBox { Text = "CAD model", AutoSize = true, Checked = true, Margin = new Padding(3, 6, 8, 0) };
            chkShowCells = new CheckBox { Text = "Buňky sítě", AutoSize = true, Margin = new Padding(3, 6, 8, 0), Enabled = false };
            chkShowBlocks = new CheckBox { Text = "Bloky", AutoSize = true, Margin = new Padding(3, 6, 8, 0), Enabled = false };
            cbMeshCut = Combo("bez řezu", "řez kolmý na x", "řez kolmý na y", "řez kolmý na z");
            cbMeshCut.Width = 120;
            cbMeshCut.Enabled = false;
            tbMeshCut = new TrackBar { Width = 200, Minimum = 0, Maximum = 0, TickStyle = TickStyle.None, Height = 30, Enabled = false };
            lblMeshCut = new Label { AutoSize = true, ForeColor = Color.DimGray, Margin = new Padding(3, 6, 3, 0) };
            chkMeshCutFlip = new CheckBox { Text = "obrátit", AutoSize = true, Margin = new Padding(3, 6, 8, 0), Enabled = false };
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
            meshBar.Controls.AddRange(new Control[] { new Label { Text = "Zobrazit:", AutoSize = true, Margin = new Padding(3, 6, 3, 0) }, chkShowModel, chkShowCells, chkShowBlocks, cbMeshCut, tbMeshCut, chkMeshCutFlip, lblMeshCut });

            var right = new Panel { Dock = DockStyle.Fill, Size = new Size(900, 800) };
            right.Controls.Add(geomView);
            right.Controls.Add(hint);
            right.Controls.Add(meshBar);
            right.Controls.Add(geomLog);
            page.Controls.Add(right);
            page.Controls.Add(column);
            return page;
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
            geom.Margin = D(tbMargin, "Okraj");
            geom.Dx = D(tbDx, "Velikost buňky");
            geom.BlockSize = int.Parse((string)cbBlock.SelectedItem);
            geom.Inlet = new PortSpec { Enabled = chkInlet.Checked, Face = (Face)cbInletFace.SelectedIndex, U = D(tbInletU, "Vtok – souřadnice 1"), V = D(tbInletV, "Vtok – souřadnice 2"), Radius = D(tbInletR, "Vtok – poloměr"), Value = D(tbInletVel, "Vtok – rychlost") };
            geom.Outlet = new PortSpec { Enabled = chkOutlet.Checked, Face = (Face)cbOutletFace.SelectedIndex, U = D(tbOutletU, "Výtok – souřadnice 1"), V = D(tbOutletV, "Výtok – souřadnice 2"), Radius = D(tbOutletR, "Výtok – poloměr"), Value = D(tbOutletP, "Výtok – tlak") };
            geom.Gravity = new Vec3(D(tbGx, "gx"), D(tbGy, "gy"), D(tbGz, "gz"));
            if (geom.Dx <= 0)
            {
                throw new FormatException("Velikost buňky musí být kladná.");
            }
        }

        void BrowseGeometry()
        {
            if (geom == null)
            {
                MessageBox.Show(this, "Nejprve vytvořte nebo vyberte geometrii vlevo.", Text);
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
                    LogGeom("Převod STEP → STL …");
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
                LogGeom("Načteno " + rawMesh.TriangleCount + " trojúhelníků.");
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
                LogGeom("Chyba: " + ex.Message);
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
                lblGeomInfo.Text = "Žádný model – načtěte soubor STEP nebo STL.";
                UpdateGeometryPreview();
                return;
            }
            mesh = new StlMesh { V = (double[])rawMesh.V.Clone() };
            mesh.Scale(UnitScales[Math.Max(0, cbUnits.SelectedIndex)]);
            Vec3 s = mesh.Max - mesh.Min;
            lblGeomInfo.Text = string.Format("Rozměry: {0} × {1} × {2} m, objem {3} m³, {4} trojúhelníků",
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
                lblMeshEstimate.Text = string.Format("Mřížka přibližně {0} × {1} × {2} = {3:N0} buněk",
                    Math.Ceiling(s.X / dx + 2), Math.Ceiling(s.Y / dx + 2), Math.Ceiling(s.Z / dx + 2), cells);
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
                MessageBox.Show(this, "Nejprve načtěte geometrii.", Text);
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
                    LogGeom("Chyba: " + msg);
                    MessageBox.Show(this, msg, Text, MessageBoxButtons.OK, MessageBoxIcon.Error);
                    return;
                }
                LogGeom("Síť je hotová. Pokračujte na záložce Simulace.");
                SetStatus("Síť: " + g.BlockCount + " bloků, " + g.FluidCells + " buněk tekutiny");
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
                    LogGeom("Soubor mesh/celltypes.bin neodpovídá síti, vytvořte síť znovu.");
                    meshTypes = null;
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
                    cellsModel = MeshDisplay.Cells(geom, meshTypes, axis, tbMeshCut.Value, keepUpper, out faces);
                }
                finally
                {
                    Cursor = Cursors.Default;
                }
                if (cellsModel == null)
                {
                    LogGeom(string.Format("Povrch sítě má {0:N0} stěn buněk, zobrazit lze nejvýše {1:N0}. Použijte řez nebo hrubší síť.", faces, MeshDisplay.MaxFaces));
                }
            }
            if (ready && chkShowBlocks.Checked && blocksModel == null)
            {
                blocksModel = MeshDisplay.BlockEdges(geom, meshTypes, Colors.DarkSlateGray);
            }
            geomView.SetLayer("cells", ready && chkShowCells.Checked ? cellsModel : null);
            geomView.SetLayer("blocks", ready && chkShowBlocks.Checked ? blocksModel : null);
        }
    }
}
