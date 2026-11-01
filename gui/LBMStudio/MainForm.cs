using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Windows.Forms;

namespace LBMStudio
{
    public partial class MainForm : Form
    {
        AppSettings settings;
        GeometryDef geom;
        StlMesh rawMesh;   // mesh in the units of the source file
        StlMesh mesh;      // mesh in metres
        readonly SolverRunner runner = new SolverRunner();

        ListBox lbGeometries;
        TabControl tabs;
        ToolStripLabel lblWorkspace;
        StatusStrip status;
        ToolStripStatusLabel statusText;

        public MainForm() : this(AppSettings.Load())
        {
        }

        public MainForm(AppSettings settings)
        {
            this.settings = settings;
            Text = "LBM Studio – lattice Boltzmann simulace proudění";
            Width = 1500;
            Height = 950;
            StartPosition = FormStartPosition.CenterScreen;
            Font = new Font("Segoe UI", 9f);

            var tool = new ToolStrip { GripStyle = ToolStripGripStyle.Hidden, Padding = new Padding(6, 2, 6, 2) };
            var btnSettings = new ToolStripButton("Nastavení…");
            btnSettings.Click += (s, e) => ShowSettings();
            var btnOpenDir = new ToolStripButton("Otevřít pracovní složku");
            btnOpenDir.Click += (s, e) => OpenFolder(settings.Workspace);
            lblWorkspace = new ToolStripLabel();
            tool.Items.Add(btnSettings);
            tool.Items.Add(btnOpenDir);
            tool.Items.Add(new ToolStripSeparator());
            tool.Items.Add(lblWorkspace);

            status = new StatusStrip();
            statusText = new ToolStripStatusLabel("Připraveno") { Spring = true, TextAlign = ContentAlignment.MiddleLeft };
            status.Items.Add(statusText);

            var split = new SplitContainer { Dock = DockStyle.Fill, FixedPanel = FixedPanel.Panel1, SplitterDistance = 220 };
            var leftPanel = new Panel { Dock = DockStyle.Fill, Padding = new Padding(6) };
            var lblG = new Label { Text = "Geometrie", Dock = DockStyle.Top, Height = 22, Font = new Font(Font, FontStyle.Bold) };
            lbGeometries = new ListBox { Dock = DockStyle.Fill, IntegralHeight = false };
            lbGeometries.SelectedIndexChanged += (s, e) => SelectGeometry(lbGeometries.SelectedItem as string);
            var btnNewGeom = new Button { Text = "Nová geometrie…", Dock = DockStyle.Bottom, Height = 30 };
            btnNewGeom.Click += (s, e) => NewGeometry();
            leftPanel.Controls.Add(lbGeometries);
            leftPanel.Controls.Add(lblG);
            leftPanel.Controls.Add(btnNewGeom);
            split.Panel1.Controls.Add(leftPanel);

            tabs = new TabControl { Dock = DockStyle.Fill, Padding = new Point(14, 5) };
            tabs.TabPages.Add(BuildGeometryTab());
            tabs.TabPages.Add(BuildSimulationTab());
            tabs.TabPages.Add(BuildResultsTab());
            tabs.SelectedIndexChanged += (s, e) =>
            {
                if (tabs.SelectedIndex == 2)
                {
                    RefreshResultsSimulations();
                }
            };
            split.Panel2.Controls.Add(tabs);

            Controls.Add(split);
            Controls.Add(tool);
            Controls.Add(status);

            runner.Line += line => BeginInvokeSafe(() => AppendSimLog(line));
            runner.Progress += p => BeginInvokeSafe(() => OnSolverProgress(p));
            runner.Exited += code => BeginInvokeSafe(() => OnSolverExited(code));
            FormClosing += (s, e) =>
            {
                if (runner.Running && MessageBox.Show(this, "Výpočet stále běží. Ukončit aplikaci a výpočet zastavit?", Text, MessageBoxButtons.YesNo, MessageBoxIcon.Question) != DialogResult.Yes)
                {
                    e.Cancel = true;
                    return;
                }
                runner.Kill();
            };
            Load += (s, e) => RefreshGeometries(null);
        }

        void BeginInvokeSafe(Action a)
        {
            if (IsDisposed || !IsHandleCreated)
            {
                return;
            }
            try
            {
                BeginInvoke(a);
            }
            catch (InvalidOperationException)
            {
            }
        }

        void SetStatus(string text)
        {
            statusText.Text = text;
        }

        static void OpenFolder(string dir)
        {
            if (Directory.Exists(dir))
            {
                Process.Start("explorer.exe", "\"" + dir + "\"");
            }
        }

        // ---------------------------------------------------------------- geometry list

        void RefreshGeometries(string select)
        {
            Directory.CreateDirectory(settings.Workspace);
            lblWorkspace.Text = "Pracovní složka: " + settings.Workspace;
            var names = Directory.GetDirectories(settings.Workspace)
                .Where(d => File.Exists(Path.Combine(d, "geometry.ini")))
                .Select(Path.GetFileName).OrderBy(n => n).ToList();
            lbGeometries.BeginUpdate();
            lbGeometries.Items.Clear();
            foreach (var n in names)
            {
                lbGeometries.Items.Add(n);
            }
            lbGeometries.EndUpdate();
            if (select != null && names.Contains(select))
            {
                lbGeometries.SelectedItem = select;
            }
            else if (names.Count > 0 && lbGeometries.SelectedIndex < 0)
            {
                lbGeometries.SelectedIndex = 0;
            }
        }

        void NewGeometry()
        {
            string name = Prompt("Nová geometrie", "Název geometrie (bez mezer a diakritiky):", "geometrie1");
            if (string.IsNullOrEmpty(name))
            {
                return;
            }
            if (name.IndexOfAny(Path.GetInvalidFileNameChars()) >= 0 || name.Contains(" "))
            {
                MessageBox.Show(this, "Neplatný název.", Text);
                return;
            }
            if (Directory.Exists(Path.Combine(settings.Workspace, name)))
            {
                MessageBox.Show(this, "Geometrie s tímto názvem už existuje.", Text);
                return;
            }
            var g = new GeometryDef { Name = name };
            g.Save(settings.Workspace);
            RefreshGeometries(name);
        }

        void SelectGeometry(string name)
        {
            if (name == null)
            {
                return;
            }
            geom = GeometryDef.Load(Path.Combine(settings.Workspace, name));
            rawMesh = null;
            mesh = null;
            string src = Path.Combine(geom.Dir(settings.Workspace), "source.stl");
            if (File.Exists(src))
            {
                try
                {
                    rawMesh = StlMesh.Read(src);
                }
                catch (Exception ex)
                {
                    LogGeom("Chyba čtení geometrie: " + ex.Message);
                }
            }
            LoadGeometryToUi();
            RefreshSimulationList(null);
            SetStatus("Geometrie " + name + (geom.MeshReady ? " – síť: " + geom.BlockCount + " bloků, " + geom.FluidCells + " buněk tekutiny" : " – bez sítě"));
        }

        // ---------------------------------------------------------------- dialogs

        public static string Prompt(string title, string text, string value)
        {
            using (var f = new Form { Text = title, Width = 420, Height = 160, FormBorderStyle = FormBorderStyle.FixedDialog, StartPosition = FormStartPosition.CenterParent, MinimizeBox = false, MaximizeBox = false })
            {
                var l = new Label { Text = text, Left = 12, Top = 12, Width = 380 };
                var t = new TextBox { Text = value, Left = 12, Top = 36, Width = 380 };
                var ok = new Button { Text = "OK", DialogResult = DialogResult.OK, Left = 236, Top = 72, Width = 75 };
                var cancel = new Button { Text = "Zrušit", DialogResult = DialogResult.Cancel, Left = 317, Top = 72, Width = 75 };
                f.Controls.AddRange(new Control[] { l, t, ok, cancel });
                f.AcceptButton = ok;
                f.CancelButton = cancel;
                return f.ShowDialog() == DialogResult.OK ? t.Text.Trim() : null;
            }
        }

        void ShowSettings()
        {
            using (var f = new Form { Text = "Nastavení", Width = 720, Height = 300, FormBorderStyle = FormBorderStyle.FixedDialog, StartPosition = FormStartPosition.CenterParent, MinimizeBox = false, MaximizeBox = false })
            {
                var table = new TableLayoutPanel { Dock = DockStyle.Top, ColumnCount = 3, Height = 200, Padding = new Padding(8) };
                table.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 170));
                table.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100));
                table.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 40));
                var rows = new List<Tuple<string, TextBox, bool>>();
                Func<string, string, bool, TextBox> add = (label, value, folder) =>
                {
                    var tb = new TextBox { Text = value, Dock = DockStyle.Fill };
                    var b = new Button { Text = "…", Width = 30 };
                    b.Click += (s, e) =>
                    {
                        if (folder)
                        {
                            using (var d = new FolderBrowserDialog { SelectedPath = tb.Text })
                            {
                                if (d.ShowDialog() == DialogResult.OK)
                                {
                                    tb.Text = d.SelectedPath;
                                }
                            }
                        }
                        else
                        {
                            using (var d = new OpenFileDialog { Filter = "Programy (*.exe)|*.exe" })
                            {
                                if (d.ShowDialog() == DialogResult.OK)
                                {
                                    tb.Text = d.FileName;
                                }
                            }
                        }
                    };
                    table.Controls.Add(new Label { Text = label, AutoSize = true, Anchor = AnchorStyles.Left });
                    table.Controls.Add(tb);
                    table.Controls.Add(b);
                    return tb;
                };
                var tWork = add("Pracovní složka", settings.Workspace, true);
                var tSolver = add("Řešič cLBMnet.exe", settings.SolverExe, false);
                var tGmsh = add("gmsh.exe (STEP)", settings.GmshExe, false);
                var tFc = add("FreeCADCmd.exe (STEP)", settings.FreeCadExe, false);
                var tPv = add("paraview.exe", settings.ParaViewExe, false);
                var info = new Label { Text = "Převod STEP: " + CadConverter.Describe(settings), Left = 12, Top = 205, Width = 500 };
                var ok = new Button { Text = "Uložit", DialogResult = DialogResult.OK, Left = 520, Top = 225, Width = 80 };
                var cancel = new Button { Text = "Zrušit", DialogResult = DialogResult.Cancel, Left = 610, Top = 225, Width = 80 };
                f.Controls.AddRange(new Control[] { table, info, ok, cancel });
                f.AcceptButton = ok;
                if (f.ShowDialog(this) == DialogResult.OK)
                {
                    settings.Workspace = tWork.Text.Trim();
                    settings.SolverExe = tSolver.Text.Trim();
                    settings.GmshExe = tGmsh.Text.Trim();
                    settings.FreeCadExe = tFc.Text.Trim();
                    settings.ParaViewExe = tPv.Text.Trim();
                    settings.Save();
                    RefreshGeometries(null);
                }
            }
        }

        // ---------------------------------------------------------------- layout helpers

        static TableLayoutPanel FieldTable()
        {
            var t = new TableLayoutPanel { ColumnCount = 3, AutoSize = true, AutoSizeMode = AutoSizeMode.GrowAndShrink, Dock = DockStyle.Top, Padding = new Padding(4, 2, 4, 4) };
            t.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 150));
            t.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 130));
            t.ColumnStyles.Add(new ColumnStyle(SizeType.AutoSize));
            return t;
        }

        static T AddRow<T>(TableLayoutPanel t, string label, T control, string unit = "") where T : Control
        {
            t.Controls.Add(new Label { Text = label, AutoSize = true, Anchor = AnchorStyles.Left, Margin = new Padding(3, 6, 3, 3) });
            control.Dock = DockStyle.Fill;
            control.Margin = new Padding(3, 3, 3, 3);
            t.Controls.Add(control);
            t.Controls.Add(new Label { Text = unit, AutoSize = true, Anchor = AnchorStyles.Left, ForeColor = Color.DimGray, Margin = new Padding(0, 6, 3, 3) });
            return control;
        }

        static TextBox Num(string value)
        {
            return new TextBox { Text = value };
        }

        static ComboBox Combo(params string[] items)
        {
            var c = new ComboBox { DropDownStyle = ComboBoxStyle.DropDownList };
            c.Items.AddRange(items);
            if (items.Length > 0)
            {
                c.SelectedIndex = 0;
            }
            return c;
        }

        static GroupBox Group(string title, Control content)
        {
            var g = new GroupBox { Text = title, Dock = DockStyle.Top, AutoSize = true, AutoSizeMode = AutoSizeMode.GrowAndShrink, Padding = new Padding(6, 4, 6, 6) };
            g.Controls.Add(content);
            return g;
        }

        /// <summary>Adds groups to a scrollable column in the given order (top to bottom).</summary>
        static Panel Column(int width, params Control[] items)
        {
            var p = new Panel { Dock = DockStyle.Left, Width = width, AutoScroll = true, Padding = new Padding(4) };
            for (int i = items.Length - 1; i >= 0; i--)
            {
                p.Controls.Add(items[i]);
            }
            return p;
        }

        static double D(TextBox t, string what)
        {
            return Fmt.Parse(t.Text, what);
        }
    }
}
