using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace LBMStudio
{
    public partial class MainForm
    {
        ListBox lbSims;
        TextBox tbSimName, tbDensity, tbVisc, tbSigma, tbAngle, tbFluidParams;
        TextBox tbSGx, tbSGy, tbSGz, tbSimInletVel, tbSimOutletP;
        TextBox tbSmag, tbLatticeU, tbCharU, tbEndTime, tbSaveInt, tbThreads, tbFillLevel, tbTauMax;
        ComboBox cbProblem, cbFluid, cbViscModel, cbCollision, cbFill;
        Label lblDerived;
        Button btnSaveSim, btnRunSim, btnStopSim;
        ProgressBar simProgress;
        TextBox simLog;
        Chart simChart;
        string runningSimName;
        int runningMaxIter;
        bool loadingSim;

        // name, density, kinematic viscosity, surface tension
        static readonly object[][] FluidPresets =
        {
            new object[] { "Water 20 °C", 998.0, 1.004e-6, 0.0728 },
            new object[] { "Aluminium (melt 700 °C)", 2375.0, 5.4e-7, 0.87 },
            new object[] { "Cast iron (melt)", 7000.0, 8.0e-7, 1.5 },
            new object[] { "Engine oil", 880.0, 1.0e-4, 0.031 },
            new object[] { "Glycerol", 1261.0, 1.18e-3, 0.063 },
            new object[] { "Air (single phase)", 1.2, 1.5e-5, 0.0 },
            new object[] { "Custom", 1000.0, 1e-6, 0.0 }
        };

        TabPage BuildSimulationTab()
        {
            var page = new TabPage("2. Simulation");

            var listPanel = new Panel { Dock = DockStyle.Left, Width = 190, Padding = new Padding(6) };
            lbSims = new ListBox { Dock = DockStyle.Fill, IntegralHeight = false };
            lbSims.SelectedIndexChanged += (s, e) => LoadSimulationToUi(lbSims.SelectedItem as string);
            var btnNew = new Button { Text = "New simulation", Dock = DockStyle.Bottom, Height = 28 };
            btnNew.Click += (s, e) => NewSimulation(false);
            var btnCopy = new Button { Text = "Copy selected", Dock = DockStyle.Bottom, Height = 28 };
            btnCopy.Click += (s, e) => NewSimulation(true);
            listPanel.Controls.Add(lbSims);
            listPanel.Controls.Add(new Label { Text = "Simulations", Dock = DockStyle.Top, Height = 20, Font = new Font(Font, FontStyle.Bold) });
            listPanel.Controls.Add(btnNew);
            listPanel.Controls.Add(btnCopy);

            var tBase = FieldTable();
            tbSimName = AddRow(tBase, "Name", new TextBox());
            cbProblem = AddRow(tBase, "Problem", Combo("Free surface (filling)", "Single phase flow"));

            var tFluid = FieldTable();
            cbFluid = AddRow(tFluid, "Fluid", Combo(FluidPresets.Select(p => (string)p[0]).ToArray()));
            tbDensity = AddRow(tFluid, "Density", Num("998"), "kg/m³");
            tbVisc = AddRow(tFluid, "Kinematic viscosity", Num("1.004e-6"), "m²/s");
            tbSigma = AddRow(tFluid, "Surface tension", Num("0.0728"), "N/m");
            tbAngle = AddRow(tFluid, "Contact angle (0 = off)", Num("0"), "°");
            cbViscModel = AddRow(tFluid, "Rheology", Combo("Newtonian", "Power-law: ν = b0·γ̇^(b1−1)", "Bingham: ν = b0 + b1·γ̇^(b2−1)", "Cross: ν = b0 + (b1−b0)/(1+b2·γ̇^b3)", "Polynomial: ν = b0 + b1·γ̇ + …"));
            tbFluidParams = AddRow(tFluid, "Parameters b0, b1, …", new TextBox(), "SI");
            cbFluid.SelectedIndexChanged += (s, e) =>
            {
                if (loadingSim)
                {
                    return;
                }
                var p = FluidPresets[cbFluid.SelectedIndex];
                if ((string)p[0] != "Custom")
                {
                    tbDensity.Text = Fmt.G((double)p[1], 6);
                    tbVisc.Text = Fmt.G((double)p[2], 6);
                    tbSigma.Text = Fmt.G((double)p[3], 6);
                }
            };
            cbViscModel.SelectedIndexChanged += (s, e) => tbFluidParams.Enabled = cbViscModel.SelectedIndex > 0;

            var tLoad = FieldTable();
            tbSGx = AddRow(tLoad, "gx", Num("0"), "m/s²");
            tbSGy = AddRow(tLoad, "gy", Num("0"), "m/s²");
            tbSGz = AddRow(tLoad, "gz", Num("-9.81"), "m/s²");
            tbSimInletVel = AddRow(tLoad, "Inlet velocity", Num("0.5"), "m/s");
            tbSimOutletP = AddRow(tLoad, "Outlet gauge pressure", Num("0"), "Pa");
            cbFill = AddRow(tLoad, "Initial state", Combo("Empty (filled by the inlet)", "Completely filled", "Filled up to a level"));
            tbFillLevel = AddRow(tLoad, "Fill level", Num("0"), "m");

            var tNum = FieldTable();
            cbCollision = AddRow(tNum, "Collision model", Combo("MRT (recommended)", "TRT", "BGK"));
            tbSmag = AddRow(tNum, "LES Smagorinsky C", Num("0.15"), "0 = laminar");
            tbLatticeU = AddRow(tNum, "Lattice velocity", Num("0.1"), "≤ 0.15");
            tbCharU = AddRow(tNum, "Character. velocity", Num("0"), "m/s, 0 = auto");
            tbTauMax = AddRow(tNum, "Max. τ (rheology)", Num("1.5"), "");
            tbEndTime = AddRow(tNum, "Simulated time", Num("1"), "s");
            tbSaveInt = AddRow(tNum, "Save interval", Num("0.05"), "s");
            tbThreads = AddRow(tNum, "Threads", Num(Environment.ProcessorCount.ToString()), "");

            lblDerived = new Label { AutoSize = true, MaximumSize = new Size(400, 0), Padding = new Padding(6), ForeColor = Color.FromArgb(40, 40, 40) };
            var derivedBox = Group("Derived quantities", lblDerived);

            foreach (var tb in new[] { tbDensity, tbVisc, tbSGx, tbSGy, tbSGz, tbSimInletVel, tbSmag, tbLatticeU, tbCharU, tbEndTime, tbSaveInt })
            {
                tb.TextChanged += (s, e) => UpdateDerived();
            }
            foreach (var cb in new[] { cbProblem, cbFill })
            {
                cb.SelectedIndexChanged += (s, e) => UpdateDerived();
            }

            var buttons = new FlowLayoutPanel { Dock = DockStyle.Top, Height = 44, Padding = new Padding(2, 6, 2, 2) };
            btnSaveSim = new Button { Text = "Save", Width = 90, Height = 30 };
            btnSaveSim.Click += (s, e) => SaveSimulation();
            btnRunSim = new Button { Text = "▶ Run", Width = 110, Height = 30, Font = new Font(Font, FontStyle.Bold) };
            btnRunSim.Click += (s, e) => RunSimulation();
            btnStopSim = new Button { Text = "■ Stop", Width = 100, Height = 30, Enabled = false };
            btnStopSim.Click += (s, e) =>
            {
                runner.RequestStop();
                AppendSimLog("Stop requested – the solver ends at the next save and writes the results.");
            };
            buttons.Controls.AddRange(new Control[] { btnSaveSim, btnRunSim, btnStopSim });

            var column = Column(440,
                Group("Simulation", tBase),
                Group("Fluid", tFluid),
                Group("Loads and initial state", tLoad),
                Group("Numerics", tNum),
                derivedBox,
                buttons);

            simChart = new Chart { Dock = DockStyle.Fill, BackColor = Color.White, MinimumSize = new Size(50, 50) };
            var area = new ChartArea("main");
            area.AxisX.Title = "time [s]";
            area.AxisY.Title = "max. velocity [m/s]";
            area.AxisY2.Title = "mass [kg]";
            area.AxisY2.Enabled = AxisEnabled.True;
            area.AxisX.MajorGrid.LineColor = Color.Gainsboro;
            area.AxisY.MajorGrid.LineColor = Color.Gainsboro;
            area.AxisY2.MajorGrid.Enabled = false;
            area.AxisX.LabelStyle.Format = "G4";
            area.AxisX.Minimum = 0;
            area.AxisY.LabelStyle.Format = "G3";
            area.AxisY2.LabelStyle.Format = "G4";
            area.AxisY2.IsStartedFromZero = false;
            simChart.ChartAreas.Add(area);
            simChart.Legends.Add(new Legend { Docking = Docking.Top });
            simChart.Series.Add(new Series("max. velocity") { ChartType = SeriesChartType.Line, BorderWidth = 2, Color = Color.FromArgb(31, 119, 180) });
            simChart.Series.Add(new Series("mass") { ChartType = SeriesChartType.Line, BorderWidth = 2, Color = Color.FromArgb(214, 39, 40), YAxisType = AxisType.Secondary });

            simProgress = new ProgressBar { Dock = DockStyle.Top, Height = 18 };
            simLog = new TextBox { Dock = DockStyle.Bottom, Height = 220, Multiline = true, ScrollBars = ScrollBars.Vertical, ReadOnly = true, Font = new Font("Consolas", 8.5f), WordWrap = false };
            var right = new Panel { Dock = DockStyle.Fill, Padding = new Padding(4), Size = new Size(900, 800) };
            right.Controls.Add(simChart);
            right.Controls.Add(simProgress);
            right.Controls.Add(simLog);

            page.Controls.Add(right);
            page.Controls.Add(column);
            page.Controls.Add(listPanel);
            return page;
        }

        void AppendSimLog(string line)
        {
            if (simLog.TextLength > 400000)
            {
                simLog.Text = simLog.Text.Substring(200000);
            }
            simLog.AppendText(line + Environment.NewLine);
        }

        void RefreshSimulationList(string select)
        {
            lbSims.Items.Clear();
            if (geom == null)
            {
                return;
            }
            string dir = geom.Dir(settings.Workspace);
            if (Directory.Exists(dir))
            {
                foreach (var d in Directory.GetDirectories(dir).Where(d => File.Exists(Path.Combine(d, "simulation.ini"))).OrderBy(d => d))
                {
                    lbSims.Items.Add(Path.GetFileName(d));
                }
            }
            if (select != null && lbSims.Items.Contains(select))
            {
                lbSims.SelectedItem = select;
            }
            else if (lbSims.Items.Count > 0)
            {
                lbSims.SelectedIndex = 0;
            }
            else
            {
                var s = DefaultSimulation("sim1");
                SimulationToUi(s);
            }
        }

        SimulationDef DefaultSimulation(string name)
        {
            var s = new SimulationDef { Name = name };
            if (geom != null)
            {
                s.Gravity = geom.Gravity;
                s.InletVelocity = geom.Inlet.Value;
                s.OutletPressure = geom.Outlet.Value;
            }
            return s;
        }

        void NewSimulation(bool copy)
        {
            if (geom == null)
            {
                return;
            }
            int n = 1;
            while (Directory.Exists(Path.Combine(geom.Dir(settings.Workspace), "sim" + n)))
            {
                n++;
            }
            string name = Prompt("New simulation", "Simulation name:", "sim" + n);
            if (string.IsNullOrEmpty(name))
            {
                return;
            }
            if (name.IndexOfAny(Path.GetInvalidFileNameChars()) >= 0 || name.Contains(" ") || name == "mesh")
            {
                MessageBox.Show(this, "Invalid name.", Text);
                return;
            }
            SimulationDef s;
            if (copy)
            {
                try
                {
                    s = SimulationFromUi();
                }
                catch (FormatException ex)
                {
                    MessageBox.Show(this, ex.Message, Text);
                    return;
                }
                s.Name = name;
            }
            else
            {
                s = DefaultSimulation(name);
            }
            s.Save(s.Dir(geom, settings.Workspace));
            RefreshSimulationList(name);
        }

        void LoadSimulationToUi(string name)
        {
            if (geom == null || name == null)
            {
                return;
            }
            var s = SimulationDef.Load(Path.Combine(geom.Dir(settings.Workspace), name));
            SimulationToUi(s);
            LoadHistory(Path.Combine(geom.Dir(settings.Workspace), name));
        }

        void SimulationToUi(SimulationDef s)
        {
            loadingSim = true;
            tbSimName.Text = s.Name;
            cbProblem.SelectedIndex = s.ProblemId == 0 ? 1 : 0;
            cbFluid.SelectedIndex = FluidPresets.Length - 1;
            for (int i = 0; i < FluidPresets.Length - 1; i++)
            {
                if (Math.Abs((double)FluidPresets[i][1] - s.Density) < 1e-9 && Math.Abs((double)FluidPresets[i][2] - s.Viscosity) < 1e-15)
                {
                    cbFluid.SelectedIndex = i;
                }
            }
            tbDensity.Text = Fmt.D(s.Density);
            tbVisc.Text = Fmt.D(s.Viscosity);
            tbSigma.Text = Fmt.D(s.SurfaceTension);
            tbAngle.Text = Fmt.D(s.ContactAngle);
            cbViscModel.SelectedIndex = Math.Min(Math.Max(s.ViscosityModel, 0), cbViscModel.Items.Count - 1);
            tbFluidParams.Text = s.FluidParameters;
            tbFluidParams.Enabled = s.ViscosityModel > 0;
            tbSGx.Text = Fmt.D(s.Gravity.X);
            tbSGy.Text = Fmt.D(s.Gravity.Y);
            tbSGz.Text = Fmt.D(s.Gravity.Z);
            tbSimInletVel.Text = Fmt.D(s.InletVelocity);
            tbSimOutletP.Text = Fmt.D(s.OutletPressure);
            cbFill.SelectedIndex = (int)s.Fill;
            tbFillLevel.Text = Fmt.D(s.FillLevel);
            cbCollision.SelectedIndex = s.CollisionModel == 1 ? 0 : (s.CollisionModel == 2 ? 1 : 2);
            tbSmag.Text = Fmt.D(s.Smagorinsky);
            tbLatticeU.Text = Fmt.D(s.LatticeVelocity);
            tbCharU.Text = Fmt.D(s.CharVelocity);
            tbTauMax.Text = Fmt.D(s.TauMax);
            tbEndTime.Text = Fmt.D(s.EndTime);
            tbSaveInt.Text = Fmt.D(s.SaveInterval);
            tbThreads.Text = s.Threads.ToString();
            loadingSim = false;
            UpdateDerived();
        }

        SimulationDef SimulationFromUi()
        {
            var s = new SimulationDef
            {
                Name = tbSimName.Text.Trim(),
                ProblemId = cbProblem.SelectedIndex == 1 ? 0 : 1,
                Density = D(tbDensity, "Density"),
                Viscosity = D(tbVisc, "Viscosity"),
                SurfaceTension = D(tbSigma, "Surface tension"),
                ContactAngle = D(tbAngle, "Contact angle"),
                ViscosityModel = cbViscModel.SelectedIndex,
                FluidParameters = tbFluidParams.Text.Trim(),
                Gravity = new Vec3(D(tbSGx, "gx"), D(tbSGy, "gy"), D(tbSGz, "gz")),
                InletVelocity = D(tbSimInletVel, "Inlet velocity"),
                OutletPressure = D(tbSimOutletP, "Outlet gauge pressure"),
                Fill = (FillMode)cbFill.SelectedIndex,
                FillLevel = D(tbFillLevel, "Fill level"),
                CollisionModel = cbCollision.SelectedIndex == 0 ? 1 : (cbCollision.SelectedIndex == 1 ? 2 : 0),
                Smagorinsky = D(tbSmag, "Smagorinsky"),
                LatticeVelocity = D(tbLatticeU, "Lattice velocity"),
                CharVelocity = D(tbCharU, "Characteristic velocity"),
                TauMax = D(tbTauMax, "Max. tau"),
                EndTime = D(tbEndTime, "Simulated time"),
                SaveInterval = D(tbSaveInt, "Save interval"),
                Threads = (int)D(tbThreads, "Threads")
            };
            if (s.Viscosity <= 0 || s.Density <= 0 || s.EndTime <= 0 || s.SaveInterval <= 0 || s.LatticeVelocity <= 0)
            {
                throw new FormatException("Density, viscosity, times and lattice velocity must be positive.");
            }
            if (s.ViscosityModel > 0 && s.FluidParameters.Length == 0)
            {
                throw new FormatException("Enter the parameters of the rheology model.");
            }
            return s;
        }

        void UpdateDerived()
        {
            if (loadingSim || geom == null)
            {
                return;
            }
            if (!geom.MeshReady)
            {
                lblDerived.Text = "The geometry has no mesh (tab 1).";
                return;
            }
            try
            {
                var s = SimulationFromUi();
                var d = Derived.Compute(geom, s);
                string text = string.Format(
                    "Time step dt = {0} s\nτ = {1} (lattice ν = {2})\nIterations: {3:N0}, saved every {4:N0}\nCharacteristic velocity {5} m/s, Ma = {6}\nCells: {7:N0}, memory ≈ {8} GB",
                    Fmt.G(d.Dt, 4), Fmt.G(d.Tau, 5), Fmt.G((d.Tau - 0.5) / 3, 3), d.MaxIter, d.SaveRate, Fmt.G(d.Uc, 3), Fmt.G(d.Mach, 3), d.Cells, Fmt.G(d.MemoryGB, 3));
                if (d.Warnings.Count > 0)
                {
                    text += "\n\n⚠ " + string.Join("\n⚠ ", d.Warnings);
                }
                lblDerived.Text = text;
                lblDerived.ForeColor = d.Warnings.Count > 0 ? Color.FromArgb(150, 60, 0) : Color.FromArgb(40, 40, 40);
            }
            catch (FormatException ex)
            {
                lblDerived.Text = ex.Message;
            }
        }

        /// <summary>Saves the simulation and writes the solver input, returns null on failure.</summary>
        SimulationDef SaveSimulation()
        {
            if (geom == null)
            {
                return null;
            }
            try
            {
                var s = SimulationFromUi();
                if (string.IsNullOrEmpty(s.Name) || s.Name.IndexOfAny(Path.GetInvalidFileNameChars()) >= 0 || s.Name.Contains(" "))
                {
                    throw new FormatException("Invalid simulation name.");
                }
                SimulationWriter.Write(geom, s, settings.Workspace);
                RefreshSimulationList(s.Name);
                SetStatus("Simulation " + s.Name + " saved.");
                return s;
            }
            catch (Exception ex)
            {
                MessageBox.Show(this, ex.Message, Text, MessageBoxButtons.OK, MessageBoxIcon.Warning);
                return null;
            }
        }

        void RunSimulation()
        {
            if (runner.Running)
            {
                MessageBox.Show(this, "A computation is already running for simulation " + runningSimName + ".", Text);
                return;
            }
            var s = SaveSimulation();
            if (s == null)
            {
                return;
            }
            var d = Derived.Compute(geom, s);
            if (d.Tau < 0.501 && s.Smagorinsky <= 0 &&
                MessageBox.Show(this, "The relaxation time τ is very close to 0.5 and the computation will probably be unstable. Run anyway?", Text, MessageBoxButtons.YesNo, MessageBoxIcon.Warning) != DialogResult.Yes)
            {
                return;
            }
            string simDir = s.Dir(geom, settings.Workspace);
            string resultDir = Path.Combine(simDir, "result");
            if (Directory.Exists(resultDir) && Directory.GetFiles(resultDir).Length > 0)
            {
                if (MessageBox.Show(this, "The simulation already has results, a new run overwrites them. Continue?", Text, MessageBoxButtons.YesNo, MessageBoxIcon.Question) != DialogResult.Yes)
                {
                    return;
                }
                foreach (var f in Directory.GetFiles(resultDir))
                {
                    File.Delete(f);
                }
            }
            simLog.Clear();
            foreach (var ser in simChart.Series)
            {
                ser.Points.Clear();
            }
            try
            {
                runningSimName = s.Name;
                runningMaxIter = d.MaxIter;
                simProgress.Value = 0;
                runner.Start(settings.SolverExe, geom.Dir(settings.Workspace), simDir);
                btnRunSim.Enabled = false;
                btnStopSim.Enabled = true;
                SetStatus("Running simulation " + geom.Name + "/" + s.Name + " …");
            }
            catch (Exception ex)
            {
                MessageBox.Show(this, ex.Message, Text, MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
        }

        void OnSolverProgress(ProgressInfo p)
        {
            if (runningMaxIter > 0)
            {
                simProgress.Value = Math.Min(100, (int)(100.0 * p.Iteration / runningMaxIter));
            }
            if ((lbSims.SelectedItem as string) != runningSimName)
            {
                return;
            }
            if (!double.IsNaN(p.UMax) && !double.IsInfinity(p.UMax))
            {
                simChart.Series[0].Points.AddXY(p.Time, p.UMax);
            }
            if (!double.IsNaN(p.Mass))
            {
                simChart.Series[1].Points.AddXY(p.Time, p.Mass);
            }
            SetStatus(string.Format("Simulation {0}: iteration {1:N0} / {2:N0}, t = {3} s, {4} MLUPS", runningSimName, p.Iteration, runningMaxIter, Fmt.G(p.Time, 4), Fmt.G(p.Mlups, 3)));
        }

        void OnSolverExited(int code)
        {
            btnRunSim.Enabled = true;
            btnStopSim.Enabled = false;
            string msg = code == 0 ? "Computation " + runningSimName + " finished." : (code == 3 ? "Computation " + runningSimName + " diverged (NaN) – lower the lattice velocity, enable LES or refine the mesh." : "Computation " + runningSimName + " failed (code " + code + ").");
            AppendSimLog(msg);
            SetStatus(msg);
            if (code == 0)
            {
                simProgress.Value = 100;
            }
        }

        void LoadHistory(string simDir)
        {
            if (runner.Running && (lbSims.SelectedItem as string) == runningSimName)
            {
                return;
            }
            foreach (var ser in simChart.Series)
            {
                ser.Points.Clear();
            }
            double dt = 0;
            try
            {
                var s = SimulationDef.Load(simDir);
                dt = Derived.Compute(geom, s).Dt;
            }
            catch (Exception)
            {
            }
            foreach (var row in SolverRunner.ReadResiduum(simDir))
            {
                double t = row[0] * dt;
                if (!double.IsNaN(row[1]))
                {
                    simChart.Series[0].Points.AddXY(t, row[1]);
                }
                if (row.Length > 2)
                {
                    simChart.Series[1].Points.AddXY(t, row[2]);
                }
            }
        }

        /// <summary>Starts the solver for the given simulation and waits for it (self test).</summary>
        public int RunSimulationSync(string name)
        {
            lbSims.SelectedItem = name;
            var s = SimulationDef.Load(Path.Combine(geom.Dir(settings.Workspace), name));
            SimulationWriter.Write(geom, s, settings.Workspace);
            int exitCode = -1;
            bool done = false;
            Action<int> handler = c => { exitCode = c; done = true; };
            runner.Exited += handler;
            runningSimName = name;
            runningMaxIter = Derived.Compute(geom, s).MaxIter;
            runner.Start(settings.SolverExe, geom.Dir(settings.Workspace), s.Dir(geom, settings.Workspace));
            while (!done)
            {
                Application.DoEvents();
                System.Threading.Thread.Sleep(50);
            }
            runner.Exited -= handler;
            Application.DoEvents();
            return exitCode;
        }
    }
}
