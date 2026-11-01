using System;
using System.Drawing;
using System.IO;
using System.Windows.Forms;

namespace LBMStudio
{
    static class Program
    {
        [STAThread]
        static int Main(string[] args)
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            if (args.Length >= 3 && args[0] == "--selftest")
            {
                return SelfTest.Run(args[1], args[2], args.Length > 3 ? args[3] : null);
            }
            Application.ThreadException += (s, e) => MessageBox.Show(e.Exception.Message, "LBM Studio", MessageBoxButtons.OK, MessageBoxIcon.Error);
            Application.Run(new MainForm());
            return 0;
        }
    }

    /// <summary>
    /// Headless check of the whole chain: STL/STEP import, mesh, simulation, solver run, results and screenshots.
    /// LBMStudio.exe --selftest &lt;cad file&gt; &lt;output dir&gt; [solver exe]
    /// </summary>
    static class SelfTest
    {
        public static int Run(string cadFile, string outDir, string solver)
        {
            Directory.CreateDirectory(outDir);
            var log = new StreamWriter(Path.Combine(outDir, "selftest.log")) { AutoFlush = true };
            try
            {
                // the test settings are not saved, settings.ini of the user stays untouched
                var settings = AppSettings.Load();
                settings.Workspace = Path.Combine(outDir, "workspace");
                if (solver != null)
                {
                    settings.SolverExe = solver;
                }
                var form = new MainForm(settings) { StartPosition = FormStartPosition.Manual, Location = new Point(0, 0), Width = 1500, Height = 950 };
                form.Show();
                Application.DoEvents();
                string geomName = "selftest";
                string geomDir = Path.Combine(settings.Workspace, geomName);
                if (Directory.Exists(geomDir))
                {
                    Directory.Delete(geomDir, true);
                }
                new GeometryDef { Name = geomName }.Save(settings.Workspace);
                form.SelectGeometryByName(geomName);
                form.ImportGeometryFile(cadFile);
                Application.DoEvents();
                form.GenerateMeshSync();
                log.WriteLine("mesh generated");
                form.Screenshot(Path.Combine(outDir, "tab_geometry.png"), 0);
                form.CreateTestSimulation("sim1");
                form.Screenshot(Path.Combine(outDir, "tab_simulation_before.png"), 1);
                int code = form.RunSimulationSync("sim1");
                log.WriteLine("solver exit code " + code);
                form.Screenshot(Path.Combine(outDir, "tab_simulation.png"), 1);
                form.ShowResultsAndSave(Path.Combine(outDir, "slice.png"), Path.Combine(outDir, "view3d.png"));
                form.Screenshot(Path.Combine(outDir, "tab_results.png"), 2);
                form.Close();
                log.WriteLine("OK");
                return code;
            }
            catch (Exception ex)
            {
                log.WriteLine("ERROR " + ex);
                return 1;
            }
            finally
            {
                log.Close();
            }
        }
    }

    public partial class MainForm
    {
        public void SelectGeometryByName(string name)
        {
            RefreshGeometries(name);
            lbGeometries.SelectedItem = name;
            Application.DoEvents();
        }

        /// <summary>Free surface filling with water through the inlet, short simulated time.</summary>
        public void CreateTestSimulation(string name)
        {
            var s = DefaultSimulation(name);
            Vec3 size = mesh.Max - mesh.Min;
            s.EndTime = 0.15;
            s.SaveInterval = 0.03;
            s.InletVelocity = 0.5;
            s.Threads = Math.Max(1, Environment.ProcessorCount / 2);
            s.Save(s.Dir(geom, settings.Workspace));
            RefreshSimulationList(name);
            tabs.SelectedIndex = 1;
            Application.DoEvents();
        }

        public void Screenshot(string file, int tab)
        {
            tabs.SelectedIndex = tab;
            Application.DoEvents();
            using (var bmp = new Bitmap(Width, Height))
            {
                DrawToBitmap(bmp, new Rectangle(0, 0, Width, Height));
                bmp.Save(file, System.Drawing.Imaging.ImageFormat.Png);
            }
        }
    }
}
