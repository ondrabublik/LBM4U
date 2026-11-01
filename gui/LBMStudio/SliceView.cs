using System;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using System.Windows.Forms;

namespace LBMStudio
{
    /// <summary>2D colour map of a slice through the result field.</summary>
    public class SliceView : Control
    {
        float[,] values;          // [a, b], NaN = no fluid
        byte[,] kind;             // 0 fluid, 1 wall, 2 gas
        float[,] vu, vv;          // optional in-plane vectors
        double extentA, extentB;  // physical size
        double minV, maxV;
        string title = "", unit = "", axisA = "", axisB = "";
        Bitmap image;
        Rectangle plotRect;
        readonly ToolTip tip = new ToolTip();
        string lastTip = "";

        public bool FixedRange;
        public double RangeMin, RangeMax;

        public SliceView()
        {
            DoubleBuffered = true;
            ResizeRedraw = true;
            BackColor = Color.White;
            MouseMove += OnHover;
        }

        public void SetData(float[,] values, byte[,] kind, double extentA, double extentB, string axisA, string axisB, string title, string unit, float[,] vu = null, float[,] vv = null)
        {
            this.values = values;
            this.kind = kind;
            this.extentA = extentA;
            this.extentB = extentB;
            this.axisA = axisA;
            this.axisB = axisB;
            this.title = title;
            this.unit = unit;
            this.vu = vu;
            this.vv = vv;
            minV = double.MaxValue;
            maxV = double.MinValue;
            foreach (float v in values)
            {
                if (!float.IsNaN(v))
                {
                    minV = Math.Min(minV, v);
                    maxV = Math.Max(maxV, v);
                }
            }
            if (FixedRange)
            {
                minV = RangeMin;
                maxV = RangeMax;
            }
            if (minV > maxV)
            {
                minV = 0;
                maxV = 1;
            }
            if (maxV - minV < 1e-12)
            {
                maxV = minV + 1e-12;
            }
            BuildImage();
            Invalidate();
        }

        public double DataMin { get { return minV; } }
        public double DataMax { get { return maxV; } }

        /// <summary>Viridis-like colour map.</summary>
        public static Color Map(double t)
        {
            t = Math.Max(0, Math.Min(1, t));
            double[,] c = {
                { 68, 1, 84 }, { 72, 40, 120 }, { 62, 74, 137 }, { 49, 104, 142 }, { 38, 130, 142 },
                { 31, 158, 137 }, { 53, 183, 121 }, { 109, 205, 89 }, { 180, 222, 44 }, { 253, 231, 37 } };
            double x = t * 9;
            int i = Math.Min(8, (int)x);
            double f = x - i;
            return Color.FromArgb((int)(c[i, 0] + f * (c[i + 1, 0] - c[i, 0])), (int)(c[i, 1] + f * (c[i + 1, 1] - c[i, 1])), (int)(c[i, 2] + f * (c[i + 1, 2] - c[i, 2])));
        }

        void BuildImage()
        {
            image?.Dispose();
            int na = values.GetLength(0), nb = values.GetLength(1);
            image = new Bitmap(na, nb, PixelFormat.Format32bppArgb);
            var data = image.LockBits(new Rectangle(0, 0, na, nb), ImageLockMode.WriteOnly, PixelFormat.Format32bppArgb);
            var px = new int[na * nb];
            for (int b = 0; b < nb; b++)
            {
                for (int a = 0; a < na; a++)
                {
                    Color col;
                    float v = values[a, b];
                    byte k = kind != null ? kind[a, b] : (byte)0;
                    if (k == 1)
                    {
                        col = Color.FromArgb(150, 150, 155);
                    }
                    else if (k == 2 || float.IsNaN(v))
                    {
                        col = k == 2 ? Color.FromArgb(235, 240, 248) : Color.FromArgb(150, 150, 155);
                    }
                    else
                    {
                        col = Map((v - minV) / (maxV - minV));
                    }
                    px[(nb - 1 - b) * na + a] = col.ToArgb(); // b axis upwards
                }
            }
            System.Runtime.InteropServices.Marshal.Copy(px, 0, data.Scan0, px.Length);
            image.UnlockBits(data);
        }

        protected override void OnPaint(PaintEventArgs e)
        {
            var g = e.Graphics;
            g.Clear(BackColor);
            if (image == null)
            {
                TextRenderer.DrawText(g, "Žádná data", Font, ClientRectangle, Color.Gray, TextFormatFlags.HorizontalCenter | TextFormatFlags.VerticalCenter);
                return;
            }
            int legendW = 80, margin = 40;
            var avail = new Rectangle(margin, 28, Math.Max(10, Width - legendW - margin - 10), Math.Max(10, Height - 28 - margin));
            double scale = Math.Min(avail.Width / extentA, avail.Height / extentB);
            int w = (int)(extentA * scale), h = (int)(extentB * scale);
            plotRect = new Rectangle(avail.X + (avail.Width - w) / 2, avail.Y + (avail.Height - h) / 2, Math.Max(w, 1), Math.Max(h, 1));
            g.InterpolationMode = InterpolationMode.NearestNeighbor;
            g.PixelOffsetMode = PixelOffsetMode.Half;
            g.DrawImage(image, plotRect);
            g.PixelOffsetMode = PixelOffsetMode.Default;
            g.SmoothingMode = SmoothingMode.AntiAlias;
            g.DrawRectangle(Pens.DimGray, plotRect);

            if (vu != null && vv != null)
            {
                DrawVectors(g);
            }

            using (var bold = new Font(Font, FontStyle.Bold))
            {
                g.DrawString(title, bold, Brushes.Black, margin, 6);
            }
            g.DrawString(axisA + " [m]: 0 … " + Fmt.G(extentA, 3), Font, Brushes.DimGray, plotRect.Left, plotRect.Bottom + 4);
            var sf = new StringFormat { FormatFlags = StringFormatFlags.DirectionVertical };
            g.DrawString(axisB + " [m]: 0 … " + Fmt.G(extentB, 3), Font, Brushes.DimGray, plotRect.Left - 20, plotRect.Top, sf);

            // colour bar
            var bar = new Rectangle(Width - legendW + 6, plotRect.Top, 16, Math.Max(plotRect.Height, 40));
            for (int y = 0; y < bar.Height; y++)
            {
                using (var pen = new Pen(Map(1 - (double)y / bar.Height)))
                {
                    g.DrawLine(pen, bar.Left, bar.Top + y, bar.Right, bar.Top + y);
                }
            }
            g.DrawRectangle(Pens.DimGray, bar);
            for (int t = 0; t <= 4; t++)
            {
                double val = maxV - (maxV - minV) * t / 4;
                int y = bar.Top + bar.Height * t / 4;
                g.DrawString(Fmt.G(val, 3), Font, Brushes.Black, bar.Right + 3, y - 7);
            }
            g.DrawString(unit, Font, Brushes.DimGray, bar.Left - 4, bar.Top - 18);
        }

        void DrawVectors(Graphics g)
        {
            int na = values.GetLength(0), nb = values.GetLength(1);
            double cell = (double)plotRect.Width / na;
            int step = Math.Max(1, (int)Math.Ceiling(14 / cell));
            double vmax = 0;
            for (int a = 0; a < na; a++)
                for (int b = 0; b < nb; b++)
                {
                    if (!float.IsNaN(vu[a, b]) && (kind == null || kind[a, b] == 0))
                    {
                        vmax = Math.Max(vmax, Math.Sqrt(vu[a, b] * vu[a, b] + vv[a, b] * vv[a, b]));
                    }
                }
            if (vmax <= 0)
            {
                return;
            }
            double len = step * cell * 0.9 / vmax;
            using (var pen = new Pen(Color.FromArgb(200, 20, 20, 20), 1f))
            {
                pen.CustomEndCap = new AdjustableArrowCap(3, 3);
                for (int a = step / 2; a < na; a += step)
                {
                    for (int b = step / 2; b < nb; b += step)
                    {
                        if (float.IsNaN(vu[a, b]) || (kind != null && kind[a, b] != 0))
                        {
                            continue;
                        }
                        float x = (float)(plotRect.Left + (a + 0.5) * cell);
                        float y = (float)(plotRect.Bottom - (b + 0.5) * cell);
                        float dx = (float)(vu[a, b] * len), dy = (float)(-vv[a, b] * len);
                        if (Math.Abs(dx) + Math.Abs(dy) > 1)
                        {
                            g.DrawLine(pen, x, y, x + dx, y + dy);
                        }
                    }
                }
            }
        }

        void OnHover(object sender, MouseEventArgs e)
        {
            if (values == null || !plotRect.Contains(e.Location))
            {
                return;
            }
            int na = values.GetLength(0), nb = values.GetLength(1);
            int a = (int)((e.X - plotRect.Left) * na / (double)plotRect.Width);
            int b = (int)((plotRect.Bottom - e.Y) * nb / (double)plotRect.Height);
            if (a < 0 || b < 0 || a >= na || b >= nb)
            {
                return;
            }
            string text = float.IsNaN(values[a, b]) ? "stěna" : Fmt.G(values[a, b], 5) + " " + unit;
            text = axisA + " = " + Fmt.G((a + 0.5) * extentA / na, 4) + " m, " + axisB + " = " + Fmt.G((b + 0.5) * extentB / nb, 4) + " m: " + text;
            if (text != lastTip)
            {
                lastTip = text;
                tip.SetToolTip(this, text);
            }
        }

        public void SaveImage(string path)
        {
            using (var bmp = new Bitmap(Math.Max(Width, 10), Math.Max(Height, 10)))
            {
                DrawToBitmap(bmp, new Rectangle(0, 0, bmp.Width, bmp.Height));
                bmp.Save(path, ImageFormat.Png);
            }
        }
    }
}
