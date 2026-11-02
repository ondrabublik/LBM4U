# LBM Studio – graphical interface of the cLBMnet solver

WinForms application (C#, .NET Framework 4.8, 3D view via WPF) for preparing the geometry and the computational mesh, running simulations and viewing the results.

## Building

```
gui\build_gui.bat          # C# compiler of Visual Studio 2022, output gui\bin\LBMStudio.exe
```

Alternatively with the .NET SDK installed: `dotnet build gui/LBMStudio/LBMStudio.csproj` (net8.0-windows).
The solver `cLBMnet.exe` is searched next to the application, in the repository root and in `build-msvc\Release` / `build-mingw`; the path can be set in **Settings**.

## Workflow

1. **Geometry and mesh**
   * *New geometry…* creates the folder `<workspace>/<name>`.
   * *Load STEP / STL…* – STL is loaded directly, STEP is converted to STL with **gmsh** (`gmsh.exe` or the Python module `gmsh`) or **FreeCAD** (`FreeCADCmd.exe`). Without them, export STL from your CAD.
   * File units (STEP from the converter is in mm).
   * Fluid region: inside the body (mold cavity) or around the body (external flow, box with a margin).
   * Cell size `dx` and block size (16³ is a good compromise).
   * *Refinement by geometry (octree)*: `dx` is the cell size of level 0, the max. level `L` gives the smallest cell `dx / 2^L`.
     * *Cells across a gap* – a ray is cast from the STL surface against the normal (into the fluid) to the opposite wall. Blocks along a gap of width `h` get the level with `dx_L ≤ h / N`. Gaps wider than `N · dx` are not refined.
     * *Cells per curvature radius* (0 = off) – the radius is estimated from the angle of neighboring triangles. Sharp edges (angle above 60°) are not taken as curvature.
     * *Level at all walls* – minimum level of the blocks containing the body surface.
     * The mesh is balanced 2:1 and each block is voxelized on its own level (a thin gap that level 0 cannot resolve is not lost). Blocks without fluid are skipped. The *Blocks* view is colored by level (0 gray, 1 blue, 2 green, 3 orange, 4 red), *Mesh cells* show each block in its own resolution.
   * Inlet and outlet: face of the bounding box, center and radius (0 = whole face). The cells in front of the fluid surface are marked.
   * Gravity (default value for the simulations).
   * *Generate mesh* writes `mesh/blocks.txt`, `mesh/boundary/*.txt` and `mesh/celltypes.bin` (a refined mesh also writes `mesh/octree.bin`, `celltypes.bin` is then in the level 0 resolution for the preview only).
2. **Simulation** – any number of simulations on one mesh (each in its own subfolder).
   * Problem: free surface (filling) or single phase flow.
   * Fluid from presets or custom, rheology models (power-law, Bingham, Cross, polynomial).
   * Inlet velocity, outlet gauge pressure, gravity, initial state (empty / completely filled / filled up to a level).
   * Numerics: MRT + LES (C = 0.15) are the defaults and stable even for water; the *Derived quantities* panel shows dt, τ, the number of iterations, the Mach number, memory and warnings.
   * *Run* / *Stop* (the solver saves the results and ends), maximum velocity and mass in a chart.
3. **Results** – for a refined mesh the slices are assembled on the finest level at which the grid does not exceed 64 million cells. Slice along any axis (|u|, velocity components, pressure, VOF, cell state), optional vectors, time line with playback, 3D free surface, image export and opening in ParaView (`result/output.pvd`, `outputFS.pvd`).

## Automatic test

```
gui\bin\LBMStudio.exe --selftest <file.stl|step> <output folder> [cLBMnet.exe] [--refine]
```

Runs the whole chain (import, mesh, water filling simulation, results) and saves screenshots of all tabs. With `--refine` it creates a refined mesh: 2 levels, 10 cells across a gap, inlet above the first fifth of the body in x.
