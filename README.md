# LBM4U

```
L       BBBBB   MM    MM     44   U     U
L       B    B  M M  M M    4 4   U     U
L       BBBBB   M  MM  M   4  4   U     U
L       B    B  M      M  444444  U     U
LLLLLL  BBBBB   M      M      4    UUUUU
```

**LBM4U** je 3D řešič proudění tekutin metodou lattice Boltzmann (LBM) v C++ s OpenMP, zaměřený na **plnění forem s volnou hladinou**. Součástí je grafické rozhraní **LBM Studio**, které provede celým postupem od CAD modelu (STEP / STL) přes výpočetní síť a nastavení simulace až po prohlížení výsledků.

![LBM Studio – výsledky plnění formy](docs/images/gui_results.png)

## Vlastnosti

**Řešič (`cLBMnet`)**

* Mřížky D3Q19, D3Q15 a D2Q9, bloková struktura sítě paralelizovaná přes OpenMP.
* Kolizní operátory **BGK, TRT a MRT**, Smagorinského **LES** model, Guo forcing.
* **Volná hladina** (VOF/mass-tracking podle Thüreyho): výměna hmoty mezi buňkami rozhraní, redistribuce přebytku hmoty, povrchové napětí a kontaktní úhel.
* **Ne-newtonské kapaliny**: power-law, Bingham, Cross, polynom; smyková rychlost z nerovnovážných momentů.
* Okrajové podmínky: bounce-back stěna, rychlostní vtok, tlakový výtok (anti-bounce-back), pohyblivá stěna.
* **Statické zahušťování sítě (octree)** podle geometrie: vyvážení 2:1, sub-cycling jemných úrovní, objemová metoda rozhraní úrovní (Rohde et al.), zachování hmoty na strojovou přesnost – i pro volnou hladinu.
* Částice a vlákna: Lagrangeovské částice, vlákna a Eulerovský popis koncentrace (FVM / LBM).
* Výstupy ve formátu **VTK XML** (`.vtu`, `.vtp`, `.pvd`) s kompresí zlib – přímo do ParaView, plus kompaktní `fields_<iter>.bin` pro GUI.
* Volitelné populace v jednoduché přesnosti (`-DLBM_FP32=ON`) – poloviční paměť.

**Grafické rozhraní (LBM Studio)**

* Import STL a STEP (převod přes gmsh nebo FreeCAD), 3D náhled modelu.
* Generování sítě včetně octree zahušťování podle nejužších mezer a zakřivení.
* Vtok a výtok zadané graficky na stěně ohraničujícího kvádru.
* Předvolby tekutin, reologické modely, odvozené veličiny (dt, τ, Machovo číslo, paměť) s varováním.
* Spuštění a zastavení výpočtu, průběh max. rychlosti a hmotnosti v grafu.
* Řezy libovolnou rovinou (|u|, složky rychlosti, tlak, VOF, stav buňky), vektory, časová osa s přehráváním, 3D hladina, otevření v ParaView.

## LBM Studio – ukázka

Příklad `test5`: plnění dutiny formy motorovým olejem, 3 úrovně octree sítě (252 bloků, 334 333 buněk tekutiny).

**1. Geometrie a síť** – CAD model, vtok (červeně), gravitace (modrá šipka) a bloky sítě obarvené podle úrovně zahuštění.

![Geometrie a síť](docs/images/gui_geometry.png)

**2. Simulace** – tekutina, zatížení, numerika a průběh výpočtu.

![Simulace](docs/images/gui_simulation.png)

**3. Výsledky** – řez rychlostním polem a 3D tvar hladiny v čase t ≈ 0,95 s.

![Výsledky](docs/images/gui_results.png)

## Sestavení

Požadavky: CMake ≥ 3.16, kompilátor C++17 s OpenMP (Visual Studio 2022 nebo MSYS2 MinGW g++). GUI potřebuje .NET Framework 4.8 (kompilátor C# z Visual Studia) nebo .NET SDK.

```bat
build.bat            :: Visual Studio 2022 (MSVC, OpenMP /openmp:llvm)
build.bat mingw      :: MSYS2 g++ + Ninja, staticky linkovaný exe
gui\build_gui.bat    :: LBM Studio -> gui\bin\LBMStudio.exe
```

Na Linuxu nebo ručně:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release   # volitelně -DLBM_FP32=ON
cmake --build build
```

Vznikne řešič `cLBMnet` a unit testy `lbm_tests`.

## Spuštění řešiče

Nejjednodušší je pracovat přes LBM Studio, které připraví síť i parametry a řešič spouští samo. Z příkazové řádky:

```bat
cLBMnet                                :: čte "geometrie simulace" z args.txt
cLBMnet <geometrie> <simulace>         :: simulations/<geometrie>/<simulace>/
cLBMnet --paths <adresářGeometrie> <adresářSimulace>
```

Struktura úlohy:

```
simulations/
└── <geometrie>/
    ├── geometry.ini, geometry.stl     # zadání geometrie (GUI)
    ├── mesh/                          # blocks.txt, boundary/*.txt, celltypes.bin, octree.bin
    └── <simulace>/
        ├── parameters.txt             # parametry výpočtu
        ├── init/                      # počáteční stav
        └── result/                    # output.pvd, outputFS.pvd, *.vtu, *.vtp, fields_<iter>.bin
```

* Výpočet lze předčasně ukončit vytvořením souboru `STOP` v adresáři simulace – výsledky se uloží.
* Návratový kód 3 znamená divergenci (NaN).
* Výsledky se otevřou v ParaView přes `result/output.pvd` (objem) a `result/outputFS.pvd` (hladina).

## Testy

```bat
cd build-msvc\Release
lbm_tests.exe              :: všechny testy
lbm_tests.exe refinement   :: jen testy, jejichž název obsahuje "refinement"
```

Testy pokrývají kolizní operátory, Poiseuilleovo proudění (přesnost TRT/MRT), volnou hladinu (zachování hmoty, protržení hráze, klidová hladina) a zahušťování sítě.

GUI má automatický test celého řetězce (import, síť, simulace, výsledky, snímky obrazovky):

```bat
gui\bin\LBMStudio.exe --selftest <soubor.stl|step> <výstupní složka> [cLBMnet.exe] [--refine]
```

## Struktura repozitáře

| Cesta | Obsah |
|---|---|
| `main.cpp` | hlavní smyčka řešiče |
| `Mesh.*` | bloková síť, octree, výměna mezi bloky, výstupy |
| `Block*.*` | bloky: jednofázové proudění, volná hladina, částice a vlákna |
| `Collision.*` | kolizní operátory BGK / TRT / MRT, LES |
| `Lagrange*.*`, `MarchingCubes.*` | částice, vlákna, rekonstrukce hladiny |
| `Deflate.*` | vlastní kompresor pro VTK a `fields_*.bin` (bez závislosti na zlib) |
| `tests/` | unit testy |
| `gui/` | LBM Studio (C#, WinForms + WPF 3D) – [popis](gui/README.md) |
| `matlab/` | starší MATLAB skripty pro přípravu úloh |
| `ZMENY.md` | [přehled změn](ZMENY.md), numerické detaily a formát `fields_<iter>.bin` |

## Autor

Ondřej Bublík
