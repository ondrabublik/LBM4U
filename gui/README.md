# LBM Studio – grafické rozhraní řešiče cLBMnet

WinForms aplikace (C#, .NET Framework 4.8, 3D náhled přes WPF) pro přípravu geometrie, výpočetní sítě, spouštění simulací a prohlížení výsledků.

## Sestavení

```
gui\build_gui.bat          # kompilátor C# z Visual Studia 2022, výstup gui\bin\LBMStudio.exe
```

Alternativně s nainstalovaným .NET SDK: `dotnet build gui/LBMStudio/LBMStudio.csproj` (net8.0-windows).
Řešič `cLBMnet.exe` se hledá vedle aplikace, v kořeni repozitáře a v `build-msvc\Release` / `build-mingw`; cestu lze nastavit v **Nastavení**.

## Postup

1. **Geometrie a síť**
   * *Nová geometrie…* vytvoří složku `<pracovní složka>/<název>`.
   * *Načíst STEP / STL…* – STL se načte přímo, STEP se převede na STL přes **gmsh** (`gmsh.exe` nebo Python modul `gmsh`) nebo **FreeCAD** (`FreeCADCmd.exe`). Bez nich je potřeba exportovat STL z CAD.
   * Jednotky souboru (STEP z převodníku je v mm).
   * Oblast tekutiny: uvnitř tělesa (dutina formy) nebo kolem tělesa (obtékání, kvádr s okrajem).
   * Velikost buňky `dx` a bloku (16³ je dobrý kompromis).
   * *Zahušťování podle geometrie (octree)*: `dx` je velikost buňky úrovně 0, nejvyšší úroveň `L` dává nejmenší buňku `dx / 2^L`.
     * *Buněk přes nejužší mezeru* – z povrchu STL vede paprsek proti normále (do tekutiny) k protější stěně. Bloky podél mezery šířky `h` dostanou úroveň s `dx_L ≤ h / N`. Mezery širší než `N · dx` se nezahušťují.
     * *Buněk na poloměr zakřivení* (0 = vypnuto) – poloměr se odhadne z úhlu sousedních trojúhelníků. Ostré hrany (úhel nad 60°) se neberou jako zakřivení.
     * *Úroveň u všech stěn* – minimální úroveň bloků s povrchem tělesa.
     * Síť se vyváží 2:1 a každý blok se voxelizuje na své úrovni (tenká mezera, kterou úroveň 0 nerozliší, se tak neztratí). Bloky bez tekutiny se vynechají. Náhled *Bloky* je obarvený podle úrovně (0 šedá, 1 modrá, 2 zelená, 3 oranžová, 4 červená), *Buňky sítě* ukazují každý blok v jeho rozlišení.
   * Vtok a výtok: stěna ohraničujícího kvádru, střed a poloměr (0 = celá stěna). Označí se buňky před povrchem tekutiny.
   * Gravitace (výchozí hodnota pro simulace).
   * *Vytvořit výpočetní síť* zapíše `mesh/blocks.txt`, `mesh/boundary/*.txt` a `mesh/celltypes.bin` (u zahuštěné sítě navíc `mesh/octree.bin`, `celltypes.bin` je pak v rozlišení úrovně 0 jen pro náhled).
2. **Simulace** – libovolný počet simulací nad jednou sítí (každá ve své podsložce).
   * Úloha: volná hladina (plnění) nebo jednofázové proudění.
   * Tekutina z předvoleb nebo vlastní, reologické modely (power-law, Bingham, Cross, polynom).
   * Rychlost vtoku, přetlak výtoku, gravitace, počáteční stav (prázdná / plná / hladina do výšky).
   * Numerika: MRT + LES (C = 0,15) jsou výchozí a stabilní i pro vodu; panel *Odvozené veličiny* ukazuje dt, τ, počet iterací, Machovo číslo, paměť a varování.
   * *Spustit* / *Zastavit* (řešič uloží výsledky a skončí), průběh max. rychlosti a hmotnosti v grafu.
3. **Výsledky** – u zahuštěné sítě se řezy skládají v nejjemnější úrovni, při které mřížka nepřesáhne 64 milionů buněk. Řez libovolnou rovinou (|u|, složky rychlosti, tlak, VOF, stav buněk), volitelně vektory, časová osa s přehráváním, 3D zobrazení hladiny, export obrázku a otevření v ParaView (`result/output.pvd`, `outputFS.pvd`).

## Automatický test

```
gui\bin\LBMStudio.exe --selftest <soubor.stl|step> <výstupní složka> [cLBMnet.exe] [--refine]
```

Projde celý řetězec (import, síť, simulace plnění vodou, výsledky) a uloží snímky obrazovky všech záložek. S `--refine` vytvoří zahuštěnou síť: 2 úrovně, 10 buněk přes mezeru, vtok nad první pětinou tělesa v ose x.
