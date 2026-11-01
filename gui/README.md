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
   * Vtok a výtok: stěna ohraničujícího kvádru, střed a poloměr (0 = celá stěna). Označí se buňky před povrchem tekutiny.
   * Gravitace (výchozí hodnota pro simulace).
   * *Vytvořit výpočetní síť* zapíše `mesh/blocks.txt`, `mesh/boundary/*.txt` a `mesh/celltypes.bin`.
2. **Simulace** – libovolný počet simulací nad jednou sítí (každá ve své podsložce).
   * Úloha: volná hladina (plnění) nebo jednofázové proudění.
   * Tekutina z předvoleb nebo vlastní, reologické modely (power-law, Bingham, Cross, polynom).
   * Rychlost vtoku, přetlak výtoku, gravitace, počáteční stav (prázdná / plná / hladina do výšky).
   * Numerika: MRT + LES (C = 0,15) jsou výchozí a stabilní i pro vodu; panel *Odvozené veličiny* ukazuje dt, τ, počet iterací, Machovo číslo, paměť a varování.
   * *Spustit* / *Zastavit* (řešič uloží výsledky a skončí), průběh max. rychlosti a hmotnosti v grafu.
3. **Výsledky** – řez libovolnou rovinou (|u|, složky rychlosti, tlak, VOF, stav buněk), volitelně vektory, časová osa s přehráváním, 3D zobrazení hladiny, export obrázku a otevření v ParaView (`result/output.pvd`, `outputFS.pvd`).

## Automatický test

```
gui\bin\LBMStudio.exe --selftest <soubor.stl|step> <výstupní složka> [cLBMnet.exe]
```

Projde celý řetězec (import, síť, simulace plnění vodou, výsledky) a uloží snímky obrazovky všech záložek.
