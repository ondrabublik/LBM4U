# cLBMnet v2 – přehled změn

Výchozí stav je první commit (kopie větve `rotace_master_new` @ 93dd6e7).

## Sestavení a testy

```
build.bat            # Visual Studio 2022 (MSVC, OpenMP /openmp:llvm)
build.bat mingw      # MSYS2 g++ + Ninja, staticky linkovaný exe
```

CMake vytvoří `cLBMnet.exe` (řešič) a `lbm_tests.exe` (26 unit testů, spouštět ve složce buildu, volitelně s filtrem názvu).

Grafické rozhraní: `gui/build_gui.bat` → `gui/bin/LBMStudio.exe`, popis v [gui/README.md](gui/README.md).

Volba `-DLBM_FP32=ON` uloží distribuční funkce v jednoduché přesnosti (poloviční paměť a datový tok).

Spuštění řešiče:

```
cLBMnet                                # jako dříve, čte args.txt
cLBMnet <geometrie> <simulace>         # simulations/<geometrie>/<simulace>/
cLBMnet --paths <adresářGeometrie> <adresářSimulace>
```

Výpočet lze předčasně ukončit vytvořením souboru `STOP` v adresáři simulace (výsledky se uloží).
Návratový kód 3 znamená divergenci (NaN).

## Opravené chyby

| Oblast | Chyba | Oprava |
|---|---|---|
| Mesh krok 8 | `makro`, `collision` a `newFsGeometry` v jedné paralelní smyčce, sousední bloky četly rozpracované `Rho/U/flagMap` | tři oddělené průchody (makro → kolize + příprava nových buněk → aplikace změn) |
| Redistribuce hmoty | souběžné `mass +=` do sousedních bloků, `massErr` se nikdy nevracel | `omp atomic`, `massErr` se rozdělí mezi buňky rozhraní bloku |
| Krok 43 | souběh na `nParticles`, `maxParticlesInCell` | sériová redukce po bloku |
| FVM částice | tok přes hranici bloku zapisoval do špatného bloku, `oppLW` měl špatné opačné směry, koncentrace se měnila během čtení sousedy | konzervativní tok `(F(A→B) − F(B→A))/2`, oddělené `computeFlux`/`applyFlux` |
| `BlockBasic` | `isBlockEmpty` neinicializováno; gravitace 3× slabší a člen `−rho0` | inicializace; Guo forcing |
| `findArtefactsFilled` | `boundary[ind]` místo `boundary[indShift]` | opraveno |
| Artefakty | současná konverze sousedních buněk rozhraní na kapalinu a plyn vytvořila kontakt kapalina–plyn a ztrátu hmoty | buňka sousedící s „filled“ artefaktem se nevyprázdní, pojistka `fixFluidGasContacts` |
| `initHydrostaticPress` | zastíněná proměnná `lo`, reference od vrcholu domény | počítá se od nejvyšší buňky kapaliny |
| Restart | `writeResultsToOutFile` bral sekvenční index bloku | správné mapování a id souborů |
| `coarseLevel` | nekonzistentní velikosti polí | zahušťování zatím není implementováno – `coarseLevel != 1` hlásí chybu |
| Ne-newtonské kapaliny | fungovaly jen s MRT; smyková rychlost z rozdílů sousedů (souběh) | libovolný kolizní model, lokální smyková rychlost z nerovnovážných momentů |
| Okrajové soubory | `line.compare("1") == 1` | čtení celočíselného typu buňky, kontrola chyb |
| Čtení souborů | `catch (int)` nezachytil výjimky | `std::exception`, srozumitelné hlášení |
| Marching cubes | `realloc` + `delete[]`, podmínka „povrch“ vždy pravdivá | `std::vector`, test min/max izohodnoty |
| `saveRate` | výchozí `1e10` v `int` | bezpečné čtení parametrů (i zápis `1e5`) |
| `particleBarrierForce` | `omp parallel` bez `for` | opraveno |
| Euler-LBM částice | nová rychlost částic používala už přepsané složky | použijí se staré hodnoty |
| Další | úniky paměti, `delete` vs `delete[]`, `getParDoubleArray`, pád D2Q9 (`nk = 1`), chybný výpis `tau`, přetečení názvů výstupů, `windows.h` | opraveno |

## Numerika

* Kolizní operátory BGK, TRT a MRT ve společné třídě `CollisionOperator` (parametr `collisionModel`, starý `MRT` funguje dál).
  * TRT s magickým parametrem `trtMagic` (výchozí 3/16 – stěna bounce-back leží přesně v polovině, Poiseuille je přesný na 1e-13).
  * MRT: liché momenty relaxují s `s_q = 8(2 − s_ν)/(8 − s_ν)` (`mrtMagicSq = 1`), poloha stěny nezávisí na viskozitě.
  * MRT matice se invertují numericky, D2Q9 má vlastní MRT matici.
* `BlockBasic`: Guo forcing (rychlost obsahuje F/2).
* Volná hladina: výstupní rychlost `u + F/2` (EDM), síly částic se dělí hustotou.
* Výměna hmoty mezi buňkami rozhraní s vahami podle typu buněk (Thürey 2007).
* Redistribuce přebytku hmoty vážená normálou.
* Křivost výškovou funkcí (`curvatureModel = 1`); výchozí zůstává divergence vyhlazených normál, která v dynamice dává menší parazitní proudy.
* Okrajové podmínky: vtok (`uInlet, vInlet, wInlet`, typ buňky 2), výtok s tlakem (`pOutlet`, typ 3, anti-bounce-back), pohyblivá stěna (`uwall, vwall, wwall`, typ 4, víko domény `movingLidTop`).
* Volitelné meze `tauMin`, `tauMax` pro proměnnou relaxaci.
* Smagorinského LES model (`smagorinsky` = C, typicky 0,1–0,17) pro jednofázový i free-surface výpočet, explicitní vztah Hou et al. (1996).
* `tauMinusMax`: horní mez relaxačního času lichých momentů (TRT τ⁻, MRT s_q). Při τ → 0,5 (voda) dává magický vztah τ⁻ ≈ 10³ a výpočet diverguje; s `tauMinusMax = 1` a LES je plnění vodou stabilní. GUI tuto hodnotu zapisuje vždy.
* Řešič upozorní, když rychlost mřížky překročí 0,3 (vysoké Machovo číslo).

## Výkon

* Odstraněno nepoužívané pole `fcOld` (1/3 paměti populací).
* Tabulky sousedů jen pro okrajové buňky bloku (místo 2·nf intů na každou buňku).
* `BlockBasic`: spojený výpočet momentů a kolize.
* Volitelné FP32 populace + posun DDF (`ddfShift = 1`).
* Vzdálenostní funkce Dijkstrovým algoritmem (stejný výsledek, O(N log N)).
* Binární VTK (`vtkBinary = 1`), kolekce `output.pvd`/`outputFS.pvd`, kompaktní `fields_<iter>.bin` pro GUI (`writeFields = 1`).
* MLUPS z `omp_get_wtime` průměrované přes interval ukládání.

## Neprovedeno

* Ghost vrstvy bloků, AA-pattern (jediné pole populací) a SoA uspořádání – větší přestavba, která má smysl spolu se zahušťováním sítě.
* Zahušťování sítě (oct-tree).

## Formát `fields_<iter>.bin`

```
char[4] "LBMF", int32 verze (1), int32 iterace, double čas, int32 problemID,
int32 nBlock, int32 ni, nj, nk, double Ch, int32 nFields, nFields × char[16] název
pro každý blok: int32 bi, bj, bk; pro každé pole float32[ni*nj*nk] (index i + ni*j + ni*nj*k)
pole: ux, uy, uz [m/s], p [Pa], type (typ buňky), u volné hladiny navíc vof, flag
```
