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
| `coarseLevel` | nekonzistentní velikosti polí | nahrazeno statickým zahušťováním (viz níže) |
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
* Výstupy ve formátu VTK XML: `output_*.vtu`, hladina `outputFS_*.vtp`, `boundary.vtu`, kolekce `output.pvd`/`outputFS.pvd` (ParaView čte v `.pvd` jen XML soubory, starší legacy `.vtk` v kolekci hlásil "Could not determine the data type"). Binární data za hlavičkou (`vtkBinary = 1`, jinak ASCII), kompaktní `fields_<iter>.bin` pro GUI (`writeFields = 1`).
* Komprese výsledků (`compressResults = 1`, výchozí): pole VTK jsou komprimovaná zlibem po blocích 1 MB (`vtkZLibDataCompressor`, ParaView je čte přímo), `fields_<iter>.bin` má verzi 3 s daty bloků v jednom proudu deflate. Kompresor je vlastní (`Deflate.cpp`, bez závislosti na zlib), bloky se komprimují paralelně. Do `.vtu` se zapisují jen body použité buňkami.
* MLUPS z `omp_get_wtime` průměrované přes interval ukládání.

## Zahušťování sítě (octree)

Statické zahušťování podle geometrie: síť vznikne jednou při generování a během výpočtu se nemění (hladinu nesleduje).

* Blok úrovně `L` má stále `ni × nj × nk` buněk velikosti `H / ni / 2^L`. Sousední bloky (i přes hranu a roh) se liší nejvýše o jednu úroveň (2:1), řešič nevyvážené sítě odmítne.
* `mesh/blocks.txt`: řádek `bi bj bk` (úroveň 0, beze změny) nebo `bi bj bk L` se souřadnicemi bloku na jeho úrovni. Soubory bloků úrovně 0 mají dál číslo `nBi·nBj·bk + nBi·bj + bi`, jemnější bloky značku `L<L>_<bi>_<bj>_<bk>` (`blockBoundaryL2_5_3_1.txt`, `blockInitL2_5_3_1.txt`).
* Akustické škálování: `Ct` i `Ch` se na každé úrovni půlí, `Cu` a tlak jsou společné. Jemnější úroveň dělá dva kroky na jeden krok hrubší (sub-cycling), `τ_{L+1} = ½ + 2(τ_L − ½)` – jemné úrovně jsou stabilnější. Iterace a čas výstupů jsou kroky úrovně 0.
* Rozhraní úrovní: objemová metoda Rohde et al. (2006), stejně jako waLBerla. Hrubé post-kolizní populace se homogenně rozdělí do ghost vrstvy jemné úrovně, kterou během dvou jemných kroků proudí. Hrubá buňka dostane průměr populací, které do jejích 8 jemných buněk přitekly z jemné oblasti.
* Ghost bloky: kolem bloků každé úrovně stojí bloky stejné třídy (`Block::ghost`), naplněné stavem jiné úrovně. Proto stávající kód sousedů (`nbBlock`, `nbLocal`, `local2globalShift`) funguje beze změn. Propagace a bilance probíhají jen v buňkách do vzdálenosti 2 od počítaných buněk.
* Hmota se zachovává přesně (strojová přesnost). Rozdíl mezi tím, co hrubá úroveň odevzdá, a tím, co jemná převezme (a naopak), se vrací do hrubé buňky (`massCorr`). Rozdíl vzniká hlavně tam, kde rozhraní úrovní protíná stěnu.
* Přeškálování nerovnovážné části populací (Dupuis–Chopard) bylo vyzkoušeno a odstraněno. V objemovém schématu porušuje zachování hybnosti po jednotlivých populacích a v Poiseuillově testu dávalo chybu +21 % / −11 %.
* Volná hladina přes rozhraní: výměnu hmoty mezi rozhraními úrovní počítá jemná strana a hrubá buňka dostane opačnou hodnotu (`massAcc`). Konverze buněk (plyn → rozhraní, kapalina → rozhraní) se na ghost buňkách provedou hned a pošlou úrovni, která buňku vlastní. Totéž platí pro hmotu redistribuovanou do ghost buněk. Vazba kapalné buňky na plynovou ghost buňku (konverze dorazí o krok později) nepřenáší hmotu. Ghost bloky počítají lokální veličiny (typy rozhraní, hladiny, vyhlazení VOF, normály, křivost). Prázdné bloky se aktivují i přes ghost bloky.
* Omezení: jen `problemID` 0 a 1, ne D2Q9, sudé `ni, nj, nk`. Kde rozhraní úrovní protíná stěnu kolmo na proud, dává Poiseuilleův kanál o výšce 32 hrubých buněk chybu průtoku 1–3,5 % (rovnoběžné rozhraní ~1e-3). Chyba je v rozích rozhraní a stěny, ne ve hmotě.
* Výstupy: VTK obsahuje bloky všech úrovní v souřadnicích buněk úrovně 0, `fields_<iter>.bin` má u zahuštěné sítě verzi 2 s úrovní bloku, MLUPS počítá `2^L` kroků jemných úrovní.
* Testy: `refinement_*` (klidový stav, odmítnutí nevyvážené sítě, zachování hmoty, Poiseuille s rozhraním rovnoběžným i kolmým na proud, volná hladina – protržení hráze přes rozhraní a klidová hladina).

## Neprovedeno

* AA-pattern (jediné pole populací) a SoA uspořádání.
* Zahušťování pro částice a vlákna (`problemID` 2, 21, 3, 4).

## Formát `fields_<iter>.bin`

```
char[4] "LBMF", int32 verze (1, zahuštěná síť 2, komprimovaný 3), int32 iterace, double čas, int32 problemID,
int32 nBlock, int32 ni, nj, nk, double Ch (úroveň 0), int32 nFields, nFields × char[16] název
pro každý blok: int32 bi, bj, bk, [verze 2, 3: int32 úroveň]; pro každé pole float32[ni*nj*nk] (index i + ni*j + ni*nj*k)
verze 3: vše za hlavičkou (od prvního bloku) je jeden raw deflate proud (RFC 1951, bez zlib hlavičky)
pole: ux, uy, uz [m/s], p [Pa], type (typ buňky), u volné hladiny navíc vof, flag
```
