# D17 — Audit di compatibilità di esempi, test e generatori di input

**Repository:** `demichie/IMEX_SfloW2D_dev` · **branch:** `feature/d17-source-refactor` · **HEAD:** `c7d81d36f6ea08ce2849cb0eac9c554b0384a91a` · **data:** 25 ottobre 2026.

> **Esito:** audit statico di inventario generale e analisi puntuale delle superfici D17. **Non è una certificazione di compatibilità completa**: emergono errori di compilazione certi in test storici, mentre la validazione dinamica di tutti gli esempi e l'ispezione approfondita di tutti gli script non sono ancora state eseguite. Non segnare un test come PASS perché la CI D17 non lo esegue.

## Avanzamento bonifica — 25 ottobre 2026

È stato preparato un primo change set unitario che aggiorna i test `TEST_IMEX_STAGES/test_stage_count.f90`, `TEST_SPATIAL_OPERATOR/test_spatial_operator.f90`, `TEST_CURVATURE_SOURCE/test_curvature_source.f90`, `TEST_MASS_EXCHANGE/test_mass_exchange.f90` e `TEST_RHEOLOGY_FROUDE/test_rheology_froude.f90` ai nomi e alle firme attuali. Il runner `TEST_RHEOLOGY_CELL_LOCAL/run_test.sh` deriva ora i sorgenti da `src/Makefile.am` e include le dipendenze C/NetCDF.

**Stato delle verifiche:** modifiche di codice proposte e pubblicate, **non compilate né eseguite** in questo ambiente. Il cambiamento del runner non rende automaticamente compatibili `probe_sources.f90` e `test_grouped_implicit_stage.f90`, che richiedono ancora la migrazione del contratto a gruppi nominali e Newton unificato; `test_imex_stages.f90` conserva ancora vecchie costanti di catalogo. I test runtime possono inoltre richiedere esplicitamente un `source_group_type` con `source_names` per esercitare un rate attivo. Non dichiarare superato il gate P0 prima del risanamento dei file residui e di un run reale delle suite.

---

## 1. Perimetro e strategia

Nel tree di HEAD sono presenti **507 file versionati**, di cui **198 sotto TESTS** e **122 sotto EXAMPLES**. Le directory trasversali `tools/`, `UTILS/`, `.github/workflows/` sono rilevanti per generazione di input, orchestrazione, build e verifica. Sono state identificate **21 unità di input/template in EXAMPLES** e **38 sorgenti Fortran di test**. Tutti i 21 input/template principali sono stati letti per i campi D17 pertinenti; i 38 sorgenti Fortran sono stati ispezionati per riferimenti a API abolite o spostate. È stata inoltre ispezionata una selezione dei runner e generatori Python/Shell centrali. Rimangono da esaminare sistematicamente tutti gli altri script e tutti gli input delle sottodirectory TESTS e da eseguire gli smoke test.

**Criteri di classificazione:** **ERRORE CERTO** = simbolo importato o campo eliminato, tale da impedire compilazione/uso del percorso; **DA VERIFICARE** = percorso compatibile per forma ma non eseguito; **NESSUN BLOCCO STATICO D17 RILEVATO** = non trovato il pattern specifico analizzato, non equivale a PASS.

## 2. Incompatibilità certe nei test Fortran (priorità P0)

| Suite | File | Incompatibilità rilevata | Intervento richiesto |
|---|---|---|---|
| TEST_IMEX_STAGES | `test_imex_stages.f90` | `rheology_source_terms`, `SOURCE_TERM_THERMAL_RADIATIVE_COOLING`, `SOURCE_TERM_THERMAL_CONVECTIVE_COOLING`, `SOURCE_TERM_PORE_PRESSURE_DISSIPATION`, `SOURCE_TERM_GEOMETRY_CURVATURE`; import iniezione da `equation_terms_2d` | Nomi fisici canonici e `source_terms_2d`; mantenere intenti e confronti del test |
| TEST_IMEX_STAGES | `test_stage_count.f90` | `SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS`, `SOURCE_TERM_BOTTOM_CIRCULAR_INJECTION`, `%term_ids` | Usare `source_names` e verificare la stessa cadenza, identità e ordine |
| TEST_MASS_EXCHANGE | `test_mass_exchange.f90` | `USE equation_terms_2d, ONLY : eval_mass_exchange_terms`, `source_solver_type`; vecchio descrittore con `no_solvers` | Richiamare la nuova API rate tramite `source_terms_2d` o coordinatore fisico per test unitario; rimuovere sotto-solver, senza aggirare il registro |
| TEST_RHEOLOGY_CELL_LOCAL | `probe_sources.f90` | `eval_implicit_terms`, `eval_semi_implicit_terms`, `eval_mass_exchange_terms` dal vecchio `equation_terms_2d`; ID `SOURCE_TERM_*` | Migrare agli assemblatori in `source_terms_2d` e ai gruppi nominali |
| TEST_RHEOLOGY_CELL_LOCAL | `test_grouped_implicit_stage.f90` | `source_solver_type`, `filter_solvers`, `term_ids`, ID termici separati e `eval_implicit_terms` | Ricostruire i test con un solo Newton accoppiato per gruppo e sorgenti unificate |
| TEST_RHEOLOGY_CELL_LOCAL | `run_test.sh` | Elenco hardcoded con `physics/mass_exchange_laws_2d` inesistente nel layout attuale | Usare `src/Makefile.am` quale fonte ordinata dei sorgenti o aggiornare l'elenco con `physics/mass_exchange/mass_exchange_laws_2d.f90` e relative dipendenze |
| TEST_CURVATURE_SOURCE | `test_curvature_source.f90` | Import da `equation_terms_2d` di curvatura/iniezione/vecchio `eval_expl_terms`; `%term_ids` | Migrare a `source_terms_2d::eval_source_explicit` e `source_names` |
| TEST_RHEOLOGY_FROUDE | `test_rheology_froude.f90` | Import `eval_semi_implicit_terms` / `eval_implicit_terms` da `equation_terms_2d` | Usare `eval_source_semi_implicit` / `eval_source_implicit`, preservando il confronto REAL/COMPLEX |
| TEST_SPATIAL_OPERATOR | `test_spatial_operator.f90` | `eval_expl_terms`, `SOURCE_TERM_GEOMETRY_CURVATURE`, `%term_ids` | Migrare ad assemblatore e gruppi canonici |

**I test di questa tabella non risultano compatibili con il codice D17 attuale nelle chiamate indicate.** Non sostituire simboli con no-op né riscrivere il riferimento scientifico per far passare un controllo.

## 3. Suite con nessun blocco D17 evidente nella porzione ispezionata

`TEST_SOURCE_REGISTRY` (suite CI corrente); `TEST_LAPACK_INTERFACES`, `TEST_EQUATION_METADATA`, `TEST_MODEL_LAYOUT`, `TEST_GEOMETRY_1D`, `TEST_HP_RECONSTRUCTION`, `TEST_HYDROSTATIC_PATH`, `TEST_DIAGNOSTICS`, `TEST_NONCONSERVATIVE_PATH`, `TEST_OUTPUT_PATHS`, `TEST_PCCU`, `TEST_STATE_CONVERSION`, `TEST_STATE_INITIALIZATION`, `TEST_STOCHASTIC_EVOLUTION`, `TEST_STOCHASTIC_RANDOM`, `TEST_TOPOGRAPHY_FACES`, `TEST_WORKSPACE_LIFECYCLE` e test locali specifici lava/pressione. I quattro test Fortran numerici `TESTS/ACCEPTANCE/test_d17_*.f90` non mostrano i pattern di API obsolete controllati; verificare sempre la propria Action e il profilo.

Il file `TESTS/TEST_HYDROSTATIC_ROUNDOFF/reference/hyperbolic_2d.f90` importa ancora alcuni operatori idrodinamici da `equation_terms_2d`: questa dipendenza va verificata simbolo per simbolo rispetto al modulo attuale; il semplice nome del modulo non basta a considerarla un errore.

## 4. Esempi, input e template

**21 input/template sotto EXAMPLES esaminati** nei percorsi `EXAMPLE_2D`, `EXAMPLE_2D_PP`, `EXAMPLE_BOX`, `EXAMPLE_BUMP` (quattro regimi), `EXAMPLE_BW`, `EXAMPLE_BW2D`, `EXAMPLE_CF`, `EXAMPLE_DAMBREAK`, `EXAMPLE_ENGWELL2016`, `EXAMPLE_ENSEMBLE_ETNA/templatedir`, `EXAMPLE_ETNA`, `EXAMPLE_ETNA_LAVA`, `EXAMPLE_INCLINED_EXCAVATION_2D`, `EXAMPLE_LATERAL_SOURCE`, `EXAMPLE_SLOPE`, `EXAMPLE_TAAL`, `EXAMPLE_TAAL_BOTTOM`.

- **RHEOLOGY_MODEL**: le configurazioni esaminate usano stringhe come `'quadratic'`, `'voellmy'`, `'mui'`, `'lava'`, `'none'` e non ID numerici di reologia. Per i casi con `RHEOLOGY_FLAG=F`, la stringa impostata ma non attiva non è di per sé un errore.
- **SOURCE_GROUPS**: non risulta necessario introdurlo automaticamente nei template legacy; il fallback di configurazione D17 è deliberato. Non migrare silenziosamente i benchmark ai gruppi nativi, perché cambierebbe la semantica di confronto.
- **THERMAL_COOLING_FLAG**: non compare nelle configurazioni esaminate. Verificare particolarmente `EXAMPLE_ETNA_LAVA`: l'assenza del flag deve ricadere sulla compatibilità legacy; aggiungere il flag soltanto dove l'aspettativa scientifica lo richiede e confrontare il raffreddamento.
- **RUN_NAME, restart e placeholder**: i template generati richiedono l'esecuzione degli script di generazione prima del run; non trattare un template con placeholder come input immediatamente eseguibile.
- **Mass exchange**: esempi con `SETTLING_FLAG=T` devono verificare deposizione/perdita dipendente e limitatore comune nel nuovo gruppo STEP, mantenendo gli inventari e la topografia.

**Nessuna equivalenza dinamica degli esempi è stata dimostrata da questa lettura statica.**

## 5. Script e generatori: verifiche e problemi

| Percorso | Osservazione | Giudizio |
|---|---|---|
| `EXAMPLES/EXAMPLE_2D/create_example.py` | Produce `IMEX_SfloW2D.inp` sostituendo placeholder di run, restart e griglia | Schema di generazione compatibile per forma; run da verificare |
| `EXAMPLES/EXAMPLE_DAMBREAK/create_dambreak.py` | Legge il template e genera input/restart per dam break | Verificare file generato e lancio |
| `EXAMPLES/EXAMPLE_ENSEMBLE_ETNA/create_inputfiles.py` | Riscrive file del template a partire da campi tabellari | Verificare che le sostituzioni non alterino nomi/opzioni D17 e non lascino placeholder |
| `EXAMPLES/EXAMPLE_ETNA_LAVA/prepare.py` | Preserva input e DEM esistenti, senza sovrascrivere parametri | Verificare compatibilità raffreddamento/fallback e output |
| `TESTS/run_tests.sh` | Esegue casi e manipola file input per verifiche di naming e sicurezza | Testare input generati; non assumere che copra le suite Fortran non eseguite |
| `TESTS/TEST_IMEX_STAGES/run_test.sh`, `TESTS/TEST_MASS_EXCHANGE/run_test.sh` | Estraggono sorgenti di build da `src/Makefile.am` | La lista delle dipendenze di produzione viene aggiornata, ma i programmi di test importano API obsolete |
| `TESTS/TEST_RHEOLOGY_CELL_LOCAL/run_test.sh` | Elenco hardcoded dei sorgenti con percorso precedente di mass exchange | **ERRORE CERTO** di build |
| `TESTS/ACCEPTANCE/run_acceptance.py` | Harness di compilazione/esecuzione; non sostituisce i singoli test | Valutare suite richieste e toolchain prima della certificazione |
| Workflow D17 | La suite corrente valida catalogo, interfacce e nuovi test quantitativi | **PASS nelle Actions confermate**, ma non certifica le suite storiche sopra elencate |

## 6. Piano di bonifica senza perdita di contenuto scientifico

1. **P0 — Riparazione test non compilabili:** migrare in blocco le suite `TEST_IMEX_STAGES`, `TEST_MASS_EXCHANGE`, `TEST_RHEOLOGY_CELL_LOCAL`, `TEST_CURVATURE_SOURCE`, `TEST_RHEOLOGY_FROUDE`, `TEST_SPATIAL_OPERATOR`; per ciascuna preservare asserzioni, casi, tolleranze e semantica, rifiutando adattatori permanenti alle API soppresse.
2. **P1 — Audit generazione input:** esaminare tutti gli script Shell/Python e tutti gli input in `TESTS` e `EXAMPLES`, renderli riproducibili, convalidare i placeholder e i fallback reologici/termici/STEP. Questa parte **non è ancora completamente svolta** nel presente audit.
3. **P2 — Esecuzione selettiva:** eseguire la suite risanata e un smoke test per ciascuna famiglia di esempi con input generato; annotare exit code, errori di parser, output, restart e metadati. Nessun nuovo benchmark scientifico indiscriminato.
4. **P3 — Tracciabilità:** mappare ogni `run_test.sh`, caso e generator al suo stato `PASS / FAIL / NON ESEGUITO`, commit e toolchain; aggiornare roadmap e `docs/acceptance` senza dichiarare completato ciò che non è stato verificato.

**Decisione:** la compatibilità D17 del repository nel suo complesso **non è ancora certificabile**. Prima di chiudere la manutenzione delle suite storiche occorre almeno correggere i blocchi P0 e completare la ricognizione di P1. Non modificare il comportamento fisico né i default scientifici come scorciatoia per ottenere test verdi.
