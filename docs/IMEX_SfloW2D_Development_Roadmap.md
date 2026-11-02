# IMEX_SfloW2D — Roadmap generale di sviluppo

**Roadmap operativa unica** · 2 novembre 2026 · Repository `demichie/IMEX_SfloW2D_dev` · Branch `master`.

> Questo è il solo documento operativo da aggiornare con priorità, risultati, lavori aperti e criteri di accettazione. I precedenti tracker D17–D20 e la roadmap storica sono conservati integralmente in `docs/archive/roadmaps/`, ma non sono piani concorrenti. Gli schemi IMEX-RK, le cadenze STAGE/STEP e i percorsi di compatibilità precedenti sono superati dal Lie splitting.

## 1. Obiettivo scientifico

Il codice simula flussi geofisici depth-averaged 2D, multifase e polidispersi, con topografia, reologie configurabili, scambi conservativi, temperatura e pressione interstiziale. Gli obiettivi sono: robustezza e conservazione, HP-PCCU e well-balancing, modularità fisica, parallelismo, riproducibilità, I/O affidabile, estensioni layer-aware, HPC e chiusure subgrid/verticali anche ML.

## 2. Architettura attuale, vincolante

### Lie e transazioni

L'integratore temporale globale è esclusivamente il **Lie splitting del primo ordine**. Il **gruppo 0**, obbligatorio EXPLICIT, valuta l'operatore spaziale sullo stato iniziale `q_n`, ne applica l'incremento e verifica l'ammissibilità **senza limiter dei rate spaziali**. I rate delle sorgenti del gruppo 0 sono valutati sullo stesso `q_n`; i loro incrementi sono limitati rispetto al candidato post-spaziale e agli inventari disponibili. I gruppi seguenti operano in ordine sullo stato aggiornato, con metodo EXPLICIT, ANALYTIC o NEWTON conforme alle capacità della sorgente. Un gruppo NEWTON risolve congiuntamente le equazioni interessate. Il tempo avanza una sola volta; in caso di fallimento, rollback dell'intero passo di conservative, primitive, inventari, fondo e geometria derivata.

### Responsabilità

- `numerics/`: integratore Lie, orchestrazione delle sorgenti, solver locale, estrazione dei dati per cella.
- `spatial/`: flussi, ricostruzione, HP-PCCU, termini spaziali; nessun limiter fisico sui rate spaziali.
- `physics/`: formule e rate dei processi, assemblaggio e limiter delle sole sorgenti fisiche; le routine ricevono **valori locali**, mai indici spaziali `j,k` o campi globali per cella.
- `state/`: `q`, `qp`, deposito, erosione, lava, inventari e aggregati; validità e conversioni di stato.
- `geometry/`: `B_vertex`, topografia e grandezze geometriche derivate, aggiornabili anche con fondo evolutivo.

Ogni sorgente conserva **quattro interfacce** ordinate EXPLICIT, SEMI_IMPLICIT, IMPLICIT, ANALYTIC; le operazioni non disponibili falliscono esplicitamente. Esplicito, semi-implicito e analitico ricevono stato reale; implicito stato complesso. Negli espliciti reologici è lecito convertire localmente lo stato per comporre semi-implicito e implicito senza duplicare la fisica. Evitare wrapper di mero inoltro. Le sorgenti specifiche calcolano anche gli effetti sulle altre equazioni conservative; limitatori condivisi carrier/packing preservano la coerenza dei trasferimenti. L'iniezione di bordo resta distinta dall'iniezione volumetrica.

### Parallelismo

L'ordine dei gruppi Lie è seriale; le celle indipendenti al loro interno sono parallelizzabili con OpenMP. I temporanei devono essere locali e thread-safe, con input immutabile durante la valutazione e commit/rollback coerenti. Il percorso futuro a MPI+OpenMP richiede domini locali, halo e tile di celle attive.

## 3. Stato verificabile al 2 novembre 2026

**Legenda:** *IMPLEMENTATO* = codice e struttura presenti; *VERIFICATO* = gate comunicato o documentato; *PARZIALE* = restano verifiche o interventi; *FUTURO* = lavoro pianificato. Non equiparare compilazione, regressione e accettazione scientifica.

| Ambito | Stato | Evidenza / limite |
| --- | --- | --- |
| Migrazione a Lie, gruppi sorgente e confini architetturali | IMPLEMENTATO | Audit statico D19 del 1 novembre; non equivale a validazione scientifica totale. |
| Ownership di stato/inventari e topografia | IMPLEMENTATO | Riordino di stato e geometria, transazioni e rollback nei percorsi Lie. |
| `physics/` cell-local | VERIFICATO staticamente | Audit dei 34 moduli fisici: nessun indice spaziale `j,k` in ingresso. Residuo da esaminare: scalari globali `pi_g` e `bottom_source_area_ratio` nelle injection. |
| Reologia, scambio di massa, injection, termica, Darcy | IMPLEMENTATO / PARZIALE | Moduli e trattamenti riorganizzati; non estendere alle future revisioni la validità di test precedenti. |
| Compilazione `strict-debug` | VERIFICATO dall'utente | Verde dopo la correzione di `source_registry_2d.f90` comunicata il 2 novembre. |
| **Lie OpenMP validation (D20)** | VERIFICATO dall'utente | Action comunicata verde il 2 novembre; non dimostra assenza universale di race o scaling ottimale. |
| Regressioni numeriche e confronto scientifico complessivo | PARZIALE | Restano da raccogliere e associare a SHA tutte le prove quantitative richieste. |
| Formattazione e Doxygen dei 64 sorgenti | PARZIALE, IN PAUSA | Audit del 2 novembre ha trovato violazioni certe; nessuna modifica ai sorgenti in questa riorganizzazione documentale. |

**Debito di formattazione registrato, temporaneamente non prioritario:** 2 righe di codice oltre 120 colonne, 18 `&` non in colonna 120, 2 righe con TAB, 11 righe con indentazione non multipla di quattro; nessuna riga Doxygen oltre 100 e nessuno spazio finale. Rimangono audit sintattico dei livelli annidati, righe vuote, compattazione argomenti e validità Doxygen. Il checker `tools/docs/check_fortran_style.py` contiene inoltre una sequenza letterale `\\n` da correggere. Verificare le eccezioni prima di intervenire.

## 4. Backlog ordinato per priorità

| Priorità | ID | Attività | Chiusura richiesta |
| --- | --- | --- | --- |
| P0 | V01 | Congelare baseline Lie verificata | Associare strict-debug e D20 al commit effettivamente testato, conservare log. |
| P0 | V02 | Matrice di regressioni scientifiche | Confronti per idrodinamica, wet/dry, well-balancing, sorgenti, massa, energia, pressione, NEWTON, analytiche e topografia. |
| P0 | V03 | Consolidare le prove D20 | Evidenze di Newton effettivo, fondo variato, rollback dell'intero passo dopo fallimento e confronti 1/2/4 thread, incluse le revisioni dei test. |
| P1 | A01 | Chiudere contratti e ownership della fisica | Nessun accesso a celle globali, quattro interfacce coerenti, inventari e limiter verificati, senza wrapper ridondanti. |
| P1 | A02 | Test fisici locali e complex-step | Tipi, rank, unità, direzione degli scambi e aggiornamento cinematico dopo variazioni di velocità. |
| P1 | F01 | Completare conformità di formattazione, **in pausa** | Correggere checker; controllare insieme indentazione reale, spazi, larghezza 120/100, continuazioni e argomenti; garantire invarianti sintattiche e ricontrollare strict-debug/D20. |
| P1 | N01 | Benchmark OpenMP e costo computazionale | Action manuale `lie-openmp-benchmark.yml` e script di timing **implementati**; risultati ancora da eseguire e verificare per 1/2/4/8 thread, con misure su macchina dedicata quando necessario. |
| P2 | S01 | Consolidamento HP-PCCU x/y | Well-balancing, topografia Q1 evolutiva, curvatura, bordi e fronti wet/dry con test quantitativi. |
| P2 | T01 | Controller rifiuto/retry del passo | Limiti di ammissibilità, dt minimo, rollback e diagnostica senza clipping conservativo artificiale. |
| P2 | C01 | Build/CI/I-O/restart | Gate proporzionati al rischio, artefatti, distribuzione e output per run. |
| P3 | L01 | Stato e storage layer-aware | Preservare `N_LAYERS=1`, predisporre halo/I-O/restart multistrato. |
| P3 | H01 | MPI+OpenMP, parallel I/O e GPU readiness | Decomposizione del dominio e scaling verificati. |
| P3 | X01 | Chiusure verticali depth-averaged | Contratti conservativi, riferimento da teoria/esperimenti/modelli verticalmente risolti. |
| P3 | X02 | Modelli subgrid topografici fine/coarse | Simulazioni ad alta e bassa risoluzione **con lo stesso codice**; distinguere effetti dissipativi da deviazioni/redistribuzioni di flussi e correzioni delle interfacce. |
| P3 | X03 | Backend ML/NN opzionale | Inferenza batch, features/versioni, domini di validità, fallback, test su terreni indipendenti. |
| P3 | L02 | Estensione fisica a due strati | Formulazione degli scambi interstrato e gate dedicati. |

## 5. Sequenza operativa

1. Stabilire SHA e prove della baseline corrente (V01–V03), quindi affrontare la validazione numerica/scientifica.
2. Chiudere le sole incoerenze architetturali reali (A01–A02), senza reintrodurre la vecchia struttura IMEX-RK.
3. Riprendere F01 quando prioritario, usando un checker eseguibile e invarianti di preservazione del codice.
4. Consolidare discretizzazione, retry, performance, I/O e CI (S01, T01, N01, C01).
5. Passare alle estensioni P3 dopo aver congelato interfacce e baseline.

Un grande refactoring può avere commit intermedi non compilabili; la verifica avviene al termine del change set completo. Le modifiche di **sola formattazione** invece devono preservare token eseguibili, ordine delle espressioni e direttive OpenMP.

## 6. Regole di manutenzione

Aggiornare **soltanto questo documento** per piano generale e stato: ogni voce deve riportare stato, SHA/evidenza, rischi aperti e prossimo passo. Non aggiungere registri cronologici duplicati. Una Action verde su un commit non certifica commit successivi. Le specifiche archiviate servono per le motivazioni storiche e l'audit, mentre codice e test correnti restano la fonte di verità.

## 7. Archivio tecnico-storico

Conservati integralmente in `docs/archive/roadmaps/`:
- `IMEX_SfloW2D_Development_Roadmap_pre_Lie_2026-10-25.md`: vecchia roadmap dettagliata, comprese sezioni IMEX-RK superate.
- `D17_SOURCE_REFACTOR_HANDOFF.md`: registro e decisioni D17.
- `D18_IMEX_RK_SOURCE_ARCHITECTURE_PRINCIPLES.md`: specifica di transizione precedente.
- `D19_LIE_SPLITTING_MIGRATION_TRACKER.md`: dettaglio della migrazione e audit statici.
- `D20_LIE_OPENMP_PARALLELIZATION_PLAN.md`: sviluppo e gate OpenMP.
- `physical_source_refactoring_plan.md`: dettagli storici del refactoring fisico.

I report di accettazione specifici possono rimanere in `docs/acceptance/` e le suite in `TESTS/`; non sono roadmap alternative.
