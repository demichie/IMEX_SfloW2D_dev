# IMEX_SfloW2D — Development Roadmap

**Piano unico, integrato e consolidato | 14 ottobre 2026 | repository fd43e4ca2ee6**

Un unico indice operativo e una sola scheda per ogni attività. La tabella espone lo stato corrente; tutte le specifiche e i gate sono riportati nella scheda a cui si riferiscono. Le attività implementate ma non ancora accettate sono distinte dalle milestone completate.

# 1. Visione scientifica e funzionale del modello a sviluppo completato

Scopo della sezione. Questa presentazione descrive il sistema scientifico e numerico che IMEX_SfloW2D dovrà diventare a conclusione del programma. È rivolta anche a chi non conosce il sorgente. Le capacità obiettivo non vanno confuse con funzionalità già implementate: la roadmap e i gate successivi distinguono baseline verificata, refactoring, estensioni scientifiche e futuri backend HPC.

### A. Problemi fisici, applicazioni e filosofia di modellazione

IMEX_SfloW2D è una piattaforma per la simulazione di flussi geofisici e vulcanici su topografie naturali, mediante equazioni multifase depth-averaged. Il sistema è pensato per flussi di materiale granulare e miscele fluide particellate, potenzialmente polidisperse, dalle correnti piroclastiche a diversa concentrazione alle colate laviche e ad altri flussi gravitativi rappresentabili dalle formulazioni costitutive disponibili. Le applicazioni comprendono la propagazione, la canalizzazione e la deposizione, la valutazione delle aree invase e l’analisi di scenari su modelli digitali del terreno.

La formulazione depth-averaged riduce il costo rispetto ai modelli tridimensionali verticalmente risolti, ma introduce problemi di chiusura: le distribuzioni verticali delle componenti, la microtopografia non risolta e diversi meccanismi di scambio non sono ricostruibili direttamente dalla sola media di profondità. L’obiettivo della nuova architettura è esplicitare tali limiti, integrare chiusure verificabili e mantenere separati i meccanismi fisici dagli algoritmi numerici e dall’infrastruttura di calcolo.

Il sistema finale dovrà poter eseguire studi scientifici controllati, simulazioni su DEM reali, confronti fra differenti modelli reologici e campagne di calcolo sistematiche, documentando la configurazione effettiva, le ipotesi costitutive, la qualità della soluzione e l’origine dei dati di riferimento.

### B. Sistema multifase, polidispersione e variabili conservative

Il modello dovrà descrivere miscele composte da fasi vettore e componenti solide trasportate, incluse classi granulometriche distinte quando richieste dal modello. Le classi trasportate mantengono identità e bilanci di massa coerenti, così che granulometria, densità, concentrazione e proprietà termodinamiche contribuiscano ai flussi e alle sorgenti. Le variabili conservative costituiscono la base del bilancio numerico; le primitive di cella sono ricostruite mediante una trasformazione fisicamente consistente, con trattamento esplicito delle celle secche e quasi secche.

La termodinamica delle miscele, la conversione fra frazioni massiche e volumetriche, la densità del vettore e le proprietà efficaci devono essere calcolate da un’unica interpretazione dello stato, preservando i kernel condivisi REAL/COMPLEX e il percorso di derivazione complex-step dei Jacobiani. L’obiettivo non è assumere un modello costitutivo unico per tutti i flussi, bensì permettere formulazioni diverse entro uno stesso quadro di conservazione e validazione.

### C. Reologie alternative e proprietà costitutive

La piattaforma dovrà consentire di selezionare reologie diverse senza modificare il ciclo IMEX o il motore generale di Newton. I modelli reologici esistenti verranno estratti in moduli autonomi, selezionati tramite identificatori leggibili e registrati all’inizializzazione. Le implementazioni dovranno conservare le specificità fisiche del modello, inclusi eventuali regimi limite, dissipazione, frenamento e dipendenze termiche, senza imporre una correzione della quantità di moto a sorgenti che non rappresentano esclusivamente attrito.

L’architettura dovrà garantire che l’aggiunta o il confronto di una reologia abbia un perimetro locale e sia accompagnata da caratterizzazione scientifica REAL/COMPLEX, controlli sulle derivate, test di stati degeneri e confronto rispetto alla baseline. Questa fase di modularizzazione riguarda per ora i soli modelli reologici: non autorizza un refactoring generalizzato e simultaneo di tutti i processi costitutivi.

### D. Uno o due strati e processi di interazione

Il sistema dovrà evolvere da una formulazione produttiva a strato singolo verso una scelta runtime fra uno e due strati fisici, preservando inizialmente il comportamento N_LAYERS=1. Il futuro modello a due strati dovrà distinguere le variabili e i processi di ciascuno strato e rappresentare i trasferimenti interstrato di massa, quantità di moto ed energia, oltre agli accoppiamenti idrostatici e non conservativi richiesti dalle equazioni. La modellazione dei due strati è una vera estensione fisica, non una duplicazione meccanica del codice a strato singolo.

La struttura delle variabili e delle equazioni dovrà essere descritta da model_layout_type, con indici, dimensioni e partizioni definite una sola volta. Lo storage sarà gradualmente reso layer-aware prima dell’attivazione delle nuove equazioni. I kernel comuni resteranno condivisi; i trasferimenti interstrato saranno distinti sia dalle sorgenti esterne sia dai prodotti spaziali non conservativi.

### E. Discretizzazione spaziale, topografia, wet/dry e well-balancing

La dinamica su terreno naturale richiede una discretizzazione capace di trattare variazioni anche forti di quota, velocità caratteristiche direzionali, gradienti geometrici, fronti di avanzamento, arresto del flusso e celle alternativamente asciutte o bagnate. L’infrastruttura HP-PCCU accoppia flussi convettivi, ricostruzione coerente, contributi di percorso idrostatici alle celle e alle interfacce, e propagazione numerica mediante velocità caratteristiche locali. L’obiettivo è preservare gli equilibri appropriati senza perdere robustezza nella dinamica di propagazione.

La topografia nodale B_vertex rappresenta la geometria autoritativa; da essa si ottengono quote al centro e alle facce coerenti e le derivate richieste dal modello. La geometria non deve diventare permanentemente statica: la stessa infrastruttura dovrà permettere l’aggiornamento di un fondo che evolve, con refresh delle quantità derivate. Il motore numerico dei prodotti non conservativi rimane generico, mentre i coefficienti fisici e i provider di equilibrio dipendono dal modello.

Il consolidamento direzionale x/y ridurrà duplicazioni e differenze accidentali tra i due assi, preservando la semantica dell’operatore e i gate di regressione. L’interfaccia dell’operatore spaziale dovrà rimanere indipendente dagli specifici stencils e dalle tecniche di parallelizzazione.

### F. Integrazione temporale robusta e stati fisicamente ammissibili

Il modello conserva gli schemi IMEX Runge–Kutta già approvati; non è previsto in questa roadmap un programma di sostituzione o ottimizzazione dei tableaux. Il controllo della rigidità e delle sorgenti rimane affidato alla struttura esplicita/implicita esistente e ai miglioramenti architetturali del solve locale. Per aumentare la robustezza introdurremo invece un meccanismo generico di time-step rejection and retry: un tentativo di passo è accettato soltanto se i suoi stadi e lo stato finale superano verifiche numeriche e fisiche esplicitamente definite.

In caso di fallimento recuperabile il tempo fisico non avanza, si ripristina lo stato di inizio passo, si riduce Δt e si ripete la stessa integrazione. Se il nuovo Δt scende sotto un minimo configurabile, oppure sono esauriti i retry, il run termina con una diagnostica completa. I limiti inferiori e superiori delle variabili primitive saranno definiti da default generali, eventuali default specifici del modello fisico/reologico e override puntuali di input, con precedenza esplicita e validazione all’avvio.

Il controllo delle primitive verrà eseguito soltanto sugli stati delle celle effettivamente risolte DOPO ciascuno stadio RK completato e DOPO l’assemblaggio finale, senza inserirlo nelle conversioni conservative→primitive usate sulle facce di ricostruzione. I fallimenti numerici che impediscono di terminare uno stadio, come il mancato successo del solve implicito, devono invece propagare tempestivamente uno stato recuperabile. Nessun limite viene applicato mediante clipping artificiale dello stato conservativo.

### G. Chiusure della struttura verticale nel modello depth-averaged

Un modello depth-averaged non risolve esplicitamente le distribuzioni lungo la verticale. Tuttavia, concentrazioni solide e distribuzioni granulometriche possono essere fortemente stratificate; anche velocità, temperatura e altre grandezze possono variare con la quota. I termini integrati di flusso, trasporto e scambio dipendono spesso da prodotti non lineari di questi campi, per cui sostituire indiscriminatamente il prodotto delle medie alla media del prodotto può introdurre errori sistematici.

La nuova famiglia di chiusure verticali dovrà rappresentare soprattutto i profili delle classi solide e le loro interazioni con altri campi verticalmente distribuiti, restituendo alle equazioni depth-averaged quantità integrali consistenti. Sono possibili formulazioni analitiche, fisiche e semiempiriche oppure chiusure ML/NN, ma tutte devono esporre lo stesso contratto fisico e rispettare vincoli di massa, normalizzazione e ammissibilità.

Per questo sviluppo non è sufficiente raffinare la griglia orizzontale di IMEX_SfloW2D: il codice non risolve la coordinata verticale e non può generare da solo il riferimento dettagliato che occorre apprendere o parametrizzare. Le fonti di validazione dovranno quindi includere teoria, osservazioni/esperimenti e, quando disponibili, simulazioni di modelli verticalmente risolti. Questi ultimi sono fonti offline di dati scientifici, non simulatori esterni da accoppiare durante il run.

### H. Modelli subgrid topografici basati sul confronto fine/coarse

Per “subgrid” in questa roadmap si intende prioritariamente l’effetto di una descrizione del terreno molto più dettagliata della griglia impiegata nella simulazione produttiva. Piccoli canali, creste, ostacoli, discontinuità morfologiche e variazioni di pendenza possono influire su avanzamento, canalizzazione, blocco locale, direzione del flusso e distribuzione spaziale degli spessori anche quando non vengono risolti dalla griglia grossolana.

A differenza delle chiusure verticali, qui il riferimento può essere ottenuto con IMEX_SfloW2D stesso: eseguire lo stesso problema fisico su DEM e griglie fini, confrontarlo con il corrispondente problema grossolano e identificare gli effetti sistematici delle strutture topografiche eliminate dal coarse-graining. Occorrono insiemi di casi con terreno, condizioni iniziali, reologie, granulometrie e risoluzioni diverse; il confronto deve esplicitare come aggregare campi, topografia, flussi e diagnostiche tra le due discretizzazioni.

Le chiusure potranno essere costruite in forma fisica, semiempirica o ML/NN e utilizzate nel solver coarse. Non si assume a priori che l’effetto subgrid sia un semplice attrito aggiuntivo: potranno emergere correzioni direzionali al trasporto, blocco parziale, effetti sulle interfacce o altri contributi, che dovranno essere identificati e verificati prima di scegliere dove inserirli nelle equazioni. La validazione richiede dati e terreni indipendenti da quelli usati per la calibrazione o l’addestramento, insieme a controlli su conservazione, simmetria, stabilità e generalizzazione.

### I. ML/NN come backend interno delle chiusure

Il machine learning non costituisce un simulatore separato e non sostituisce automaticamente le equazioni di bilancio. Verrà introdotto come backend opzionale per valutare chiusure verticali o topografiche subgrid attraverso inferenza durante l’esecuzione di IMEX_SfloW2D. Le medesime API dovranno consentire alternative analitiche, semiempiriche e neurali; la scelta della tecnologia non deve propagarsi nell’integratore IMEX, nel calcolo dei flussi o nell’infrastruttura MPI.

Il contratto ML comprende definizione/versionamento delle feature, unità e normalizzazioni, intervalli di validità, pesi e metadati del modello, riproducibilità della valutazione, gestione di output non ammissibili e fallback espliciti. Per la futura esecuzione MPI/GPU l’inferenza dovrà poter avvenire in batch, evitando allocazioni per cella e assicurando indipendenza dal partizionamento quando richiesto. L’accoppiamento operativo con modelli numerici esterni a IMEX_SfloW2D è escluso dal perimetro di questa roadmap.

### J. Esecuzione parallela, I/O scientifico e riproducibilità

Il modello finale dovrà essere utilizzabile in modalità seriale, OpenMP e MPI+OpenMP, con decomposizione in domini locali, halo espliciti, tile di celle attive e kernel progressivamente predisposti a GPU. Lo storage e la comunicazione dovranno conoscere il layout numerico e le dimensioni, non il significato fisico delle variabili. L’ottimizzazione del carico e l’inferenza neurale in batch verranno dopo la stabilizzazione della correttezza e delle interfacce.

Le simulazioni dovranno produrre campi scientifici e diagnostiche in formati documentati, con gestione coerente dei restart, degli output per run, delle metainformazioni e del NetCDF seriale o parallelo. Lo stato necessario a ripetere un passo deve essere distinguibile dal workspace temporaneo: questa distinzione è essenziale tanto per il restart deterministico quanto per il rollback dopo un tentativo di timestep rifiutato. Il sistema di build, i test, la CI, Doxygen, MkDocs e la documentazione verificata degli esempi fanno parte dell’infrastruttura di ricerca, non sono accessori opzionali.

### K. Risultato atteso e criteri generali di credibilità scientifica

A sviluppo completato, IMEX_SfloW2D sarà una piattaforma depth-averaged multifase, polidispersa, selezionabile a uno o due strati, dotata di reologie modulari, operatore spaziale HP-PCCU generico e validato, robustezza temporale con controllo di ammissibilità, chiusure verticali e topografiche subgrid opzionali anche ML/NN, supporto a topografia evolutiva, I/O scientifico riproducibile e calcolo distribuito predisposto ad acceleratori.

Il principio trasversale resta la separazione fra (i) formulazione e validazione fisica, (ii) discretizzazione spaziale e temporale, (iii) gestione dello stato e ammissibilità, (iv) infrastrutture HPC, dati e documentazione. Ogni nuovo sviluppo deve specificare che cosa cambia scientificamente, rispetto a quale baseline viene confrontato e quale evidenza permette di considerarlo accettato.

# 2. Principi e obiettivi dello sviluppo

Sviluppare IMEX_SfloW2D mantenendo la baseline scientifica HP-PCCU, l’equivalenza numerica dei refactoring e la possibilità di estendere progressivamente reologie, metodi impliciti, storage layer-aware, elaborazione MPI/HPC, chiusure verticali e subgrid e backend ML/NN. Non sono previsti accoppiamenti runtime con modelli esterni.

Stati: COMPLETATO = verifiche e gate accettati; IMPLEMENTATO = codice presente ma gate non ancora approvato; PARZIALE = attività iniziata o copertura incompleta; DA FARE = non implementato; GATE = verifica/approvazione da effettuare. Lo stato non equivale alla descrizione storica di una vecchia roadmap.

# 3. Master Development Plan — ordine e stato corrente

## Fase A — Fondamenti numerici e baseline scientifica

|     | Vecchia label | Attività                                             | Stato      |
|-----|---------------|------------------------------------------------------|------------|
| A01 | N0            | Baseline, test di riferimento e manifest complessivo | PARZIALE   |
| A02 | N1            | Primitive canoniche e mass fractions                 | COMPLETATO |
| A03 | N2            | Topografia Q1 e geometria B_vertex                   | COMPLETATO |
| A04 | N3            | Ricostruzione HP                                     | COMPLETATO |
| A05 | N4            | Operatore HP-PCCU, G=1                               | COMPLETATO |
| A06 | N5            | Slope correction                                     | COMPLETATO |
| A07 | N6            | Curvature source                                     | COMPLETATO |
| A08 | N7            | Validazione idrodinamica 1D e 2D                     | COMPLETATO |
| A09 | D-N7-CFL      | Due policy CFL                                       | COMPLETATO |
| A10 | N8            | Topografia evolutiva e scambio di massa              | COMPLETATO |
| A11 | N9            | Cutover e baseline scientifica                       | APPROVATO  |

## Fase B — Build e consolidamento dell’operatore x/y

|     | Vecchia label | Attività                                       | Stato                     |
|-----|---------------|------------------------------------------------|---------------------------|
| B01 | M12c          | Build, CI e portabilita                        | COMPLETATO                |
| B02 | P1            | Mappa dipendenze e relocation path-only        | COMPLETATO                |
| B03 | M12d          | Stocastico OU, RNG e restart seriale           | IMPLEMENTATO              |
| B04 | M13           | Model layout, N_LAYERS=1                       | IMPLEMENTATO              |
| B05 | XY-1          | Velocita caratteristiche comuni                | IMPLEMENTATO              |
| B06 | XY-2          | Ricostruzione HP sulle linee                   | IMPLEMENTATO              |
| B07 | XY-3          | Flusso PCCU di interfaccia                     | IMPLEMENTATO              |
| B08 | XY-4          | Velocita faccia e path idrostatici di cella    | IMPLEMENTATO, GATE APERTO |
| B09 | XY-5          | Condizioni al contorno e ricostruzione lineare | IMPLEMENTATO, GATE APERTO |
| B10 | XY-6          | Derivate e halo geometrici                     | IMPLEMENTATO, GATE APERTO |
| B11 | XY-7          | Riduzione CFL direzionale                      | IMPLEMENTATO, GATE APERTO |
| B12 | XY-GATE       | Accettazione complessiva x/y                   | GATE                      |

Evidenza B05–B11 raccolta al 15 ottobre 2026: le fattorizzazioni XY-1…XY-7 sono implementate con le evidenze elencate nelle rispettive schede. Sono passate le regressioni spaziali e dinamiche strict/optimized, incluse prove a 1/4 thread, composizione variabile, fixture ruotate e confronti binari pre/post refactoring sui casi esplicitamente elencati. B12 resta aperto per l’accettazione complessiva.

## Fase C — Output e distribuzione

|     | Vecchia label | Attività                                | Stato   |
|-----|---------------|-----------------------------------------|---------|
| C01 | IO-1          | Output segregati per run                | DA FARE |
| C02 | IO-2          | Contratto di configurazione e sicurezza | DA FARE |
| C03 | IO-GATE       | Accettazione I/O per run                | GATE    |
| C04 | DIST-1        | Autotools policy B                      | DA FARE |
| C05 | DIST-2        | Doppio archivio di distribuzione        | DA FARE |
| C06 | DIST-GATE     | Validazione distributiva                | GATE    |

## Fase D — Reologie e solver IMEX

|     | Vecchia label | Attività                                                    | Stato    |
|-----|---------------|-------------------------------------------------------------|----------|
| D01 | M12e-A        | Inventario formule e rami REAL/COMPLEX                      | PARZIALE |
| D02 | M12e-B        | Caratterizzazione numerica cell-local                       | DA FARE  |
| D03 | —             | Audit e ownership di parametri, variabili e stato condiviso | DA FARE  |
| D04 | P2            | Eliminazione reologia 8                                     | DA FARE  |
| D05 | P3            | Estrazione moduli reologici                                 | DA FARE  |
| D06 | P3            | Sorgenti termiche e pore pressure                           | DA FARE  |
| D07 | —             | Migrazione progressiva dei parametri non reologici          | DA FARE  |
| D08 | P4            | RHEOLOGY_MODEL da numeri a stringhe                         | DA FARE  |
| D09 | P5            | Correzione del solo attrito                                 | DA FARE  |
| D10 | P6            | EX_terms / IM_terms                                         | DA FARE  |
| D11 | P7            | IM_group e metadati                                         | DA FARE  |
| D12 | P7            | Newton per gruppo e Jacobiano ridotto                       | DA FARE  |
| D13 | P8            | Algebriche ANALYTIC per stage DIRK                          | DA FARE  |
| D14 | P8            | Gruppi misti e validazione                                  | DA FARE  |
| D15 | P9            | Acceptance della nuova architettura                         | GATE     |

## Fase E — Architettura e robustezza del passo

|     | Vecchia label | Attività                              | Stato   |
|-----|---------------|---------------------------------------|---------|
| E01 | ARCH-1        | Audit dipendenze e ownership stato    | DA FARE |
| E02 | ARCH-2        | Architecture fitness checks           | DA FARE |
| E03 | ARCH-GATE     | Architecture readiness                | GATE    |
| E04 | TIME-1        | Policy timestep e input               | DA FARE |
| E05 | TIME-2        | Bounds primitive e override           | DA FARE |
| E06 | TIME-3        | Controlli di ammissibilita cell-local | DA FARE |
| E07 | TIME-4        | Error propagation, rollback e retry   | DA FARE |
| E08 | TIME-5        | Test del timestep controller          | DA FARE |
| E09 | P-TIME        | Gate timestep                         | GATE    |

## Fase F — Interfacce per chiusure future

|     | Vecchia label | Attività                       | Stato   |
|-----|---------------|--------------------------------|---------|
| F01 | V-API         | Interfaccia closure verticali  | DA FARE |
| F02 | S-API         | Interfaccia subgrid terreno    | DA FARE |
| F03 | ML-API        | Backend ML/NN opzionale        | DA FARE |
| F04 | EXT-TEST      | Test di neutralita API         | DA FARE |
| F05 | EXT-GATE      | Architecture closure readiness | GATE    |

## Fase G — Migrazione layer-aware

|     | Vecchia label | Attività                         | Stato   |
|-----|---------------|----------------------------------|---------|
| G01 | M14           | Storage rank-4 con N_LAYERS=1    | DA FARE |
| G02 | M15           | Activity/restart/I-O layer-aware | DA FARE |
| G03 | M16           | Architecture freeze              | GATE    |

## Fase H — MPI/HPC e prestazioni

|     | Vecchia label | Attività                        | Stato   |
|-----|---------------|---------------------------------|---------|
| H01 | M17           | Local domain e halo abstraction | DA FARE |
| H02 | M18           | Gerarchia active tiles          | DA FARE |
| H03 | M19           | Decomposizione MPI              | DA FARE |
| H04 | M20           | MPI+OpenMP e Parallel NetCDF    | DA FARE |
| H05 | M21           | Tuning e GPU readiness          | DA FARE |

## Fase I — Estensioni scientifiche

|     | Vecchia label | Attività                             | Stato  |
|-----|---------------|--------------------------------------|--------|
| I01 | LAYERS-2      | Formulazione fisica due strati       | FUTURO |
| I02 | LAYERS-EX     | Scambi interstrato                   | FUTURO |
| I03 | LAYERS-GATE   | Validazione due strati               | GATE   |
| I04 | V-1           | Closure verticali fisiche/analitiche | FUTURO |
| I05 | V-GATE        | Validazione verticale                | GATE   |
| I06 | S-1           | Subgrid topografia fine/coarse       | FUTURO |
| I07 | S-GATE        | Validazione subgrid                  | GATE   |
| I08 | ML-1          | Inferenza ML/NN delle closure        | FUTURO |
| I09 | ML-GATE       | Validazione ML/NN                    | GATE   |

# 4. Detailed Development Plan — schede operative

## Fase A — Fondamenti numerici e baseline scientifica

### A01 — Baseline, test di riferimento e manifest complessivo

Stato ufficiale: PARZIALE

test/reference presenti; manifest di accettazione complessivo da chiudere

**Attività previste:** Congelare commit HEAD, input, checksum, toolchain e tolleranze; collegare reference Python 1D/2D al solver HP-PCCU; identificare comparatori legacy storici senza ricreare backend rimossi; predisporre test ripetibili.

**File interessati:** TESTS/TEST_SPATIAL_OPERATOR/\*; TESTS/ACCEPTANCE/\*; docs/acceptance/\*; TESTS/run_tests.sh; README.md.

**Criteri di accettazione:** Reference e manifest riproducibili; sorgenti Python e Fortran verificati; risultati numerici tracciabili e test lanciati con comando documentato.

#### Pre-HP structural/reference gate

Strict/debug and optimized builds of 93314f81.

Existing deterministic and stochastic restart regressions.

Representative KT cases with one and four OpenMP threads.

Frozen clean Python HP-PCCU references and checksums.

Focused state-conversion/path/geometry tests available before production edits.

### A02 — Primitive canoniche e mass fractions

Stato ufficiale: COMPLETATO

### A03 — Topografia Q1 e geometria B_vertex

Stato ufficiale: COMPLETATO

### A04 — Ricostruzione HP

Stato ufficiale: COMPLETATO

### A05 — Operatore HP-PCCU, G=1

Stato ufficiale: COMPLETATO

### A06 — Slope correction

Stato ufficiale: COMPLETATO

### A07 — Curvature source

Stato ufficiale: COMPLETATO

### A08 — Validazione idrodinamica 1D e 2D

Stato ufficiale: COMPLETATO

#### A02 - Canonical mass-fraction primitive state

Attività e verifiche operative: Usare Y (frazioni massiche) come variabili primitive; verificare catena q↔qp, T, c_p,mix e closure; conversione sorgenti/BC, flussi di componenti con segno positivo/negativo e comportamento a secco.

File interessati: src/state_conversion_2d.f90; src/state_conversion_phys_var.inc; src/mixture_from_mass.inc; src/reconstruction_2d.f90; src/inpout_2d.f90; TESTS/TEST_STATE_CONVERSION/\*.

Criteri di accettazione: Round-trip e bilanci di massa per composizioni variabili; trattamento wet/dry; test sorgenti e BC; nessuna ricomparsa alpha_flag come switch della ricostruzione.

Change qp composition slots to independently transported mass fractions Y.

Use the conservative chain q -\> Y -\> c_p,mix -\> T -\> rho/alpha/gprime -\> h; do not independently reconstruct alpha, density, Gamma or G.

Update qp_to_qc to rebuild conservative component masses from M\*Y after the local face-state closure guard.

Reconstruct each Y component with the existing monotone/TVD scalar limiter machinery.

Generalize the current post-KT solid-flux safeguard to all independently transported component masses and make the signed correction symmetric for positive/negative total mass flux.

Migrate BC/source semantics and examples; only after all call sites are converted, remove alpha_flag and rename misleading alpha-specific structural indices where practical.

Keep maximum solid packing as a separate constitutive condition after the Y,T -\> alpha map.

#### A03 - Continuous Q1 topography and geometry refresh

Attività e verifiche operative: B_vertex autoritativo; derivare B_cent e quote condivise alle facce da Q1; rigenerare LS pendenze/Hessiane e G dopo variazioni topografiche; non congelare geometria nel tempo.

File interessati: src/geometry_2d.f90; src/domain_2d.f90; src/init_2d.f90; TESTS/TEST_TOPOGRAPHY_FACES/\*.

Criteri di accettazione: Uguaglianza esatta delle facce condivise; geometria coerente dopo restart/scambio; pendenze/Hessiane corrette su letti noti.

Introduce B_vertex with dimensions (Nx+1,Ny+1) as the authoritative bed.

Derive one shared x-face and y-face elevation from vertex averages; derive B_cent from the same four vertices.

Replace independent limited topography-face reconstruction in HP mode with direct access to the shared Q1 face arrays.

Retain the current five-point quadratic LS fit on B_cent for Bx, By, Bxx, Bxy and Byy; use the resulting slopes to build G_center and neighbor-averaged shared G_face.

Migrate analytic/DEM initialization and subtract_init_flag so they modify the authoritative nodal geometry and then call one refresh_topography_geometry routine.

Do not permanently cache geometry in a way that prevents bed evolution.

#### A04 - Port the HP reconstruction only

Attività e verifiche operative: Ricostruzione h/eta con positività e blending; velocity interval e momento normale compatibile; velocità tangenziale 2D limitata; workset attivi senza accessi esterni ai bounds.

File interessati: src/hp_reconstruction_2d.f90; src/reconstruction_2d.f90; TESTS/TEST_HP_RECONSTRUCTION/\*.

Criteri di accettazione: Non-negatività h, test wet/dry, recovery da composizione differenziata e agreement 1D/2D con reference.

Port the clean-core direct-h candidate and positivity-preserving eta=h+B candidate using the shared Q1 face bed.

Port the parameter-free continuity-defect blend and shoreline priority.

Port compatible normal volumetric-momentum reconstruction and final velocity-interval admissibility.

Rebuild final thermodynamics after h, normal momentum, T and Y are fixed; overwrite stored physical face velocity from the final conservative face state.

In 2-D, reconstruct tangential velocity separately and form tangential momentum with the final face mass.

Keep the old flux/source operator temporarily so reconstruction differences can be isolated.

#### A05 - HP-PCCU spatial operator with G=1

Attività e verifiche operative: Operator HP-PCCU con flussi convettivi e percorso idrostatico; mantenere separazione fra parte alla faccia e percorso interno; evitare doppio conteggio della sorgente idrostatica; percorso termico nullo.

File interessati: src/spatial_operator_2d.f90; src/pccu_2d.f90; src/nonconservative_2d.f90; src/hyperbolic_2d.f90; TESTS/TEST_PCCU/\*; TESTS/TEST_HYDROSTATIC_PATH/\*.

Criteri di accettazione: Lake-at-rest roundoff in G=1; test dam break/escarpment e semidiscrete operator; test di riferimento direttamente sull’operatore spaziale.

Implement the inertial conservative flux, hydrostatic path and oriented PCCU left/right pair from the clean baseline.

Connect face and cell-internal path contributions to the production semidiscrete residual.

Disable the legacy hydrostatic source contribution in the HP path so it is not double-counted.

Keep the fourth thermal equation purely advective in the retained model: hydrostatic path component 4 is zero.

Use the final reconstructed thermodynamic endpoints consistently for Gamma, wave speeds and path data.

#### A06 - Validated slope correction

Attività e verifiche operative: G_center dai gradienti LS del letto Q1; G_face unico per faccia; velocità d’onda e quadratura del percorso coerenti con G; preservare lake at rest.

File interessati: src/geometry_2d.f90; src/pccu_2d.f90; src/hyperbolic_2d.f90; TESTS/TEST_PCCU_LAKE_AT_REST/\*.

Criteri di accettazione: Test G=1 invariato, lake at rest roundoff con forti pendii, residui/CFL confrontabili con referenza.

Compute G_center = 1/(1+Bx<sup>2+By</sup>2) from the retained LS slopes of the current B_cent.

Use one shared G_face at each internal interface.

Use c = sqrt(G_face*gprime*h) consistently in interface wave speeds.

Use the validated G-weighted hydrostatic one-form; explicit dG terms cancel in the combined pressure/bed contribution.

For the cell-internal path, vary G linearly between the two shared face values and integrate with the same quadrature as the other path variables.

Do not let the LS derivative model redefine the shared Q1 bed elevation.

#### A07 - Curvature as a separate velocity-dependent source

Attività e verifiche operative: Sorgente di curvatura esterna al percorso idrostatico: a_curv=Bxx\*u²+2\*Bxy\*u\*v+Byy\*v², con fattore slope coerente; campionare il termine misto Bxy uv.

File interessati: src/equation_terms_2d.f90; src/geometry_2d.f90; TESTS/TEST_CURVATURE_SOURCE/\*.

Criteri di accettazione: Sorgente nulla a riposo, test misti 2D, confronto algebraico e non regressione G=1.

Use the retained LS Hessian Bxx, Bxy, Byy.

Implement a_curv = Bxx*u^2 + 2*Bxy*u*v + Byy\*v^2 and the existing model-form curvature source.

Keep curvature outside the hydrostatic path.

Verify the source vanishes identically at lake at rest and that mixed Bxy/uv terms are exercised in a dynamic 2-D test.

#### A08 - Full two-dimensional hydrodynamic validation

A/B/C/D SODDISFATTE nei test definiti; resta solo A09 (implementazione delle due opzioni e decisione documentata).

Attivita gia completate: verifiche 1D/2D, composizione, equilibrio, dry/wet, topografia Q1 inclinata/ruvida, Gate-H full-field, bilanci, simmetria, riproducibilita OpenMP e restart. Commit: e4333a5c, aed31540, 691a3e21. Il pacchetto A08-C/D registra 138 casi dinamici, 24 confronti Gate-H e audit 41/41. Resta il solo contratto A09 descritto di seguito.

File interessati: TESTS/TEST_PCCU_GATE_H/\*; TESTS/TEST_PCCU_INCLINED_EXCAVATION/\*; TESTS/TEST_SPATIAL_OPERATOR/\*; TESTS/ACCEPTANCE/\*; docs/acceptance/N7_N8_CLOSURE_2026_10_10.md.

Criteri di accettazione: tutte le verifiche A08-A/B/C/D restano chiuse. A08 diventa soddisfatta quando i due criteri CFL sono disponibili con default retrocompatibile, con test direzionali e combinati, e il contratto di scelta/stabilita IMEX e documentato e approvato. Nessun altro test e aggiunto senza una decisione esplicita di modifica del perimetro.

Run the complete 1-D and 2-D reference hierarchy with G and curvature combinations.

Include the four-sided continuous one-cell inclined excavation, smooth/rough continuous beds, random continuous Q1 beds, wet/dry fronts and variable composition.

Check semidiscrete residuals, CFL, finite-time solutions, mass/component conservation, symmetry and absence of uphill leakage in the validated continuous-wall representation.

Keep deterministic and stochastic restart regressions active around the new hydrodynamic operator.

#### A09 — scelta funzionale IMPLEMENTATA e VALIDATA (8fc540eba791)

Formula directional_min: dt = min(MAX_DT, CFL\*dx/lambda_x se lambda_x\>0, CFL\*dy/lambda_y se lambda_y\>0). Formula combined: dt = min(MAX_DT, CFL/(lambda_x/dx + lambda_y/dy)) quando il denominatore e positivo. Per la direzione inattiva o senza propagazione, il contributo e zero; se entrambi nulli si restituisce MAX_DT. Per CFL=-1 permane il sentinel che restituisce MAX_DT senza valutare velocita di faccia. Conservare il cap MAX_DT e il trattamento dei casi 1D x/y.

Implementazione circoscritta: src/parameters_2d.f90 (parametro/default), src/inpout_2d.f90 (lettura, validazione e documentazione input), src/spatial_operator_2d.f90 (combinazione dei limiti), templates/examples e TESTS (unitari/integrazione). La scelta deve rifiutare stringhe non valide con errore chiaro. I due metodi vanno applicati alle medesime velocita caratteristiche ricostruite; evitare di duplicare i kernel numerici. Correggere il percorso di validazione end-to-end del sentinel CFL=-1: oggi l'operatore spaziale lo riconosce, ma la lettura input rifiuta i CFL negativi; permettere precisamente -1 e continuare a rifiutare altri valori negativi. Distinguere il passo suggerito dall'operatore dal timestep EFFETTIVAMENTE ACCETTATO dall'integratore: restano applicabili MAX_DT, times di output, tempo finale, limitatore di crescita e gli altri vincoli temporali attivi. Non modificare implicitamente limiti massimi CFL o retry/adattivita.

Gate di accettazione: SUPERATO. Report docs/acceptance/CFL_POLICY_COMPLETION_2026_10_11.md e risultati results_cfl_policy_2026_10_11.json: test 1D/2D, entrambe le politiche, CFL=-1, MAX_DT, N_RK 2/3/4, thread 1/4, strict/debug e optimized. Non affermare stabilita IMEX generale ne equivalenza universale con il riferimento Python autononomo.

#### A10 - Vertex-first evolving topography and mass exchange

A/A10-B/D-A10-AREA completate. Non riaprire automaticamente A10 per modifiche successive.

Attività e verifiche operative: Exchange limitato per cella prima dell’assemblaggio nodale; topo_rate_cell unico; gather Q1 A_C/4 verso Bdot_vertex, contributor dry zero nel denominatore; restart esatto B_vertex.

File interessati: src/mass_exchange_2d.f90; src/geometry_2d.f90; TESTS/TEST_MASS_EXCHANGE/\*; docs/acceptance/N8_A_CELL_EXCHANGE_2026_10_10.md; docs/acceptance/N8_B_PROJECTION_2026_10_10.md.

Criteri di accettazione: Nessun inventario negativo; global bed-volume exchange a roundoff; geometria ricostruita coerente; regressioni con restart e OpenMP.

Keep erosion/deposition physics finite-volume and cell local. Use the cell-centered physical state as a flat one-sided state at each of the cell four vertices; evaluate the constitutive law once per cell.

Apply solid availability, erodible inventory, continuous-phase loss, source masks and other admissibility limits before any nodal averaging.

Update q, deposit, erosion and erodible from that same cell-local limited source, exactly preserving local inventory semantics.

Store the final limited scalar topo_rate_cell for each physical cell.

At every vertex assemble Bdot_vertex with Q1 mass-lumped weights A_C/4. On a uniform interior grid this is the arithmetic mean of the four adjacent cell proposals; dry cells contribute zero but remain in the denominator.

Update only B_vertex and then refresh face/center/LS/G geometry once. Remove the independent B_cent += dt\*topo_term update in this mode.

Retain diagnostics of the local geometric mismatch while requiring exact global bed-volume exchange to floating-point roundoff.

#### A11 - Production cutover and numerical baseline freeze

A08-CFL gia completati. La baseline A11 precede B01 e il refactoring.

Gate A11: CHIUSO per lo scope autorizzato. I casi ETNA 200 s ed ETNA_LAVA 3600 s passano la matrice full/split/restart in strict/debug e optimized e 1/4 thread. Per ETNA_LAVA il bilancio di massa e valutato solo fino a 3180 s (prima del contatto al bordo evidenziato nei salvataggi), mentre confronti campo/restart continuano a 3600 s. D15 approvera separatamente la nuova architettura rispetto a questa baseline; non riaprire A11.

Attivita e verifiche operative: COMPLETE per il perimetro A11 approvato, con checksum, profili, input, output, policy dei comparatori e limite del budget lava archiviati. Eventuali nuove campagne saranno regressioni dei cambiamenti successivi e non una ripetizione obbligatoria del gate scientifico A11.

File interessati: src/spatial_operator_2d.f90; TESTS/ACCEPTANCE/\*; docs/acceptance/\*; EXAMPLES/\*; README.md.

Criteri di accettazione: SODDISFATTI per il perimetro scientifico autorizzato e la revisione congelata. D15 resta un gate distinto, successivo a relocation, reologie e solver IMEX per gruppi, e deve trattare separatamente equivalenza strutturale e cambiamenti numerici autorizzati.

Make HP-PCCU the canonical production spatial operator after A02-A10 gates pass.

Nota storica, non lavoro corrente: gli output KT appartengono al porting e a confronti legacy archiviati. HP-PCCU e l'unico backend di produzione; non pianificare il mantenimento di percorsi KT eseguibili nel solver.

Update examples, default inputs, restart versioning, NetCDF/diagnostic metadata and documentation to match mass-fraction/continuous-bed semantics.

Freeze a reproducible set of Fortran regression outputs and checksums for subsequent structural/HPC work.

Re-run strict/debug and optimized builds, 1-thread/4-thread OpenMP checks and deterministic/stochastic restart regressions.

14\. Sviluppi scientifici futuri preservati nel piano

14\. Sviluppi scientifici futuri preservati nel piano

Questi sviluppi sono indipendenti dai nove passi prioritari del refactoring reologico/IMEX. Sono mantenuti come milestone scientifiche esplicite e con proprio gate. La loro descrizione di dettaglio non era presente nei documenti roadmap del 18/24 settembre disponibili al controllo: specifiche quantitative di training, modello multilayer e closure sub-grid richiedono un successivo documento di progettazione scientifica. Non si sostituiscono tali scelte con assunzioni arbitrarie.

S0 — Progettazione delle interfacce di estensione prima del freeze G03

S0 — Progettazione delle interfacce di estensione prima del freeze G03

necessario per la compatibilità architetturale, non per l’avvio del path-only move

**Attività previste:** Definire contratti di input/output per vertical profiles, closure ML e sub-grid; variabili depth-averaged richieste, condizioni di ammissibilità, fallback deterministico, geometria temporale, capacità REAL/COMPLEX e derivazione della source; definire la distinzione fra chiusure fisiche e infrastruttura di inferenza. Evitare di incorporare tensori ML direttamente in time_integration_2d o spatial_operator_2d.

**File interessati / nuovi:** docs/ARCHITECTURE.md; docs/CODE_MAP.md; src/model_layout_2d.f90; nuove interfacce in src/vertical/, src/subgrid/, src/ml/ solo se giustificate; TESTS/ nuovi contract tests.

**Criteri di accettazione:** Interfacce documentate, dipendenze senza cicli, nessun cambiamento dei risultati standard, fallback esplicito e indipendenza dal backend/numero di layer.

**Dipendenze e collocazione:** Dopo A11, B01, reologie/IMEX; prima del gate G03. Implementazione scientifica S1–S4 successiva e indipendente dal primo refactoring.

S1 — Chiusure per profili verticali e strategie basate su ML

S1 — Chiusure per profili verticali e strategie basate su ML

**Attività previste:** Formalizzare ricostruzione verticale di velocità, concentrazione solida polidispersa, temperatura e delle variabili pertinenti al modello; distinguere profilo analitico/fisico, closure data-driven e inferenza ML. Imporre coerenza degli integrali in profondità e dei flussi (masse, quantità di moto, energia), limiti fisici, disponibilità di training/reference 3D, trattamento di regimi e fuori distribuzione. Scegliere separatamente predizione dei profili e correzioni dei termini depth-averaged, evitando doppio conteggio.

**File interessati / nuovi:** src/vertical/vertical_closure_2d.f90 \[FUTURE\]; vertical_profiles_2d.f90; vertical_quadrature_2d.f90; src/ml/ml_inference_2d.f90 e ml_feature_scaling_2d.f90 \[FUTURE\]; src/physics/equation_terms_2d.f90; TESTS/TEST_VERTICAL_CLOSURES/\* \[FUTURE\]; assets ML versionati fuori src.

**Criteri di accettazione:** Profili coerenti con medie depth-averaged a quadratura controllata, positività/frazioni, limiti wet/dry, consistenza di flussi, serie di benchmark indipendenti, accuratezza in distribuzione e robustezza OOD, prova seriale/OpenMP/MPI e fallback senza ML.

**Dipendenze e collocazione:** S0 e baseline fisica A11; può iniziare scientificamente durante G01–H04 ma non blocca MPI. Interfaccia layer-agnostic; algoritmo operativo dopo validazione dedicata.

S2 — Topografia ad alta risoluzione e parametrizzazione sub-grid

S2 — Topografia ad alta risoluzione e parametrizzazione sub-grid

**Attività previste:** Definire rapporto DEM fine/griglia di calcolo; estrarre statistiche di microtopografia, distribuzione delle pendenze/asperità, connettività e anisotropia; sviluppare correzioni fisicamente interpretabili su drag/frizione, propagazione e/o conservazione compatibili con HP-PCCU. Separare quota Q1 risolta B_vertex dalle statistiche irrisolte e specificare aggiornamento quando erosione/deposizione modifica il fondo. Verificare impatto su CFL, bilanci, nonconservatività e lake-at-rest; distinguere offline preprocessing e campi dinamici.

**File interessati / nuovi:** src/subgrid/terrain_statistics_2d.f90 \[FUTURE\]; subgrid_topography_2d.f90 \[FUTURE\]; subgrid_closure_2d.f90 \[FUTURE\]; src/geometry/geometry_2d.f90; src/spatial/spatial_operator_2d.f90; src/physics/mass_exchange_2d.f90; src/io/\*; TESTS/TEST_SUBGRID_TOPOGRAPHY/\* \[FUTURE\].

**Criteri di accettazione:** Confronti sistematici risoluzione fine/coarse, casi limite topografia uniforme (correzione nulla), lake-at-rest/positività conservati, scambi di massa e restart riproducibili, dipendenza da orientazione e scala documentata; gestione esplicita della topografia evolutiva.

**Dipendenze e collocazione:** S0, A11 e geometria Q1, dopo gate di compatibilità. Può essere studiata in parallelo a MPI; esecuzione completa non è prerequisito per H01–H04.

S3 — Fisica multilayer realmente attiva

S3 — Fisica multilayer realmente attiva

**Attività previste:** Progettare equazioni per ogni layer, flussi verticali, sorgenti interfaccia, chiusure reologiche e mass/heat/momentum exchange; definire topografia condivisa e accoppiamento alla layer inferiore, positività e wet/dry. Distinguere discretizzazione in layer da profili verticali diagnostici/ML; validare che il limite N_LAYERS=1 restituisca il depth-averaged standard. La sola dimensione rank-4 non abilita N_LAYERS\>1.

**File interessati / nuovi:** src/layers/layer_layout_2d.f90 \[FUTURE\], layer_exchange_2d.f90, multilayer_sources_2d.f90; src/core/model_layout_2d.f90; src/state/state_2d.f90; src/numerics/time_integration_2d.f90; src/physics/equation_terms_2d.f90; src/io/\*; TESTS/TEST_MULTILAYER/\* \[FUTURE\].

**Criteri di accettazione:** Conservazione massa/energia/momento fra layer, robustezza per layer sottili/asciutti, riduzione esatta o entro tolleranza al caso a un layer, campi e restart dimensionati correttamente, parallel scaling senza halo verticali inconsistenti.

**Dipendenze e collocazione:** Richiede G01/G02 e S0; fisica da approvare separatamente dopo freeze delle API G03; non confondere con attuale milestone G01.

S4 — Validazione dati, ML e workflow riproducibile

S4 — Validazione dati, ML e workflow riproducibile

**Attività previste:** Catalogare dataset training/validation/test, simulazioni 3D o esperimenti e provenance di DEM; versionare pesi ML, scaler, set di feature e versioni schema; costruire benchmark predizione vs simulazione risolta; misurare costi, fallback e comportamento a parametri non visti; mettere dataset pesanti fuori da src e dalla build standard. Separare addestramento (pipeline Python e ambiente dipendenze) da inferenza Fortran in produzione.

**File interessati / nuovi:** ml/ \[root, FUTURE\] oppure tools/ml/ per training; src/ml/ per backend inferenza; docs/models/\*; docs/acceptance/\*; TESTS/TEST_ML_INFERENCE/\* \[FUTURE\]; EXAMPLES/ dedicati.

**Criteri di accettazione:** Replay deterministico e metadati versionati, benchmarks indipendenti, assenza di dipendenza ML nella build standard, corretta diagnostica di modelli incompatibili, comportamenti per missing model/fallback documentati.

**Dipendenze e collocazione:** Supporta S1 e le eventuali varianti ML di S2; milestone indipendente dall’implementazione MPI e da H05 tuning.

Ulteriori funzionalità da preservare

Ulteriori funzionalità da preservare

Restano esplicitamente nel piano anche: sorgenti fessurali e accoppiamenti ingresso/uscita; entrainment/erosion/deposition e letto Q1 dinamico; stochastic OU, reproducible RNG e restart; modello termico e pore pressure; stato polidisperso e frazioni massiche; possibilità di nuovi termini nonconservativi tramite infrastruttura generica; I/O/versioned restart, GIS/DEM, NetCDF seriale e parallelo; OpenMP, local-domain/halos, tiles, MPI, GPU readiness; test scientifici e mantenimento REAL/COMPLEX per complex-step. Queste sono funzionalità già implementate o obiettivi documentati, non nuovi algoritmi ML dichiarati esistenti.

### A09 — Due policy CFL

Stato ufficiale: COMPLETATO

### A10 — Topografia evolutiva e scambio di massa

Stato ufficiale: COMPLETATO

### A11 — Cutover e baseline scientifica

Stato ufficiale: APPROVATO

#### HP-PCCU scientific gate

Mass-fraction q \<-\> qp round trip and closure/admissibility tests.

Exact shared Q1 face/center geometry identities.

HP reconstruction positivity, cell-mean and velocity-interval diagnostics.

G=1 lake-at-rest residual at roundoff.

Slope-corrected lake-at-rest at roundoff for smooth and one-cell 1-D/2-D topographies.

Curvature source zero at rest and analytical agreement in moving 2-D tests.

Semidiscrete residual/CFL comparison with the clean Python reference.

1-D/2-D wet/dry, rough/random continuous topography and excavation regressions.

Component-mass conservation and signed flux-admissibility tests.

Evolving-topography availability limits and global bed-volume identity at roundoff.

Exact B_vertex restart when evolving topography is enabled.

## Fase B — Build e consolidamento dell’operatore x/y

### B01 — Build, CI e portabilita

Stato ufficiale: COMPLETATO

commit 05989f5e7a del 2026-10-11. Verifiche Linux GNU Fortran optimized/strict-debug (1,4 thread), Docker optimized/strict-debug (1,4 thread) e Intel oneAPI ifx 2025.1 optimized/diagnostic (1,4 thread). Evidenze e scope: docs/acceptance/M12C_BUILD_SYSTEM_PROGRESS_2026_10_11.md e M12C_REPRODUCIBILITY_MANIFEST_2026_10_11.json.

Attivita eseguite: compilatore configurabile tramite FC, FCFLAGS utente preservati, flag OpenMP/NetCDF/LAPACK distinti, dipendenze moduli compatibili con build parallela/out-of-tree, make check collegato ai nove test di integrazione, suite di restart e Gate-H in CI; Docker e Intel ifx testati. Le correzioni diagnostiche durante la campagna ifx comprendono guardie sugli array allocabili della namelist opzionale e lettura DEM con implied DO senza array temporanei. Non interpretare questi interventi come refactoring numerico autorizzato.

**File interessati:** configure.ac; src/Makefile.am; Makefile.am; Makefile.in generati; TESTS/run_tests.sh; TESTS/\*/run_test.sh; .github/workflows/\*; Dockerfile; README.md; docs/acceptance/\*; template.doxy.

Criteri di accettazione: SODDISFATTI nel perimetro documentato. GitHub Actions GNU (run 38154444911), Docker (run 38158816015), Intel ifx (run 38168602274); test make check, restart binario/runout, restart stocastico e Gate-H. Il manifest limita il supporto a GNU e ifx su Linux; Docker usa ubuntu:24.04 senza digest immutabile e dipendenze non interamente pinnate. I run storici verificati si riferiscono a SHA specifici, non dichiarano che ogni futuro commit sia automaticamente verificato.

#### Verifiche trasversali B01 (riepilogo; requisiti completi nella scheda B01)

No active batch-blocking READ(*,*) path.

Typed LAPACK wrappers remain the only solver interface to SGESV/DGESV.

Compiler selection is honored; compile and link flags are separated.

Strict/debug and optimized/release modes are documented and automated.

make check or equivalent runs focused tests; CI executes representative solver/restart/HP regressions.

Current warning baseline is measured and correctness/portability warnings are fixed or justified.

### B02 — Mappa dipendenze e relocation path-only

Stato ufficiale: COMPLETATO

### B03 — Stocastico OU, RNG e restart seriale

Stato ufficiale: IMPLEMENTATO

**Attività previste:** Preservare modello OU, flussi RNG indipendenti, seed riproducibili, memoria latente nelle celle asciutte e continuità del restart; mantenere test esistenti durante refactoring e definire policy di correlazione tra rank in H03.

**File interessati:** src/stochastic_module.f90; src/stochastic_random_2d.f90; src/runtime_2d.f90; src/inpout_2d.f90; TESTS/TEST_STOCHASTIC_EVOLUTION/\*; TESTS/TEST_STOCHASTIC_RESTART/\*; TESTS/TEST_STOCHASTIC_RANDOM/\*.

**Criteri di accettazione:** Test seriale/OpenMP e restart stocastico invariati; verifica statistica della correlazione; nessun consumo RNG introdotto dal nuovo solver o cambio di numerazione delle celle.

### B04 — Model layout, N_LAYERS=1

Stato ufficiale: IMPLEMENTATO

**Attività previste:** Conservare model_layout_type quale unica autorità su n_layers, variabili ed indici; rifiutare N_LAYERS\>1 finché la fisica multilayer non esiste; usare metadati per G01 anziché crearne una seconda copia.

**File interessati:** src/model_layout_2d.f90; src/parameters_2d.f90; src/init_2d.f90; src/simulation_2d.f90; TESTS/TEST_MODEL_LAYOUT/\*.

**Criteri di accettazione:** Layout N_LAYERS=1 coerente con dimensioni e indici; rifiuto esplicito configurazioni non supportate; nessuna variazione dei risultati.

### B05 — Velocita caratteristiche comuni

Stato ufficiale: IMPLEMENTATO

Obiettivo e vincolo: condividere il codice che realizza lo stesso algoritmo nelle due direzioni cartesiane, mantenendo separati storage orientati, cicli sulle interfacce e ownership OpenMP. La direction determina solo l’orientamento; composizione, layer, closure verticale, reologia e sub-grid restano gestiti da interfacce fisiche distinte. Vietato introdurre array temporanei dinamici per faccia, transpose obbligatori o dipendenze circolari.

File interessati: src/equation_terms_2d.f90; src/hyperbolic_2d.f90.

Attività: Sostituire eval_local_speeds_x/y con un kernel direzionale condiviso, selezione idx_u/idx_v, stessa chiamata mixt_var e identico ramo per gravità ridotta negativa. Valutare piccoli wrapper compatibili temporanei; non modificare la formula delle velocità caratteristiche.

Criteri di accettazione: Stesse velocità min/max per entrambi gli assi, anche con dry/wet e gravità ridotta non positiva; confronto unitario su x/y e regressioni CFL/HP.

Evidenza: implementato nel commit 00400001. La suite TEST_SPATIAL_OPERATOR è passata in profilo strict e optimized con verifiche CFL direzionali e team OpenMP effettivi da 1/4 thread; TEST_N7_GEOMETRY_DYNAMIC è passata con 69 casi per profilo, inclusi stati dry/wet e casi x/y ruotati. Queste regressioni supportano lo step; il confronto esplicito isolato delle velocità per il caso di gravità ridotta non positiva resta da includere nella matrice B12.

### B06 — Ricostruzione HP sulle linee

Stato ufficiale: IMPLEMENTATO

File interessati: src/reconstruction_2d.f90; src/hp_reconstruction_2d.f90.

Attività: Consolidare evaluate_hp_cell_x/y in un kernel che prepari gli stencil con direzione esplicita e chiami reconstruct_hp_line invariata. Conservare scelta delle facce W/E/S/N, normal momentum, mask e hp_scratch thread-private; preservare la barriera prima del commit.

Criteri di accettazione: Tracce HP, mean momentum, positività, quiet-water e controlli x↔y invariati; output deterministici 1/4 thread bitwise dove previsto.

Stato aggiornato: implementato nel commit 00400001. TEST_HP_RECONSTRUCTION passa i controlli di quiet-water, positività, media del momento e indipendenza direzionale; TEST_N7_GEOMETRY_DYNAMIC passa con 69 casi per profilo e output binari identici tra 1 e 4 thread. La precedente dicitura «DA IMPLEMENTARE» descriveva uno stato superato e non è applicabile.

### B07 — Flusso PCCU di interfaccia

Stato ufficiale: IMPLEMENTATO

File interessati: src/hyperbolic_2d.f90; src/pccu_2d.f90.

Attività: Estrarre la sequenza comune eval_inertial_flux, central-upwind, limit_component_mass_flux, exact-rest scalar guard, eval_hydrostatic_coefficient, regularize_dry_gamma_pair, eval_hydrostatic_path ed eval_oriented_pccu_pair. Lasciare separati i due DO OpenMP e le matrici di faccia.

Criteri di accettazione: Flussi, path e residual identici rispetto alla baseline per gli stessi input; conservare l’ordine delle operazioni numeriche e il guard di riposo esatto introdotto nei commit recenti.

Stato aggiornato: kernel PCCU comune implementato nel commit f2065730. In entrambi i profili strict e optimized, i payload spaziali completi sono risultati bitwise identici alla baseline pre-refactoring per dynamic_rough, dynamic_drainage e lake_ramp, a 1 e 4 thread. TEST_PCCU, TEST_SPATIAL_OPERATOR e TEST_N7_GEOMETRY_DYNAMIC (69 casi per profilo) sono passati. Il test storico TEST_HYDROSTATIC_ROUNDOFF non è compatibile con l’interfaccia post-XY-1 perché il suo modulo di riferimento importa i wrapper x/y rimossi; l’equivalenza è stata verificata con un confronto diretto HEAD/candidato. B12 resta aperto per l’accettazione complessiva.

### B08 — Velocita faccia e path idrostatici di cella

Stato ufficiale: IMPLEMENTATO, GATE APERTO

File interessati: src/hyperbolic_2d.f90.

Attività: Centralizzare i due stati di faccia, riduzione min/max, percorso idrostatico dentro cella, convenzioni di orientamento e segni. Evitare la fusione dei loop OpenMP e rispettare accessi allineati a memoria.

Criteri di accettazione: Equivalenza x/y su casi 1D ruotati, 2D e composizione variabile; assenza di differenze su CFL direzionali.

Stato aggiornato: velocità di faccia e cell hydrostatic path consolidati nel commit fd43e4ca. I payload spaziali completi sono risultati bitwise identici alla baseline pre-XY-4 per dynamic_rough, dynamic_drainage e lake_ramp, in strict e optimized, con 1/4 thread. TEST_N7_GEOMETRY_DYNAMIC è passata (69 casi per profilo); TEST_N7_COMPOSITION è passata (162 casi per profilo, composizione variabile e confronti a 1/4 thread); TEST_SPATIAL_OPERATOR ha passato controlli di tracce, residui e CFL direzionale, inclusi casi x/y ruotati. B08 rimane IMPLEMENTATO, GATE APERTO fino all’accettazione formale in B12.

### B09 — Condizioni al contorno e ricostruzione lineare

Stato ufficiale: IMPLEMENTATO, GATE APERTO

File interessati: src/spatial/reconstruction_2d.f90; src/geometry/geometry_2d.f90; src/physics/equation_terms_2d.f90.

Attività: Estrarre logica riutilizzabile dei casi Dirichlet/Neumann/interpolazione per lato senza cancellare semantica dei boundary source, radial source o lateral source; nessuna modifica simultanea allo schema.

Criteri di accettazione: Tests mirati W/E/S/N e sorgenti, incluso nx=1/ny=1; output identici sulle configurazioni supportate.

Stato aggiornato al 15 ottobre 2026: in `reconstruction_2d.f90` è stato estratto `boundary_reconstruction_slope`, kernel comune per le pendenze Dirichlet, Neumann e secante sui lati W/E/S/N. Il test in `TEST_IMEX_STAGES` verifica in strict e optimized i valori analitici delle tracce per le tre modalità su ciascun lato; verifica lateral-source W/E, riflessioni radial-source W/E/S/N e griglie 5x5, 1x5, 5x1 e 1x1. Per tutti i casi, i payload completi di ricostruzione sono bitwise identici alla baseline pre-XY-5 (sorgente `reconstruction_2d.f90` del commit c648e9c5356c65d68d0624369d408cec781f0816, SHA-256 90d5d50923069798cd6135f18a940cdaa423a7d26644044ca9989e483d6bfd43), sia in strict sia in optimized, a 1 e 4 thread; i payload sono inoltre identici fra i due conteggi di thread. `TEST_HP_RECONSTRUCTION`, `TEST_SPATIAL_OPERATOR` e `TEST_IMEX_STAGES` sono passati. B09 è IMPLEMENTATO; il gate formale di accettazione complessiva resta in B12.

| Revisione               | Modifica / evidenza                                  | Interpretazione per roadmap                                                |
|-------------------------|------------------------------------------------------|----------------------------------------------------------------------------|
| 752b47fe85              | Doxygen comments in Fortran sources                  | Documentazione dei sorgenti: implementata; qualità da mantenere.           |
| 793edd3482 / 05989f5e7a | N9 approval / M12c accepted                          | Baseline N9 pinned 3f7329c85d; M12c chiusa nel suo scope.                  |
| 2bb0efff87              | Relocate Fortran sources into functional directories | Relocation path-only completata; matrice solver/Docker/ifx verde.          |
| 44329ad06e              | Standardize EXAMPLES README.md                       | 20 casi disponibili per il catalogo; confronto integrale storico pendente. |
| ec7a11324a              | MkDocs, Doxygen and auto Examples publishing         | Sito pubblicato tramite workflow di documentazione.                        |
| 1d28ca6bd4              | Documentation theme contrast fix                     | Palette Ocean Blue & Teal aggiornata; publishing passata.                  |

CI e pubblicazione: relocation: https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/38215124368 (solver), /38215124389 (Docker), /38215124380 (ifx). Documentazione: /38220579536 e /38221400795 (riuscite). Le run solver per gli ultimi commit del sito risultavano ancora in corso al momento della ricognizione: non attribuire loro un esito prima della conclusione.

### B10 — Derivate e halo geometrici

Stato ufficiale: IMPLEMENTATO, GATE APERTO

File interessati: src/geometry/geometry_2d.f90.

Attività: Valutare line kernels per derivate prime e seconde e utility di copia bordi senza variare ordine di accumulazione, formule Q1/LS/Hessian, trattamento delle dimensioni degeneri o mixed derivative. Non precomputare definitivamente la geometria: topografia futura time-dependent.

Criteri di accettazione: TEST_GEOMETRY_1D per x e y e fixture 2D, ricostruzione dopo aggiornamenti B_vertex, confronto 1/4 thread e output di topografia invariati.

Stato aggiornato al 15 ottobre 2026: le copie dei bordi dei campi derivati sono state estratte in utility private x/y; i loop 2D, le formule e le somme miste non sono stati modificati. L'estrazione del kernel aritmetico 1D è stata scartata perché il confronto bitwise ha rilevato differenze in optimized; mantenere i loop scalari originali conserva l'ordine numerico. `TEST_GEOMETRY_1D` è passato in strict e optimized (112 casi per profilo, fixture 1D x/y e 2D, cache avvelenate e snapshot identici a 1/4 thread); `TEST_TOPOGRAPHY_FACES` è passato in entrambi i profili. Il confronto A/B con il sorgente pre-B10 ha confermato payload identici bitwise su griglie 5x1, 9x1, 17x1, 5x5 e 7x9, tutti i quattro profili di topografia e slope-correction OFF/ON, a 1 e 4 thread, in strict e optimized. B10 è IMPLEMENTATO; il gate formale resta in B12.

### B11 — Riduzione CFL direzionale

Stato ufficiale: IMPLEMENTATO, GATE APERTO

File interessati: src/spatial/spatial_operator_2d.f90.

Attivita: dopo A11, fattorizzare ove utile la raccolta dei massimi direzionali senza alterare i criteri scelti da CFL_POLICY. La funzionalita duale directional_min/combined e un lavoro circoscritto del gate A09, PRIMA di A11; B11 e soltanto refactoring interno successivo.

Criteri di accettazione: stesso risultato numerico per ciascuna delle due policy prima e dopo il refactoring, conservando comportamento CFL=-1, MAX_DT, griglie 1D x/y e ordine di operazioni appropriato. Nessuna variazione di default.

Stato aggiornato al 15 ottobre 2026: la raccolta dei quattro segnali di velocità di faccia per cella è stata estratta in una utility privata; le due riduzioni `MAX` mantengono lo stesso ordine di operandi e la riduzione OpenMP, mentre le formule `directional_min`/`combined`, il default e i loro rami 1D non cambiano. `TEST_SPATIAL_OPERATOR` è passato in strict e optimized: verifica le due formule, griglie x/y/point, propagazione nulla, cap `MAX_DT`, sentinel `CFL=-1` e team effettivi 1/4. `TEST_CFL_POLICY` è passato con solver strict e optimized, 18 traiettorie per profilo a 4 thread (N_RK=2,3,4; policy omessa, directional_min, combined; CFL adattivo e sentinel), confrontando output canonici, checkpoint e tempi accettati; sono passati anche i tre controlli di input non valido. B11 è IMPLEMENTATO; il gate complessivo resta in B12.

### B12 — Accettazione complessiva x/y

Stato ufficiale: GATE

#### Sequenziamento e disciplina di regressione

Gate di ingresso: A11 \[APPROVATO\] e B01 \[ACCETTATO\]. Prima eseguire la relocation path-only con manifest sorgenti e dipendenze, senza modifica di istruzioni o ordine floating-point. Successivamente sviluppare in commit separati le fattorizzazioni x/y e la riorganizzazione reologica/IMEX; eseguire CI GNU, Docker e Intel, oltre ai confronti A11 pertinenti. Non confondere la riuscita dei test di compilazione con l'accettazione scientifica D15.

## Fase C — Output e distribuzione

### C01 — Output segregati per run

Stato ufficiale: DA FARE

#### Struttura di directory concordata

> \<OUTPUT_BASE_DIRECTORY\>/\<RUN_NAME\>/
>
> metadata/ effective input backup, output manifest
>
> fields/ RUN_NAME.nc (NetCDF multidimensionale)
>
> snapshots/ \<indice\>\_t\_\<tempo\>/ per ASCII ESRI e q_2d
>
> diagnostics/ runout.csv, mass_center.csv, maximum_thickness.asc
>
> probes/ CSV delle sonde
>
> vulnerability/ RUN_NAME_VT.txt, RUN_NAME_VT\_\*.asc
>
> terrain/ dem_esri\*.asc
>
> restart/ restart.bin (nome fisso)

### C02 — Contratto di configurazione e sicurezza

Stato ufficiale: DA FARE

#### Contratto di configurazione e sicurezza

RUN_NAME già esistente identifica univocamente il run; validare caratteri ammessi ed escludere nomi vuoti, separatori, “.”, “..” e traversals.

OUTPUT_BASE_DIRECTORY (nuovo parametro proposto) è configurabile, di default OUTPUT; sono ammessi percorsi relativi o assoluti, inclusi filesystem scratch HPC.

OUTPUT_OVERWRITE (nuovo parametro proposto) è FALSE per default: un nuovo run non deve sovrascrivere output di una directory esistente. Un’autorizzazione TRUE è esplicita e va confinata rigorosamente al run selezionato, con difese contro symlink, percorsi malformati o riferimenti esterni.

Un restart ordinario riutilizza la directory dello stesso run, senza richiedere overwrite, ma verifica compatibilità checkpoint e integrità degli output preesistenti. La lettura esplicita di checkpoint esterni resta supportata.

Nome del checkpoint di scrittura FISSO: restart.bin in \<OUTPUT_BASE_DIRECTORY\>/\<RUN_NAME\>/restart/restart.bin. Nessun parametro utente per rinominarlo. Il chiamante può continuare a passare il nome letterale alla routine di I/O; il path viene costruito centralmente.

Snapshot ASCII organizzati in sottodirectory con indice univoco e tempo fisico leggibile; registrare il tempo esatto in un manifest. Il NetCDF resta un unico file di serie temporale in fields/, non duplicato negli snapshot.

Vulnerability è categoria autonoma di primo livello, NON una sottodirectory di diagnostics/. Preservare inizialmente i nomi dei file prodotti dal solver.

Tutti i file scritti direttamente dal solver, inclusi backup di input, DEM esportati, raster, diagnostiche, serie, massimi e checkpoint, devono rimanere sotto il run. File di input originali e prodotti di script esterni non vengono spostati automaticamente.

Centralizzare la costruzione/validazione dei percorsi di output in un componente I/O dedicato. Definire eventuale modalità legacy temporanea solo per il periodo di validazione e migrazione degli script di post-processing.

### C03 — Accettazione I/O per run

Stato ufficiale: GATE

#### Test, documentazione e gate di accettazione

Inventario delle operazioni di creazione/scrittura file in tutto il codice e classificazione completa degli output; test che nessun output solver sia creato nella working directory.

Test multi-run, directory già esistente, overwrite esplicito, percorso personalizzato, validazione RUN_NAME e protezione da cancellazioni fuori dal run.

Test restart v5 in place, da checkpoint esterno, append NetCDF/CSV, interruzione durante snapshot e ripartenza senza sovrascrivere risultati validi.

Confronto quantitativo/bitwise, quando applicabile, con baseline A11 per verificare che la modifica dei path non alteri risultati numerici né tempi di output.

Aggiornamento di EXAMPLES e script in UTILS, regressioni e documentazione GitHub Pages: nuova pagina “Output files and formats”, naming, manifest, semantica overwrite e restart.

Solo dopo accettazione della riorganizzazione dell’I/O procedere con politica Autotools B, make distcheck e distribuzione di archivi distinti codice/examples dallo stesso tag.

### C04 — Autotools policy B

Stato ufficiale: DA FARE

#### Milestone post-XY — Politica Autotools B e distribuzione con due archivi

Sequenziamento: iniziare soltanto dopo il completamento e l’accettazione del refactoring x/y, della riorganizzazione dell’I/O e della relocation dei sorgenti nelle directory funzionali. L’attività è di packaging/build e non deve essere confusa con l’accettazione scientifica D15. Stato al 12 ottobre 2026: decisione di indirizzo presa; implementazione e test non ancora eseguiti.

Politica Autotools B: in Git versionare configure.ac, Makefile.am, src/Makefile.am e le eventuali macro sorgente; non versionare configure, Makefile.in (anche in src/), aclocal.m4, depcomp, install-sh, missing e altri output rigenerabili. Aggiornare .gitignore e le istruzioni per compilare un checkout tramite autoreconf -fi; assicurare compatibilità GNU, ifx e Docker.

Due artefatti della medesima release/tag, con versione e commit coerenti: (1) IMEX_SfloW2D-\<version\>.tar.gz, archivio sorgente ufficiale creato da Autotools, contenente Fortran, include, configure e file ausiliari generati, documentazione minima, licenza e test necessari; (2) IMEX_SfloW2D-\<version\>-examples.tar.gz, contenente EXAMPLES/ con input, script, DEM/dati necessari e tutti i README.md standardizzati, mantenendo i percorsi relativi. I nomi sono contratti proposti, con versione derivata dal tag e non fissata in anticipo.

Audit preliminare obbligatorio: verificare make distcheck e i contenuti EXTRA_DIST, l’inclusione di tutti i .inc e dei sorgenti nelle nuove sottodirectory, l’eseguibilità dei test da un archivio estratto e l’assenza di dipendenze implicite dal checkout Git. Controllare la dimensione e le licenze dei dati Examples, gli archivi DEM, le dipendenze tra esempi e il percorso dell’eseguibile dopo installazione. Non dichiarare gli archivi autonomi finché questi controlli non passano.

Automazione release: creare un workflow GitHub Actions che generi e testi entrambi gli archivi dalla stessa revisione, confronti un manifest dei file e checksum, poi li alleghi alla release ufficiale soltanto dopo esito positivo. Chiarire nel README e nel sito che i download automatici «Source code (zip/tar.gz)» di GitHub NON equivalgono all’archivio Autotools preparato e possono non contenere configure.

Gate di accettazione: (a) build/test da clone Git con autoreconf -fi; (b) configure + make + make check senza Autotools dall’archivio sorgente; (c) make distcheck; (d) riproducibilità degli Examples estratti con l’archivio del codice corrispondente, senza dipendenze non dichiarate; (e) validazione su GNU/ifx/Docker secondo matrice supportata; (f) documentazione aggiornata su GitHub Pages e CONTRIBUTING.md. Non eliminare i vecchi README.txt prima dell’audit integrale già pianificato.

### C05 — Doppio archivio di distribuzione

Stato ufficiale: DA FARE

### C06 — Validazione distributiva

Stato ufficiale: GATE

#### Regola di rilascio per documentazione pubblica

Pubblicazione unificata: solo il workflow publish-documentation.yml deve distribuire contenuti su GitHub Pages; non creare un secondo workflow Doxygen che possa sostituire la root MkDocs. Eseguire manualmente o da tag/release secondo il workflow tracciato. La homepage resta a /, il catalogo in /examples/ e l’HTML Doxygen in /doxygen/. La documentazione deve essere costruita dalla medesima revisione del sorgente e il catalogo legge direttamente EXAMPLES/\*\*/README.md. Testare tema chiaro, tema scuro, stampa e risoluzione dei link quando si modifica il CSS.

## Fase D — Reologie e solver IMEX

### D01 — Inventario formule e rami REAL/COMPLEX

Stato ufficiale: PARZIALE

#### M12e — scheda operativa

**Attività previste:** Caratterizzare sorgenti per ciascuna reologia, temperatura e pressione interstiziale su stati REAL/COMPLEX; confrontare Jacobiano complex-step e finite difference controllata; coprire U→0, h→0, soglie, ABS/MAX e branch non olomorfi; conservare formule condivise e test modello-per-modello prima dell’estrazione.

**File interessati:** src/equation_terms_2d.f90; src/state_conversion_2d.f90; src/complexify.f90; src/state/\*.inc dopo A2; src/rheology/\*.f90 futuri; TESTS/TEST_RHEOLOGY_FROUDE/\* e nuovi TEST_RHEOLOGY\_\*.

**Criteri di accettazione:** Parità numerica reale/complessa per punti regolari; differenze delle derivate spiegate nei punti non lisci; nessuna regressione scientifica alla migrazione delle formule.

#### Cosa rimane da verificare prima delle modifiche

Ci sono ancora tre aspetti che considero opportuno chiarire direttamente dal codice:

Identificazione del modello reologico 7: prima di assegnargli un nome e un file, occorre stabilire esattamente quale formulazione rappresenti e se sia effettivamente supportato.

Classificazione delle reologie analitiche: dobbiamo verificare le formule disponibili per ogni modello e distinguere le integrazioni dell’ODE dalle soluzioni algebriche dello stage DIRK.

Dipendenze delle sorgenti: durante l’estrazione delle reologie occorre preservare l’unica implementazione delle formule REAL/COMPLEX, necessaria per il complex-step.

Questi sono approfondimenti implementativi, non motivi per cambiare la struttura proposta.

### D02 — Caratterizzazione numerica cell-local

Stato ufficiale: DA FARE

### D03 — Audit e ownership di parametri, variabili e stato condiviso

Stato ufficiale: DA FARE

#### Obiettivo e posizione nella sequenza

Definire un owner univoco per ogni dichiarazione di modulo in constitutive_parameters_2d.f90 e parameters_2d.f90 prima dei trasferimenti. Conservare i risultati della caratterizzazione REAL/COMPLEX D01-D02. Non duplicare il dato nei moduli consumatori: owner significa il modulo responsabile della definizione e del ciclo di vita, non ogni utilizzatore.

#### Inventario richiesto e deliverable

Produrre un manifest machine-readable con una riga PER SIMBOLO (non solo per famiglia): dichiarazione e tipo/kind/rank, inizializzazione da namelist o derivazione, tutti gli USE/ONLY e riferimenti indiretti, lettori e scrittori, eventuale mutabilità per time-step/stage, local/shared/allocated lifetime, restart/I-O, OpenMP/MPI, owner attuale, owner proposto e gate richiesto. La matrice seguente e preliminare: i consumatori/scrittori definitivi richiedono ricerca code-wide, non sono dedotti dalla sola dichiarazione.

#### Regole vincolanti

Distinguere input immutabili, grandezze derivate, stato evolutivo, diagnostiche e workspace. Evitare namespace globali generici, duplicazioni e dipendenze circolari; usare PUBLIC selettivi, USE ONLY, configurazioni validate e interfacce esplicite. Non modificare formule, algoritmi, ordine floating-point, serializzazione, semantica COMPLEX o stato persistente solo per spostare simboli. Preservare compatibilità namelist, restart, configurazione N_LAYERS=1 e test/CI. Ogni MODULE, elemento e procedura modificati sono soggetti al contratto Doxygen completo.

#### Criterio di accettazione

Manifest completo e revisionato, dependency graph aciclico, piano di migrazione per dominio fisico, baseline registrate, nessuna proprietà senza owner, nessun workspace globale non motivato. Questo gate precede D05 (estrazione reologie), ma la rimozione isolata della reologia 8 in D04 puo procedere senza attendere la migrazione fisica se le regressioni sono preservate.

| Origine      | Variabili / famiglie                                                                                                                                      | Owner candidato                                 | Verifiche prima del trasferimento                                         |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------|---------------------------------------------------------------------------|
| constitutive | mu, xi, xi_temp, mu_0, mu_inf, Fr_0, U_w                                                                                                                  | Rheology / each law                             | xi_temp: verify local workspace vs shared configuration; stochastic paths |
| constitutive | I_0, mu_2, mu_s, muI_inf, friction_factor, tau                                                                                                            | Rheology / each law                             | Cross-model parameters and read/write/validation paths                    |
| constitutive | T_ref, nu_ref, visc_par, tau0, emme, enne                                                                                                                 | Lava rheology / thermophysical definition       | Avoid duplicating shared thermal parameters                               |
| constitutive | alpha1_coeff, beta1, alpha2, beta2, Kappa, n_td, n_td2                                                                                                    | Lahar rheology                                  | Initialization/order of derived n_td2                                     |
| constitutive | SBconst, emissivity, convective_term_coeff, radiative_term_coeff, atm_heat_transf_coeff, exp_area_fract                                                   | Thermal exchange                                | Separate environmental configuration from derived rates                   |
| constitutive | T_env, T_ground, thermal_conductivity, c_p                                                                                                                | Thermal / ambient or phase-property owner       | Resolve whether lava-only or general temperature                          |
| constitutive | grav, inv_grav, pres, inv_pres, T_ambient, rho_a_amb                                                                                                      | Physical constants / ambient conditions         | Immutable runtime configuration vs derived reciprocals                    |
| constitutive | rho_l, inv_rho_l, kin_visc_l, sp_heat_l, gamma_steam                                                                                                      | Liquid / vaporization                           | Check vaporization and composition clients                                |
| constitutive | sp_heat_a, sp_gas_const_a, kin_visc_a, sp_heat_g, sp_gas_const_g                                                                                          | Gas and phase properties                        | Allocated gas vectors, n_add_gas and restart                              |
| constitutive | sp_heat_c, rho_c_sub, kin_visc_c, Tref_Suth, muRef_Suth, S_mu                                                                                             | Carrier gas / liquid / Sutherland               | Differentiate raw phase properties vs model-derived                       |
| constitutive | rho_s, inv_rho_s, diam_s, sphericity_s, sp_heat_s                                                                                                         | Particle/solid properties                       | Allocation/lifetime, inverse density, PSD initialization                  |
| constitutive | A_drag, B_drag, collective_settling_flag, settling_flag, alphastot_min, vonK                                                                              | Particle drag / settling                        | Drag may also belong to common closure; validate consumers                |
| constitutive | entrainment_flag, erosion_coeff, loss_rate, erodible_fract, erodible_porosity, coeff_porosity, T_erodible                                                 | Mass exchange / erodible substrate              | Source and deposition coupling; ownership of allocated data               |
| constitutive | hydraulic_permeability, dynamic_permeability_flag, maximum_solid_packing, alpha_trans, N_inh, f_inhibit_mode                                              | Pore pressure / inhibition                      | Complex-step semantics and branch/regularization                          |
| constitutive | pascal_coeff(:), pascal_coeff_precomputed, precompute_pascal_coefficient                                                                                  | Pore-pressure inhibition coefficients           | Derived precompute state, thread safety and restart                       |
| parameters   | sp, dp, wp, xinf, four_thirds, neg_four_thirds, dry_thickness_tolerance                                                                                   | Numerical kinds and numerical constants         | Retain low-level dependency direction; no cycles                          |
| parameters   | dt0, max_dt, cfl, cfl_policy, n_RK, t_start, t_end                                                                                                        | Time integrator configuration                   | Temporal config vs mutable runtime clock                                  |
| parameters   | max_nl_iter, tol_abs, tol_rel, eps_sing, eps_sing4                                                                                                        | Nonlinear solver / regularization               | Confirm which solvers and spatial kernels consume                         |
| parameters   | reconstr_coeff, limiter(:), theta                                                                                                                         | Spatial reconstruction configuration            | No numerically visible changes during relocation                          |
| parameters   | n_vars, n_eqns, n_layers, n_solid, n_add_gas, n_stoch_vars, n_pore_vars, idx\_\*                                                                          | Model layout / equation metadata                | Index initialization, physical state and no circular USE                  |
| parameters   | rheology_flag, rheology_model, liquid_flag, gas_flag, pore_pressure_flag, sutherland_flag                                                                 | Model physics configuration / rheology selector | Rheology selection migration is separate D08                              |
| parameters   | slope_correction_flag, curvature_term_flag, topo_change_flag, water_level                                                                                 | Geometry / spatial physics configuration        | Time-dependent geometry remains updateable                                |
| parameters   | stochastic_flag, stoch_transport_flag, length_spatial_corr                                                                                                | Stochastic model configuration                  | RNG/OU state and restart distinct from input settings                     |
| parameters   | x_source, y_source, r_source, r2_source, azimuth_source, arc_width_source, mfr_source, vel_source, T_source, h_source, alphas_source(:), alphag_source(:) | Source configuration (radial/initial)           | Inspect activation flags, units and source initialization                 |
| parameters   | lateral_source_flag, source_side, x1_source, x2_source, y1_source, y2_source, bottom\_\*\_source_flag, n_fissures, \*\_fissures                           | Source configuration (lateral/fissural)         | Time histories, fixed arrays and source validation                        |
| parameters   | collapsing_volume_flag, x_collapse, y_collapse, r_collapse, T_collapse, h_collapse, alphas_collapse(:)                                                    | Initial/collapse source configuration           | Shared initialization and startup state                                   |
| parameters   | particle_distribution_flag, fractal_dim, diam_min, diam_max, rho_s_min, rho_s_max, sigmoid_k, sigmoid_d0                                                  | Particle-size distribution input                | Separate PSD parameters from derived per-class arrays                     |
| parameters   | bc type, bcW(:), bcE(:), bcS(:), bcN(:)                                                                                                                   | Boundary-condition descriptors                  | Ownership and allocation/restart semantics                                |
| parameters   | t_output, dt_output, t_runout, t_steady, t_probes, verbose_level, n\_\*\_levels, thickness_levels(:), dyn_pres_levels(:)                                  | Runtime scheduler / output and diagnostics      | Separate event scheduling from output preferences                         |

### D04 — Eliminazione reologia 8

Stato ufficiale: DA FARE

Questa decisione è definitiva.

Il modello 8 è un percorso sperimentale che integra analiticamente l’ODE della frizione quadratica, anziché risolvere l’equazione algebrica dello stage DIRK.

La verifica del repository non ha individuato nessun input di esempio o test che selezioni RHEOLOGY_MODEL = 8.

Si dovranno eliminare:

Il ramo speciale in nonlinear_solver_2d.f90.

La procedura integrate_friction_term in equation_terms_2d.f90.

Il trattamento speciale in time_integration_2d.f90.

La selezione del modello 8 nella lettura e validazione dell’input.

Il ramo specifico del modello 8 per shearVel in inpout_2d.f90.

Eventuali altri riferimenti esclusivamente associati a questo modello.

La futura soluzione analitica della frizione quadratica apparterrà a rheology_quadratic.f90 e risolverà correttamente lo stage DIRK.

Non verrà mantenuto un metodo speciale di integrazione dell’ODE esterno alla nuova architettura.

#### Caso \#8 ed identificatori testuali

Rimuovere rheology_model==8 da nonlinear_solver_2d, time_integration_2d, inpout_2d (incluso ramo shearVel) e l’unica routine integrate_friction_term; eliminare dipendenze e documentazione *solo* del modello 8. Il modello 6 conserva la propria fisica e ottiene un’opzione algebrica ANALYTIC. Tutti gli identificatori esterni sono stringhe corrispondenti al suffisso file rheology\_\<nome\>.f90; none non richiede un modulo. La selezione va risolta una volta in rheology_2d.f90, senza SELECT CASE su stringa nelle celle. Modello 7 da auditare prima di assegnare nome e file. Migrare anche gli script generatori e respingere esplicitamente i vecchi input numerici.

Decisione approvata: eliminare completamente il modello reologico 8, mantenere il modello 6 per la legge Bursik–Woods; la soluzione algebrica ANALYTIC dello stage DIRK verrà introdotta separatamente in D13. L’implementazione nel repository non è ancora stata rimossa.

### D05 — Estrazione moduli reologici

Prerequisito: audit D03. Durante l’estrazione dei moduli reologici, trasferire agli owner di legge mu/xi, parametri lava, lahar, plastic, quadratic, Froude, Lucas e mu(I); non trascinare proprietà fisiche non reologiche in un generico modulo reologico. Registrare read/write e mantenere la numerica REAL/COMPLEX.

Stato ufficiale: DA FARE

Questa è la parte operativa: cosa deve cambiare rispetto all’implementazione attuale.

| File                           | Modifiche richieste                                                                                                                                     |
|--------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------|
| parameters_2d.f90              | Cambiare rheology_model da intero a stringa; aggiornare le definizioni e le opzioni di configurazione.                                                  |
| inpout_2d.f90                  | Leggere i nuovi nomi, eliminare i confronti numerici, delegare la validazione dei parametri reologici; aggiornare gli output che dipendono dal modello. |
| init_2d.f90                    | Inizializzare il modello reologico, la partizione implicita e i gruppi.                                                                                 |
| constitutive_parameters_2d.f90 | Conservare i parametri fisici; riorganizzarli soltanto se necessario per eliminare dipendenze improprie.                                                |
| equation_metadata_2d.f90       | Aggiungere IM_group, IM_group_method, n_IM_groups, mapping delle equazioni e controlli di validità.                                                     |
| time_integration_2d.f90        | Introdurre EX_terms e IM_terms; spostare fuori dall’integratore la sequenza SI+Newton e le correzioni reologiche.                                       |
| nonlinear_solver_2d.f90        | Implementare il coordinamento dei gruppi, Newton con Jacobiano ridotto e la selezione ANALYTIC/NEWTON.                                                  |
| equation_terms_2d.f90          | Estrarre le formule specifiche delle reologie, della sovrapressione e del raffreddamento; mantenere l’assemblaggio generale delle sorgenti.             |
| rheology_2d.f90                | Nuovo dispatcher e procedura comune correct_friction_momentum.                                                                                          |
| rheology\_\*.f90               | Implementare le singole leggi fisiche, i loro contributi SI/IM e le formule chiuse disponibili.                                                         |
| pore_pressure_2d.f90           | Concentrare le sorgenti di sovrapressione interstiziale e la relativa soluzione analitica.                                                              |
| thermal_sources_2d.f90         | Concentrare le sorgenti termiche e le formule analitiche supportate.                                                                                    |
| simulation_2d.f90              | Aggiornare workspace e inizializzazione dei nuovi metadati.                                                                                             |
| Makefile.am                    | Aggiornare sorgenti e dipendenze Fortran; rigenerare i file Automake derivati.                                                                          |

I moduli spaziali HP-PCCU, le conversioni di stato e la geometria dovrebbero richiedere prevalentemente aggiornamenti dei percorsi e delle dipendenze, non modifiche degli algoritmi.

#### Struttura prevista della nuova versione

I file esistenti vengono spostati senza modificarne inutilmente i nomi. I nuovi moduli sono indicati con \[NEW\].

src/  
├── IMEX_SfloW2D.f90  
├── Makefile.am \[Automake\]  
├── core/ simulation_2d, runtime_2d, init_2d, parameters_2d,  
│ constitutive_parameters_2d, model_layout_2d,  
│ diagnostics_2d, complexify  
├── numerics/ time_integration_2d, nonlinear_solver_2d,  
│ equation_metadata_2d, lapack_interfaces_2d  
├── spatial/ spatial_operator_2d, hyperbolic_2d,  
│ hp_reconstruction_2d, reconstruction_2d,  
│ pccu_2d, nonconservative_2d  
├── state/ state_2d, state_conversion_2d, associated \*.inc  
├── geometry/ geometry_2d, domain_2d  
├── physics/ equation_terms_2d, mass_exchange_2d,  
│ pore_pressure_2d \[NEW\], thermal_sources_2d \[NEW\]  
├── rheology/ rheology_2d \[NEW\], rheology_voellmy \[NEW\],  
│ rheology_plastic \[NEW\], rheology_lava \[NEW\],  
│ rheology_lahar \[NEW\], rheology_empirical \[NEW\],  
│ rheology_quadratic \[NEW\], rheology_froude \[NEW\],  
│ rheology_lucas \[NEW\], rheology_mui \[NEW\],  
│ rheology_mui_voellmy \[NEW\] \[\*.f90\]  
├── stochastic/ stochastic_module, stochastic_random_2d  
├── io/ inpout_2d; future input, restart,  
│ serial/parallel NetCDF and diagnostics modules  
├── layers/ layer_layout_2d, layer_exchange_2d, \[FUTURE\]  
│ multilayer_sources_2d \[FUTURE\]  
├── vertical/ vertical_closure_2d, vertical_profiles_2d,  
│ vertical_quadrature_2d \[FUTURE\]  
├── subgrid/ subgrid_topography_2d, terrain_statistics_2d,  
│ subgrid_closure_2d \[FUTURE\]  
├── ml/ ml_inference_2d, ml_model_registry_2d,  
│ ml_feature_scaling_2d \[FUTURE\]  
├── parallel/ local_domain_2d, halo_exchange_2d,  
│ tiles_2d, mpi_backend_2d \[FUTURE\]  
└── coupling/ external_forcing_2d, data_interfaces_2d  
\[FUTURE/OPTIONAL\]

La struttura prevede un solo livello di directory sotto src/. Non introdurrei ulteriori sottodirectory per i metodi analitici, gli algoritmi Newton o le singole sorgenti.

Le directory EXAMPLES/, TESTS/ e docs/ rimangono nella posizione attuale, aggiornandone i contenuti quando necessario.

Compatibilità con gli sviluppi scientifici futuri

Compatibilità con gli sviluppi scientifici futuri

Le directory layers/, vertical/, subgrid/, ml/, parallel/ e coupling/ costituiscono confini di responsabilità previsti, NON moduli attualmente implementati. Durante il refactoring iniziale si possono predisporre API o documentazione senza creare file vuoti. I nomi elencati sono proposte di collocazione e non affermazioni sull’esistenza del codice.

La dinamica risolta resta nel solver fisico; vertical/ calcola profili o quadrature condizionati dallo stato depth-averaged; ml/ fornisce un backend opzionale con versionamento di pesi e feature; subgrid/ fornisce correzioni o statistiche di terreno senza modificare l’autorità di B_vertex. layers/ gestisce le interazioni verticali solo quando il modello fisico multilayer sarà approvato. geometry/ è proprietario delle quote e derivate risolte, spatial/ dell’operatore HP-PCCU, io/ dei contratti di restart e NetCDF, parallel/ degli halo; nessuna chiamata MPI o I/O nei kernel reologici/ML.

Le API proposte devono accettare geometria e stato cell-local, essere testabili senza allocazioni per cella, supportare N_LAYERS=1 come baseline ed essere esplicitamente disattivabili per la regressione. La descrizione dell’interfaccia viene fissata prima di G03, non l’implementazione scientifica completa degli algoritmi futuri.

#### Struttura architetturale finale

time_integration_2d

Stage RK, spatial_terms, EX_terms, IM_terms

nonlinear_solver_2d

Gruppi impliciti, metodo ANALYTIC/NEWTON, Jacobiano ridotto

equation_terms_2d

Assemblaggio delle sorgenti fisiche

Rheology

rheology_2d e modelli individuali

Thermal

thermal_sources_2d

Pore pressure

pore_pressure_2d

Schema delle responsabilità e delle interfacce concettuali, non una sequenza obbligatoria di chiamate per la soluzione analitica.

La caratteristica fondamentale della nuova struttura sarà la separazione fra il modello fisico, il trattamento numerico delle sorgenti e l’integratore IMEX.

In particolare, aggiungere una nuova reologia dovrà richiedere principalmente un nuovo file rheology\_\<nome\>.f90 e la sua registrazione nel dispatcher, senza modificare il ciclo RK o l’algoritmo Newton.

Questo è il piano consolidato che utilizzerei come riferimento per l’implementazione e per il successivo aggiornamento della Development Roadmap.

### D06 — Sorgenti termiche e pore pressure

Stato ufficiale: DA FARE

### D07 — Migrazione progressiva dei parametri non reologici

Stato ufficiale: DA FARE

#### Ambito e dipendenze

Eseguire per famiglie, DOPO il gate D03 e coordinandosi con D05/D06 (reologie e sorgenti). Non spostare ripetutamente le variabili da un contenitore globale a un altro. La riallocazione deve portare ownership semantica, API esplicite e minor accoppiamento verificato; E01 restera un audit architetturale indipendente post-refactoring.

#### Sequenza operativa

\(1\) Parametri ambientali, proprieta di gas/liquidi/solidi e loro reciprocali; (2) trasporto/deposizione e mass exchange; (3) pressione interstiziale e coefficienti di inibizione; (4) configurazione sorgenti, condizioni al contorno, output e runtime da parameters_2d; (5) metadati layout, solver e geometria. La migrazione termica va sincronizzata con D06 e non duplicata. Ogni tranche deve identificare gli USE, le inizializzazioni, le allocazioni e la proprieta dello stato.

#### Policy dello stato

Parametri letti da namelist convalidati una volta; grandezze derivate calcolate al momento corretto e aggiornabili se la fisica richiede variabilita; tempo corrente, diagnostica e workspace di stage non devono diventare configurazioni globali. Non assumere geometria/topografia statica: il DEM puo evolvere nel tempo. Gli array allocabili devono avere ownership esplicita, senza copie divergenti.

#### Acceptance per tranche

Compilazione GNU e ifx, OpenMP 1/4 thread, test di input/restart, confronto numerico strict/optimized e casi fisici mirati. Registrare symbol manifest before/after, dependency graph, confronto di default e prestazioni; ogni modifica alla semantica delle fonti/risultati e un task scientifico separato. Aggiornare Doxygen per ogni simbolo e procedura alterati.

### D08 — RHEOLOGY_MODEL da numeri a stringhe

Stato ufficiale: DA FARE

#### Identificazione delle reologie

La convenzione definitiva sarà:

La stringa specificata nell’input corrisponde al suffisso del file sorgente.

| Modello numerico attuale | Nuovo nome                  | File                     |
|--------------------------|-----------------------------|--------------------------|
| 0                        | none                        | Nessuna reologia         |
| 1                        | voellmy                     | rheology_voellmy.f90     |
| 2                        | plastic                     | rheology_plastic.f90     |
| 3                        | lava                        | rheology_lava.f90        |
| 4                        | lahar                       | rheology_lahar.f90       |
| 5                        | empirical                   | rheology_empirical.f90   |
| 6                        | quadratic                   | rheology_quadratic.f90   |
| 7                        | Da verificare nel dettaglio | Da definire              |
| 8                        | Eliminato                   | Nessuno                  |
| 9                        | froude                      | rheology_froude.f90      |
| 10                       | lucas                       | rheology_lucas.f90       |
| 11                       | mui                         | rheology_mui.f90         |
| 12                       | mui_voellmy                 | rheology_mui_voellmy.f90 |

La tabella costituisce la mappatura di progetto. I nomi descrittivi e, in particolare, il modello 7 devono essere verificati rispetto alle formulazioni realmente implementate prima della migrazione definitiva.

Ad esempio, l’attuale input:

&RHEOLOGY_PARAMETERS  
RHEOLOGY_MODEL = 6  
FRICTION_FACTOR = 1.0D-3  
/

diventerà:

&RHEOLOGY_PARAMETERS  
RHEOLOGY_MODEL = 'quadratic'  
FRICTION_FACTOR = 1.0D-3  
/

Il dispatcher comune sarà in rheology_2d.f90. La selezione della reologia avverrà all’inizializzazione; nei cicli numerici non dovranno essere eseguiti ripetutamente confronti tra stringhe.

#### Modifiche agli input, esempi e test

La trasformazione da numeri a stringhe richiede di aggiornare tutti i file presenti in EXAMPLES/, TESTS/ e nella documentazione.

In particolare:

Sostituire tutti i valori numerici di RHEOLOGY_MODEL con i nomi testuali.

Aggiornare gli script che generano automaticamente gli input.

Verificare che gli output e i messaggi diagnostici riportino il nome del modello.

Introdurre test per modelli sconosciuti e configurazioni dei gruppi non valide.

Aggiungere test specifici per ANALYTIC, NEWTON e gruppi misti.

Verificare il comportamento con numero diverso di incognite e in esecuzione OpenMP.

Per i vecchi file di input numerici prevederei un errore esplicito con indicazione del nuovo formato, anziché conservare due sistemi paralleli di identificazione.

### D09 — Correzione del solo attrito

Stato ufficiale: DA FARE

La nuova procedura:

correct_friction_momentum(...)

sarà contenuta in:

src/rheology/rheology_2d.f90

Sostituirà i due blocchi di correzione attualmente presenti dopo SI e dopo Newton.

La correzione deve preservare la direzione e non aumentare il modulo della quantita di moto SOLO per il sottopasso in cui agisce esclusivamente la frizione. La stessa proprieta non si applica allo stato complessivo di uno stage DIRK o di un gruppo in cui siano presenti gravita, termini termici, pore pressure o altre sorgenti che possono ruotare o accelerare il momento.

Pertanto, per un aggiornamento esclusivamente reologico:

La correzione non deve essere applicata a contributi spaziali, gravitativi o ad altre forze che possono modificare la direzione del moto.

Inoltre, occorre eliminare l’incoerenza attuale fra OR e AND nei controlli di inversione.

La scelta precisa dell’algoritmo di proiezione e arresto sarà validata mediante test specifici, evitando modifiche numeriche non controllate durante il semplice spostamento del codice.

#### Correzione esclusivamente reologica

correct_friction_momentum vive in rheology_2d.f90 ed opera solo sul contributo della frizione isolato: m_after = alpha\*m_before, 0 \<= alpha \<= 1, con arresto ammissibile e nessuna rotazione/inversione. L'interfaccia deve rendere espliciti momentum_before_friction e momentum_after_friction e impedire di applicare questo vincolo allo stato globale derivante da altre sorgenti o dal solve accoppiato. SI, NEWTON e ANALYTIC possono usare il correttore esclusivamente nel sottopasso di frizione definito; i gruppi misti richiedono una prova separata del contributo isolato.

### D10 — EX_terms / IM_terms

Stato ufficiale: DA FARE

L’integratore dovrà conoscere soltanto:

spatial_terms  
EX_terms  
IM_terms

L’attuale separazione NH e SI_NH non sarà più esposta all’integratore, poiché i due contributi utilizzano lo stesso tableau implicito.

Lo stage sarà scritto concettualmente come:

Il risolutore restituirà lo stato aggiornato e il contributo implicito effettivo:

quando $a_{ii} \neq 0$.

Questo contributo deve comprendere anche le eventuali correzioni apportate durante la soluzione.

Il trattamento SI rimane disponibile internamente quando richiesto dalle reologie, ma non deve più essere gestito direttamente dal ciclo RK.

Il caso predefinito con un unico gruppo NEWTON dovrà riprodurre il comportamento SI+Newton attuale.

#### SI legacy, correzioni e incremento implicito effettivo

Nel default un unico gruppo NEWTON deve riprodurre in sequenza la parte SI, lo stato intermedio q_si, la soluzione Newton e le correzioni attuali, salvo la successiva correzione scientifica esplicitamente testata del criterio di arresto. Non confondere una sorgente fisica valutata dopo clipping con il contributo effettivo RK. Dopo le correzioni ricostruire:

Con tau=0, lo stage non effettua un solve implicito: evitare divisioni per zero e definire il contributo/stage coerentemente con il tableau e la sua necessità nel passo. I tableau IMEX e la distinzione sorgenti esplicite/operatore spaziale non vengono modificati.

### D11 — IM_group e metadati

Stato ufficiale: DA FARE

#### Metadati

Manteniamo:

implicit(:)  
implicit_map(:)  
n_implicit

e aggiungiamo:

IM_group(:)  
IM_group_method(:)  
n_IM_groups

IM_group(i) identifica il gruppo dell’equazione i; il valore zero significa che l’equazione non viene trattata implicitamente.

IM_group_method(g) identifica il metodo utilizzato per il gruppo g.

Le due possibilità sono:

ANALYTIC  
NEWTON

Non introduciamo un terzo metodo SI_IM.

La partizione implicita globale resta distinta dalla dimensione dei singoli sistemi Newton.

#### Esempio

Supponiamo di avere:

IM_group = \[0, 3, 3, 1, 2\]  
  
IM_group_method = \[ANALYTIC, ANALYTIC, NEWTON\]

L’ordine è:

Gruppo 1 — Energia

ANALYTIC

Gruppo 2 — Sovrapressione

ANALYTIC

Gruppo 3 — Momento x/y

NEWTON, Jacobiano 2×2

Ogni gruppo utilizza i valori più aggiornati disponibili delle variabili esterne. Non introduciamo ulteriori maschere di lagging.

### D12 — Newton per gruppo e Jacobiano ridotto

Stato ufficiale: DA FARE

#### Jacobiano Newton

Attualmente la dimensione dipende da n_implicit.

Nella nuova versione dipenderà invece dal singolo gruppo:

Il Jacobiano sarà:

Le derivate continueranno a essere calcolate mediante complex-step, perturbando esclusivamente le incognite del gruppo.

Le mappe degli indici saranno costruite all’inizializzazione, evitando allocazioni ripetute durante la soluzione delle celle.

#### Inizializzazione e validazione dei gruppi

equation_partition_type conserva n_equations, n_explicit, n_implicit, implicit(:), implicit_map(:) e explicit_map(:). Aggiunge n_IM_groups, IM_group(:) per equazione, IM_group_method(:) per gruppo e le mappe compatte precomputate. La relazione IM_group(i)=0 se e solo se .NOT. implicit(i) è un’invariante. I gruppi devono essere numerati da 1 a n_IM_groups senza vuoti, ogni metodo deve essere valido, ogni richiesta ANALYTIC deve avere una formula chiusa disponibile per **tutte** le sorgenti implicite agenti sulle equazioni del gruppo. Il mapping globale non può essere sovrascritto con indici locali temporanei. implicit(:) deve avere una sola fonte di verità, evitando configurazioni indipendenti incoerenti.

#### Contratto algebrico dei gruppi

Sia tau = dt\*a_diag e q_known lo stato noto dopo l’assemblaggio dei termini RK dagli stage precedenti. Un gruppo g risolve soltanto le proprie componenti, sempre con il termine noto dello stage q_known(i); le sorgenti possono leggere lo stato aggiornato dai gruppi precedenti e il valore ancora noto dai gruppi successivi. Questa strategia è una singola sweep block Gauss–Seidel, non identica alla soluzione simultanea se esiste coupling intergruppo. Vanno testati ordine di convergenza temporale, stabilità e sensibilità all’ordine dei gruppi. Nessuna maschera separata di lagging. Il gruppo unico NEWTON costituisce un gate di equivalenza con il solver monolitico preesistente, con confronto di residuali, stati e convergenza agli stage secondo le tolleranze scientifiche approvate. I gruppi multipli ordinati e la singola sweep Gauss-Seidel costituiscono invece una MODIFICA NUMERICA, non un refactoring equivalente: devono superare un distinto gate di accuratezza, convergenza rispetto al coupling e dipendenza dall'ordine dei gruppi, anche quando l'esito non coincide bitwise col solve simultaneo.

Per NEWTON, costruire il residuo locale R_g(x_g) e il Jacobiano J_g = dR_g/dx_g con dimensione n_g × n_g, n_g = COUNT(IM_group==g). Complex-step perturba solo colonne del gruppo e mantiene le altre variabili fisse. Riutilizzare gli ottimizzati solver lineari 2×2/3×3 quando appropriato, calcolare gradienti e line search nello spazio locale e garantire nessuna modifica a componenti esterne. Preallocare workspace per thread, dimensionati per il massimo gruppo Newton, **non** allocare per cella.

### D13 — Algebriche ANALYTIC per stage DIRK

Stato ufficiale: DA FARE

#### Limitazioni di ANALYTIC

Per la prima implementazione:

Una sola equazione scalare, se è disponibile una formula chiusa.

Oppure le due equazioni del momento, se esiste una soluzione vettoriale chiusa.

Nessuna soluzione analitica generica per gruppi con più equazioni accoppiate.

Questi vincoli sono controllati durante l’inizializzazione.

#### Contratto ANALYTIC e soluzione quadratic drag

ANALYTIC ammette esclusivamente (a) una singola equazione scalare oppure (b) esattamente entrambe le componenti della quantità di moto, quando esiste la formula congiunta. Si risolve **l’equazione algebrica dello stage DIRK**, non l’ODE su dt. Per la frizione quadratica isotropa, con k=f/h≥0 e modulo noto U\*, la relazione è:

La radice non negativa, in forma stabile, è:

Ricostruire il momento mantenendo il verso di riferimento, con casi robusti U\*=0, k=0, tau=0. Per raffreddamento e pore-pressure la formula è disponibile soltanto se rappresenta l’insieme delle sorgenti implicite attive per l’equazione. Un’eventuale formula quartica radiativa deve essere verificata quanto a radici, selezione fisica e stabilità; non utilizzare Newton implicitamente dentro ANALYTIC.

### D14 — Gruppi misti e validazione

Stato ufficiale: DA FARE

### D15 — Acceptance della nuova architettura

Stato ufficiale: GATE

#### Manifest, ripartenza e gate

Ogni commit funzionale porta evidenza di compilazione pulita, test associati, SHA e descrizione della tolleranza. Gli input e l’output devono registrare stringa reologica e gruppo/metodo/ordine; i restart non devono ripartire silenziosamente con un diverso modello/gruppo. Le ottimizzazioni OpenMP non possono introdurre workspace Newton condivisi. Conservare i .inc di conversione REAL/COMPLEX e i controlli sulla non analiticità di ABS, MAX, segni e soglie nel complex-step. Prima di G01 eseguire il gate di compatibilità del default, regressioni scientifiche HP-PCCU, storage/restart e accuratezza per gruppi.

## Fase E — Architettura e robustezza del passo

### E01 — Audit dipendenze e ownership stato

Questo audit resta un gate INDIPENDENTE dopo D15: verificare che D03 e D07 abbiano effettivamente risolto ownership e dipendenze e che non vi siano stati globali o temporanei introdotti dai refactoring. Non usare E01 come primo inventario dei parametri.

Stato ufficiale: DA FARE

DECISIONE DI PIANIFICAZIONE: la conclusione del refactoring e il superamento dei test scientifici D15 sono condizioni necessarie, ma non sufficienti, per avviare i nuovi sviluppi. Prima di aprire nuove linee di sviluppo scientifico, deve essere completato un audit indipendente di architettura, manutenibilità ed estensibilità, con evidenze riproducibili e una decisione esplicita di accettazione (E03).

#### Obiettivo: architettura funzionale agli sviluppi futuri

L’architettura deve consentire l’aggiunta e la modifica di reologie, equazioni e gruppi impliciti, output, livelli verticali, geometria evolutiva e successive strategie di parallelismo senza ripetere modifiche pervasive nei moduli non pertinenti. “Moderno” significa responsabilità chiare, contratti verificabili, prestazioni preservate e workflow riproducibili, non l’introduzione indiscriminata di astrazioni o classi.

#### Sequenziamento e vincolo di ingresso

Relocation sorgenti (completata) → refactoring x/y → riorganizzazione I/O e OUTPUT/\<RUN_NAME\>/ → politica Autotools B e packaging con due archivi → separazione reologie e riorganizzazione solver IMEX → D15 (accettazione scientifica) → E03 (accettazione architetturale) → nuovi sviluppi scientifici pianificati (G01–H05 e altre estensioni).

Le attività preparatorie e i test esplorativi possono proseguire prima di E03, ma non devono imporre nuove dipendenze strutturali o essere assunti come architettura definitiva. L’ingresso nelle milestone funzionali successive richiede la chiusura esplicita di E03.

#### Matrice di verifica post-D15

| **Area**                | **Criterio verificabile**                                                                                                                                                                  |
|-------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Dipendenze tra moduli   | Mappa delle dipendenze USE/call, assenza di cicli non giustificati; responsabilità e direzione delle dipendenze documentate.                                                               |
| Stato e configurazione  | Inventario di variabili di modulo mutabili e side effect; stato fisico, configurazione, geometria, runtime e workspace distinti; trasferimento esplicito dove utile.                       |
| Operatori x/y           | Kernel comuni senza duplicazioni sostanziali; differenze geometriche e di bordo rese esplicite, proprietà numeriche preservate.                                                            |
| Interfacce interne      | Contratti per operatori spaziali, reologie, IMEX, I/O, stato e geometria: input/output, proprietà, ownership ed errori definiti e documentati.                                             |
| Estensibilità reologica | Aggiungere o modificare una reologia richiede interventi localizzati; nessun branching diffuso nell’integratore o in moduli non correlati.                                                 |
| Estensibilità IMEX      | Gruppi di incognite, metadata, Jacobiani e trattamenti numerici estendibili senza modificare in modo pervasivo il ciclo temporale.                                                         |
| Geometria e livelli     | Geometria/topografia ricalcolabili quando necessario (anche in futuro nel tempo); niente precomputazioni statiche assunte come vincolo; predisposizione a N_LAYERS=1 e futuro multi-layer. |
| Memoria e prestazioni   | Ownership di allocazioni e workspace; assenza di temporanei/allocazioni indesiderate negli hot path; confronto prestazioni e memoria con baseline fissata.                                 |
| Error handling          | Gestione coerente degli errori di input, numerici, I/O e restart; messaggi contestualizzati; nessuna cancellazione accidentale di risultati.                                               |
| Test e scienza          | Test unitari/operatoriali, regressioni riproducibili e benchmark scientifici indipendenti; conservazione, positività, well-balancing, 1/4 thread, restart e profili build.                 |
| Build e distribuzione   | Build GNU/ifx/Docker, make distcheck e archivi code/examples riproducibili e compatibili con README e dipendenze degli esempi.                                                             |
| Documentazione          | Doxygen, MkDocs, Examples, output formats, parametri e contratti interni coerenti con il codice; riferimenti scientifici mantenuti.                                                        |

### E02 — Architecture fitness checks

Stato ufficiale: DA FARE

#### Prove di estensione (architecture fitness checks)

Prima del via libera, svolgere almeno tre analisi di impatto su casi realistici: (1) aggiunta di una nuova reologia o variazione significativa di una esistente; (2) estensione dei gruppi impliciti/equazioni con metadata e Jacobiano; (3) introduzione di uno strato addizionale o di aggiornamenti geometrici nel tempo. In ciascun caso documentare moduli interessati, interfacce da modificare, dipendenze introdotte, necessità di riallocazione e test da estendere. Sono prove di fattibilità architetturale, non obblighi di implementare anticipatamente le nuove funzionalità.

Per il percorso verso G01–H05 verificare inoltre che l’organizzazione dello stato e della geometria non ostacoli domain decomposition, halo exchange, tile/local domain, MPI, Parallel NetCDF e futuri acceleratori. Non si richiede adesso di implementare MPI o GPU.

### E03 — Architecture readiness

Stato ufficiale: GATE

#### Evidenze, esiti e gestione delle non conformità

Produrre docs/acceptance/ARCHITECTURE_READINESS_AUDIT.md e un manifest machine-readable con revisione/commit, matrice PASS/FAIL/WAIVER, link ai workflow, risultati dei benchmark e tracciabilità dei requisiti. Documentare ogni non conformità con impatto su correttezza, estensibilità, prestazioni o manutenzione, responsabile e azione correttiva.

E03 può essere approvato soltanto con D15 già accettato, evidenze complete e nessun problema aperto ad alta criticità che renda onerosi o rischiosi gli sviluppi previsti. Eventuali eccezioni a bassa criticità devono essere motivate, circoscritte e accompagnate da test/mitigazioni. Se il gate fallisce, pianificare correzioni mirate e ripetere le verifiche; evitare una riscrittura generale senza benefici misurabili.

Criterio di uscita: “Il codice e il repository sono sufficientemente organizzati, validati, manutenibili ed estendibili da costituire una piattaforma stabile per i nuovi sviluppi scientifici”. Solo dopo tale approvazione si sblocca il percorso di nuove funzionalità.

### E04 — Policy timestep e input

Stato ufficiale: DA FARE

Introdurre enable_step_retry, retry_dt_factor (default candidato 0.5, 0\<factor\<1), max_step_retries (default candidato 10), min_time_step \[s\] configurabile e parametri di diagnostica. min_time_step è il limite assoluto di timestep: nessun tentativo può essere lanciato con Δt inferiore al minimo; se la prossima riduzione lo violasse, la simulazione termina con un errore controllato. Validare valori finiti, consistenza con il timestep iniziale e compatibilità con il controllo CFL; non accettare un timestep che ecceda il limite temporale di output/fine run per effetto del controllore. I default numerici proposti sono da finalizzare durante l’implementazione, non da presentare come soglie fisiche universali.

La policy per il passo successivo rimane quella del calcolo CFL esistente. Un retry modifica solamente Δt del tentativo corrente: non è un controller adattativo di accuratezza basato su stime embedded. Registrare accepted_dt, proposed_dt, rejected_attempts e minimum_accepted_dt nelle diagnostiche; definire i contatori in modo coerente con il restart.

### E05 — Bounds primitive e override

Stato ufficiale: DA FARE

Definire un descrittore centralizzato dei limiti per ciascuna variabile primitiva fisicamente significativa del model_layout; ciascun minimo/massimo deve poter essere attivo o disattivo in modo indipendente. L’input può sovrascrivere soltanto alcune variabili e soltanto uno dei lati dell’intervallo. La precedenza è: limiti generali ampi → override/default del modello reologico o della particolare chiusura attiva → override esplicito dell’utente. Validare i bounds al setup, inclusa l’intersezione fra vincoli di modelli concorrenti: un intervallo vuoto è un errore di configurazione, non un fallimento recuperabile del timestep.

Distinguere condizioni di dominio matematico (p.es. frazioni e proprietà termodinamiche valide) e intervalli prudenziali di controllo fisico/numerico. Consentire tolleranze solo ove compatibili con la definizione del vincolo e le unità della variabile; non usare clipping dello stato come sostituto del rigetto. I default non devono aggiungere vincoli ingiustificati a velocità, densità, pressioni o temperature. Il riferimento ambientale della legge di raffreddamento della lava può essere proposto dalla reologia/chiusura termica come default condizionato alle ipotesi fisiche, non come limite universale di ogni stato termico.

Lo sviluppo comprende un inventario completo delle primitive effettivamente esposte dal codice, con indici e significati ricavati da model_layout e state_conversion_2d; il controllo deve trattare gli stati asciutti secondo una semantica esplicita, senza segnalare falsi errori per variabili convenzionali non definite nel secco. Le componenti non finite (NaN/Inf) devono produrre un fallimento recuperabile o fatale appropriato senza corrompere le diagnostiche.

### E06 — Controlli di ammissibilita cell-local

Stato ufficiale: DA FARE

Contratto concettuale: CALL check_primitive_bounds(state, valid, failure_info). La routine non ha effetti collaterali: legge unicamente lo stato primitivo di cella e i bounds risolti, restituisce il flag e un record diagnostico con variabile/indice, (j,k), valore, limite, tipo di violazione e identificatore di stadio o finale. Non stampa, non clippa e non invoca fatal_error. Il chiamante determina il retry e la messaggistica.

Punti di chiamata ESCLUSIVI: (a) dopo il completamento di ciascuno stadio RK, prima dell’uso dello stato di stage nelle successive valutazioni esplicite/spaziali; (b) dopo l’assemblaggio e le conversioni/proiezioni finali, prima dell’accettazione del passo. Eseguire la scansione soltanto su l=1,domain%solve_cells, usando j_cent(l),k_cent(l); non effettuare controlli in ogni conversione conservative→primitive e non controllare automaticamente gli stati ricostruiti alle facce.

time_integration_2d.f90 conserva this%q_rk per ogni cella risolta dopo lo stage. this%qp_rk è aggiornato oggi solo se need_explicit_stage=.TRUE.: la conversione qc_to_qp per lo stato di stage dovrà essere resa disponibile anche negli altri stadi, senza duplicarla quando già necessaria all’operatore esplicito. qp ordinario, nel loop implicito, può contenere la conversione del noto di stage e non va assunto uguale al risultato finale; dopo l’assemblaggio conclusivo qp è invece aggiornato ed è la scelta naturale per il controllo finale. Verificare che conversioni di diagnostica in verbose non alterino le conclusioni. Conservare la semantica halo/dry esistente.

### E07 — Error propagation, rollback e retry

Stato ufficiale: DA FARE

Separare esplicitamente i fallimenti recuperabili (Newton non convergente, solve lineare fallito, stato di cella non ammissibile, NaN/Inf riconducibile al tentativo) dagli errori strutturali (input incoerente, file mancanti, allocazione e incompatibilità di configurazione/restart). Le routine coinvolte dovranno restituire status/failure_info invece di arrestare incondizionatamente il processo per le sole condizioni recuperabili. I fallimenti dei thread OpenMP vanno raccolti senza race e senza I/O concorrente non controllato; interrompere il passo solo dopo una sincronizzazione sicura, evitando sottopassi successivi su workspace parziali.

All’inizio del timestep salvare lo stato persistente necessario a ripetere il tentativo: conservative e primitive accettate, variabili dinamiche del fondo e geometria aggiornata ove modificate dal passo, maschere e liste attive, campi e memoria stocastica, stato RNG, contatori persistenti, eventuali accumulatori diagnostici e stato dei modelli con memoria. Non è obbligatorio copiare workspace ricalcolabili, ma è obbligatorio re-inizializzarli prima del retry. Il tempo, gli indici del passo accettato, i file scientifici e le scadenze di output non avanzano durante un tentativo respinto. I metadati del retry possono essere diagnosticati separatamente senza alterare il percorso di restart riproducibile.

Se un tentativo fallisce: registrare failure_info; ripristinare lo snapshot; calcolare Δt_new=factor·Δt_old; validare min_time_step e limite retry; ripetere tutti gli stadi IMEX dall’inizio. Nessuna modifica alle equazioni di bilancio o ai pesi RK. Il comportamento deve essere deterministico rispetto al tentativo effettivamente accettato, entro la policy di RNG documentata per l’attivazione stocastica.

### E08 — Test del timestep controller

Stato ufficiale: DA FARE

| **Test**                  | **Requisito di accettazione**                                                                                                           |
|---------------------------|-----------------------------------------------------------------------------------------------------------------------------------------|
| T0 - Passo valido         | Quando il primo tentativo è accettato, risultati e sequenza di output coincidono con baseline senza retry; overhead misurato.           |
| T1 - Violazione di stadio | Bounds violati dopo un RK stage rifiutano l’intero passo; nessun passaggio di tempo o output contaminato.                               |
| T2 - Violazione finale    | Controllo qp dopo assembly, ripristino e riduzione Δt; validare stage/finale separatamente.                                             |
| T3 - Newton               | Non convergenza o solve fallito propagati come errore recuperabile, senza deadlock OpenMP.                                              |
| T4 - Soglia minima        | Se Δt_new \< min_time_step, errore finale con Δt, tempo, stadio, cella, variabile e ultimo limite violato.                              |
| T5 - Riproducibilità      | Tentativo fallito + retry confrontato con run partito direttamente dallo stesso stato con Δt accettato; includere 1/4 thread e restart. |
| T6 - Fisica Etna          | Riprodurre la cella termica (77,68), T iniziale circa 859.62 K e Δt 0.67107 s; bounds coerenti con modello attivo.                      |
| T7 - Wet/dry Taal         | Masse/frazioni e fronti wet/dry; confrontare con timestep ridotto, senza clamping di massa introdotto dal controller.                   |
| T8 - Stochastic           | Assenza di consumo RNG permanente nei tentativi respinti, regressione restart e coerenza con policy M12d.                               |

Per l’implementazione: creare un gate E09 indipendente dalle ottimizzazioni delle matrici IMEX; registrare SHA, input, compiler, variabili monitorate, massimo numero di tentativi, dt minimo effettivo, differenze scientifiche e condizioni di arresto. E09 deve precedere la migrazione layer-aware e la parallelizzazione MPI, così da definire una semantica locale di rollback riutilizzabile.

### E09 — Gate timestep

Stato ufficiale: GATE

## Fase F — Interfacce per chiusure future

### F01 — Interfaccia closure verticali

Stato ufficiale: DA FARE

### F02 — Interfaccia subgrid terreno

Stato ufficiale: DA FARE

### F03 — Backend ML/NN opzionale

Stato ufficiale: DA FARE

### F04 — Test di neutralita API

Stato ufficiale: DA FARE

### F05 — Architecture closure readiness

Stato ufficiale: GATE

## Fase G — Migrazione layer-aware

### G01 — Storage rank-4 con N_LAYERS=1

G01 non parte senza (1) baseline scientifica HP-PCCU A11 approvata e congelata e (2) accettazione D15 della nuova architettura reologica/IMEX. D15 include due gate distinti: equivalenza del gruppo unico NEWTON rispetto alla baseline e accuratezza dei gruppi multipli ordinati, con effetti numerici documentati. G02 adatta restart e I/O allo stato layer-aware senza cambiare la fisica per N_LAYERS=1.

Stato ufficiale: DA FARE

Attività e verifiche operative: Migrare q/qp e workspace a (var,layer,x_storage,y_storage), N_LAYERS=1; conservare prima dimensione contigua; letto Q1 condiviso; nessuna fisica multilayer.

File interessati: src/state_2d.f90; src/domain_2d.f90; src/simulation_2d.f90; src/spatial_operator_2d.f90; src/time_integration_2d.f90; src/mass_exchange_2d.f90; TESTS/TEST_MODEL_LAYOUT/\*.

Criteri di accettazione: Regressioni A11 dentro tolleranze, memoria e shape verificate; n_layers diverso da 1 rifiutato.

Migrate q and qp to (variable,layer,x,y), keeping variable first so q(:,layer,j,k) is contiguous.

Migrate reconstruction, HP-PCCU hyperbolic/path workspace, IMEX stage/residual storage, relevant mass-exchange interfaces and stochastic transported storage where model semantics require a layer extent.

Keep common per-layer physics single-sourced with an outer layer loop; N_LAYERS=1 remains the only executable physics mode.

B_vertex and derived bed geometry remain shared geometric fields, not replicated per layer, unless a future multilayer model explicitly requires layer-specific beds.

Use the frozen A11 regression as the G01 reference.

#### Runtime layer count

N_LAYERS is already introduced structurally. G01 changes storage rank, not physics. Unsupported layer counts must continue to fail explicitly until a separate multi-layer scientific workstream is approved.

#### State layout

The variable index remains first for contiguous cell vectors. After A02, primitive composition slots mean mass fractions. The layer dimension is infrastructure; it does not reactivate alpha/h\*alpha coordinate branches.

#### Geometry layout

Only B_vertex is authoritative for the evolving bed. Face, center, derivative and G fields are derived. Their distributed storage may include halos/cache, but ownership and restart semantics must never imply independent geometric states.

### G02 — Activity/restart/I-O layer-aware

Stato ufficiale: DA FARE

Attività e verifiche operative: Maschere di attività spaziali non per-layer; restart binario e NetCDF seriale layer-aware; serializzare B_vertex autoritativo e rigenerare geometry; test multistrato sintetici senza fisica.

File interessati: src/inpout_2d.f90; src/state_2d.f90; src/geometry_2d.f90; TESTS/TEST_BINARY_RESTART/\*; TESTS/TEST_STOCHASTIC_RESTART/\*.

Criteri di accettazione: Restart identico per N_LAYERS=1, versioning/schema documentato; roundtrip sintetici e NetCDF coerente.

Define active spatial-cell semantics through model/layout metadata rather than fixed q(1) assumptions where practical.

Make binary restart and NetCDF internals layer-aware while preserving exact B_vertex restart when evolving topography is enabled.

Regenerate B_cent/face/LS/G geometry from restarted B_vertex instead of storing inconsistent independent geometry state.

Decide explicitly whether the external N_LAYERS=1 NetCDF schema preserves legacy shapes or exposes a length-one layer dimension.

Add synthetic multi-layer pack/unpack and serialization tests without enabling multi-layer physics.

Continue the incremental inpout_2d split along model-input, restart, scientific output and diagnostics boundaries.

#### Active-domain semantics

The active spatial workset remains distinct from the layer dimension and from vertex ownership. With N_LAYERS=1 it must reproduce current compact lists during G01/G02. Wet/dry activity affects mass-exchange proposals, but physical dry cells remain part of the geometric mass-lumped vertex stencil with zero proposal.

### G03 — Architecture freeze

Stato ufficiale: GATE

#### Layer migration gate

N_LAYERS omitted and N_LAYERS=1 explicitly selected reproduce the frozen HP-PCCU baseline.

Unsupported layer counts fail clearly.

All layer-dependent core state/workspaces carry a layer extent.

Layer-aware restart/serial NetCDF tests pass, including B_vertex restart.

Synthetic multi-layer storage can be packed/unpacked and serialized without enabling new physics.

Attività e verifiche operative: Congelare ownership di model_layout, equation_partition, spatial_operator, physical sources, io e topography; documentare API e assenza di dipendenze circolari.

File interessati: src/model_layout_2d.f90; src/equation_metadata_2d.f90; src/spatial_operator_2d.f90; src/equation_terms_2d.f90; docs/ARCHITECTURE.md; docs/CODE_MAP.md.

Criteri di accettazione: A11, B01, refactoring reologie/IMEX, G01/G02 e regressioni passate; API documentate e compatibili con MPI.

HP-PCCU A11 baseline and all mandatory numerical regressions are frozen.

State/workspace storage contracts are layer-aware and stable.

B_vertex ownership, restart semantics and geometry refresh API are stable.

Model/layout and equation-partition ownership is explicit and immutable after configuration.

Serial restart/NetCDF paths are layer-aware for supported modes.

Strict/debug and optimized build/test paths are automated enough to diagnose distributed regressions.

No batch-blocking interactive input remains.

## Fase H — MPI/HPC e prestazioni

### H01 — Local domain e halo abstraction

Stato ufficiale: DA FARE

Attività e verifiche operative: Introdurre local_domain (global/local/storage), halo espliciti a larghezza adeguata alla derivata LS (±2); local==global per primo; separare boundary fill fisico da exchange.

File interessati: src/domain_2d.f90; src/geometry_2d.f90; src/spatial_operator_2d.f90; src/reconstruction_2d.f90; test locali sintetici nuovi.

Criteri di accettazione: Equivalenza seriale local==global, test halo LS e topografia nodale, nessun MPI nei kernel.

Introduce global extents, local owned extents, storage extents, halo width and global/local transforms.

First instantiate local==global so all HP-PCCU kernels run through local extents without communication.

Derive halo width from the maximum active stencil reach. The retained five-point LS derivative fit requires two cell-center neighbors on each side; this must be represented explicitly rather than hidden in boundary copying.

Define nodal B_vertex storage/ownership consistently with cell ownership. A rank must be able to reconstruct its owned cell face geometry and LS derivatives from synchronized local/halo data.

Keep physical-boundary filling separate from future inter-rank halo exchange.

#### Halo-width contract

The halo width must be derived from active stencils. HP face reconstruction is local, but the retained five-point quadratic LS slope/Hessian operator uses j-2…j+2 or k-2…k+2. The local-domain design therefore needs at least the corresponding two-cell derivative reach unless the operator is reorganized with an equivalent communication/reconstruction strategy that passes the same numerical regressions.

### H02 — Gerarchia active tiles

Stato ufficiale: DA FARE

#### Local-domain and tile gate

local==global reproduces G03 exactly or to the documented tolerance.

Owned interior, cell halos and vertex halos/extents are correct for 1-D-like and fully 2-D grids.

Two-cell LS derivative reach is valid at local/rank boundaries.

Physical boundary filling and halo storage do not double-update interior state.

Tile traversal selects the same active cell/interface work and produces the same vertex bed update as the compact-list reference.

Different tile sizes produce equivalent scientific results.

Attività e verifiche operative: Tile rettangolari owned e scheduling locality; conservare compact active-list; verificare vicini cross tile e gather ai vertici.

File interessati: src/domain_2d.f90; src/hp_reconstruction_2d.f90; src/spatial_operator_2d.f90; src/mass_exchange_2d.f90; nuovi test tile.

Criteri di accettazione: Risultati identici per tile sizes multipli, gestione stencil e rispetto bounds, OpenMP stabile.

Build rectangular tile descriptors over the owned local domain while retaining compact active lists as the correctness reference.

Tile halo/stencil reach must cover HP reconstruction, two-cell LS derivative geometry and curvature/source needs.

Treat tiles only as scheduling/locality units; they must not change wet/dry activity or vertex mass-exchange semantics.

Add tile-vs-compact-list regressions before removing any old traversal path.

### H03 — Decomposizione MPI

Stato ufficiale: DA FARE

Attività e verifiche operative: Decomposizione cartesiana MPI e halo exchange generico; B_vertex ownership condivisa e gather nodale across rank; continuità correlazioni stocastiche ai rank boundaries.

File interessati: Nuovo communication/MPI module; src/domain_2d.f90; src/geometry_2d.f90; src/mass_exchange_2d.f90; src/stochastic_module.f90; nuovi test MPI.

Criteri di accettazione: Confronto np=1/2/4 con seriale, bilancio globale e restart coerenti; nessun MPI interno al solver per cella.

Create a 2-D Cartesian communicator and map ranks to local_domain ownership.

Implement generic exchange for all layer/state components without encoding physical meanings in MPI code.

Synchronize authoritative B_vertex across rank boundaries or define an equivalent unique owner/reduction protocol so shared face geometry is bitwise/sufficiently identical on both sides.

For evolving topography, ensure vertex assembly combines cell proposals from all cells sharing a vertex even when those cells are on different ranks; do not renormalize over rank-local contributors.

Exchange/recompute the data required by the LS slope/Hessian stencil before evaluating G or curvature near rank boundaries.

Use global reductions only for true global quantities such as CFL extrema, conservation diagnostics and termination conditions.

Validate np=1 through the MPI path, then 2- and 4-rank square and non-square decompositions.

### H04 — MPI+OpenMP e Parallel NetCDF

Stato ufficiale: DA FARE

#### MPI and Parallel NetCDF gate

MPI np=1 matches the serial/local-domain HP-PCCU reference.

2- and 4-rank Cartesian decompositions reproduce accepted fields within the defined floating-point tolerance.

Rank-boundary HP reconstruction/path, LS slopes/Hessian, G, curvature and evolving-topography assembly match equivalent interior operations.

Global CFL, conservation and bed-volume diagnostics are correct.

MPI+OpenMP combinations reproduce the accepted solution.

Parallel NetCDF owned-interior/owned-vertex writes reconstruct the same global scientific state as serial output.

Distributed restart reproduces supported continuation semantics, including B_vertex and stochastic state.

Dipendenze fra milestone scientifiche ed HPC

Dipendenze fra milestone scientifiche ed HPC

Percorso critico aggiornato: A08/A10/A09 \[CHIUSI\] -\> A11 baseline \[APPROVATA, 3f7329c85d2e\] -\> B01 \[ACCETTATA, 05989f5e7a\] -\> dipendenze e relocation path-only \[DA AVVIARE\] -\> consolidamento x/y e reologie/solver implicito \[PIANIFICATI\] -\> D15 \[DA ESEGUIRE\] -\> G01/G02 e sviluppi multilayer/MPI successivi.

Attività e verifiche operative: Combinazioni rank/thread; Parallel NetCDF single-file collective hyperslabs owned, esclusione halos e vertices duplicati; restart distribuito.

File interessati: src/inpout_2d.f90 e futuri io/\*; src/domain_2d.f90; nuovi test MPI+OpenMP/Parallel NetCDF e restart.

Criteri di accettazione: Confronto campi scientifici seriale/distribuito, nessuna duplice scrittura, completezza metadata e ripartenza stocastica.

Validate representative MPI rank/thread combinations against the frozen A11/G03 reference.

Implement one-file Parallel NetCDF using owned-interior hyperslab writes; halos are never written as owned scientific state.

Define the global nodal B_vertex NetCDF/restart ownership convention so shared rank-boundary vertices are written exactly once or collectively consistently.

Preserve serial NetCDF as a decoded-field reference and compare serial/parallel metadata and fields.

Define distributed restart ownership and validate continuation, including evolving-topography state.

### H05 — Tuning e GPU readiness

Stato ufficiale: DA FARE

Attività e verifiche operative: Profilare hot spots/NUMA/bandwidth/SIMD e tile sizes; documentare percorso GPU-ready senza introdurre dipendenze sul target acceleratore; non precalcolare topografia assunta statica.

File interessati: src/spatial_operator_2d.f90; src/hp_reconstruction_2d.f90; src/nonlinear_solver_2d.f90; src/geometry_2d.f90; benchmark e performance logs.

Criteri di accettazione: Baseline CPU/MPI validata intatta; prestazioni misurate su workload rappresentativi con throughput, scaling e memoria.

Profile only after distributed correctness is established.

Prioritize memory bandwidth, halo overlap, NUMA placement, vectorization and tile size using representative HP-PCCU workloads.

Treat derived geometry recomputation versus communication as a measured optimization decision, not an early architectural assumption.

Use tiles as the future accelerator unit instead of introducing a second GPU-specific state representation.

Preserve CPU/MPI reference execution and the same model/layout/local-domain contracts.

## Fase I — Estensioni scientifiche

### I01 — Formulazione fisica due strati

Stato ufficiale: FUTURO

### I02 — Scambi interstrato

Stato ufficiale: FUTURO

### I03 — Validazione due strati

Stato ufficiale: GATE

### I04 — Closure verticali fisiche/analitiche

Stato ufficiale: FUTURO

Definire un’interfaccia vertical_closure capace di ricevere stato depth-averaged e parametri di modello e restituire momenti verticali, fattori correttivi e termini di accoppiamento previsti dalle equazioni, con normalizzazione e unità documentate. La prima priorità scientifica è la distribuzione verticale delle classi particellari, con effetti sui flussi e prodotti non lineari, collegata a distribuzioni di velocità e temperatura quando rilevanti. Preservare la distinzione tra quantità mediate e medie di prodotti; non sostituire una struttura verticale inesistente con semplici campi diagnostici a profilo arbitrario.

Implementazioni alternative: closure analitica, fisica, semiempirica e ML/NN sotto la stessa API. Dati di verifica esterni alla sola simulazione depth-averaged: teoria, esperimenti e calcoli verticalmente risolti utilizzati offline. Requisiti: limiti e positività di concentrazioni, coerenza del bilancio depth-integrated, conservazione, differenziabilità necessaria ai metodi numerici e prestazioni misurate; non assumere che la vecchia vertical_profiles, già rimossa, sia un riferimento da ripristinare.

### I05 — Validazione verticale

Stato ufficiale: GATE

### I06 — Subgrid topografia fine/coarse

Stato ufficiale: FUTURO

Definire una procedura scientifica tracciabile di simulazioni IMEX_SfloW2D su DEM dettagliati e su DEM/griglie grossolane associati: allineamento geometrico, filtraggio/aggregazione della quota, mapping dei campi, medesima fisica e condizioni iniziali, rapporti di risoluzione, trattamento dei bordi e confronto di flussi, runout, deposizione, canali, blocco e direzioni del moto. Identificare il residuo di coarse-graining che la closure dovrà riprodurre, senza fissare aprioristicamente la sua forma come attrito basale.

Separare training/calibrazione da verifica su terreni e condizioni non visti. Prevedere statistiche topografiche multiscala e feature direzionali coerenti con i vincoli di simmetria e con la geometria B_vertex; stabilire se la correzione opera su un flusso, una sorgente fisica o un provider di interfaccia soltanto dopo una caratterizzazione matematica. Ogni correzione deve superare test di conservazione, positivà/ammissibilità, riproduzione del limite fine e convergenza al diminuire del rapporto di coarse-graining.

### I07 — Validazione subgrid

Stato ufficiale: GATE

### I08 — Inferenza ML/NN delle closure

Stato ufficiale: FUTURO

Prevedere un’interfaccia uniforme tra kernel di chiusura deterministici e modelli ML/NN, comprendente registro dei modelli disponibili, metadata/versioni dei pesi, normalizzazione e dimensione delle feature, unità, dominio valido di input, forma/rank delle uscite e stato di valutazione. Evitare dipendenze dirette del ciclo RK da framework specifici; cache e batch inference devono appartenere a workspace dedicati, senza allocazione per cella o stato globale nascosto. Il percorso senza ML deve mantenere i risultati di riferimento.

Un output ML non finito o fuori dal dominio fisico deve essere intercettato da regole locali della closure e da politiche documentate di errore o fallback, senza clipping silenzioso. L’addestramento è un’attività offline; soltanto l’inferenza è destinata al run del solver. Verificare scalabilità 1/4 thread e futura esecuzione MPI/GPU, riproducibilità su uno stesso backend e compatibilità con restart quando i modelli abbiano memoria.

#### Ambito esplicitamente escluso

Non pianificare coupling runtime bidirezionale o orchestrazione numerica con altri codici di simulazione. Modelli 3D o osservazioni esterne possono alimentare dataset offline per la chiusura verticale; questa attività non crea una dipendenza esecutiva del solver. Le interfacce input, dati, meteorologia e forcing già approvate restano distinte dal coupling fra simulatori.

### I09 — Validazione ML/NN

Stato ufficiale: GATE

#### Ordine, dipendenze e verifiche delle estensioni scientifiche

Prima definire model_layout, conservazione e semantica degli operatori, completare D15 e il gate post-D15 di architecture readiness. Poi introdurre API e test synthetic-only per vertical/subgrid/ML senza alterare il run di produzione. In seguito sviluppare i dataset e le chiusure fisiche con gate indipendenti; soltanto dopo integrare i backend ML. Le API devono essere preparate prima di bloccare la struttura layer-aware e MPI, mentre la validazione di nuove chiusure può proseguire successivamente.

| **Gate**  | **Evidenza richiesta**                                                                                                                                      |
|-----------|-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| V-GATE    | Profilo verticale integrato riproduce medie e flussi di riferimento; classi polidisperse e vincoli di massa; confronto offline teoria/esperimenti/3D.       |
| S-GATE    | Confronti fine/coarse sullo stesso solver, split indipendente per DEM e reologia, identificazione contributi non risolti, conservazione e generalizzazione. |
| ML-GATE   | Equivalenza del contratto delle closure con backend neurale e non neurale, versionamento modello/feature, fallback documentato, test OpenMP e performance.  |
| ARCH-GATE | Disattivazione delle nuove closure riproduce la baseline; nessun nuovo accoppiamento runtime con simulatori esterni.                                        |

# 5. Requisiti trasversali di completamento

Le prescrizioni seguenti sono obbligatorie per tutte le fasi del piano. I vincoli specifici di una attività sono riportati nella relativa scheda. La pubblicazione del sito non equivale al completamento del gate di copertura Doxygen.

**Regola trasversale obbligatoria — documentazione Doxygen e commenti nel sorgente**

REQUISITO DI REPOSITORY (applicabile a tutto src/, ai nuovi moduli e a ogni file modificato): ogni elemento del codice deve essere documentato in modo che un nuovo sviluppatore possa comprenderne scopo, interfaccia, assunzioni fisiche/numeriche, effetti collaterali e vincoli d’uso senza ricostruirli da altre parti del sorgente. Questa regola vale per tutte le milestone, incluse le attività XY, reologie, IMEX, ML, multilayer e HPC.

**Copertura minima e contratto per le procedure**

Ogni MODULE, TYPE, componente/variabile o parametro significativo (inclusi campi dei tipi, costanti, configurazioni, workspace, array di stato e variabili locali dal ruolo non immediato) riceve un commento Doxygen (!\> o !!) che ne descriva significato, unità fisiche ove applicabili, intervalli/valori speciali, dimensioni e ownership quando rilevanti. Sono documentate anche le interfacce PUBLIC/PRIVATE e le procedure TYPE-BOUND; evitare descrizioni tautologiche.

Ogni FUNCTION e SUBROUTINE, comprese quelle PRIVATE, interne, PURE, ELEMENTAL e legate a tipi, DEVE avere immediatamente prima della definizione un blocco Doxygen con almeno: @brief; @details quando necessari; una riga @param\[in\], @param\[out\] oppure @param\[in,out\] PER CIASCUN ARGOMENTO DUMMY, senza eccezioni per this/self; @return per FUNCTION (o @retval quando utile); unità, forma/rank, condizioni ammesse, precondizioni, effetti sullo stato, thread-safety e failure mode quando pertinenti. I tag @param devono corrispondere esattamente ai nomi nella signature e alle INTENT effettive. Per procedure con argomenti senza INTENT, esplicitarne la semantica nella documentazione e aggiungere INTENT dove legalmente possibile.

**Commenti esplicativi all’interno delle procedure**

Commentare i blocchi logici del corpo: inizializzazione, preparazione stencil, condizioni al bordo, conversioni di stato, ricostruzione, flussi, sorgenti, limiti di positività, accumulazione residui, aggiornamenti di massa/energia, sincronizzazione OpenMP/MPI, I/O e cleanup. I commenti devono spiegare COSA si sta facendo e, per i passaggi non ovvi, PERCHÉ: invarianti, scelte numeriche, ordine delle operazioni floating-point, casi degeneri, unità e conservazione. Non aggiungere commenti per ogni istruzione elementare né limitarsi a ripetere il codice.

**Modello Doxygen obbligatorio per Fortran**

!\> @brief Evaluate a cell contribution to the explicit spatial residual.  
!\> @details Uses reconstructed HP-PCCU traces; preserves the exact-rest guard.  
!\> @param\[in\] this Spatial-operator object and configuration.  
!\> @param\[in\] q Conservative cell state \[model units\].  
!\> @param\[in\] direction Normal Cartesian direction (1=x, 2=y).  
!\> @param\[out\] residual Conservative residual \[state units/s\].  
SUBROUTINE evaluate_cell(this, q, direction, residual)  
...  
END SUBROUTINE evaluate_cell

**Verifica automatica, code review e applicazione graduale**

Per ogni nuovo file e per ogni funzione/subroutine nuova o modificata, la copertura Doxygen e una condizione di accettazione della relativa attivita. Il controllo statico dei blocchi Doxygen e dei tag @param, nonche la generazione senza warning non giustificati, sono obblighi da realizzare e verificare progressivamente durante la relocation e il refactoring post-M12c. L'accettazione M12c dell'infrastruttura build/CI NON dimostra di per se il completamento del gate Doxygen.

Criterio di completamento: (a) ogni elemento introdotto o modificato è documentato; (b) tutti i dummy arguments hanno tag direzionali corretti; (c) i blocchi non banali sono spiegati; (d) il generatore Doxygen e lo static checker passano in CI per il perimetro previsto; (e) la documentazione è aggiornata nello stesso commit che modifica l’interfaccia o il comportamento. Nelle schede XY-1…XY-7, il refactoring deve preservare la documentazione esistente e aggiornarla prima del merge.

#### Vincoli di sequenziamento, documentazione e decisione di rilascio

I criteri di accettazione dello schema HP-PCCU comprendono topografia Q1 nodale B_vertex, mass fractions Y, ricostruzione wet/dry, operatori non conservativi e path idrostatici, correzione di slope, sorgente di curvatura, cell-local exchange e vertex gather con bilancio globale, ripartenza geometrica esatta. L’ottimizzazione non deve rendere statiche quantità topografiche che sono evolutive. Stato (var,layer,x_storage,y_storage) con N_LAYERS=1 è il prossimo passo, non una funzionalità già acquisita.

**Risultato finale richiesto:** codice leggibile e rintracciabile, identificazione testuale univoca delle reologie, nessun modello 8, solver RK generico e sorgenti fisiche modulari, gruppi misti Newton/analitici, test scientifici e build reproducibili, e una sola architettura congelata per il futuro lavoro MPI.

#### Attività confermate, non autorizzate e criteri di non-perdita

Decisioni confermate: gerarchia directory funzionali già avviata (con eventuale completamento e riordino mirato, non una seconda relocation indiscriminata), pulizia dei USE e delle dipendenze, documentazione Doxygen capillare, consolidamento x/y, modularizzazione reologica, organizzazione del solve implicito a gruppi nel perimetro già deliberato, retry/bounds, vertical closure e subgrid topografia anche ML/NN, MPI/OpenMP, Parallel NetCDF e GPU readiness. La semantica complex-step e i kernel REAL/COMPLEX non devono essere riscritti durante mere relocations.

Decisioni escluse: ricerca di nuovi tableaux RK/IMEX e accoppiamento runtime con simulatori esterni. Non estendere la modularizzazione reologica a un refactoring generale di tutte le sorgenti costitutive per il solo fatto che alcuni file architetturali futuri sono elencati nel piano. Dove il testo storico usa nomi di file o flag non più coerenti con master, trattarli come riferimenti di migrazione da verificare, non come autorizzazione a reintrodurre rami obsoleti.

Le attività realizzate nei commit del 12 ottobre non rendono automaticamente completati D04–D15, il gate Doxygen generale, la validazione esecutiva di tutti gli esempi, il rollback, le chiusure ML, il layer-aware rank-4 o MPI. Ogni avanzamento va supportato da verifiche puntuali e manifest di accettazione.

Ogni milestone deve produrre evidenze riproducibili di regressione, documentazione Doxygen delle procedure e degli argomenti, test dei casi degeneri pertinenti, verifica della semantica REAL/COMPLEX e complex-step quando applicabile, e controllo delle build e CI GNU/ifx/Docker/OpenMP. Ogni passaggio di stato deve citare commit e gate specifici. Le linee scientifiche future possono procedere in parallelo quando non esiste una dipendenza tecnica esplicita.

Regola generale sull’ownership: ogni parametro, configurazione, quantità derivata e stato persistente possiede un solo modulo responsabile. Il consumatore importa tramite API esplicite, senza duplicazione; le variabili temporanee hanno scope locale o workspace esplicito, con lifetime e thread-safety documentati. L’audit D03 e la migrazione D07 applicano questa regola e la verificano tramite regressioni; il gate E01 ne controlla indipendentemente l’efficacia.

# 6. Corrispondenza con gli identificativi precedenti

Le sigle storiche non sono identificativi operativi e sono conservate esclusivamente qui, per consultare versioni precedenti della roadmap e riferimenti nei commit. Ogni attività usa unicamente il codice della fase (A01–I09).

| Nuovo ID | Identificativo precedente | Attivita                                             |
|----------|---------------------------|------------------------------------------------------|
| A01      | N0                        | Baseline, test di riferimento e manifest complessivo |
| A02      | N1                        | Primitive canoniche e mass fractions                 |
| A03      | N2                        | Topografia Q1 e geometria B_vertex                   |
| A04      | N3                        | Ricostruzione HP                                     |
| A05      | N4                        | Operatore HP-PCCU, G=1                               |
| A06      | N5                        | Slope correction                                     |
| A07      | N6                        | Curvature source                                     |
| A08      | N7                        | Validazione idrodinamica 1D e 2D                     |
| A09      | D-N7-CFL                  | Due policy CFL                                       |
| A10      | N8                        | Topografia evolutiva e scambio di massa              |
| A11      | N9                        | Cutover e baseline scientifica                       |
| B01      | M12c                      | Build, CI e portabilita                              |
| B02      | P1                        | Mappa dipendenze e relocation path-only              |
| B03      | M12d                      | Stocastico OU, RNG e restart seriale                 |
| B04      | M13                       | Model layout, N_LAYERS=1                             |
| B05      | XY-1                      | Velocita caratteristiche comuni                      |
| B06      | XY-2                      | Ricostruzione HP sulle linee                         |
| B07      | XY-3                      | Flusso PCCU di interfaccia                           |
| B08      | XY-4                      | Velocita faccia e path idrostatici di cella          |
| B09      | XY-5                      | Condizioni al contorno e ricostruzione lineare       |
| B10      | XY-6                      | Derivate e halo geometrici                           |
| B11      | XY-7                      | Riduzione CFL direzionale                            |
| B12      | XY-GATE                   | Accettazione complessiva x/y                         |
| C01      | IO-1                      | Output segregati per run                             |
| C02      | IO-2                      | Contratto di configurazione e sicurezza              |
| C03      | IO-GATE                   | Accettazione I/O per run                             |
| C04      | DIST-1                    | Autotools policy B                                   |
| C05      | DIST-2                    | Doppio archivio di distribuzione                     |
| C06      | DIST-GATE                 | Validazione distributiva                             |
| D01      | M12e-A                    | Inventario formule e rami REAL/COMPLEX               |
| D02      | M12e-B                    | Caratterizzazione numerica cell-local                |
| D04      | P2                        | Eliminazione reologia 8                              |
| D05      | P3                        | Estrazione moduli reologici                          |
| D06      | P3                        | Sorgenti termiche e pore pressure                    |
| D08      | P4                        | RHEOLOGY_MODEL da numeri a stringhe                  |
| D09      | P5                        | Correzione del solo attrito                          |
| D10      | P6                        | EX_terms / IM_terms                                  |
| D11      | P7                        | IM_group e metadati                                  |
| D12      | P7                        | Newton per gruppo e Jacobiano ridotto                |
| D13      | P8                        | Algebriche ANALYTIC per stage DIRK                   |
| D14      | P8                        | Gruppi misti e validazione                           |
| D15      | P9                        | Acceptance della nuova architettura                  |
| E01      | ARCH-1                    | Audit dipendenze e ownership stato                   |
| E02      | ARCH-2                    | Architecture fitness checks                          |
| E03      | ARCH-GATE                 | Architecture readiness                               |
| E04      | TIME-1                    | Policy timestep e input                              |
| E05      | TIME-2                    | Bounds primitive e override                          |
| E06      | TIME-3                    | Controlli di ammissibilita cell-local                |
| E07      | TIME-4                    | Error propagation, rollback e retry                  |
| E08      | TIME-5                    | Test del timestep controller                         |
| E09      | P-TIME                    | Gate timestep                                        |
| F01      | V-API                     | Interfaccia closure verticali                        |
| F02      | S-API                     | Interfaccia subgrid terreno                          |
| F03      | ML-API                    | Backend ML/NN opzionale                              |
| F04      | EXT-TEST                  | Test di neutralita API                               |
| F05      | EXT-GATE                  | Architecture closure readiness                       |
| G01      | M14                       | Storage rank-4 con N_LAYERS=1                        |
| G02      | M15                       | Activity/restart/I-O layer-aware                     |
| G03      | M16                       | Architecture freeze                                  |
| H01      | M17                       | Local domain e halo abstraction                      |
| H02      | M18                       | Gerarchia active tiles                               |
| H03      | M19                       | Decomposizione MPI                                   |
| H04      | M20                       | MPI+OpenMP e Parallel NetCDF                         |
| H05      | M21                       | Tuning e GPU readiness                               |
| I01      | LAYERS-2                  | Formulazione fisica due strati                       |
| I02      | LAYERS-EX                 | Scambi interstrato                                   |
| I03      | LAYERS-GATE               | Validazione due strati                               |
| I04      | V-1                       | Closure verticali fisiche/analitiche                 |
| I05      | V-GATE                    | Validazione verticale                                |
| I06      | S-1                       | Subgrid topografia fine/coarse                       |
| I07      | S-GATE                    | Validazione subgrid                                  |
| I08      | ML-1                      | Inferenza ML/NN delle closure                        |
| I09      | ML-GATE                   | Validazione ML/NN                                    |
