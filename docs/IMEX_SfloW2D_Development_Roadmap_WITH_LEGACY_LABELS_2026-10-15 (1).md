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

## Regola generale vincolante — Validazione proporzionata al rischio (risk-based regression testing)

> **PRINCIPIO TRASVERSALE PER TUTTE LE FASI A–I.** Ogni modifica viene verificata in funzione del suo impatto tecnico e scientifico; **non è richiesto ripetere automaticamente l'intera suite numerica ad ogni modifica**. La validazione completa è riservata ai gate scientifici pertinenti, alle modifiche trasversali ad alto rischio e alle release. Un gate già approvato resta approvato nel proprio perimetro storico: le modifiche successive richiedono una nuova verifica d'impatto, non la riapertura retroattiva del gate.

### Classi di verifica e quando eseguirle

| Livello | Momento e finalità | Contenuto minimo | Non implica |
|---|---|---|---|
| **T1 — Mirato** | Ogni piccola modifica durante lo sviluppo | Compilazione/lint dei componenti toccati quando applicabile; test unitari, di contratto, negativi e edge case direttamente coinvolti | Riesecuzione di tutte le simulazioni scientifiche |
| **T2 — Integrazione** | Fine di una modifica funzionale e prima della revisione del relativo step | Un numero contenuto di casi end-to-end pertinenti; compatibilità delle interfacce, restart/I/O o output numerici secondo la superficie modificata | Campagna A11/N9 completa |
| **T3 — Regressione selettiva di gate** | Per accettazione di step o gate, e per modifiche trasversali | Sottoinsieme rappresentativo, congelato e riproducibile di regressioni che copre i percorsi realmente impattati; strict/optimized, 1/4 thread o toolchain aggiuntive **solo se pertinenti** | Riesecuzione automatica di ogni test storico |
| **T4 — Campagna estesa/scientifica e release** | Nuovi modelli fisici, cambiamenti numerici importanti, HPC con nuovi ordini/riduzioni, gate scientifici o release | Matrici, benchmark, conservazione, comparatori indipendenti, riproducibilità e CI secondo lo scope formalmente approvato | Approvazione scientifica oltre il perimetro effettivamente testato |

### Procedura obbligatoria per ogni change set

1. **Classificare la modifica**: solo documentazione, packaging, I/O/path, interfaccia/stato/restart, algoritmo numerico, fisica, parallelismo; includere dipendenze indirette reali. Se un file contiene più responsabilità (es. `inpout_2d.f90`), valutare *le routine e i flussi di dati modificati*, non soltanto il nome del file.
2. **Definire prima dei test la matrice minima**: T1 sempre ove applicabile, T2 quando cambia un comportamento osservabile, T3/T4 solo per rischi e gate motivati. Annotare motivazione dei test inclusi **e di quelli costosi esclusi**. Non usare l'assenza di modifiche alle equazioni come unica prova di equivalenza: percorsi I/O, checkpoint, scheduling e stato globale possono avere effetti indiretti.
3. **Usare comparatori appropriati**: bitwise se formato, ordine delle operazioni e ambiente lo consentono; confronti numerici con tolleranze dichiarate altrimenti; per I/O distinguere bytes scientifici/metadati intenzionalmente variati (path, manifest, timestamp). Conservare il riferimento alla baseline A11 con revisioni, input e checksum.
4. **Escalation guidata**: ampliare i test se un controllo fallisce, cambiano precisione/ordine delle operazioni, OpenMP/MPI, bilanci, discretizzazione, stato/restart, o emergono modifiche collaterali non previste. Non dichiarare un gate chiuso con soli smoke test.
5. **Registrare l'evidenza**: revisione Git, toolchain, profilo, conteggio thread/rank se pertinente, comando o workflow, casi, comparatore, risultato, esclusioni motivate, limiti e responsabile della decisione. Le evidenze storiche non sono automaticamente test eseguiti sul nuovo HEAD.
6. **Separare frequenze**: T1/T2 su sviluppo/PR; T3 nei gate e, ove conveniente, nella CI; T4 per release o milestone scientifiche. È lecito eseguire una suite ampia in CI a basso costo, ma ciò non è un requisito per ogni micro-commit o per ogni revisione locale.

### Prescrizione specifica per modifiche esclusivamente di I/O (Fase C)

Per C01/C02 sono richiesti test di naming, sicurezza directory, overwrite, manifest, serializzazione, output NetCDF/ASCII/CSV e almeno un caso full/split/restart con append. È richiesta almeno una regressione end-to-end breve sul contenuto scientifico, confrontando output omogenei dopo normalizzazione dei soli metadati deliberatamente cambiati. **I test specifici di ricostruzione HP, PCCU, Jacobiani, reologie e schemi RK non si ripetono tutti** se non sono stati toccati percorsi numerici/stato e i confronti end-to-end selezionati sono positivi. La chiusura C03 richiede il proprio gate I/O e la relativa evidenza; non riapre A11 o B12.

### Applicazione alle schede operative

Le voci **"Test proporzionati"** inserite sotto le schede A01–I09 sono *matrici prospettiche di selezione* e non dichiarazioni di test già eseguiti. I gate, gli stati ufficiali e gli esiti storici riportati nel testo originale restano invariati. Per gli step già completati queste voci valgono per eventuali modifiche future che li riguardino, non come richieste retroattive di nuova validazione.

### Regola generale per la chiusura degli step implementati ma non ancora accettati

Per ogni step con stato **`IMPLEMENTATO, GATE APERTO`** la scheda deve distinguere chiaramente: (i) i requisiti del contratto non ancora coperti da evidenze, (ii) le verifiche minime ancora necessarie e (iii) l'evidenza o decisione richiesta per il passaggio a **`COMPLETATO`**. Se non rimangono requisiti scoperti, il responsabile può chiudere il gate con una decisione documentata basata sulle evidenze già disponibili, senza ripetere test. Lo stato di uno step **non cambia automaticamente** perché è stato chiuso un gate più ampio: l'accettazione deve essere registrata esplicitamente nella sua scheda e nella tabella del Master Development Plan.

Le evidenze di un gate trasversale possono essere riutilizzate da più step soltanto entro lo **scope, la revisione e l'ambiente effettivamente verificati**. Un elemento non coperto richiede soltanto test mirati T1/T2 (o T3 se il rischio residuo lo giustifica), secondo la regola di validazione proporzionata al rischio; evitare campagne scientifiche generalizzate in assenza di impatto numerico.


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
| B08 | XY-4          | Velocita faccia e path idrostatici di cella    | COMPLETATO                |
| B09 | XY-5          | Condizioni al contorno e ricostruzione lineare | COMPLETATO                |
| B10 | XY-6          | Derivate e halo geometrici                     | COMPLETATO                |
| B11 | XY-7          | Riduzione CFL direzionale                      | COMPLETATO                |
| B12 | XY-GATE       | Accettazione complessiva x/y                   | COMPLETATO                |

Evidenza B05–B12 raccolta al 15 ottobre 2026: le fattorizzazioni XY-1…XY-7 e il gate complessivo x/y sono completati nel perimetro verificato. Sono passate le regressioni spaziali e dinamiche strict/optimized, incluse prove a 1/4 thread, composizione variabile, fixture ruotate e confronti binari pre/post refactoring sui casi esplicitamente elencati. L'accettazione B12 non estende l'approvazione scientifica A11 oltre la revisione e lo scope approvati, né chiude D15.

## Fase C — Output e distribuzione

|     | Vecchia label | Attività                                | Stato   |
|-----|---------------|-----------------------------------------|---------|
| C01 | IO-1          | Output segregati per run                | COMPLETATO |
| C02 | IO-2          | Contratto di configurazione e sicurezza | COMPLETATO |
| C03 | IO-GATE       | Accettazione I/O per run                | COMPLETATO |
| C04 | DIST-1        | Autotools policy B                      | COMPLETATO |
| C05 | DIST-2        | Doppio archivio di distribuzione        | COMPLETATO |
| C06 | DIST-GATE     | Validazione distributiva                | COMPLETATO |

La fase C è **COMPLETATA per il perimetro tecnico di validazione**. L'automazione per allegare gli archivi sorgente ed Examples a una release GitHub e la pubblicazione di una release ufficiale sono follow-up separati, esplicitamente differiti; il candidato `9.9.27` resta sintetico e non pubblicato.

## Fase D — Reologie e solver IMEX

|     | Vecchia label | Attività                                                    | Stato    |
|-----|---------------|-------------------------------------------------------------|----------|
| D01 | M12e-A        | Inventario formule e rami REAL/COMPLEX                      | COMPLETATO (audit statico) |
| D02 | M12e-B        | Caratterizzazione numerica cell-local                       | IMPLEMENTATO, GATE APERTO |
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

**Test proporzionati (nuove modifiche a questo step):** T1 manifest/checksum, harness e reference; T2 riproduzione di un caso di riferimento; T3 audit della baseline solo se cambia il contratto di riferimento. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 round-trip q/qp, mass fractions e wet/dry; T2 bilanci di componenti con sorgenti/BC; T3 reference e composizione variabile solo per cambi al mapping. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO

### A03 — Topografia Q1 e geometria B_vertex

**Test proporzionati (nuove modifiche a questo step):** T1 identità Q1 facce/centri e refresh; T2 fondo variabile/restart; T3 regressioni geometriche e 1/4 thread se cambia accesso o riduzione. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO

### A04 — Ricostruzione HP

**Test proporzionati (nuove modifiche a questo step):** T1 positività, velocity interval, momenti e bordi; T2 caso 1D/2D wet/dry; T3 operator/reference se varia la ricostruzione. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO

### A05 — Operatore HP-PCCU, G=1

**Test proporzionati (nuove modifiche a questo step):** T1 flussi/path/residui e quiet water; T2 dam-break e topografia continua; T3 confronto con reference HP-PCCU per modifiche all’operatore. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO

### A06 — Slope correction

**Test proporzionati (nuove modifiche a questo step):** T1 G_center/G_face e path corretto; T2 lake-at-rest su pendii; T3 confronto residui/CFL se cambiano le formule. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO

### A07 — Curvature source

**Test proporzionati (nuove modifiche a questo step):** T1 sorgente a riposo, termine misto Bxy uv; T2 caso dinamico curvo; T3 regressione fisica se cambia la sorgente. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO

### A08 — Validazione idrodinamica 1D e 2D

**Test proporzionati (nuove modifiche a questo step):** T1 comparatori e controlli di bilancio; T2 casi rappresentativi 1D/2D; T3/T4 matrice scientifica N7 solo quando si modifica il suo scope. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 due formule CFL, sentinel e input non valido; T2 traiettorie N_RK 2/3/4 rappresentative; T3 profili/thread selezionati se cambia scheduling. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO

### A10 — Topografia evolutiva e scambio di massa

**Test proporzionati (nuove modifiche a questo step):** T1 rate locale, disponibilità e assemblaggio nodale; T2 fondo evolutivo e restart B_vertex; T3 bilanci/1-4 thread se cambiano le riduzioni. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO

### A11 — Cutover e baseline scientifica

**Test proporzionati (nuove modifiche a questo step):** T1 integrità comparatori, manifest e checksum; T2 casi ETNA/ETNA_LAVA selezionati; T4 solo per nuova decisione scientifica esplicita, senza riaprire il gate storico. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 configure/build e make check mirato; T2 build GNU nel profilo modificato; T3 GNU/ifx/Docker soltanto quando cambia portabilità/build/release. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 inventario dei sorgenti, riferimenti ai path e dependency graph; T2 build pulita e out-of-tree; T3 matrice di toolchain se il relocation cambia build. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO

### B03 — Stocastico OU, RNG e restart seriale

**Test proporzionati (nuove modifiche a questo step):** T1 sequenza RNG/OU e celle asciutte; T2 full/split restart stocastico; T3 1/4 thread se cambia ordine o ownership RNG. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: IMPLEMENTATO

**Attività previste:** Preservare modello OU, flussi RNG indipendenti, seed riproducibili, memoria latente nelle celle asciutte e continuità del restart; mantenere test esistenti durante refactoring e definire policy di correlazione tra rank in H03.

**File interessati:** src/stochastic_module.f90; src/stochastic_random_2d.f90; src/runtime_2d.f90; src/inpout_2d.f90; TESTS/TEST_STOCHASTIC_EVOLUTION/\*; TESTS/TEST_STOCHASTIC_RESTART/\*; TESTS/TEST_STOCHASTIC_RANDOM/\*.

**Criteri di accettazione:** Test seriale/OpenMP e restart stocastico invariati; verifica statistica della correlazione; nessun consumo RNG introdotto dal nuovo solver o cambio di numerazione delle celle.

### B04 — Model layout, N_LAYERS=1

**Test proporzionati (nuove modifiche a questo step):** T1 layout, indici e rifiuto N_LAYERS>1; T2 smoke del solver; T3 regressione di stato se cambia memoria o serializzazione. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: IMPLEMENTATO

**Attività previste:** Conservare model_layout_type quale unica autorità su n_layers, variabili ed indici; rifiutare N_LAYERS\>1 finché la fisica multilayer non esiste; usare metadati per G01 anziché crearne una seconda copia.

**File interessati:** src/model_layout_2d.f90; src/parameters_2d.f90; src/init_2d.f90; src/simulation_2d.f90; TESTS/TEST_MODEL_LAYOUT/\*.

**Criteri di accettazione:** Layout N_LAYERS=1 coerente con dimensioni e indici; rifiuto esplicito configurazioni non supportate; nessuna variazione dei risultati.

### B05 — Velocita caratteristiche comuni

**Test proporzionati (nuove modifiche a questo step):** T1 velocità x/y, dry/wet e g’ non positiva; T2 residui/CFL ruotati; T3 bitwise strict/optimized se si modifica nuovamente il kernel. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: IMPLEMENTATO

Obiettivo e vincolo: condividere il codice che realizza lo stesso algoritmo nelle due direzioni cartesiane, mantenendo separati storage orientati, cicli sulle interfacce e ownership OpenMP. La direction determina solo l’orientamento; composizione, layer, closure verticale, reologia e sub-grid restano gestiti da interfacce fisiche distinte. Vietato introdurre array temporanei dinamici per faccia, transpose obbligatori o dipendenze circolari.

File interessati: src/equation_terms_2d.f90; src/hyperbolic_2d.f90.

Attività: Sostituire eval_local_speeds_x/y con un kernel direzionale condiviso, selezione idx_u/idx_v, stessa chiamata mixt_var e identico ramo per gravità ridotta negativa. Valutare piccoli wrapper compatibili temporanei; non modificare la formula delle velocità caratteristiche.

Criteri di accettazione: Stesse velocità min/max per entrambi gli assi, anche con dry/wet e gravità ridotta non positiva; confronto unitario su x/y e regressioni CFL/HP.

Evidenza: implementato nel commit 00400001. La suite TEST_SPATIAL_OPERATOR è passata in profilo strict e optimized con verifiche CFL direzionali e team OpenMP effettivi da 1/4 thread; TEST_N7_GEOMETRY_DYNAMIC è passata con 69 casi per profilo, inclusi stati dry/wet e casi x/y ruotati. Queste regressioni supportano lo step; il confronto esplicito isolato delle velocità per il caso di gravità ridotta non positiva resta da includere nella matrice B12.

### B06 — Ricostruzione HP sulle linee

**Test proporzionati (nuove modifiche a questo step):** T1 tracce HP, quiet-water, positività e indici x/y; T2 caso dinamico ruotato; T3 payload bitwise se cambia la ricostruzione. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: IMPLEMENTATO

File interessati: src/reconstruction_2d.f90; src/hp_reconstruction_2d.f90.

Attività: Consolidare evaluate_hp_cell_x/y in un kernel che prepari gli stencil con direzione esplicita e chiami reconstruct_hp_line invariata. Conservare scelta delle facce W/E/S/N, normal momentum, mask e hp_scratch thread-private; preservare la barriera prima del commit.

Criteri di accettazione: Tracce HP, mean momentum, positività, quiet-water e controlli x↔y invariati; output deterministici 1/4 thread bitwise dove previsto.

Stato aggiornato: implementato nel commit 00400001. TEST_HP_RECONSTRUCTION passa i controlli di quiet-water, positività, media del momento e indipendenza direzionale; TEST_N7_GEOMETRY_DYNAMIC passa con 69 casi per profilo e output binari identici tra 1 e 4 thread. La precedente dicitura «DA IMPLEMENTARE» descriveva uno stato superato e non è applicabile.

### B07 — Flusso PCCU di interfaccia

**Test proporzionati (nuove modifiche a questo step):** T1 PCCU/path e guard exact rest; T2 casi rough/drainage/lake; T3 confronto pre/post e thread se cambiano i kernel. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: IMPLEMENTATO

File interessati: src/hyperbolic_2d.f90; src/pccu_2d.f90.

Attività: Estrarre la sequenza comune eval_inertial_flux, central-upwind, limit_component_mass_flux, exact-rest scalar guard, eval_hydrostatic_coefficient, regularize_dry_gamma_pair, eval_hydrostatic_path ed eval_oriented_pccu_pair. Lasciare separati i due DO OpenMP e le matrici di faccia.

Criteri di accettazione: Flussi, path e residual identici rispetto alla baseline per gli stessi input; conservare l’ordine delle operazioni numeriche e il guard di riposo esatto introdotto nei commit recenti.

Stato aggiornato: kernel PCCU comune implementato nel commit f2065730. In entrambi i profili strict e optimized, i payload spaziali completi sono risultati bitwise identici alla baseline pre-refactoring per dynamic_rough, dynamic_drainage e lake_ramp, a 1 e 4 thread. TEST_PCCU, TEST_SPATIAL_OPERATOR e TEST_N7_GEOMETRY_DYNAMIC (69 casi per profilo) sono passati. Il test storico TEST_HYDROSTATIC_ROUNDOFF non è compatibile con l’interfaccia post-XY-1 perché il suo modulo di riferimento importa i wrapper x/y rimossi; l’equivalenza è stata verificata con un confronto diretto HEAD/candidato. B12 resta aperto per l’accettazione complessiva.

### B08 — Velocita faccia e path idrostatici di cella

**Test proporzionati (nuove modifiche a questo step):** T1 segnali di faccia, orientazioni e path di cella; T2 griglie ruotate/composizione; T3 confronti bitwise solo per nuovi cambi ai kernel; B12 storico invariato. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: IMPLEMENTATO, GATE APERTO

File interessati: src/hyperbolic_2d.f90.

Attività: Centralizzare i due stati di faccia, riduzione min/max, percorso idrostatico dentro cella, convenzioni di orientamento e segni. Evitare la fusione dei loop OpenMP e rispettare accessi allineati a memoria.

Criteri di accettazione: Equivalenza x/y su casi 1D ruotati, 2D e composizione variabile; assenza di differenze su CFL direzionali.

Stato aggiornato: velocità di faccia e cell hydrostatic path consolidati nel commit fd43e4ca. I payload spaziali completi sono risultati bitwise identici alla baseline pre-XY-4 per dynamic_rough, dynamic_drainage e lake_ramp, in strict e optimized, con 1/4 thread. TEST_N7_GEOMETRY_DYNAMIC è passata (69 casi per profilo); TEST_N7_COMPOSITION è passata (162 casi per profilo, composizione variabile e confronti a 1/4 thread); TEST_SPATIAL_OPERATOR ha passato controlli di tracce, residui e CFL direzionale, inclusi casi x/y ruotati. B08 è COMPLETATO con la chiusura del gate B12.

### B09 — Condizioni al contorno e ricostruzione lineare

**Test proporzionati (nuove modifiche a questo step):** T1 Dirichlet/Neumann/secante W/E/S/N e dimensioni 1xN; T2 lateral/radial source; T3 payload bitwise se si interviene di nuovo sulle BC. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: IMPLEMENTATO, GATE APERTO

File interessati: src/spatial/reconstruction_2d.f90; src/geometry/geometry_2d.f90; src/physics/equation_terms_2d.f90.

Attività: Estrarre logica riutilizzabile dei casi Dirichlet/Neumann/interpolazione per lato senza cancellare semantica dei boundary source, radial source o lateral source; nessuna modifica simultanea allo schema.

Criteri di accettazione: Tests mirati W/E/S/N e sorgenti, incluso nx=1/ny=1; output identici sulle configurazioni supportate.

Stato aggiornato al 15 ottobre 2026: in `reconstruction_2d.f90` è stato estratto `boundary_reconstruction_slope`, kernel comune per le pendenze Dirichlet, Neumann e secante sui lati W/E/S/N. Il test in `TEST_IMEX_STAGES` verifica in strict e optimized i valori analitici delle tracce per le tre modalità su ciascun lato; verifica lateral-source W/E, riflessioni radial-source W/E/S/N e griglie 5x5, 1x5, 5x1 e 1x1. Per tutti i casi, i payload completi di ricostruzione sono bitwise identici alla baseline pre-XY-5 (sorgente `reconstruction_2d.f90` del commit c648e9c5356c65d68d0624369d408cec781f0816, SHA-256 90d5d50923069798cd6135f18a940cdaa423a7d26644044ca9989e483d6bfd43), sia in strict sia in optimized, a 1 e 4 thread; i payload sono inoltre identici fra i due conteggi di thread. `TEST_HP_RECONSTRUCTION`, `TEST_SPATIAL_OPERATOR` e `TEST_IMEX_STAGES` sono passati. B09 è COMPLETATO con la chiusura del gate B12.

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

**Test proporzionati (nuove modifiche a questo step):** T1 geometria 1D x/y, 2D e halo; T2 topo dinamica; T3 bitwise 1/4 thread se cambiano ordine, copie o somme. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: IMPLEMENTATO, GATE APERTO

File interessati: src/geometry/geometry_2d.f90.

Attività: Valutare line kernels per derivate prime e seconde e utility di copia bordi senza variare ordine di accumulazione, formule Q1/LS/Hessian, trattamento delle dimensioni degeneri o mixed derivative. Non precomputare definitivamente la geometria: topografia futura time-dependent.

Criteri di accettazione: TEST_GEOMETRY_1D per x e y e fixture 2D, ricostruzione dopo aggiornamenti B_vertex, confronto 1/4 thread e output di topografia invariati.

Stato aggiornato al 15 ottobre 2026: le copie dei bordi dei campi derivati sono state estratte in utility private x/y; i loop 2D, le formule e le somme miste non sono stati modificati. L'estrazione del kernel aritmetico 1D è stata scartata perché il confronto bitwise ha rilevato differenze in optimized; mantenere i loop scalari originali conserva l'ordine numerico. `TEST_GEOMETRY_1D` è passato in strict e optimized (112 casi per profilo, fixture 1D x/y e 2D, cache avvelenate e snapshot identici a 1/4 thread); `TEST_TOPOGRAPHY_FACES` è passato in entrambi i profili. Il confronto A/B con il sorgente pre-B10 ha confermato payload identici bitwise su griglie 5x1, 9x1, 17x1, 5x5 e 7x9, tutti i quattro profili di topografia e slope-correction OFF/ON, a 1 e 4 thread, in strict e optimized. B10 è COMPLETATO con la chiusura del gate B12.

### B11 — Riduzione CFL direzionale

**Test proporzionati (nuove modifiche a questo step):** T1 directional_min/combined, CFL=-1, MAX_DT; T2 N_RK e output temporali; T3 thread/profili se si modifica la riduzione. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: IMPLEMENTATO, GATE APERTO

File interessati: src/spatial/spatial_operator_2d.f90.

Attivita: dopo A11, fattorizzare ove utile la raccolta dei massimi direzionali senza alterare i criteri scelti da CFL_POLICY. La funzionalita duale directional_min/combined e un lavoro circoscritto del gate A09, PRIMA di A11; B11 e soltanto refactoring interno successivo.

Criteri di accettazione: stesso risultato numerico per ciascuna delle due policy prima e dopo il refactoring, conservando comportamento CFL=-1, MAX_DT, griglie 1D x/y e ordine di operazioni appropriato. Nessuna variazione di default.

Stato aggiornato al 15 ottobre 2026: la raccolta dei quattro segnali di velocità di faccia per cella è stata estratta in una utility privata; le due riduzioni `MAX` mantengono lo stesso ordine di operandi e la riduzione OpenMP, mentre le formule `directional_min`/`combined`, il default e i loro rami 1D non cambiano. `TEST_SPATIAL_OPERATOR` è passato in strict e optimized: verifica le due formule, griglie x/y/point, propagazione nulla, cap `MAX_DT`, sentinel `CFL=-1` e team effettivi 1/4. `TEST_CFL_POLICY` è passato con solver strict e optimized, 18 traiettorie per profilo a 4 thread (N_RK=2,3,4; policy omessa, directional_min, combined; CFL adattivo e sentinel), confrontando output canonici, checkpoint e tempi accettati; sono passati anche i tre controlli di input non valido. B11 è COMPLETATO con la chiusura del gate B12.

### B12 — Accettazione complessiva x/y

**Test proporzionati (nuove modifiche a questo step):** T3 matrice x/y e comparatori strict/optimized nel perimetro del gate; T4 non implicito: non riaprire B12 per modifiche I/O o documentazione. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO

#### Sequenziamento e disciplina di regressione

Gate di ingresso: A11 \[APPROVATO\] e B01 \[ACCETTATO\]. Prima eseguire la relocation path-only con manifest sorgenti e dipendenze, senza modifica di istruzioni o ordine floating-point. Successivamente sviluppare in commit separati le fattorizzazioni x/y e la riorganizzazione reologica/IMEX; eseguire CI GNU, Docker e Intel, oltre ai confronti A11 pertinenti. Non confondere la riuscita dei test di compilazione con l'accettazione scientifica D15.

Evidenza aggiornata al 15 ottobre 2026: l'audit completo sul candidato basato sul commit `f3f6240ede126c8c577b7e3f2da3ca638bec0352` ha superato 45/45 test, con build `strict_debug` e `historical_optimized` e 8/8 confronti di thread. Il patch candidato di codice/test verificato (prima di questa annotazione in roadmap) è identificato da SHA-256 `9e7fac3ed441b77dd3ed5afe2f2e7facac7193df22c2b6bb35ac55d46720f091`; include solo documentazione e adattamenti dei test/audit, senza variazioni numeriche del solver. Sul commit B12 `307818fb7be2730768bf437c646a77cc2d707e11` sono passate le workflow Solver validation (run 38704914944), Docker validation (run 38704915006) e Intel ifx portability (run 38704914990). B08–B11 e B12 sono quindi COMPLETATI nel perimetro x/y verificato. L'audit mantiene esplicitamente `n9_decision=open` e non concede sign-off: A11 resta limitato al proprio scope approvato e D15 resta separato.

## Fase C — Output e distribuzione

### C01 — Output segregati per run

**Test proporzionati (nuove modifiche a questo step):** T1 layout di directory, naming, indice/tempo snapshot e manifest; T2 output NetCDF/ASCII/CSV, multi-run e assenza di file solver nella cwd; T3 confronto scientifico breve pre/post paths. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: PARZIALE

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

#### Inventario iniziale dei file prodotti dal solver

Ricognizione del sorgente al 15 ottobre 2026: le scritture di output scientifico e diagnostico sono centralizzate in `src/io/inpout_2d.f90`; il checkpoint viene invocato da `src/IMEX_SfloW2D.f90`. I percorsi correnti sono relativi alla working directory e vanno ricondotti alla directory del run senza cambiare i contenuti dei file:

| Categoria | File attuali | Destinazione C01 |
|-----------|--------------|------------------|
| Backup parametri effettivi | `<RUN_NAME>.bak` | `metadata/` |
| Serie NetCDF | `<RUN_NAME>.nc` | `fields/` |
| Snapshot conservativi | `<RUN_NAME>_<idx>.q_2d`, `_erosion_`, `_erodible_`, `_deposit_`, `_B_` | `snapshots/<idx>_t_<time>/` |
| Snapshot raster | `<RUN_NAME>_<idx>.asc`, `_T_`, `_dep_`, `_depTot_`, `_ers_`, `_ersTot_` | `snapshots/<idx>_t_<time>/` |
| Diagnostica massima e runout | `<RUN_NAME>_max.asc`, `<RUN_NAME>_runout.csv`, `<RUN_NAME>_mass_center.csv` | `diagnostics/` |
| Probe | `<RUN_NAME>_prb_<probe>.csv` | `probes/` |
| Vulnerability table | `<RUN_NAME>_VT.txt`, `<RUN_NAME>_VT_<idx>.asc` | `vulnerability/` |
| Raster del terreno computazionale | `dem_esri.asc`, `dem_esri_nodata.asc`, opzionalmente `dem_esri_water.asc` | `terrain/` |
| Restart scritto | `restart.bin` | `restart/` |

I file input `topography_file`, `erodible_file` e i checkpoint indicati esplicitamente per la lettura sono aperti in lettura e non vanno spostati. L'inventario è il riferimento per l'implementazione e per l'audit completo delle scritture, inclusi i percorsi di restart e le routine di post-processing.

#### Implementazione C01/IO-1 — 15 ottobre 2026

È stato aggiunto il modulo `src/io/output_paths.f90` con l'helper POSIX `src/io/output_dirs.c`; tutti gli output del solver censiti sono instradati sotto `<OUTPUT_BASE_DIRECTORY>/<RUN_NAME>/`. Il valore predefinito della base è `OUTPUT`; `OUTPUT_OVERWRITE` è falso per default e consente il riuso esplicito senza cancellare file preesistenti. I restart binari in-place possono riutilizzare la propria directory senza overwrite e risolvono il checkpoint scritto `restart/restart.bin`; i checkpoint esterni restano input esterni. I percorsi rifiutano componenti symlink e `..`; prima del riuso vengono inoltre rifiutati symlink, hard link o elementi non regolari già presenti nell'albero del run. `RUN_NAME` ammette solo lettere ASCII, numeri, `_` e `-`. Ogni snapshot ha una cartella indice/tempo e il manifest in `metadata/output_manifest.csv` registra l'indice e il tempo esatto.

La documentazione introduttiva e la pagina website “Output files and formats” descrivono il layout e il contratto. Verifiche passate sul build GNU locale: `TESTS/run_tests.sh` (casi legacy, struttura, manifest, collisione senza overwrite, riuso esplicito, `RUN_NAME` e base con traversal rifiutati, basi symlink, file symlink e hard link rifiutati, basi personalizzate relative e assolute con spazi); `TESTS/TEST_BINARY_RESTART/run_test.sh`; `TESTS/TEST_STOCHASTIC_RESTART/run_test.sh`; `TESTS/TEST_PCCU_GATE_H/run_test.sh`; `TESTS/TEST_PCCU_LAKE_AT_REST/run_test.sh`; `TESTS/TEST_PCCU_INCLINED_EXCAVATION/run_test.sh`; `TESTS/TEST_CFL_POLICY/run_test.sh`. L'harness Acceptance/N9 è stato aggiornato per i nuovi percorsi, ma non è stato eseguito integralmente. Con l'accettazione di C03, C01 è COMPLETATO; l'inventario delle scritture e le verifiche mirate di percorso, restart e integrazione post-processing sono registrati nella sezione C03.

### C02 — Contratto di configurazione e sicurezza

**Test proporzionati (nuove modifiche a questo step):** T1 input RUN_NAME, path relativi/assoluti, traversal, symlink/hardlink, overwrite; T2 existing-run/restart esterno, directory HPC lunghe; T3 stress/race pertinenti al gate C03. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO

#### Contratto di configurazione e sicurezza

RUN_NAME già esistente identifica univocamente il run; validare caratteri ammessi ed escludere nomi vuoti, separatori, “.”, “..” e traversals.

OUTPUT_BASE_DIRECTORY è configurabile, di default OUTPUT; sono ammessi percorsi relativi o assoluti, inclusi filesystem scratch HPC, purché non contengano componenti `..` e non attraversino symlink.

OUTPUT_OVERWRITE è FALSE per default: un nuovo run non deve sovrascrivere output di una directory esistente. Un’autorizzazione TRUE è esplicita e resta confinata al run selezionato; non rimuove i file preesistenti. La creazione usa componenti directory senza seguire symlink.

Un restart ordinario riutilizza la directory dello stesso run, senza richiedere overwrite, ma verifica compatibilità checkpoint e integrità degli output preesistenti. La lettura esplicita di checkpoint esterni resta supportata.

Nome del checkpoint di scrittura FISSO: restart.bin in \<OUTPUT_BASE_DIRECTORY\>/\<RUN_NAME\>/restart/restart.bin. Nessun parametro utente per rinominarlo. Il chiamante può continuare a passare il nome letterale alla routine di I/O; il path viene costruito centralmente.

Snapshot ASCII organizzati in sottodirectory con indice univoco e tempo fisico leggibile; registrare il tempo esatto in un manifest. Il NetCDF resta un unico file di serie temporale in fields/, non duplicato negli snapshot.

Vulnerability è categoria autonoma di primo livello, NON una sottodirectory di diagnostics/. Preservare inizialmente i nomi dei file prodotti dal solver.

Tutti i file scritti direttamente dal solver, inclusi backup di input, DEM esportati, raster, diagnostiche, serie, massimi e checkpoint, devono rimanere sotto il run. File di input originali e prodotti di script esterni non vengono spostati automaticamente.

Centralizzare la costruzione/validazione dei percorsi di output in un componente I/O dedicato. Definire eventuale modalità legacy temporanea solo per il periodo di validazione e migrazione degli script di post-processing.

#### Criterio esplicito di chiusura del gate C02 — verifica documentale delle evidenze

**Decisione formale del 16 ottobre 2026:** dopo la verifica della matrice requisito-evidenza, della suite locale e della CI GNU/ifx/Docker, il responsabile approva la chiusura di C02 per la revisione sorgente `bd4a87d8a8a7f9876950a7aa3024b41e8c7612d6`. La decisione accetta esplicitamente come limite il mancato test su un filesystem scratch HPC reale, già escluso dal gate C03; non modifica il contratto dei percorsi né le verifiche locali di path relativi/assoluti, traversal e symlink. C02 è COMPLETATO; C03 resta COMPLETATO.

**Copertura da verificare e tracciare:** associare ogni requisito C02 al test o all'evidenza C03 che effettivamente lo copre, indicando revisione Git e ambito verificato. La verifica include: validazione di `RUN_NAME`; percorsi relativi/assoluti di `OUTPUT_BASE_DIRECTORY` e rifiuto di `..`/symlink; semantica di `OUTPUT_OVERWRITE` e isolamento del singolo run; collisioni e directory preesistenti; restart interno e checkpoint esterno con compatibilità e protezione delle scritture; assenza di output solver fuori dal run; naming degli snapshot/manifest e destinazioni dei file; controlli su symlink/hard link e preservazione degli output in caso di input o storia di restart non validi.

**Evidenze già riportate in C03:** i test mirati su sicurezza dei percorsi, multi-run, overwrite, restart, manifest, fault injection e confronti di output sono descritti nella scheda C03; possono essere richiamati **senza rieseguirli**, purché coprano il requisito specifico e la revisione pertinente. La decisione scientifica N9 associata a C03 non sostituisce il controllo puntuale del contratto C02.

#### Matrice requisito-evidenza C02

| Requisito | Evidenza e risultato | Revisione e ambito | Lacuna o limite |
|---|---|---|---|
| Validazione di `RUN_NAME`, rifiuto del traversal e assenza di scritture su input non valido | `TESTS/run_tests.sh`: dopo un run valido, rifiuta `../escape` e confronta file, checksum e directory dell'albero `OUTPUT`; verifica l'errore atteso e l'assenza di output esterno. Suite GNU locale e CI: PASS. | Test committato in `bd4a87d8a8a7f9876950a7aa3024b41e8c7612d6`; il solver GNU locale è costruito da `432c2bbac6f09c6e94ad1bf9e95bac9ba16802b7`. Solver run `38877762627`, ifx `38877778700`, Docker `38877794457`: tutti i quattro job di ciascuna matrice PASS sul commit `bd4a87d8`. | Evidenza sufficiente per il requisito nel perimetro accettato. |
| `OUTPUT_BASE_DIRECTORY` relativo/assoluto, traversal, symlink e hard link | `TESTS/run_tests.sh` e verifiche C03 di percorsi personalizzati, symlink/hard link e sentinelle: PASS come riportato nella sezione C03. | Evidenze locali C03; il change set I/O/restart registra base `96b4cd3a` e build GNU derivata da `0a965c67`. Percorsi relativi/assoluti con spazi verificati. | Il responsabile accetta esplicitamente di non eseguire il test su filesystem scratch HPC reale; i controlli del contratto dei percorsi sono stati verificati localmente. |
| `OUTPUT_OVERWRITE`, collisioni e isolamento tra run | `TESTS/run_tests.sh` e `TESTS/TEST_OUTPUT_PATHS/run_test.sh`: collisione senza overwrite, riuso esplicito e due run distinti nella stessa base; PASS nelle evidenze locali C03. | Evidenze C03 e suite GNU locale; le verifiche di layout includono i casi documentati nella sezione C03. | Nessuna lacuna funzionale identificata nel perimetro provato. |
| Restart in-place/esterno, compatibilità e preservazione degli output validi | `TESTS/TEST_BINARY_RESTART/run_test.sh` e `test_runout_restart.py`: restart full/split, checkpoint esterno, manifest e output incoerenti rifiutati senza alterare i file validi; fault injection 29/29 PASS per la campagna registrata in C03. | Evidenze C03 su build GNU/historical optimized e revisione candidata indicate nella relativa sezione; non sono una nuova esecuzione sulla revisione attuale. | L'eventuale eccezione di filesystem resta quella già dichiarata; non è richiesta una nuova campagna scientifica per C02. |
| Destinazioni degli output, snapshot/manifest, indici e assenza di output fuori dal run | Audit delle scritture e `TESTS/TEST_OUTPUT_PATHS/run_test.sh`: output sotto il run, manifest e indici 0/9999/10000/10001; PASS nelle evidenze C03. | Evidenze locali C03, inclusi i test T1/T2 riportati per la revisione `96b4cd3a` e la build derivata da `0a965c67`. | Nessuna lacuna funzionale identificata nel perimetro provato. |

La suite `TESTS/run_tests.sh` è stata rieseguita localmente con il solver GNU indicato nella prima riga: tutti i test legacy e di sicurezza dei percorsi sono PASS. La regressione è stata poi verificata nella CI sul commit `bd4a87d8`: tutti i job Solver, ifx e Docker sono passati.

**Valutazione delle lacune:** per i requisiti del contratto C02 documentati nella matrice non è emerso un ulteriore test T1/T2 indispensabile; T3 non è motivato da un rischio residuo specifico. Il limite relativo al filesystem scratch HPC è stato esplicitamente accettato dal responsabile nella decisione di chiusura registrata sopra.

**Procedura minima per chiudere C02:**

1. Compilare una matrice requisito C02 → evidenza/test C03 → risultato → revisione/scope → eventuale lacuna; non basta dichiarare genericamente che C03 è completo.
2. Identificare in modo esplicito gli eventuali requisiti non verificati o esclusi dallo scope. Per ciascuno distinguere fra test mirato ancora necessario e limitazione consapevolmente accettata. Non introdurre nuovi requisiti retroattivamente.
3. Eseguire **solo** eventuali test T1/T2 indispensabili a colmare lacune reali; T3 soltanto in presenza di un rischio residuo motivato. Nessuna nuova campagna di regressione HP-PCCU/IMEX o ripetizione della matrice A11/N9 è richiesta per una chiusura puramente documentale.
4. Registrare la decisione del responsabile, con data, revisione Git, evidenze riutilizzate, eccezioni e limiti; **solo dopo** aggiornare sia la scheda C02 sia la riga C02 nella tabella generale a `COMPLETATO`.

**Condizione di uscita:** tutti i requisiti del contratto C02 risultano coperti da evidenza appropriata oppure da eccezioni esplicite approvate; non rimangono controlli obbligatori senza decisione. **Condizione soddisfatta; stato C02: `COMPLETATO`.**


### C03 — Accettazione I/O per run

**Rapporto tra C03 e C02:** la chiusura di C03 documenta evidenze utilizzabili anche per C02, limitatamente ai requisiti effettivamente verificati. **Non determina automaticamente lo stato di C02**: la chiusura di quest'ultimo richiede la matrice di copertura e la decisione specifica descritte nella scheda C02. Non ripetere la campagna C03 per produrre una seconda evidenza identica.

**Test proporzionati (nuove modifiche a questo step):** T1 tutti i negativi di path e collisione; T2 full/split/external restart, append NetCDF/CSV, snapshot interrotto; T3 2-3 scenari numerici rappresentativi con baseline A11; nessuna campagna HP-PCCU completa automatica. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO

#### Test, documentazione e gate di accettazione

Inventario delle operazioni di creazione/scrittura file in tutto il codice e classificazione completa degli output; test che nessun output solver sia creato nella working directory.

Test multi-run, directory già esistente, overwrite esplicito, percorso personalizzato, validazione RUN_NAME e protezione da cancellazioni fuori dal run.

Test restart v5 in place, da checkpoint esterno, append NetCDF/CSV, interruzione durante snapshot e ripartenza senza sovrascrivere risultati validi.

Confronto quantitativo/bitwise, quando applicabile, con baseline A11 per verificare che la modifica dei path non alteri risultati numerici né tempi di output.

Aggiornamento di EXAMPLES e script in UTILS, regressioni e documentazione GitHub Pages: nuova pagina “Output files and formats”, naming, manifest, semantica overwrite e restart.

Solo dopo accettazione della riorganizzazione dell’I/O procedere con politica Autotools B, make distcheck e distribuzione di archivi distinti codice/examples dallo stesso tag.

#### Verifiche integrative della code review — 15 ottobre 2026 (C03-R1–R4)

Le verifiche seguenti **integrano** i criteri C03 già elencati, senza sostituirli né riaprire le attività B08–B12. Erano osservazioni di revisione statica sui file locali forniti (`IMEX_changes.patch` e quattro nuovi file), non regressioni osservate né test eseguiti durante quella revisione. Le relative evidenze riproducibili sono registrate di seguito.

**C03-R1 — Indici degli snapshot oltre quattro cifre.** In `src/io/output_paths.f90`, `begin_snapshot` usava `CHARACTER(LEN=4) :: index_string` e `WRITE(index_string, '(I4.4)') index`. La correzione estende l'etichetta; l'audit ha inoltre trovato e corretto il troncamento a quattro caratteri in `output_esri`, che poteva riutilizzare nomi raster dopo l'indice 9999. Il test integrato ora verifica etichette e directory 0/9999/10000/10001 e produce il raster ESRI `TEST_2D_10000.asc` durante un restart ad alto indice.

**C03-R2 — Lunghezza dei nomi di checkpoint.** `src/io/inpout_2d.f90` conserva `CHARACTER(LEN=40) :: restart_file`, a fronte del supporto a percorsi assoluti e basi HPC più lunghi. La dichiarazione non prova da sola un guasto, ma richiede di tracciare parsing, array `restart_files`, copie intermedie e chiamate di apertura. Provare checkpoint esterni con nome assoluto oltre 40 caratteri, con spazi e con `OUTPUT_BASE_DIRECTORY` lungo; correggere eventuali troncamenti.

**C03-R3 — Restart esterno e autorizzazione al riuso.** Verificare end-to-end quando il restart binario, riconosciuto anche tramite estensione `.bin`, abilita il riuso di una directory esistente. Distinguere: checkpoint interno dello stesso run; checkpoint esterno con nuovo output assente; checkpoint esterno con output già esistente e `OUTPUT_OVERWRITE=FALSE`; medesimo caso con overwrite esplicito; checkpoint mancante/incompatibile. L'uso di un checkpoint esterno non deve autorizzare implicitamente la scrittura su output preesistenti senza la semantica di restart in-place verificata. Controllare anche l'assenza di nuovi output nella working directory.

**C03-R4 — Integrità di output, riuso e interruzioni.** I controlli su tipo di file, symlink e hard link in `src/io/output_dirs.c` non attestano da soli la compatibilità dei contenuti di un run esistente. Definire i controlli obbligatori su checkpoint, manifest, snapshot e append NetCDF/CSV; simulare interruzioni durante creazione snapshot, scrittura manifest e append, quindi riprendere verificando indici univoci, tempi coerenti/monotoni, dati validi preservati e nessuna sovrascrittura indebita.

#### Evidenze locali aggiuntive

Le correzioni C03 estendono `snapshot_index_label` oltre quattro cifre, usano la stessa etichetta per directory e file snapshot, e rimuovono il troncamento dell'indice raster in `output_esri`; i buffer dei checkpoint sono stati allungati e il riuso automatico è limitato al restart binario interno al run. `TESTS/TEST_OUTPUT_PATHS/run_test.sh` verifica le etichette e directory distinte per gli indici 0, 9999, 10000 e 10001; ora esegue anche due run con `RUN_NAME` distinti nella stessa base e verifica alberi di output separati. Il test end-to-end ad alto indice riprende un checkpoint con `output_idx=10000`, verifica manifest monotono senza indici duplicati e controlla sia lo snapshot conservativo sia il raster ESRI con suffisso `10000`.

Per questo change set di test I/O/restart (T1 mirato e T2 integrazione), la base Git è `96b4cd3a` e il solver usato è il build out-of-tree GNU Fortran derivato da `0a965c67`, configurato con NetCDF Homebrew e OpenMP (`-fopenmp`). Sono passati `make -j2 -C /tmp/imex-c03-build.SjtHPP`, `OMP_NUM_THREADS=1 make -C /tmp/imex-c03-build.SjtHPP check` e `OMP_NUM_THREADS=1 TESTS/TEST_BINARY_RESTART/run_test.sh /tmp/imex-c03-build.SjtHPP/src/IMEX_SfloW2D`. La validazione del manifest avviene subito dopo il caricamento del checkpoint binario e prima delle scritture di output; anche il backup dell'input effettivo è differito fino al superamento del preflight. Le righe del manifest devono avere indici contigui e tempi finiti e crescenti. I test verificano il restart full/split bitwise, l'append all'indice 10000, il rifiuto di un manifest interrotto senza modificare alcun file preesistente del run (backup, snapshot, CSV e checkpoint inclusi), e la ripresa dopo la riparazione; verificano inoltre il rifiuto byte-preserving di un checkpoint esterno incompatibile anche con `OUTPUT_OVERWRITE=T`, e il riuso con un checkpoint esterno compatibile. In quella selezione storica non erano inclusi l'interruzione snapshot e l'append NetCDF interrotto; tali verifiche sono state aggiunte e registrate nella sezione fault-injection qui sotto. A11, N9, multi-toolchain/thread e GitHub Actions/CI appartengono a verifiche distinte; non sono rieseguiti automaticamente a ogni modifica I/O.

Il confronto con la baseline A11 è stato eseguito sul commit `3f7329c85d2e30e33cf68348721329f7726b01c1`, con GNU Fortran 15.2 e profilo optimized. Per `TEST_2D` e `TEST_BOX` gli snapshot conservativi e `restart.bin` sono identici byte per byte; il backup input di `TEST_BOX` differisce solo per i default registrati e la dimensione del buffer dei checkpoint. Per ETNA a 200 s sono identici byte per byte NetCDF, checkpoint, 115 raster ASCII e i due CSV scientifici; il manifest del candidato contiene 41 record da 0 a 200 s ogni 5 s. Le due esecuzioni hanno usato quattro thread. Hash NetCDF ETNA: `a603f6b69c841e0d7127022c8ee192ea3bc00372654d6e43bb157e40659e4463`; hash checkpoint: `36b7dbca8d5bc6c33c1dbbfe0a1af556d02dd9313959c34bb848850cd30a6333`.

Acceptance r4 è stata eseguita sul candidato working-tree basato su `6b7ff8dd56acbe9b1eaa518c1313fad0f8e15aad` (snapshot di 414 file, SHA-256 della patch `87bdfd72edf3fe067e5bdac6de051733dd4d588acdeeba12bf5ef1f0b2124669`), su macOS arm64, con profili `strict_debug` e `historical_optimized`: build riuscite, 45/45 esecuzioni passate e 8 confronti bitwise tra thread passati. Manifest locale: `/tmp/imex-c03-acceptance-6b7ff8dd-20261015-r4/manifest.json` (SHA-256 `ae8ca4bcfb6be5ba03c91cca5fe28f9e7aa83f89244cd5a9337738c4496c1b06`). Il risultato riporta tuttavia `n9_accepted=false` e `n9_decision=open`: i casi N9 long non sono stati rieseguiti, perché l’harness rifiuta l’evidenza v5 con hash sorgente non corrispondenti al candidato. Non trasferire quindi al candidato la precedente decisione/evidenza N9. I test Python `TESTS/ACCEPTANCE` risultano 23/24; l’unico fallimento è il controllo storico dell’hash del fixture `TESTS/TEST_GEOMETRY_1D/run_test.sh`, che non è stato modificato.

#### Fault injection end-to-end — 15 ottobre 2026

`TESTS/TEST_BINARY_RESTART/run_test.sh` è stata eseguita sul build audit `historical_optimized` con un thread effettivo: **29/29 controlli PASS**. Il helper d'interposizione è compilabile per macOS/Linux; questa esecuzione è stata validata su macOS arm64. Per snapshot, manifest e CSV interrompe una vera chiamata `write` sul file aperto dopo aver scritto metà del buffer; per NetCDF si arresta dopo la sincronizzazione del record `time`, prima dei campi del record. Il restart rifiuta i file incompleti senza ulteriori modifiche; dopo la riparazione esplicita dei soli artefatti incoerenti, checkpoint, manifest, snapshot, CSV e dati NetCDF recuperati corrispondono alla traiettoria di controllo. Evidenza: `/var/folders/gh/4whp2yg162q27q7q47ng3s200000gn/T/imex-runout-restart-40gdlxoa/evidence.json` (SHA-256 `965840b1451a1e24b9bf3e769736f486614e7713d2928f013053d0aa4b5d345b`); eseguibile SHA-256 `d9c01b6ee47b565df6022757de6d60b5cb68ceac4891476f9d5c67ba68afc722`. Questo risultato è limitato al build, al sistema operativo e al thread indicati; non sostituisce la matrice multi-toolchain/thread né chiude C03.

#### Verifiche mirate del worktree C03

Sul build GNU out-of-tree aggiornato (`/tmp/imex-c03-targeted.T3fjtY/src/IMEX_SfloW2D`), `TESTS/TEST_BINARY_RESTART/run_test.sh` ha superato **29/29 controlli**, includendo l'asserzione sul raster ESRI con indice 10000 e l'inventario dei file fuori da `OUTPUT/<RUN_NAME>/`. Evidenza: `/var/folders/gh/4whp2yg162q27q7q47ng3s200000gn/T/imex-runout-restart-1sazwac5/evidence.json` (SHA-256 `e0ad2daf083708c749d0f878bb33567cdc899d98609a6624f7ab75460ca7fa09`); eseguibile SHA-256 `be95244a8dc666010a762586b217146b688fca9c8cf624a36ff4f9a1f63e7dde`. `TESTS/TEST_OUTPUT_PATHS/run_test.sh`, eseguito via `sh` sullo stesso build, ha superato entrambi i controlli: indici 0/9999/10000/10001 e due `RUN_NAME` distinti nella stessa base, con manifest e directory separate.

L'audit dei consumer ha individuato l'uso di `delim_whitespace=True` di pandas 3 anche in `create_hazard_maps.py`; il parser ESRI di entrambi gli script ETNA è ora compatibile tramite `sep=r'\s+'`. Il test su raster reale del solver per `create_prob_maps.calculate_map` ha prodotto `result.001/TEST_2D_max.asc`. Anche `create_hazard_maps.read_asc_file` ha letto lo stesso raster (68 x 100); per questo controllo mirato le dipendenze top-level non necessarie alla funzione (`netCDF4` e GDAL/`osgeo`) sono state sostituite con stub temporanei, poiché non sono installate nell'ambiente locale. Non è stata eseguita l'intera pipeline hazard-map né una campagna completa degli esempi. Gli altri usi di `delim_whitespace` individuati sono nei generatori/preparatori di input, non in questi consumer di output, e restano fuori da questo intervento.

Queste verifiche chiudono l'audit locale dei writer, il test multi-run e il controllo mirato del post-processing. Il test su filesystem scratch HPC è escluso dal gate C03, secondo la decisione dell'utente, perché non è l'ambiente di deployment previsto. Nessuna campagna CI o suite generale è stata avviata.

#### Acceptance/N9 diagnostica r4 — sorgenti correnti

La campagna diagnostica N9 è stata completata sul build audit r4 e sulla revisione sorgente `0a46c555bdaa8a9337213a8dfa0231dedaf522a9`, in modalità `current_build_audit_diagnostic`. Risultato: **8/8 continuazioni fisiche PASS**, nessun fallimento, nessuna anomalia diagnostica e **4/4 confronti di thread PASS** per ETNA ed ETNA_LAVA nei profili `historical_optimized` e `strict_debug`. Evidenza: `/tmp/imex-c03-n9-diagnostic-0a46c555-r4/evidence.json` (SHA-256 `350a7e5540089cc075d06ab12b4ba2504dbf4e384545d548841daf42b347fd49`). È una verifica diagnostica sui sorgenti correnti, non un confronto scientifico approvativo con la baseline: `n9_decision=open` e `N9_scientific_approval=false` restano invariati.

L'evidenza r4 sopra è storica e non rappresenta la revisione più recente. Il test su filesystem scratch HPC è escluso dal gate C03, secondo la decisione dell'utente, perché non è l'ambiente di deployment previsto.

#### Acceptance/N9 diagnostica — commit C03 `48a92a785460`

Sul commit `48a92a785460be30f958fcc665cc29a440079316` sono state create da `git archive` due build pulite out-of-tree GNU Fortran 15.2, `historical_optimized` (`-Ofast -funroll-all-loops`) e `strict_debug` (`-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`). Il manifest mirato di provenienza build (non una nuova acceptance completa N7/N8) è `/tmp/imex-n9-commit-48a92a785460/manifest.json` (SHA-256 `d998746af32fcd57a5a72f03ed559861e1dfaedc0128d9fe9cc3dccfb6fcff18`).

`TESTS/ACCEPTANCE/n9_long_cases.py --diagnostic-current-source` ha completato ETNA (200 s) ed ETNA_LAVA (3600 s) con entrambi i profili e team effettivi da uno e quattro thread: **8/8 continuazioni PASS**, inclusi controlli fisici e CSV runout/mass-centre; **4/4 confronti bitwise tra thread PASS**. Evidenza: `/tmp/imex-n9-commit-48a92a785460-diagnostic/evidence.json` (SHA-256 `af25317a837190c24afbcffbc3dc0fb16a327a1947eacfc4020850e6f7a0edfe`). La modalità diagnostica non confronta i risultati con la baseline scientifica approvata; il manifest mantiene `n9_decision=open` e `N9_scientific_approval=false` perché il runner non assegna approvazioni scientifiche.

#### Decisione scientifica N9 e chiusura C03 — 16 ottobre 2026

Il responsabile scientifico ha approvato N9 per il commit `48a92a785460be30f958fcc665cc29a440079316`, dopo aver esaminato l'evidenza diagnostica sopra. Questa approvazione è la decisione del responsabile e non un confronto aggiuntivo con la baseline: il manifest del runner resta invariato e conserva correttamente `n9_decision=open` e `N9_scientific_approval=false`.

Con tale decisione e le verifiche C03 registrate, **C03 è COMPLETATO**. C01 è COMPLETATO; C04–C06 restano DA FARE/GATE e possono procedere secondo il sequenziamento della roadmap.

#### Checklist consolidata per la chiusura di C03 (inclusi i requisiti già previsti)

- [x] Audit delle scritture del solver e controllo dell'inventario fuori da `<OUTPUT_BASE_DIRECTORY>/<RUN_NAME>/` per i run di restart verificati.
- [x] Confronti quantitativi/bitwise appropriati con A11 per campi, diagnostiche, file e tempi di output, senza estendere lo scope scientifico A11.
- [x] Valutazione scientifica e decisione N9 sulla revisione corrente; approvazione del responsabile scientifico registrata il 16 ottobre 2026.
- [x] Matrice locale restart in-place/esterno, append NetCDF/CSV, interruzione e ripresa senza perdita o sovrascrittura di risultati validi; vedi fault injection e test mirati sopra.
- [x] Multi-run, collisioni, `RUN_NAME`, overwrite, traversal, symlink/hard link e percorsi relativi/assoluti verificati nell'harness locale.
- [x] Filesystem scratch HPC escluso dal gate C03 perché non è l'ambiente di deployment previsto.
- [x] Verifiche aggiuntive C03-R1–R4 e registrazione di risultati e riproducibilità.
- [x] Audit dei percorsi documentati in EXAMPLES/UTILS e GitHub Pages, più test mirato dei consumer ETNA su raster reale; la pipeline hazard-map completa non è stata eseguita (netCDF4/GDAL non disponibili localmente).
- [x] Sequenziamento rispettato: C04–C06 non sono stati anticipati rispetto all'accettazione di C03 e restano il lavoro successivo.

**Stato ed evidenze:** l'implementazione C01/C02 e i test GNU selezionati sono descritti sopra e risultano *riportati* nella roadmap; non sono stati rieseguiti in questa revisione. La diagnostica N9 mirata sul commit corrente è passata e il responsabile scientifico ha approvato N9 il 16 ottobre 2026. C03 è COMPLETATO; C04–C06 non sono stati avviati.

### C04 — Autotools policy B

**Test proporzionati (nuove modifiche a questo step):** T1 autoreconf e manifest EXTRA_DIST; T2 configure/make/check da checkout; T3 make distcheck e GNU/ifx/Docker al gate packaging, non a ogni modifica documentale. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO (verifiche C04 concluse nel gate C06)

#### Evidenza C04 — 16 ottobre 2026

I file generati Autotools non sono più versionati: `configure`, `Makefile.in` (root e `src/`), `aclocal.m4`, `compile`, `depcomp`, `install-sh` e `missing` sono stati rimossi dall'indice Git. `.gitignore` copre questi output e quelli di configurazione; `README.md` istruisce a rigenerarli con `autoreconf -fi`. Non è stato necessario modificare `src/Makefile.am`: l'inventario ha verificato che i quattro file `.inc` presenti sotto `src/` sono tutti inclusi in `EXTRA_DIST`.

Su un worktree temporaneo pulito del commit C03, dopo l'applicazione delle modifiche C04, sono passati `autoreconf -fi`, `configure --with-netcdf=/opt/homebrew/opt/netcdf-fortran`, `make -j2` e `make check`. La rigenerazione ha inoltre verificato che gli output Autotools risultino ignorati da Git. Evidenza locale conservata in `/tmp/imex-c04-git-ffbea0b9` e `/tmp/imex-c04-build-ffbea0b9`; `git diff --check` passa.

Al 16 ottobre, `make distcheck` e la matrice ifx/Docker erano differiti al gate di packaging C06 come previsto. Sono stati completati e accettati nel successivo audit C06; C04 non crea autonomamente gli archivi C05/C06.

#### Milestone post-XY — Politica Autotools B e distribuzione con due archivi

Sequenziamento: iniziare soltanto dopo il completamento e l’accettazione del refactoring x/y, della riorganizzazione dell’I/O e della relocation dei sorgenti nelle directory funzionali. L’attività è di packaging/build e non deve essere confusa con l’accettazione scientifica D15. Stato al 12 ottobre 2026: decisione di indirizzo presa; implementazione e test non ancora eseguiti.

Politica Autotools B: in Git versionare configure.ac, Makefile.am, src/Makefile.am e le eventuali macro sorgente; non versionare configure, Makefile.in (anche in src/), aclocal.m4, depcomp, install-sh, missing e altri output rigenerabili. Aggiornare .gitignore e le istruzioni per compilare un checkout tramite autoreconf -fi; assicurare compatibilità GNU, ifx e Docker.

Due artefatti della medesima release/tag, con versione e commit coerenti: (1) IMEX_SfloW2D-\<version\>.tar.gz, archivio sorgente ufficiale creato da Autotools, contenente Fortran, include, configure e file ausiliari generati, documentazione minima, licenza e test necessari; (2) IMEX_SfloW2D-\<version\>-examples.tar.gz, contenente EXAMPLES/ con input, script, DEM/dati necessari e tutti i README.md standardizzati, mantenendo i percorsi relativi. I nomi sono contratti proposti, con versione derivata dal tag e non fissata in anticipo.

Audit preliminare obbligatorio: verificare make distcheck e i contenuti EXTRA_DIST, l’inclusione di tutti i .inc e dei sorgenti nelle nuove sottodirectory, l’eseguibilità dei test da un archivio estratto e l’assenza di dipendenze implicite dal checkout Git. Controllare la dimensione e le licenze dei dati Examples, gli archivi DEM, le dipendenze tra esempi e il percorso dell’eseguibile dopo installazione. Non dichiarare gli archivi autonomi finché questi controlli non passano.

Automazione release: creare un workflow GitHub Actions che generi e testi entrambi gli archivi dalla stessa revisione, confronti un manifest dei file e checksum, poi li alleghi alla release ufficiale soltanto dopo esito positivo. Chiarire nel README e nel sito che i download automatici «Source code (zip/tar.gz)» di GitHub NON equivalgono all’archivio Autotools preparato e possono non contenere configure.

Gate di accettazione: (a) build/test da clone Git con autoreconf -fi; (b) configure + make + make check senza Autotools dall’archivio sorgente; (c) make distcheck; (d) riproducibilità degli Examples estratti con l’archivio del codice corrispondente, senza dipendenze non dichiarate; (e) validazione su GNU/ifx/Docker secondo matrice supportata; (f) documentazione aggiornata su GitHub Pages e CONTRIBUTING.md. Non eliminare i vecchi README.txt prima dell’audit integrale già pianificato.

### C05 — Doppio archivio di distribuzione

**Test proporzionati (nuove modifiche a questo step):** T1 versioni, checksum e inventario di entrambi gli archivi; T2 esempi estratti; T3 coerenza stessa revisione/tag prima del rilascio. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Implementato `tools/release/build_archives.py`: il tool risolve un tag Git, genera l'archivio sorgente con `make dist` e l'archivio Examples dall'identico commit, stampando la versione Autotools dal tag. Il manifest JSON registra revisione, checksum SHA-256 e inventario completo. `Makefile.am` distribuisce documentazione, licenza, i fixture richiesti da `make check` e i due dataset usati da tali fixture; i README storici restano inclusi nell'archivio Examples. Le fixture conservano in Git un solo input minuscolo, così da evitare collisioni di case sui filesystem macOS. Nei test, le copie temporanee vengono rinominate nel nome esatto `IMEX_SfloW2D.inp` richiesto dal solver.

T1–T3 verificati in un repository temporaneo isolato con tag sintetico `9.9.13` (revisione `1309bd62dcdecdda25b17b4d6033dcdf0c349ba8`): entrambi gli archivi sono stati generati, i checksum e la versione nel `configure` verificati indipendentemente, l'inventario Examples estratto ha conservato tutti i percorsi e i 20 `README.md`, e il manifest riporta lo stesso tag/revisione per entrambi. Non erano ancora stati eseguiti `make check` dall'archivio estratto, `make distcheck`, prove ifx/Docker o pubblicazione: sono stati affrontati nel gate C06. Provenienza/licenze dei dataset restano non verificate e non sono un blocco per la validazione tecnica; non dichiarare diritti accertati. Le dipendenze operative e la riproducibilità di alcuni esempi restano un audit distinto.

Stato ufficiale: COMPLETATO (validazione T1–T3; automazione e pubblicazione della release differite come follow-up)

### C06 — Validazione distributiva

**Test proporzionati (nuove modifiche a questo step):** T3 sorgente Autotools e Examples estratti, installazione, GNU/ifx/Docker, assenza dipendenze Git; T4 release solo dopo gate distributivo. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO (validazione distributiva; release non pubblicata)

**Evidenza locale:** in un clone temporaneo isolato è stato creato il tag sintetico `9.9.16` sulla revisione `0ee23e80a3bd93847d7199323cb558c78d2656f5`; non è un tag del repository di lavoro né una release. `build_archives.py` ha generato entrambi gli archivi e il manifest dalla stessa revisione. Il tarball sorgente è stato estratto senza `.git`; `configure`, `make distcheck` e l'installazione/suite `make check` incluse nel gate sono passati con GNU Fortran e NetCDF-Fortran su macOS. Un ulteriore smoke test ha eseguito `TEST_2D` con l'eseguibile installato dal tarball. L'archivio Examples estratto conserva tutti i 20 `README.md` e i file di input/dati verificati per `EXAMPLE_CF`.

Il primo `distcheck` ha evidenziato due requisiti di pulizia: rendere scrivibili le fixture copiate dai sorgenti di sola lettura e rimuovere i moduli Fortran `.mod` durante `clean`. Sono stati corretti `TESTS/run_tests.sh` e `src/Makefile.am`; il `distcheck` finale sull'archivio generato ha superato anche `distcleancheck`.

#### Chiusura tecnica della fase C — 21 ottobre 2026

Per decisione dell'utente, la fase C si considera completata per la validazione tecnica. Le evidenze includono la validazione degli archivi estratti, `make distcheck`, gli smoke test e la matrice GitHub Actions accettata: Solver e Intel ifx (run `39563559928` e `39563583552`) sul commit `c3b47cbf9a07162f1f7cbaf8b26230c05ef87798`; Docker, tutti e quattro i job (run `39567697006`), sul commit `cf34b8d13e1739f0b4113b30f8365a2253f547e3`. Il secondo commit aggiunge solo la dipendenza Automake di `init_2d`; l'utente ha accettato le Actions sui due commit consecutivi senza richiederne il rerun.

Non è stato creato un workflow per allegare automaticamente entrambi gli archivi a una release GitHub e non è stata pubblicata una release ufficiale. L'automazione e la pubblicazione sono differite a un follow-up distinto; `docker_create_publish.yml` pubblica solo l'immagine Docker. Il candidato `9.9.27` resta quindi un artefatto sintetico non pubblicato.

Le Actions Solver, Intel ifx e Docker sono state completate con successo sul commit `bd4a87d8a8a7f9876950a7aa3024b41e8c7612d6`: run `38877762627`, `38877778700` e `38877794457`, rispettivamente, con tutti e quattro i job di ciascuna matrice passati. La CI comprende il test di preservazione dell'albero output su `RUN_NAME` non valido e conferma che la precedente regressione di naming/input è risolta. La documentazione aggiornata è stata pubblicata con successo sul commit `320499a3270d7e09d65f29f72506cc52a4812a0e` dal workflow unificato, run `38900609895` (build e deploy riusciti). Provenienza e licenze dei dataset restano non verificate, ma non sono un blocco per la validazione tecnica C06; le dipendenze e la riproducibilità di alcuni esempi restano da verificare.

#### Audit documentale degli Examples — 16 ottobre 2026

L'inventario locale comprende 18 directory di esempio, 122 file tracciati e 20 README Markdown standardizzati. L'audit documentale non ha eseguito gli esempi e non attribuisce provenienza o licenze sulla base dei soli nomi dei file.

| Area | Evidenza locale | Esito/azione richiesta |
|---|---|---|
| Completezza e generatori | `EXAMPLE_1D_VERT` non contiene file tracciati e non è incluso nell'archivio; `EXAMPLE_BUMP` (quattro regimi), `EXAMPLE_BW`, `EXAMPLE_BOX` ed `EXAMPLE_SLOPE` hanno README, template e generatori. L'esecuzione ha rivelato che i generatori indicizzavano come scalari array di forma `(n_points, 1)`, incompatibili con NumPy 2.x quando valorizzano le medie cell-centered. | I generatori sono stati corretti con l'indicizzazione esplicita della colonna e i sette casi sono stati eseguiti con esito positivo dall'ambiente di test. `EXAMPLE_1D_VERT` non fa parte della distribuzione corrente; ripristinarlo solo se è un caso previsto. |
| Provenienza e diritti dei dati | Non sono indicate fonti originali o licenze verificabili per `EXAMPLE_CF/topo.zip`, i DEM di `EXAMPLE_TAAL` e `EXAMPLE_TAAL_BOTTOM`, `EXAMPLE_ENSEMBLE_ETNA/sectors.asc` e `Etna2014_crop.asc`. | Se reperibili, aggiungere evidenze da fonti autorevoli (titolare, URL/DOI, termini di riuso e metadati); in assenza, mantenere la provenienza e i diritti esplicitamente non verificati senza bloccare le verifiche tecniche di C06. |
| DEM Etna documentato | `Etna2014_crop.asc` non è presente come file estratto, ma è contenuto in `EXAMPLE_ETNA/DEM/Etna2014_crop.zip`; il README documenta l'estrazione. La ricerca non ha identificato la fonte originale o i termini applicabili. | La via di acquisizione dal repository è verificata; provenienza e licenza originale restano non verificate e non bloccano la validazione tecnica. |
| Dipendenze e riproducibilità | I README dei casi verificati indicano Python 3 e NumPy; Matplotlib è opzionale per i plot. `EXAMPLE_ENSEMBLE_ETNA` documenta i placeholder `ENSEMBLE_` e le fasi di preprocessing. Dall'archivio estratto sono stati completati tutti i 20 membri; aggregazione e plotting hanno prodotto 31 raster e 31 PNG. L'hazard mapping opzionale non è stato eseguito: nell'ambiente di test mancano netCDF4 e GDAL e il template disabilita l'output NetCDF. | Continuare l'audit degli esempi non ancora provati; mantenere espliciti requisiti e limiti delle fasi opzionali. |
| Coerenza documentazione/codice | L'audit dei casi residui ha verificato i README contro generatori, input e codice sorgente: `EXAMPLE_2D` usa `create_example.py`, richiede NumPy (Matplotlib solo per il plot) e abilita NetCDF/ESRI; `EXAMPLE_CF` e i due casi TAAL abilitano NetCDF ma non ESRI. In `EXAMPLE_TAAL`, sorgente laterale, `TIME_PARAM(4)` è il cut-off assoluto; l'esempio legacy con periodo 100 s e cut-off 200 s genera due impulsi, non uno. In `EXAMPLE_TAAL_BOTTOM`, sorgente bottom, `TIME_PARAM(4)` è invece uno sfasamento di fase: il valore corrente 150 s equivale a cinque periodi da 30 s e non arresta la sorgente. | README standardizzati allineati ai file e alla semantica implementata; conservare i README storici come archivio. Resta aperta la comparazione integrale dei testi legacy e l'esecuzione fresca di tutti gli esempi; non inferire interpretazioni fisiche non verificate né diritti dei dataset. |

La verifica esterna ha confermato che il repository contiene gli archivi/file distribuiti, ma non ne stabilisce l'origine o i diritti. I due raster Taal sono copie identiche; SRTM è soltanto un possibile indizio non collegato ai file da evidenze di provenienza o confronto. La licenza GPL-3.0 del solver non prova i termini applicabili ai dataset. Provenienza e diritti restano quindi non verificati: questa incertezza va riportata senza attribuire licenze o diritti accertati, ma non costituisce una priorità né un blocco per le verifiche tecniche e la chiusura di C06.

#### Validazione degli archivi estratti e della documentazione — 16 ottobre 2026

In un clone temporaneo isolato è stato creato il solo tag sintetico `9.9.17` sulla revisione `355f504591696cc045759ebeb6d43d27855562ab`; non è un tag del repository di lavoro e non è stata creata una release. Il manifest associa alla stessa revisione l'archivio sorgente (`SHA-256 1936216bc12b3e5e9bd4ef9c5f2aac53d0035faf221fa64b5ed1d8e379c3a9fd`) e l'archivio Examples (`SHA-256 520eef48a218841e19cbbefaa7d33e437c21608c84b395b334580d715c9d1b9c`).

Il sorgente è stato estratto senza `.git`: `configure --with-netcdf=/opt/homebrew/opt/netcdf-fortran`, `make -j2` e `make check` sono passati. `make distcheck` è passato con `DISTCHECK_CONFIGURE_FLAGS='--with-netcdf=/opt/homebrew/opt/netcdf-fortran'`; un primo tentativo senza propagare questa opzione al configure annidato non trovava `netcdf.mod` e non costituisce un difetto del sorgente. L'installazione in un prefisso temporaneo è riuscita. Dall'archivio Examples estratto, `EXAMPLE_2D_PP` è stato generato con `create_example.py 100 0.5 300 false` ed eseguito con l'eseguibile installato: il run è terminato e ha prodotto `OUTPUT/example2D/fields/example2D.nc`.

La documentazione è stata costruita localmente in modalità MkDocs `--strict`: il catalogo ha generato 20 Examples e zero avvisi di contenuto. Doxygen e lo styling sono riusciti e sono stati verificati `index.html`, `examples/index.html` e `doxygen/index.html`. La pubblicazione GitHub Pages è stata completata con successo dal workflow `38900609895` sul commit `320499a3270d7e09d65f29f72506cc52a4812a0e`. C06 resta GATE per l'audit ancora aperto sulla completezza e riproducibilità degli esempi; le incertezze su provenienza/licenza dei DEM non sono un blocco tecnico né una priorità.

L'estensione dello smoke test ha trovato che i generatori di `EXAMPLE_BOX`, `EXAMPLE_BW`, i quattro casi `EXAMPLE_BUMP` e `EXAMPLE_SLOPE` fallivano con NumPy 2.x durante la conversione delle medie cell-centered. I sette generatori sono stati corretti usando esplicitamente la colonna 0 degli array di topografia/altezza/velocità. I README indicano ora le dipendenze Python essenziali e il requisito opzionale per il plotting.

Dopo il commit `60474fe25d372a008059b0caf433ebea0e7fdf12` è stato creato un tag sintetico `9.9.18` soltanto in un clone temporaneo, senza aggiungere tag o release al repository principale. `build_archives.py` ha generato gli archivi sorgente e Examples e il manifest dalla stessa revisione; l'inventario dell'archivio Examples comprende tutti i file tracciati e 20 README standardizzati. SHA-256: sorgente `380302e66c3e16a1e1a3e05d78703c87f0e0e6d2dfeb1dc0f4c70cea6e908eb6`; Examples `8d1c86ce3092cb660c62c59a2e2cb0e1a26058d4e9259f945a5eff2145afd468`. L'archivio sorgente ha superato la generazione Autotools, `configure` e `make dist`; non sono stati ripetuti `make check` e `make distcheck` sul candidato perché il commit non modifica sorgenti o fixture, già validati nel gate C06.

Tutti e sette i generatori corretti sono stati eseguiti con successo, insieme al solver, direttamente dall'archivio Examples `9.9.18` estratto. Anche il preprocessore di `EXAMPLE_ETNA` ha completato la creazione di `DEM/init_volume.asc`. Per `EXAMPLE_ENSEMBLE_ETNA`, il test ha rilevato che `delim_whitespace=True` non è più accettato da pandas 3.0; i parser nei due preprocessori Etna e negli helper di campionamento sono stati aggiornati a `sep=r"\s+"`. Dall'archivio estratto sono stati generati 20 campioni e 20 directory. Tutti i 20 membri hanno raggiunto `T_END=50` e prodotto gli output ESRI/CSV e i restart; `create_prob_maps.py` ha generato 31 raster `.asc` e `plot_ensemble.py` 31 PNG. Il template ha `OUTPUT_CONS_FLAG=F`, quindi non produce NetCDF per default; l'hazard mapping opzionale non è stato eseguito perché netCDF4 e GDAL non sono disponibili nell'ambiente.

La prova dell'ensemble ha inoltre rivelato che `launch_jobs.py` terminava prima di attendere l'ultimo gruppo di processi e non propagava gli exit code dei run. La routine è stata corretta per attendere tutti i processi e segnalare i membri falliti; test isolati con otto job sintetici hanno verificato sia l'attesa dell'ultimo gruppo sia la propagazione di un exit code nonzero.

Dopo il commit `8d776f8f7065179fc920eab7bc12c28169321e51` è stato creato un nuovo candidato con tag sintetico `9.9.19` soltanto in un clone temporaneo. Il manifest associa gli archivi sorgente (`SHA-256 4c1139ed95d98f33063e29e2b20f1682c75a6f6764ba4511c4e83fcd56acc635`) ed Examples (`SHA-256 99e6cf9d2d874540dab77ee8c8cbeeb6003018be282d142203e074b18291aad1`) alla stessa revisione; l'archivio Examples include la correzione del launcher. Dall'archivio estratto, `launch_jobs.py` ha superato i test isolati di attesa dell'ultimo gruppo e propagazione dell'errore di un membro. La suite completa di 20 run non è stata ripetuta sul candidato `9.9.19`: le simulazioni e le mappe sono già state verificate dal precedente archivio `9.9.18`, e il nuovo commit cambia il coordinamento dei processi, non il solver o i dati di input. Nessun tag o release ufficiale è stato creato.

#### Prove runtime addizionali dall'estrazione Examples — 16 ottobre 2026

La sorgente sotto `src` dell'eseguibile installato è identica a quella di `07b740f7`; le prove sono state lanciate da un'estrazione temporanea pulita degli Examples a quel commit. `EXAMPLE_2D`, `EXAMPLE_BW2D`, `EXAMPLE_ENGWELL2016`, `EXAMPLE_TAAL` e `EXAMPLE_TAAL_BOTTOM` hanno raggiunto il rispettivo `T_END` e prodotto NetCDF non vuoti. `EXAMPLE_CF` non ha completato: dopo circa nove minuti il tempo stampato era fermo a 105.584 s con `dt` arrotondato a 0.00000; il run è stato interrotto senza raggiungere `T_END`.

Nelle prove del 16 ottobre, la configurazione documentata di `EXAMPLE_DAMBREAK` con frazione di pressione interstiziale 1.0 generava NaN e terminava con errore fatale per spessore negativo a circa 0.36 s. Anche 0.5 falliva a circa 0.49 s; 0.0 raggiungeva 3 s, ma non rappresentava lo stesso stato iniziale pressurizzato. Con frazione 1.0, la pressione iniziale cancellava esattamente `rho_m * g_red * h` e il codice portava lo sforzo efficace a zero. Disabilitare separatamente la permeabilità dinamica o la perdita di gas non evitava il guasto; neppure un test diagnostico con reologia 6 e coefficiente d'attrito valido. Le prove indicavano il percorso pressione/sforzo efficace come area interessata, senza isolare il meccanismo numerico preciso. Il generatore stampava inoltre la pressione idrostatica completa prima di applicare la frazione richiesta. Il run standard di `EXAMPLE_LATERAL_SOURCE` terminava con codice 0 ma calcolava gravità ridotta negativa, parametri sorgente NaN e massa nulla; una prova diagnostica con frazione massica 0.8 rendeva finiti i parametri derivati ma continuava a produrre massa nulla. Non era stata sostituita nessuna configurazione fisica.

Il generatore documentale ha prodotto le 20 pagine Examples con zero avvisi di contenuto e la build MkDocs `--strict` con la configurazione generata dal workflow è riuscita.

Esito C06 al 16 ottobre: **GATE invariato**. Le prove runtime avevano verificato altri cinque casi, ma l'audit integrale dei 19 README storici non era concluso. I problemi di riproducibilità di `EXAMPLE_CF`, `EXAMPLE_DAMBREAK` e `EXAMPLE_LATERAL_SOURCE` erano ancora aperti sulla versione allora provata.

#### Aggiornamento solver spaziale e validazione della massa — 20 ottobre 2026

I commit `e0c8c616`, `e94463c9` e `0920761e` hanno introdotto il flusso outflow-only legacy in HP-PCCU, adottato l'operatore spaziale KT e corretto la perdita di massa alla frontiera della sorgente radiale nella ricostruzione KT. Il solver al commit `0920761e` è stato ricompilato prima dei test.

Con il template standard di `EXAMPLE_DAMBREAK` (`KT_TRACE_RECONSTRUCTION='froude_blend'`, pressione interstiziale e perdita di gas attive, flag delle condizioni al contorno pari a 1), il run ha raggiunto `T_END=3 s` e prodotto il NetCDF senza il precedente runaway. La massa iniziale era circa 17.4059 kg e nei log è rimasta stabile a circa 17.41 kg. Questo esito aggiorna e supera la diagnosi del 16 ottobre per DamBreak: il fallimento osservato allora non si riproduce con il solver corrente.

`TEST_SOURCE_MASS` è passato con tolleranza relativa `1e-9` in tutti e cinque i casi: `bottom_intervals`, `bottom_window`, `fissure_map_coordinates`, `fissures` e `lateral_window`. Quest'ultimo include la frontiera della sorgente radiale. L'errore relativo massimo osservato è stato `1.30e-11` (`fissure_map_coordinates`). Il test è stato eseguito con l'eseguibile compilato dal medesimo commit e l'ambiente Python `corso_geo` (NumPy 2.0.2).

Non è stata eseguita la suite completa del repository. C06 resta **GATE**: l'audit di tutti gli esempi e dei README storici non è concluso; `EXAMPLE_CF` non è stato ritestato con il solver corrente e l'esito di `TEST_SOURCE_MASS` non dimostra da solo la riproducibilità di `EXAMPLE_LATERAL_SOURCE`. Le evidenze dei run e dei test sono state conservate in directory temporanee, non nel repository.

#### Aggiornamento prove runtime e packaging — 20 ottobre 2026

Le prove sono state ripetute con il solver compilato dal commit corrente `7f05b9b2` e un'estrazione pulita dell'archivio Examples sintetico `9.9.21`. `EXAMPLE_CF` ha raggiunto `T_END=150 s` (NetCDF non vuoto di circa 25 MB), `EXAMPLE_TAAL` `T_END=200 s` (circa 172 MB) e `EXAMPLE_TAAL_BOTTOM` `T_END=50 s` (circa 39 MB). Questi risultati aggiornano le osservazioni precedenti: il collasso del passo temporale di `EXAMPLE_CF` non si riproduce con il solver corrente. `EXAMPLE_2D` ha inoltre completato lo smoke test dall'archivio Examples estratto.

Il primo candidato `9.9.21` aveva omesso `src/python_diagnostics/kt_ritter_1d.f90`; l'inventario lo ha rilevato e la correzione locale a `EXTRA_DIST` è stata inclusa nel candidato successivo. Dopo gli aggiornamenti documentali è stato generato il candidato temporaneo `9.9.24` (revisione sintetica `40e48589c5f1db9c09c351073083cbd7c5f9063e`; nessun tag o commit aggiunto al repository principale). Il builder ha verificato gli inventari di entrambi gli archivi e tutti i 20 README. Il manifest registra SHA-256 `0d75ba496143a34d0a9458c3f103fae819b288e997ddd535e8172983fe6781e9` per l'archivio sorgente e `f47e3f18bfdf8abf51770c15be83134021c929ef3df2c40ab8aeca4f9c777215` per Examples. Il tarball sorgente estratto ha superato `configure`, `make`, `make check`, `make distcheck` e `make install` con GNU Fortran e NetCDF-Fortran su macOS.

Dall'archivio Examples `9.9.24` estratto, `EXAMPLE_2D` è stato generato e avviato con il solver installato; il run ha raggiunto 1.5 s e scritto un NetCDF non vuoto di circa 5.1 MB. Il generatore del catalogo ha prodotto 20 Examples senza avvisi di contenuto e la build MkDocs generata dal workflow è riuscita con `--strict`. Al 20 ottobre il confronto completo dei 19 README legacy restava da completare; è stato chiuso nell'aggiornamento del 21 ottobre riportato di seguito.

`EXAMPLE_LATERAL_SOURCE` continua a non dimostrare un flusso funzionante: l'invocazione documentata termina a 100 s con codice 0 ma produce gravità ridotta negativa, proprietà della sorgente non finite e massa nulla. Una prova diagnostica a 300 K rende finite le proprietà derivate (gravità ridotta circa 0.098 m/s², velocità circa 18.76 m/s, spessore circa 1794.6 m e portata richiesta 10,000,000 kg/s), ma la massa simulata resta nulla. Non è stata approvata né introdotta una configurazione fisica sostitutiva. MkDocs `--strict` e il confronto completo, sezione per sezione, dei 19 README legacy non sono ancora stati ripetuti/completati sul worktree corrente.

Esito aggiornato C06: **GATE invariato**. Sono stati chiusi i run correnti di CF e TAAL, ma restano da chiarire la mancata iniezione di massa nell'esempio laterale, rigenerare/verificare l'archivio Examples aggiornato, ripetere la build documentale strict e completare il confronto dei README legacy. Provenienza e licenze dei DEM restano non verificate e non sono un blocco tecnico.

#### Correzione della sorgente laterale e regressione — 20–21 ottobre 2026

La causa della mancata iniezione con `froude_blend` era la sostituzione dello spessore prescritto alla frontiera con una traccia Froude ricostruita nulla. La correzione preserva ora lo spessore della sorgente nei limiti delle velocità, nel flusso avvettivo e nello stato del path idrostatico. La verifica dei quattro orientamenti ha inoltre rilevato e corretto due problemi associati: la ricostruzione dei lati sud/nord non ripristinava la traccia prescritta e la geometria della sorgente nord veniva assegnata al bordo sud.

`TEST_SOURCE_MASS` è stato esteso a nove casi di bilancio e passa con tolleranza relativa `1e-9`, inclusi gli ingressi ai quattro bordi; il massimo errore relativo dei casi di bordo è `1.51e-11`. La build del solver dal worktree corrente è riuscita con GNU Fortran e NetCDF-Fortran. La simulazione diagnostica di `EXAMPLE_LATERAL_SOURCE` a 300 K ha raggiunto 100 s e scritto un NetCDF non vuoto, con massa non nulla. La ricetta documentata a 900 K resta invalida: la gravità ridotta è negativa. La temperatura diagnostica non è una sostituzione fisica approvata.

L'ispezione del percorso di inizializzazione ha chiarito il meccanismo della configurazione a 900 K: la densità della miscela sorgente scende sotto la densità dell'aria ambiente, rendendo negativa la gravità ridotta `g' = (rho_m - rho_ambient) / rho_m * g`. Per il ramo con portata massica e Richardson positivo, il codice deriva lo spessore tramite `sqrt(RI_SOURCE / g')`; prima del controllo aggiunto il 20 ottobre, la radice di un valore negativo generava proprietà non finite e il solver poteva comunque terminare con codice 0. Non è stata scelta una temperatura alternativa: la scelta fisica resta aperta.

Al momento della registrazione del fix della sorgente laterale, C06 restava **GATE** in attesa del confronto integrale dei 19 README legacy, della verifica degli esempi/dipendenze di post-processing rimanenti e dell'accettazione scientifica o correzione della configurazione documentata a 900 K. Il confronto README è stato completato nell'aggiornamento immediatamente seguente; gli altri gate restano aperti. La provenienza/licenza dei DEM resta non verificata ma non è un blocco tecnico.

#### Confronto sezionale dei README legacy — 21 ottobre 2026

Il confronto sezione per sezione di tutte le 19 coppie `README.txt`/`README.md` è completo. Sono stati verificati setup, parametri, comandi, output, riferimenti e dipendenze; i 19 file storici sono rimasti invariati. Le discrepanze con la configurazione corrente sono state riportate nelle guide Markdown e in `EXAMPLES/README_SECOND_PASS_AUDIT.md`, invece di riprodurre indicazioni legacy non valide. In particolare:

- `EXAMPLE_DAMBREAK`: il README legacy dichiara pareti laterali impermeabili e bordo destro senza flusso, ma il template usa `FLAG=1` per i momenti e il solver corrente applica condizioni outflow-only.
- `EXAMPLE_BW2D`: ripristinata la possibilità documentata di una sorgente stazionaria rimuovendo `TIME_PARAM`.
- I README sorgente descrivono esplicitamente le discrepanze verificate su numero/unità degli argomenti, flag di output e semantica dei parametri temporali; non sono stati promossi errori legacy a istruzioni correnti.

Il confronto documentale è quindi chiuso. Dopo l'aggiornamento, il generatore ha prodotto 20 pagine senza avvisi e MkDocs `--strict` è terminato con successo; il messaggio informativo su Material for MkDocs 2.0 non è bloccante. Il gate C06 resta **GATE**: la ricetta 900 K è ora rifiutata esplicitamente ma resta fisicamente invalida e richiede parametri approvati; sono inoltre aperte la validazione di esempi/dipendenze di post-processing non ancora provati e la generazione/verifica di un nuovo candidato Examples dal worktree corrente. La provenienza/licenza dei DEM resta non verificata e non blocca la validazione tecnica.

#### Validazione esplicita della sorgente a gravità ridotta non positiva — 20 ottobre 2026

Il parser ora richiede gravità ridotta positiva quando una sorgente laterale/radiale usa un Richardson prescritto. Se `RI_SOURCE >= 0` e `g' <= 0`, stampa un errore esplicito ed esegue `ERROR STOP` prima di calcolare lo spessore o iniziare il run. Il test di regressione riproduce la composizione documentata a 900 K e verifica il codice d'uscita nonzero, il messaggio diagnostico e l'assenza della directory di output.

Il solver è stato ricompilato fuori dal repository con GNU Fortran e NetCDF-Fortran; `TEST_SOURCE_MASS` passa con nove bilanci di massa e il decimo test di rifiuto input atteso. Il confronto massimo di massa resta `1.51e-11`, sotto la tolleranza `1e-9`. Non sono stati modificati temperatura o altri parametri fisici dell'esempio.

#### Correzione della composizione di `EXAMPLE_LATERAL_SOURCE`

La decisione precedente di usare 300 K con `XS_SOURCE=0.01` è superata: l'utente ha chiarito che la concentrazione originale era 1% di solidi in volume, non in massa, e ha scelto di mantenere `T_SOURCE=900 K`. Il generatore ora scrive il sesto argomento in `ALPHAS_SOURCE`, che rappresenta la frazione volumetrica; `XS_SOURCE` resta la frazione massica. Con `PRES=101300 Pa`, `RHO_S=2000 kg/m³` e aria a 900 K, il gas ha densità circa `0.392 kg/m³`; la miscela al 1% in volume ha densità circa `20.39 kg/m³`. Rispetto all'aria ambiente a 300 K (`1.176 kg/m³`), la gravità ridotta è positiva, circa `9.24 m/s²`. Il caso a 900 K con `XS_SOURCE=0.01` rimane correttamente invalido, ma descrive una composizione diversa (1% in massa).

La guida attiva è stata corretta e la suite `TEST_SOURCE_MASS` è stata estesa con un caso positivo a 900 K e `ALPHAS_SOURCE=0.01`, mantenendo anche il caso di rifiuto della composizione a 1% in massa. La build corrente del solver e il run dell'esempio sono riusciti: il solver riporta densità sorgente `20.3882 kg/m³`, gravità ridotta `9.2440 m/s²`, raggiunge `T_END=100 s` e scrive un NetCDF non vuoto; le masse finali sono circa `2.287e8 kg` per la miscela e `2.025e8 kg` per i solidi. `TEST_SOURCE_MASS` passa con dieci bilanci e il caso atteso-fallimento; il nuovo caso caldo in volume ha errore relativo `3.01e-13`. Il catalogo genera 20 pagine senza avvisi e MkDocs `--strict` passa. Le evidenze dell'archivio `9.9.25` sopra riportate precedono questa modifica e non ne verificano il contenuto; il nuovo candidato e le relative prove sono registrati nella sezione seguente.

#### Validazione aggiornata del candidato distributivo — 21 ottobre 2026

Questo aggiornamento sostituisce lo stato dei gate riportato nelle note precedenti. Gli archivi sorgente ed Examples `9.9.25` sono stati generati in un clone temporaneo dalla revisione `63db5d73341c58fd5e824d5ba9d362edefeda21c`, usando un tag sintetico non pubblicato. I checksum SHA-256 sono rispettivamente `99dc7f7087cd96b3bd6d9779389443de562f077320779ecf0802e85833e93093` e `57f590929cfcc1277c1932767d5e600062c1162272fced48141a414c4e30a242`. Dal tarball sorgente estratto sono passati `configure`, build, `make check`, `make distcheck` e `make install` con GNU Fortran e NetCDF-Fortran; l'eseguibile installato ha completato lo smoke test `EXAMPLE_2D` a 1.5 s e la suite distribuita `TEST_SOURCE_MASS` (nove bilanci a tolleranza `1e-9` e caso atteso-fallimento).

Tutti i 20 esempi dell'archivio Examples estratto hanno raggiunto il rispettivo tempo finale configurato e prodotto output non vuoto. `EXAMPLE_ETNA` ha raggiunto 200 s, `EXAMPLE_ETNA_LAVA` 3600 s e l'ensemble Etna ha completato tutti i 20 membri a 50 s per membro. `create_prob_maps.py` e `plot_ensemble.py` sono terminati con successo, producendo 93 raster di probabilità e 93 PNG in tre gruppi di risultati. Il catalogo aggiornato ha generato 20 pagine senza avvisi e MkDocs `--strict` è terminato con successo; l'avviso su Material for MkDocs 2.0 è informativo.

Al momento della verifica del candidato `9.9.25`, le Actions riuscite più recenti appartenevano a `bd4a87d8a8a7f9876950a7aa3024b41e8c7612d6`; quello stato è storico e viene superato dalla matrice registrata nella sezione `9.9.27`. Il post-processing opzionale delle mappe di superamento non è stato verificato perché nell'ambiente selezionato mancano `netCDF4` e i binding GDAL `osgeo`; questa attività è distinta dal gate obbligatorio di solver. La provenienza/licenza dei DEM resta non verificata e non è un blocco tecnico.

#### Nuovo candidato distributivo `9.9.27` — 21 ottobre 2026

Gli archivi sorgente ed Examples sono stati generati in un clone temporaneo isolato dal commit sintetico `152781f8f4f0532da6523b53d3b938ef7eec40d6`, basato sul commit del repository `6b7120b9` e comprendente le correzioni delle date. Il tag `9.9.27` esiste solo nel clone temporaneo; non sono stati creati tag o release nel repository principale. Il manifest registra SHA-256 `47cb9e5d5196c7b1b99a5ecd1dfb10817a9ade32acc5bce9d248d7fb49ed1fc5` per l'archivio sorgente e `1c7157c1cd2a46dc10340a1caffb29b58e60dbc96efa6457523d319626a0ccbb` per l'archivio Examples.

Dall'archivio sorgente estratto sono passati `configure`, compilazione, `make check`, `make distcheck` e `make install` con GNU Fortran e NetCDF-Fortran. La suite distribuita `TEST_SOURCE_MASS` ha superato i dieci bilanci di massa e il caso atteso-fallimento a tolleranza relativa `1e-9`; l'errore massimo osservato è `1.51e-11`. `TEST_PCCU_GATE_H` ha superato tutti e sei i confronti full-field e i relativi budget di massa usando l'eseguibile e gli oggetti della build estratta.

Dall'archivio Examples estratto, `EXAMPLE_2D` ha raggiunto 1.5 s e prodotto un NetCDF non vuoto. `EXAMPLE_LATERAL_SOURCE`, con 900 K e `ALPHAS_SOURCE=0.01`, ha raggiunto 100 s e prodotto un NetCDF non vuoto; le masse finali sono circa `2.287e8 kg` per la miscela e `2.025e8 kg` per i solidi. Il catalogo documentale ha generato 20 pagine senza avvisi e MkDocs `--strict` è passato. La suite completa di 20 esempi era stata eseguita sul candidato precedente `9.9.25` e non è stata ripetuta sul `9.9.27`.

La matrice GitHub Actions è completata e accettata per C06: Solver (run `39563559928`) e Intel ifx (run `39563583552`) sono passati sul commit `c3b47cbf9a07162f1f7cbaf8b26230c05ef87798`; Docker (run `39567697006`) ha superato tutti e quattro i job sul commit `cf34b8d13e1739f0b4113b30f8365a2253f547e3`. Il commit Docker aggiunge soltanto la dipendenza Automake necessaria a ordinare la compilazione del modulo usato da `init_2d`; l'utente accetta i risultati delle Actions sui due commit consecutivi, senza richiedere un nuovo run Solver/ifx. C06 è quindi **COMPLETATO** per la validazione distributiva. Il candidato `9.9.27` resta sintetico e non pubblicato; il post-processing opzionale delle mappe di superamento resta non verificato e non è un blocco tecnico.

#### Regola di rilascio per documentazione pubblica

Pubblicazione unificata: solo il workflow publish-documentation.yml deve distribuire contenuti su GitHub Pages; non creare un secondo workflow Doxygen che possa sostituire la root MkDocs. Eseguire manualmente o da tag/release secondo il workflow tracciato. La homepage resta a /, il catalogo in /examples/ e l’HTML Doxygen in /doxygen/. La documentazione deve essere costruita dalla medesima revisione del sorgente e il catalogo legge direttamente EXAMPLES/\*\*/README.md. Testare tema chiaro, tema scuro, stampa e risoluzione dei link quando si modifica il CSS.

## Fase D — Reologie e solver IMEX

### D01 — Inventario formule e rami REAL/COMPLEX

**Test proporzionati (nuove modifiche a questo step):** T1 inventario branch REAL/COMPLEX, formule e proprietari; nessun run numerico se l’intervento è solo documentale. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO (inventario statico; caratterizzazione D02 implementata con gate aperto)

#### Evidenza D01 — 21 ottobre 2026

Audit della revisione `65888eeceacbf78497bb3f7d361dc5e598d8ef2d`, documentato in
[D01_RHEOLOGY_INVENTORY_2026_10_21.md](acceptance/D01_RHEOLOGY_INVENTORY_2026_10_21.md)
e nel relativo [inventario JSON](acceptance/results_d01_inventory_2026_10_21.json).
Sono inventariati selettori 0–12, formule e guardie dei termini impliciti e
semi-impliciti, sorgenti termiche e di pressione interstiziale, kernel condivisi
di conversione REAL/COMPLEX, dipendenze di input/stato e percorso SI–Newton–RK.
Il manifest registra hash di 17 file, riferimenti a 19 procedure e dieci punti
di attenzione per D02 e le successive estrazioni. I proprietari riportati sono
quelli attuali per responsabilità: non sostituiscono il manifest per simbolo D03.

Il modello 7 non ha una legge implementata ed è rifiutato dal parser quando la
reologia è attiva; non gli viene assegnato un nome fisico o un modulo futuro.
Il modello 8 integra analiticamente un'ODE usando `dt`, non risolve l'equazione
algebrica dello stage DIRK; la rimozione D04 e il solver ANALYTIC D13 restano
attività separate. I termini SI sono solo REAL; l'API implicita REAL e quella
COMPLEX eseguono la stessa algebra COMPLEX, con conversioni di stato condivise.

Verifica T1 documentale: provenienza/hash, riferimenti alle procedure e ai rami,
copertura dei selettori e dichiarazioni dei parametri, struttura JSON, link
locali e whitespace. Nessun sorgente o input modificato, nessuna simulazione,
compilazione o CI rieseguita. La chiusura riguarda l'inventario D01 e non attesta
parità numerica, correttezza dei Jacobiani o equivalenza di futuri refactoring.

#### M12e — contratto storico condiviso fra D01 e D02

Il vecchio contratto M12e riuniva inventario e caratterizzazione numerica.
L'inventario formule/rami/dipendenze è coperto da D01; le prove di valore e
derivata seguenti appartengono a D02, e l'equivalenza durante l'estrazione ai
successivi gate D05–D15. La chiusura statica non anticipa questi esiti.

**Attività previste:** Caratterizzare sorgenti per ciascuna reologia, temperatura e pressione interstiziale su stati REAL/COMPLEX; confrontare Jacobiano complex-step e finite difference controllata; coprire U→0, h→0, soglie, ABS/MAX e branch non olomorfi; conservare formule condivise e test modello-per-modello prima dell’estrazione.

**File interessati:** src/physics/equation_terms_2d.f90; src/state/state_conversion_2d.f90; src/core/complexify.f90; src/state/\*.inc; src/io/inpout_2d.f90; src/numerics/nonlinear_solver_2d.f90; src/numerics/time_integration_2d.f90; src/rheology/\*.f90 futuri; TESTS/TEST_RHEOLOGY_FROUDE/\* e nuovi TEST_RHEOLOGY\_\*.

**Criteri di accettazione:** Parità numerica reale/complessa per punti regolari; differenze delle derivate spiegate nei punti non lisci; nessuna regressione scientifica alla migrazione delle formule.

#### Esito dei chiarimenti preliminari

I tre chiarimenti preliminari sono ora coperti dall'inventario:

Modello 7: non supportato e privo di una formulazione identificabile nel codice corrente.

Reologie analitiche: il percorso speciale 8 è un'integrazione ODE; non esiste un dispatcher ANALYTIC per gli stage DIRK. Le riduzioni scalari possibili non costituiscono implementazioni già validate.

Dipendenze REAL/COMPLEX: sono mappati algebra implicita condivisa, include di conversione e routine solo REAL. Sono registrate le proiezioni REAL che perdono sensibilità, le guardie non lisce e le scritture globali da verificare prima dell'estrazione.

I rilievi statici non autorizzano correzioni silenziose delle formule durante il refactoring.

### D02 — Caratterizzazione numerica cell-local

**Test proporzionati (nuove modifiche a questo step):** T1 test cell-local di valore/derivata e finite/complex step, casi degeneri; T2 integrazione della sorgente; T3 tolleranze e confronti indipendenti se cambiano formule. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: IMPLEMENTATO, GATE APERTO

#### Evidenza D02 — 21 ottobre 2026

Caratterizzazione della baseline D01 `5aee5d9e`, senza modificare sorgenti
Fortran produttivi, input o tolleranze. Harness riproducibile in
`TESTS/TEST_RHEOLOGY_CELL_LOCAL`, rapporto
[D02_CELL_LOCAL_CHARACTERIZATION_2026_10_21.md](acceptance/D02_CELL_LOCAL_CHARACTERIZATION_2026_10_21.md)
e [risultati sintetici](acceptance/results_d02_cell_local_2026_10_21.json).
I nuovi check sono distribuiti come suite autonoma: non vengono aggiunti a
`make check` o CI come se i fallimenti produttivi fossero già risolti.

La matrice locale GNU/macOS conta 249 fixture per profilo: **221 PASS / 28 FAIL
strict**, **213 PASS / 36 FAIL optimized con trap aritmetici**. I valori finiti
delle sorgenti concordano con gli oracoli indipendenti entro `6.53e-12`; lo
scaling delle colonne complex-step normalizzate concorda entro `1e-13`.
Passano tutti i 16 bilanci isolati di perdita di gas per profilo. Passano anche
12 prove con l'integratore produttivo (`n_RK=2`, optimized/strict-debug,
1/4 thread) per quadratic DIRK, arresto Coulomb e lava accoppiata. I check del
parser reale passano 5/10 casi: cinque input con parametri mancanti/nulli sono
accettati invece di essere respinti.

#### Requisiti ancora scoperti e chiusura del gate

Restano sei rilievi documentati nel rapporto: (1) 20 fixture con derivate
CS/finite difference discordanti, incluse colonne termiche usate da Newton;
(2) quattro casi di inizializzazione/dominio della diffusione di pressione;
(3) validazione incompleta dei parametri; (4) eccezioni nel percorso SI a riposo
ottimizzato con trap, da distinguere dai corrispondenti output finiti senza
trap; (5) `T_env` configurato sovrascritto a 300 K; (6) uscita Newton per
correzione piccola che non garantisce la soglia del residuo ricalcolato.
I conteggi includono riproduzioni isolate: quattro difetti SI ottimizzati sono
contati sia nelle fixture combinate sia nelle riproduzioni, non sono otto
difetti fisici distinti.

Per passare a COMPLETATO occorrono correzioni o contratti esplicitamente
giustificati per questi punti, nuova evidenza T1/T2 mirata, equivalenza dei
valori REAL quando si ripristinano le dipendenze COMPLEX e decisione documentata
sui contratti termici, di trapping e di arresto. Non aumentare le tolleranze né
convertire i rilievi in expected PASS. Il confronto dei diversi provider di
permeabilità resta un vincolo D06; l'ownership D03 e la rimozione D04 non sono
anticipate. Nessuna nuova validazione ifx/Linux/Docker o campagna scientifica
estesa è attestata da D02.

Matrice iniziale e rilievi da verificare: sezione «Findings and D02 worklist»
del [rapporto D01](acceptance/D01_RHEOLOGY_INVENTORY_2026_10_21.md).
Separare prove REAL/COMPLEX delle sorgenti implicite, verifiche REAL delle
leggi SI e confronto della mappa di stage corretta. Includere diffusione della
pressione, inibizione, permeabilità mono/multi-classe, stati secchi e
configurazione pore-only, oltre alle soglie e agli accoppiamenti termici.
I test preesistenti da soli non coprivano tutte le leggi e tutte le derivate;
la nuova matrice rende riproducibili le discrepanze rilevate, senza approvarle.

### D03 — Audit e ownership di parametri, variabili e stato condiviso

**Test proporzionati (nuove modifiche a questo step):** T1 audit dipendenze, ownership e accessi globali; T2 inizializzazione e serializzazione; T3 integrazione se cambia stato persistente. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 rifiuto reologia 8, riferimenti e parsing; T2 un esempio con reologia supportata; T3 regressione altre reologie solo se toccate. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 contratti di ogni modulo reologico e REAL/COMPLEX; T2 un caso per modello toccato; T3 confronto pre/post se il refactoring cambia il percorso numerico. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 sorgenti termiche/pore pressure e limiti sottili; T2 casi lava/raffreddamento; T3 energia e stiff RK se cambiano i termini. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: DA FARE

### D07 — Migrazione progressiva dei parametri non reologici

**Test proporzionati (nuove modifiche a questo step):** T1 ownership parametri, inizializzazione e default; T2 input/restart rappresentativi; T3 solo se mutano effetti numerici. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 parsing nomi e compatibilità input/errori; T2 un caso per modello attivo; T3 non necessario per mera rinomina senza cambi di dispatch. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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
| 7                        | Non supportato (audit D01)   | Nessuno                  |
| 8                        | Eliminato                   | Nessuno                  |
| 9                        | froude                      | rheology_froude.f90      |
| 10                       | lucas                       | rheology_lucas.f90       |
| 11                       | mui                         | rheology_mui.f90         |
| 12                       | mui_voellmy                 | rheology_mui_voellmy.f90 |

La tabella costituisce la mappatura di progetto, non un'interfaccia già implementata. D01 conferma le identità delle leggi presenti e l'assenza di una legge supportata per il modello 7: nessun nome o file gli viene assegnato. Il modello 8 è ancora nel codice e resta da rimuovere in D04. Valori, derivate e limiti delle formulazioni saranno caratterizzati in D02 prima della migrazione definitiva.

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

**Test proporzionati (nuove modifiche a questo step):** T1 isolamento frizione e segno dei termini; T2 casi dissipativi/non dissipativi; T3 conservazione e regressione reologica se cambia forza effettiva. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 contratto EX_terms/IM_terms e consistenza delle somme; T2 uno step IMEX con tutte le componenti; T3 casi stiff se cambia partizione. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 metadati e copertura IM_group; T2 selezioni e subset; T3 regressioni RK quando cambia il raggruppamento. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 Jacobiani ridotti, complex step e convergenza Newton; T2 problemi stiff; T3 full/split e profili se cambiano traiettorie. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 variabili ANALYTIC per stage DIRK e derivate; T2 casi con/algebriche senza; T3 N_RK 2/3/4 se cambia stage coupling. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 gruppi misti e fallback; T2 combinazioni rappresentative; T3 confronto delle traiettorie se cambia la strategia implicita. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: DA FARE

### D15 — Acceptance della nuova architettura

**Test proporzionati (nuove modifiche a questo step):** T3/T4 gate architetturale: matrice reologie e IMEX, conservazione, REAL/COMPLEX, restart, baseline A11 e variazioni numeriche autorizzate; decisione scientifica autonoma. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: GATE

#### Manifest, ripartenza e gate

Ogni commit funzionale porta evidenza di compilazione pulita, test associati, SHA e descrizione della tolleranza. Gli input e l’output devono registrare stringa reologica e gruppo/metodo/ordine; i restart non devono ripartire silenziosamente con un diverso modello/gruppo. Le ottimizzazioni OpenMP non possono introdurre workspace Newton condivisi. Conservare i .inc di conversione REAL/COMPLEX e i controlli sulla non analiticità di ABS, MAX, segni e soglie nel complex-step. Prima di G01 eseguire il gate di compatibilità del default, regressioni scientifiche HP-PCCU, storage/restart e accuratezza per gruppi.

## Fase E — Architettura e robustezza del passo

### E01 — Audit dipendenze e ownership stato

**Test proporzionati (nuove modifiche a questo step):** T1 audit import/dipendenze e ownership; T2 inizializzazione e stato condiviso; T3 solo per cambi strutturali con effetti runtime. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 check statici automatici, lint/dep graph e negative fixtures; T2 CI dei contratti; nessuna campagna scientifica solo per check statici. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: DA FARE

#### Prove di estensione (architecture fitness checks)

Prima del via libera, svolgere almeno tre analisi di impatto su casi realistici: (1) aggiunta di una nuova reologia o variazione significativa di una esistente; (2) estensione dei gruppi impliciti/equazioni con metadata e Jacobiano; (3) introduzione di uno strato addizionale o di aggiornamenti geometrici nel tempo. In ciascun caso documentare moduli interessati, interfacce da modificare, dipendenze introdotte, necessità di riallocazione e test da estendere. Sono prove di fattibilità architetturale, non obblighi di implementare anticipatamente le nuove funzionalità.

Per il percorso verso G01–H05 verificare inoltre che l’organizzazione dello stato e della geometria non ostacoli domain decomposition, halo exchange, tile/local domain, MPI, Parallel NetCDF e futuri acceleratori. Non si richiede adesso di implementare MPI o GPU.

### E03 — Architecture readiness

**Test proporzionati (nuove modifiche a questo step):** T3 gate di fitness architetturale e dipendenze; casi end-to-end minimi e decisione motivata. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: GATE

#### Evidenze, esiti e gestione delle non conformità

Produrre docs/acceptance/ARCHITECTURE_READINESS_AUDIT.md e un manifest machine-readable con revisione/commit, matrice PASS/FAIL/WAIVER, link ai workflow, risultati dei benchmark e tracciabilità dei requisiti. Documentare ogni non conformità con impatto su correttezza, estensibilità, prestazioni o manutenzione, responsabile e azione correttiva.

E03 può essere approvato soltanto con D15 già accettato, evidenze complete e nessun problema aperto ad alta criticità che renda onerosi o rischiosi gli sviluppi previsti. Eventuali eccezioni a bassa criticità devono essere motivate, circoscritte e accompagnate da test/mitigazioni. Se il gate fallisce, pianificare correzioni mirate e ripetere le verifiche; evitare una riscrittura generale senza benefici misurabili.

Criterio di uscita: “Il codice e il repository sono sufficientemente organizzati, validati, manutenibili ed estendibili da costituire una piattaforma stabile per i nuovi sviluppi scientifici”. Solo dopo tale approvazione si sblocca il percorso di nuove funzionalità.

### E04 — Policy timestep e input

**Test proporzionati (nuove modifiche a questo step):** T1 priorità e validazione policy Δt; T2 sequenza di passi e output times; T3 CFL/RK rappresentativi se cambia il controller. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: DA FARE

Introdurre enable_step_retry, retry_dt_factor (default candidato 0.5, 0\<factor\<1), max_step_retries (default candidato 10), min_time_step \[s\] configurabile e parametri di diagnostica. min_time_step è il limite assoluto di timestep: nessun tentativo può essere lanciato con Δt inferiore al minimo; se la prossima riduzione lo violasse, la simulazione termina con un errore controllato. Validare valori finiti, consistenza con il timestep iniziale e compatibilità con il controllo CFL; non accettare un timestep che ecceda il limite temporale di output/fine run per effetto del controllore. I default numerici proposti sono da finalizzare durante l’implementazione, non da presentare come soglie fisiche universali.

La policy per il passo successivo rimane quella del calcolo CFL esistente. Un retry modifica solamente Δt del tentativo corrente: non è un controller adattativo di accuratezza basato su stime embedded. Registrare accepted_dt, proposed_dt, rejected_attempts e minimum_accepted_dt nelle diagnostiche; definire i contatori in modo coerente con il restart.

### E05 — Bounds primitive e override

**Test proporzionati (nuove modifiche a questo step):** T1 limiti/precedenza configurazioni e casi invalide; T2 celle dry/wet e varianti reologiche; T3 solo se il bound modifica la dinamica. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: DA FARE

Definire un descrittore centralizzato dei limiti per ciascuna variabile primitiva fisicamente significativa del model_layout; ciascun minimo/massimo deve poter essere attivo o disattivo in modo indipendente. L’input può sovrascrivere soltanto alcune variabili e soltanto uno dei lati dell’intervallo. La precedenza è: limiti generali ampi → override/default del modello reologico o della particolare chiusura attiva → override esplicito dell’utente. Validare i bounds al setup, inclusa l’intersezione fra vincoli di modelli concorrenti: un intervallo vuoto è un errore di configurazione, non un fallimento recuperabile del timestep.

Distinguere condizioni di dominio matematico (p.es. frazioni e proprietà termodinamiche valide) e intervalli prudenziali di controllo fisico/numerico. Consentire tolleranze solo ove compatibili con la definizione del vincolo e le unità della variabile; non usare clipping dello stato come sostituto del rigetto. I default non devono aggiungere vincoli ingiustificati a velocità, densità, pressioni o temperature. Il riferimento ambientale della legge di raffreddamento della lava può essere proposto dalla reologia/chiusura termica come default condizionato alle ipotesi fisiche, non come limite universale di ogni stato termico.

Lo sviluppo comprende un inventario completo delle primitive effettivamente esposte dal codice, con indici e significati ricavati da model_layout e state_conversion_2d; il controllo deve trattare gli stati asciutti secondo una semantica esplicita, senza segnalare falsi errori per variabili convenzionali non definite nel secco. Le componenti non finite (NaN/Inf) devono produrre un fallimento recuperabile o fatale appropriato senza corrompere le diagnostiche.

### E06 — Controlli di ammissibilita cell-local

**Test proporzionati (nuove modifiche a questo step):** T1 ammissibilità dopo stadi/finale, nessun clipping; T2 caso che attiva la guardia; T3 integrazione numerica se influisce sugli step. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: DA FARE

Contratto concettuale: CALL check_primitive_bounds(state, valid, failure_info). La routine non ha effetti collaterali: legge unicamente lo stato primitivo di cella e i bounds risolti, restituisce il flag e un record diagnostico con variabile/indice, (j,k), valore, limite, tipo di violazione e identificatore di stadio o finale. Non stampa, non clippa e non invoca fatal_error. Il chiamante determina il retry e la messaggistica.

Punti di chiamata ESCLUSIVI: (a) dopo il completamento di ciascuno stadio RK, prima dell’uso dello stato di stage nelle successive valutazioni esplicite/spaziali; (b) dopo l’assemblaggio e le conversioni/proiezioni finali, prima dell’accettazione del passo. Eseguire la scansione soltanto su l=1,domain%solve_cells, usando j_cent(l),k_cent(l); non effettuare controlli in ogni conversione conservative→primitive e non controllare automaticamente gli stati ricostruiti alle facce.

time_integration_2d.f90 conserva this%q_rk per ogni cella risolta dopo lo stage. this%qp_rk è aggiornato oggi solo se need_explicit_stage=.TRUE.: la conversione qc_to_qp per lo stato di stage dovrà essere resa disponibile anche negli altri stadi, senza duplicarla quando già necessaria all’operatore esplicito. qp ordinario, nel loop implicito, può contenere la conversione del noto di stage e non va assunto uguale al risultato finale; dopo l’assemblaggio conclusivo qp è invece aggiornato ed è la scelta naturale per il controllo finale. Verificare che conversioni di diagnostica in verbose non alterino le conclusioni. Conservare la semantica halo/dry esistente.

### E07 — Error propagation, rollback e retry

**Test proporzionati (nuove modifiche a questo step):** T1 rollback completo stato/RNG e limiti retry; T2 fallimento recuperabile + ripartenza; T3 traiettorie 1/4 thread se mutano stadi. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: DA FARE

Separare esplicitamente i fallimenti recuperabili (Newton non convergente, solve lineare fallito, stato di cella non ammissibile, NaN/Inf riconducibile al tentativo) dagli errori strutturali (input incoerente, file mancanti, allocazione e incompatibilità di configurazione/restart). Le routine coinvolte dovranno restituire status/failure_info invece di arrestare incondizionatamente il processo per le sole condizioni recuperabili. I fallimenti dei thread OpenMP vanno raccolti senza race e senza I/O concorrente non controllato; interrompere il passo solo dopo una sincronizzazione sicura, evitando sottopassi successivi su workspace parziali.

All’inizio del timestep salvare lo stato persistente necessario a ripetere il tentativo: conservative e primitive accettate, variabili dinamiche del fondo e geometria aggiornata ove modificate dal passo, maschere e liste attive, campi e memoria stocastica, stato RNG, contatori persistenti, eventuali accumulatori diagnostici e stato dei modelli con memoria. Non è obbligatorio copiare workspace ricalcolabili, ma è obbligatorio re-inizializzarli prima del retry. Il tempo, gli indici del passo accettato, i file scientifici e le scadenze di output non avanzano durante un tentativo respinto. I metadati del retry possono essere diagnosticati separatamente senza alterare il percorso di restart riproducibile.

Se un tentativo fallisce: registrare failure_info; ripristinare lo snapshot; calcolare Δt_new=factor·Δt_old; validare min_time_step e limite retry; ripetere tutti gli stadi IMEX dall’inizio. Nessuna modifica alle equazioni di bilancio o ai pesi RK. Il comportamento deve essere deterministico rispetto al tentativo effettivamente accettato, entro la policy di RNG documentata per l’attivazione stocastica.

### E08 — Test del timestep controller

**Test proporzionati (nuove modifiche a questo step):** T1 unità del controller e decisioni accept/reject; T2 scenari provocati e terminal failure; T3 casi stiff rappresentativi. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T3 gate timestep con budget, conservation/restart e benchmark strict/optimized; T4 solo per nuove garanzie scientifiche. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: GATE

## Fase F — Interfacce per chiusure future

### F01 — Interfaccia closure verticali

**Test proporzionati (nuove modifiche a questo step):** T1 schema di input/output vertical closure, limiti e fallback; T2 backend stub; niente test 3D finché l’API è solo neutrale. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: DA FARE

### F02 — Interfaccia subgrid terreno

**Test proporzionati (nuove modifiche a questo step):** T1 contratto DEM fine/coarse e topografia dinamica; T2 backend neutro; T3 solo per chiusure subgrid attivate. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: DA FARE

### F03 — Backend ML/NN opzionale

**Test proporzionati (nuove modifiche a questo step):** T1 feature schema/versione, NaN/OOD e fallback; T2 inferenza stub opzionale; T3 determinismo/batch quando attivo. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: DA FARE

### F04 — Test di neutralita API

**Test proporzionati (nuove modifiche a questo step):** T1/T2 test di neutralità con backend disattivato e stati identici; T3 breve confronto numerico, non una nuova validazione scientifica. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: DA FARE

### F05 — Architecture closure readiness

**Test proporzionati (nuove modifiche a questo step):** T3 gate contract API/portabilità, senza confondere con validazione fisica di ML o closure. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: GATE

## Fase G — Migrazione layer-aware

### G01 — Storage rank-4 con N_LAYERS=1

**Test proporzionati (nuove modifiche a questo step):** T1 indicizzazione rank-4, allocation e N_LAYERS=1; T2 confronto con storage precedente; T3 restart/thread se lo stato cambia. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 halo/activity e schema di restart/I-O layer-aware; T2 full/split; T3 confronto di serializzazione e campi su 1 layer. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T3 freeze architettura, ABI/API, serializzazione e N_LAYERS=1; T4 non sostituisce la validazione fisica I01. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 indici dominio locale/halo e ghost; T2 1 vs partizioni multiple; T3 equivalenza bilanci/campi se cambia ownership. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 selezione tile, attività e bordi; T2 confronti active/full; T3 timing più conservazione se cambiano riduzioni. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 MPI exchange e partizioni; T2 1/2/4 rank; T3 conservazione e regressioni rank-invariant selezionate. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 halo ibridi e parallel NetCDF; T2 MPI+OpenMP e restart; T3 portabilità/rank/thread al gate. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 profiling e correttezza kernel accelerati; T2 benchmark rappresentativi; T3 confronti numerici e performance con hardware documentato. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

**Test proporzionati (nuove modifiche a questo step):** T1 bilanci per layer, limite 1 layer e wet/dry; T2 benchmark a due strati; T4 validazione fisica dedicata. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: FUTURO

### I02 — Scambi interstrato

**Test proporzionati (nuove modifiche a questo step):** T1 scambi interstrato massa/momento/energia; T2 limiti analitici; T4 conservazione e comparatori indipendenti. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: FUTURO

### I03 — Validazione due strati

**Test proporzionati (nuove modifiche a questo step):** T4 gate scientifico a due strati, test convergenza e bilanci; non riutilizzare automaticamente il gate A11 come prova multilayer. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: GATE

### I04 — Closure verticali fisiche/analitiche

**Test proporzionati (nuove modifiche a questo step):** T1 quadrature e ricostruzione profili, frazioni e positività; T2 confronto reference verticale; T4 se cambia closure fisica. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: FUTURO

Definire un’interfaccia vertical_closure capace di ricevere stato depth-averaged e parametri di modello e restituire momenti verticali, fattori correttivi e termini di accoppiamento previsti dalle equazioni, con normalizzazione e unità documentate. La prima priorità scientifica è la distribuzione verticale delle classi particellari, con effetti sui flussi e prodotti non lineari, collegata a distribuzioni di velocità e temperatura quando rilevanti. Preservare la distinzione tra quantità mediate e medie di prodotti; non sostituire una struttura verticale inesistente con semplici campi diagnostici a profilo arbitrario.

Implementazioni alternative: closure analitica, fisica, semiempirica e ML/NN sotto la stessa API. Dati di verifica esterni alla sola simulazione depth-averaged: teoria, esperimenti e calcoli verticalmente risolti utilizzati offline. Requisiti: limiti e positività di concentrazioni, coerenza del bilancio depth-integrated, conservazione, differenziabilità necessaria ai metodi numerici e prestazioni misurate; non assumere che la vecchia vertical_profiles, già rimossa, sia un riferimento da ripristinare.

### I05 — Validazione verticale

**Test proporzionati (nuove modifiche a questo step):** T4 gate verticale con dati indipendenti, accuratezza e limiti di applicabilità. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: GATE

### I06 — Subgrid topografia fine/coarse

**Test proporzionati (nuove modifiche a questo step):** T1 coarse-graining DEM, casi uniformi e invarianza; T2 fine/coarse; T4 validazione su DEM indipendenti. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: FUTURO

Definire una procedura scientifica tracciabile di simulazioni IMEX_SfloW2D su DEM dettagliati e su DEM/griglie grossolane associati: allineamento geometrico, filtraggio/aggregazione della quota, mapping dei campi, medesima fisica e condizioni iniziali, rapporti di risoluzione, trattamento dei bordi e confronto di flussi, runout, deposizione, canali, blocco e direzioni del moto. Identificare il residuo di coarse-graining che la closure dovrà riprodurre, senza fissare aprioristicamente la sua forma come attrito basale.

Separare training/calibrazione da verifica su terreni e condizioni non visti. Prevedere statistiche topografiche multiscala e feature direzionali coerenti con i vincoli di simmetria e con la geometria B_vertex; stabilire se la correzione opera su un flusso, una sorgente fisica o un provider di interfaccia soltanto dopo una caratterizzazione matematica. Ogni correzione deve superare test di conservazione, positivà/ammissibilità, riproduzione del limite fine e convergenza al diminuire del rapporto di coarse-graining.

### I07 — Validazione subgrid

**Test proporzionati (nuove modifiche a questo step):** T4 gate subgrid: conservazione, generalizzazione, scale/rotazioni e casi indipendenti. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: GATE

### I08 — Inferenza ML/NN delle closure

**Test proporzionati (nuove modifiche a questo step):** T1 versionamento feature/weights, fallback e OOD; T2 inferenza riproducibile; T4 benchmark indipendenti quando impiegata nella fisica. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: FUTURO

Prevedere un’interfaccia uniforme tra kernel di chiusura deterministici e modelli ML/NN, comprendente registro dei modelli disponibili, metadata/versioni dei pesi, normalizzazione e dimensione delle feature, unità, dominio valido di input, forma/rank delle uscite e stato di valutazione. Evitare dipendenze dirette del ciclo RK da framework specifici; cache e batch inference devono appartenere a workspace dedicati, senza allocazione per cella o stato globale nascosto. Il percorso senza ML deve mantenere i risultati di riferimento.

Un output ML non finito o fuori dal dominio fisico deve essere intercettato da regole locali della closure e da politiche documentate di errore o fallback, senza clipping silenzioso. L’addestramento è un’attività offline; soltanto l’inferenza è destinata al run del solver. Verificare scalabilità 1/4 thread e futura esecuzione MPI/GPU, riproducibilità su uno stesso backend e compatibilità con restart quando i modelli abbiano memoria.

#### Ambito esplicitamente escluso

Non pianificare coupling runtime bidirezionale o orchestrazione numerica con altri codici di simulazione. Modelli 3D o osservazioni esterne possono alimentare dataset offline per la chiusura verticale; questa attività non crea una dipendenza esecutiva del solver. Le interfacce input, dati, meteorologia e forcing già approvate restano distinte dal coupling fra simulatori.

### I09 — Validazione ML/NN

**Test proporzionati (nuove modifiche a questo step):** T4 gate ML su test non visti, errore, robustezza, costo e fallback; distinguere accuratezza del modello e del solver. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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
