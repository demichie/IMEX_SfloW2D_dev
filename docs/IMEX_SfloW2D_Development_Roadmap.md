# IMEX_SfloW2D — Development Roadmap

**Roadmap integrata e consolidata | 25 ottobre 2026 | baseline tecnica D17: `0e545ac9`**

Questa roadmap definisce lo stato attuale verificabile e il piano di sviluppo verso la visione scientifica finale, con le responsabilità e le interfacce architetturali vigenti. Le sezioni introduttive e il Master Development Plan costituiscono il riferimento generale; le schede delle fasi descrivono gli interventi implementativi e i relativi criteri di accettazione. Le attività implementate ma non ancora accettate sono distinte dalle milestone completate e dai gate scientifici futuri; gli stati si riferiscono al perimetro di evidenza dichiarato.

**Avanzamento STEP successivo (27 ottobre 2026, PR #3; CI da verificare):** abilitati `NEWTON/STEP` mediante il solver locale `solve_rk_step` gia in produzione, con singolo gruppo, partizione locale, `a_diag=1` e `dt` fisico; passato il campo stocastico effettivo dal contesto della simulazione. Abilitati `EXPLICIT/STEP` con il medesimo `eval_source_explicit` usato negli stage e geometria locale, usando l'inizio del passo come tempo di valutazione dell'Eulero esplicito. `ANALYTIC/STEP` continua a usare l'operatore analitico comune. Il mass exchange rimane una transazione post-RK atomica ordinata per ultima. Aggiunta verifica del registry per combinazioni miste `NEWTON/STEP`, `EXPLICIT/STEP`, `ANALYTIC/STAGE`. La presenza delle tre cadenze/metodologie non dimostra ancora l'accuratezza o la stabilita scientifica di tutte le combinazioni: necessarie compilazione CI, regressioni e test quantitativi di STEP.

**Avanzamento STEP (branch PR #3, CI da ripetere):** il registry accetta ora gruppi `ANALYTIC/STEP` per sorgenti che dichiarano la capacita analitica. L'executor post-RK attraversa i gruppi analitici nell'ordine di `GROUP`, riusa `eval_source_analytic` con intervallo `dt`, aggiorna lo stato primitivo e applica successivamente il batch mass exchange `EXPLICIT/STEP` senza dividerlo. Sono stati aggiunti test di configurazione e ordinamento per raffreddamento e pressione interstiziale. `NEWTON/STEP` e `EXPLICIT/STEP` per le altre sorgenti fisiche rimangono non implementati e vengono rifiutati; non dichiarare D18 STEP completo. La semantica di splitting post-RK non coincide con l'integrazione interna agli stage.

**D18 operational baseline — 27 ottobre 2026 (branch `d18-remove-legacy-imex-config`):** la configurazione delle sorgenti usa esclusivamente le triplette fisiche `METHOD/CADENCE/GROUP` nelle namelist delle sorgenti; `initialize_configured_sources` costruisce registry e partizione. `SOURCE_GROUPS`, i campi `IMEX_*_GROUP`, `IMEX_GROUP_METHOD`, `init_problem_param`, `initialize_legacy_registry`, `resolve_imex_group_method`, il modulo `mass_exchange_step_2d` e il wrapper `solve_implicit_stage` sono stati rimossi. Le descrizioni antecedenti nelle sezioni D17 restano **evidenze storiche, non istruzioni operative**. La validazione dell'input corrente e descritta in `docs/acceptance/D17_16_VALIDATION.md`. Le combinazioni `STEP` fisiche non ancora implementate restano esplicitamente respinte; le formule di calcolo non sono state sostituite. I commit di questa revisione documentale richiedono una nuova verifica CI.

**PRIORITA' ASSOLUTA - D18.MX (28 ottobre 2026; DA IMPLEMENTARE).** Il prossimo intervento e' la riorganizzazione completa del mass exchange. La sezione D18.MX e' la specifica architetturale corrente e **SUPERSEDE le prescrizioni organizzative incompatibili di D17 e D18**, inclusi il raggruppamento precoce delle perdite della fase continua, l'uso del limiter aggregato come contratto definitivo e il percorso separato della deposizione lavica. Le sezioni D17/D18 precedenti restano evidenze storiche e requisiti per le parti NON in conflitto (catalogo D17, quattro interfacce per legge, scheduler, contratti temporali e conservazione). Non dichiarare la riorganizzazione gia' implementata.

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

La piattaforma dovrà consentire di selezionare reologie diverse senza modificare il ciclo IMEX o il motore generale di Newton. I modelli reologici sono organizzati in moduli autonomi, selezionati tramite nomi canonici e validati all’inizializzazione. Le implementazioni dovranno conservare le specificità fisiche del modello, inclusi eventuali regimi limite, dissipazione, frenamento e dipendenze termiche, senza imporre una correzione della quantità di moto a sorgenti che non rappresentano esclusivamente attrito.

L’architettura dovrà garantire che l’aggiunta o il confronto di una reologia abbia un perimetro locale e sia accompagnata da caratterizzazione scientifica REAL/COMPLEX, controlli sulle derivate, test di stati degeneri e confronto rispetto alla baseline. L’architettura delle sorgenti fisiche è ora comune anche ai processi termici, di pressione interstiziale, iniezione e scambio di massa; future estensioni costitutive conserveranno le stesse interfacce e richiederanno validazione specifica.

### C1. Sorgenti fisiche, integrazione e trasferimenti conservativi

L’architettura unificata delle sorgenti permette di comporre reologie, termica, pressione interstiziale, iniezioni e scambi di massa attraverso nomi canonici e capacità numeriche dichiarate. Ogni legge espone quattro interfacce fisiche uniformi; la selezione dei metodi ammissibili è validata in fase di configurazione. Un gruppo Newton viene trattato come un unico sistema accoppiato, anche quando coinvolge più equazioni fisiche; le correzioni numeriche della quantità di moto sono proprietà dei blocchi di equazioni, non della scelta reologica. Il raffreddamento termico è attivabile indipendentemente dal modello reologico, rispettando i limiti energetici e la struttura depth-averaged.

Gli scambi fra flusso e letto sono valutati come una transazione esplicita post-RK coerente: i contributi di erosione, deposizione, perdite della fase vettore, entrainment e vaporizzazione condividono limiti e bilanci. La perdita indotta dalla deposizione è parte della legge di deposizione; i rate conservativi, gli inventari e il fondo nodale Q1 sono aggiornati una sola volta. La separazione fra formulazione fisica, assemblaggio, integrazione temporale e gestione dello stato costituisce l’interfaccia sulla quale innestare anche i processi scientifici futuri senza imporre dipendenze reologiche al solver generale.

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
| D02 | M12e-B        | Caratterizzazione numerica cell-local                       | COMPLETATO (gate locale GNU/macOS) |
| D03 | —             | Audit e ownership di parametri, variabili e stato condiviso | IN PAUSA |
| D04 | P2            | Eliminazione reologia 8                                     | COMPLETATO |
| D05 | P3            | Estrazione moduli reologici e selezione testuale            | COMPLETATO |
| D06 | P3            | Sorgenti termiche e pore pressure                           | COMPLETATO |
| D07 | —             | Migrazione progressiva dei parametri non reologici          | IN PAUSA |
| D08 | P4            | RHEOLOGY_MODEL da numeri a stringhe                         | ASSORBITO IN D05 |
| D09 | P5            | Correzione del solo attrito                                 | COMPLETATO (gate locale GNU/macOS) |
| D10 | P6            | EX_terms / IM_terms                                         | COMPLETATO (gate locale GNU/macOS) |
| D11 | P7            | IM_group e metadati                                  | COMPLETATO (gate locale GNU/macOS) |
| D12 | P7            | Newton per gruppo e Jacobiano ridotto                       | IN CORSO |
| D13 | P8            | Integrazione analitica esatta dei termini IM                | COMPLETATO (gate locale GNU/macOS) |
| D14 | P8            | Gruppi misti e validazione                                  | COMPLETATO (gate locale GNU/macOS) |
| D16 | —             | Consolidamento leggi fisiche, rate e assemblaggio IMEX      | COMPLETATO (gate locale + CI ifx completa) |
| D17 | —             | Sorgenti canoniche, gruppi omogenei, Newton accoppiato e mass exchange STEP | IMPLEMENTATO (CI verde; gate scientifico D15 aperto) |
| D15 | P9            | Acceptance finale della nuova architettura                  | GATE     |


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

File interessati: src/physics/mass_exchange/; src/numerics/step_source_executor_2d.f90; src/geometry_2d.f90; TESTS/TEST_MASS_EXCHANGE/\*; docs/acceptance/N8_A_CELL_EXCHANGE_2026_10_10.md; docs/acceptance/N8_B_PROJECTION_2026_10_10.md.

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

Criteri di accettazione: SODDISFATTI per il perimetro scientifico autorizzato e la revisione congelata. D15 resta un gate distinto, successivo alla relocation, alle reologie, al solver IMEX per gruppi e al consolidamento D16; deve trattare separatamente equivalenza strutturale e cambiamenti numerici autorizzati.

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

**File interessati / nuovi:** src/subgrid/terrain_statistics_2d.f90 \[FUTURE\]; subgrid_topography_2d.f90 \[FUTURE\]; subgrid_closure_2d.f90 \[FUTURE\]; src/geometry/geometry_2d.f90; src/spatial/spatial_operator_2d.f90; src/physics/mass_exchange/; src/numerics/step_source_executor_2d.f90; src/io/\*; TESTS/TEST_SUBGRID_TOPOGRAPHY/\* \[FUTURE\].

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

**File interessati:** src/physics/equation_terms_2d.f90; src/state/state_conversion_2d.f90; src/core/complexify.f90; src/state/\*.inc; src/io/inpout_2d.f90; src/numerics/nonlinear_solver_2d.f90; src/numerics/time_integration_2d.f90; src/physics/rheology_laws/\*.f90; TESTS/TEST_RHEOLOGY_FROUDE/\* e nuovi TEST_RHEOLOGY\_\*.

**Criteri di accettazione:** Parità numerica reale/complessa per punti regolari; differenze delle derivate spiegate nei punti non lisci; nessuna regressione scientifica alla migrazione delle formule.

#### Esito dei chiarimenti preliminari

I tre chiarimenti preliminari sono ora coperti dall'inventario:

Modello 7: non supportato e privo di una formulazione identificabile nel codice corrente.

Reologie analitiche: il percorso speciale 8 è un'integrazione ODE; non esiste un dispatcher ANALYTIC per gli stage DIRK. Le riduzioni scalari possibili non costituiscono implementazioni già validate.

Dipendenze REAL/COMPLEX: sono mappati algebra implicita condivisa, include di conversione e routine solo REAL. Sono registrate le proiezioni REAL che perdono sensibilità, le guardie non lisce e le scritture globali da verificare prima dell'estrazione.

I rilievi statici non autorizzano correzioni silenziose delle formule durante il refactoring.

### D02 — Caratterizzazione numerica cell-local

**Test proporzionati (nuove modifiche a questo step):** T1 test cell-local di valore/derivata e finite/complex step, casi degeneri; T2 integrazione della sorgente; T3 tolleranze e confronti indipendenti se cambiano formule. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO (caratterizzazione e gate locale GNU/macOS; validazione cross-platform non attestata)

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

#### Correzioni ed evidenza di chiusura

I sei rilievi della caratterizzazione iniziale sono stati corretti o chiusi con
contratti espliciti. Le derivate delle sorgenti sono state verificate contro
differenze finite; i domini degeneri della diffusione di pressione e i
parametri non finiti o singolari sono controllati; il percorso SI a riposo,
`T_env` configurato e il criterio di convergenza Newton rispettano i contratti
documentati nel rapporto
[D02_NUMERICAL_FIXES_2026_10_21.md](acceptance/D02_NUMERICAL_FIXES_2026_10_21.md).
Le 257 fixture per profilo, le 793 colonne di derivata, i 28 casi del parser e
i sei controlli di integrazione per profilo passano senza aumentare le
tolleranze. I valori REAL comuni restano equivalenti entro gli scarti riportati
nel follow-up; la sola differenza significativa è il rispetto intenzionale
della temperatura ambiente configurata.

La matrice regressiva completa è passata: **26/26 comandi**, 13 per profilo
ottimizzato e strict-debug. Include build, parser, stage effettivi, `make check`,
restart binario, restart stocastico, Gate H e source-mass; i test multi-thread
sono stati eseguiti a 1 e 4 thread per entrambi i profili. I risultati
riepilogativi sono in
[results_d02_numerical_fixes_2026_10_21.json](acceptance/results_d02_numerical_fixes_2026_10_21.json);
i log riproducibili sono in `build-local/d02-fixes/validation-logs/`.

La chiusura attesta solo la caratterizzazione e la matrice locale GNU/macOS:
non sono state eseguite nuove validazioni ifx/Linux/Docker né una campagna
scientifica estesa. Il confronto dei diversi provider di permeabilità resta
vincolo D06; ownership e thread-safety dello stato condiviso restano D03;
la rimozione del bypass ODE del modello 8 resta D04; l'accettazione scientifica
cross-model resta D15.

La matrice iniziale e i rilievi storici restano documentati nella sezione
«Findings and D02 worklist» del
[rapporto D01](acceptance/D01_RHEOLOGY_INVENTORY_2026_10_21.md); la loro
risoluzione è registrata nel follow-up D02 collegato sopra.

### D03 — Audit e ownership di parametri, variabili e stato condiviso

**Test proporzionati (nuove modifiche a questo step):** T1 audit dipendenze, ownership e accessi globali; T2 inizializzazione e serializzazione; T3 integrazione se cambia stato persistente. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: IN PAUSA (ON HOLD)

La caratterizzazione e le migrazioni isolate già registrate restano valide; la
revisione generale del manifest e dell'ownership è sospesa per dare priorità a
D04 e D05. Il manifest conserva le voci residue e i limiti di revisione aperti.
Prima di D05 verrà svolto soltanto il sotto-audit mirato dei parametri e dei
contratti reologici interessati dall'estrazione; non è richiesto chiudere
l'inventario condiviso completo.

#### Obiettivo e posizione nella sequenza

Definire un owner univoco per ogni dichiarazione di modulo in constitutive_parameters_2d.f90 e parameters_2d.f90 prima dei trasferimenti. Conservare i risultati della caratterizzazione REAL/COMPLEX D01-D02. Non duplicare il dato nei moduli consumatori: owner significa il modulo responsabile della definizione e del ciclo di vita, non ogni utilizzatore.

#### Avanzamento iniziale

Il workspace `xi_temp`, scritto come variabile di modulo nel kernel Voellmy,
è stato rimosso da `constitutive_parameters_2d`; `eval_implicit_terms` usa ora
lo scalare automatico `xi_effective`. Il regressore
`parallel-voellmy-workspace` esegue 100.000 valutazioni con quattro thread e
quattro valori congelati di `Zij`, confrontandole con i riferimenti seriali.
Strict-debug e optimized passano entrambi 258/258 fixture; l'errore massimo del
regressore parallelo è zero. Questo verifica il percorso interessato, non
certifica la thread-safety complessiva del solver.

È disponibile una versione aggiornata di
[`D03_SHARED_STATE_MANIFEST.json`](acceptance/D03_SHARED_STATE_MANIFEST.json):
227 voci attive per simboli e componenti di `bc`, con dichiarazioni, riferimenti
diretti, import, accessi sintattici, candidati owner e gate. Quattro timer sono
stati trasferiti da `parameters_2d` a `runtime_state_type`: `t_output`,
`t_runout` e `t_probes` sono ora stato di schedulazione runtime; `t_steady`
riusa `runtime%steady_deadline`. L'indice statico non è ancora una revisione
semantica: i riferimenti di riga nei file toccati, i riferimenti indiretti e i
call-context, la classificazione input-derivato-stato-workspace e gli accessi
OpenMP/MPI richiedono verifica. Anche il grafo presente copre solo la relazione
diretta tra i due moduli; l'audit completo e la revisione del manifest restano
aperti. `xi` resta un parametro condiviso di configurazione: il suo valore non
è stato spostato o duplicato. Le rimozioni di stato condiviso sono descritte
qui sopra.

La revisione delle scritture ha inoltre trovato `sp_heat_c` aggiornato sia
durante la lettura input sia a ogni decodifica REAL dello stato. Non risultano
lettori del dato globale nel repository: le procedure di chiusura usano e
restituiscono capacità termiche locali. La variabile, i due aggiornamenti e
l'import del parser sono stati rimossi, eliminando una scrittura condivisa
nelle chiamate `qc_to_qp` parallele senza cambiare le formule di chiusura.
La suite rheology passa 258/258 fixture in strict-debug e optimized; il test
del ciclo di vita del workspace passa. Le build production strict-debug e
optimized superano gli stage quadratic, Coulomb-stop e lava-coupled con uno e
quattro thread; il controllo input passa 28/28 casi. Queste verifiche coprono
la modifica, ma non chiudono l'audit semantico D03.

Il triage aggregato dell'inventario distingue ora 227 voci attive: 225 variabili
di modulo e due componenti di `bc`; 11 sono dichiarate immutabili e 216 restano
classificate come mutabili dall'indice lessicale. Le quattro voci trasferite
sono tracciate separatamente nel manifest. Questa classificazione non implica
che tutte richiedano una modifica. Il controllo per gruppi ha già separato le
costanti/configurazioni di setup, le grandezze geometriche derivate
(`eps_sing`, `eps_sing4`, assegnate in `init_grid`) e i contatori di
dimensionamento inizializzati dall'input (`n_dyn_pres_levels`). I riferimenti a
`n_vars`, `n_eqns` e `n_RK` in istruzioni di allocazione sono falsi positivi
come scritture. La schedulazione (`t_output`, `t_runout`, `t_probes`) è stata
trasferita in `runtime_state_type`, insieme allo stato temporale e viene
aggiornata nei percorsi seriali di output. `t_steady` è stato rimosso come alias
e il ciclo principale usa `runtime%steady_deadline`. La serializzazione mantiene
ordine e numero dei valori del record binario; la prova locale di checkpoint e
continuazione del test mass-exchange passa in strict-debug a uno e quattro
thread, con artefatti confrontati identici; passa anche in optimized a uno e
quattro thread. Il checker Python della suite è stato allineato ai record
diagnostici effettivamente emessi dal fixture, verificando copertura, valori
finiti e coerenza geometrica anziché richiedere righe assenti. L'intera
`TEST_MASS_EXCHANGE/run_test.sh` passa in strict-debug e optimized, inclusi i
controlli projection-stress. Dopo la migrazione, le build production optimized
e strict-debug passano gli stage quadratic, Coulomb-stop e lava-coupled con uno
e quattro thread; `check_inputs.py` passa 28/28 casi. La cache `pascal_coeff`
e il flag sono inizializzati durante `init_param` e poi letti dalle equazioni:
non risultano scritture nel percorso per-cella, quindi sono setup/cache
condivisa, non workspace per thread. La migrazione non modifica il formato
restart.

#### Inventario richiesto e deliverable

Produrre un manifest machine-readable con una riga PER SIMBOLO (non solo per famiglia): dichiarazione e tipo/kind/rank, inizializzazione da namelist o derivazione, tutti gli USE/ONLY e riferimenti indiretti, lettori e scrittori, eventuale mutabilità per time-step/stage, local/shared/allocated lifetime, restart/I-O, OpenMP/MPI, owner attuale, owner proposto e gate richiesto. La matrice seguente e preliminare: i consumatori/scrittori definitivi richiedono ricerca code-wide, non sono dedotti dalla sola dichiarazione.

#### Regole vincolanti

Distinguere input immutabili, grandezze derivate, stato evolutivo, diagnostiche e workspace. Evitare namespace globali generici, duplicazioni e dipendenze circolari; usare PUBLIC selettivi, USE ONLY, configurazioni validate e interfacce esplicite. Non modificare formule, algoritmi, ordine floating-point, serializzazione, semantica COMPLEX o stato persistente solo per spostare simboli. Preservare compatibilità namelist, restart, configurazione N_LAYERS=1 e test/CI. Ogni MODULE, elemento e procedura modificati sono soggetti al contratto Doxygen completo.

#### Criterio di accettazione

Manifest completo e revisionato, dependency graph aciclico, piano di migrazione per dominio fisico, baseline registrate, nessuna proprietà senza owner, nessun workspace globale non motivato. Questi criteri restano il gate per riprendere e chiudere D03, ma non bloccano D04. D05 richiederà invece un sotto-gate circoscritto: ownership, read/write, contratti REAL/COMPLEX e dipendenze dei soli parametri e delle interfacce reologiche coinvolti, senza attendere la revisione generale dell'inventario.

| Origine      | Variabili / famiglie                                                                                                                                      | Owner candidato                                 | Verifiche prima del trasferimento                                         |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------|---------------------------------------------------------------------------|
| constitutive | mu, xi, mu_0, mu_inf, Fr_0, U_w                                                                                                                  | Rheology / each law                             | xi_temp is local now; parallel Voellmy regression passes; audit xi ownership |
| constitutive | I_0, mu_2, mu_s, muI_inf, friction_factor, tau                                                                                                            | Rheology / each law                             | Cross-model parameters and read/write/validation paths                    |
| constitutive | T_ref, nu_ref, visc_par, tau0, emme, enne                                                                                                                 | Lava rheology / thermophysical definition       | Avoid duplicating shared thermal parameters                               |
| constitutive | alpha1_coeff, beta1, alpha2, beta2, Kappa, n_td, n_td2                                                                                                    | Lahar rheology                                  | Initialization/order of derived n_td2                                     |
| constitutive | SBconst, emissivity, convective_term_coeff, radiative_term_coeff, atm_heat_transf_coeff, exp_area_fract                                                   | Thermal exchange                                | Separate environmental configuration from derived rates                   |
| constitutive | T_env, T_ground, thermal_conductivity, c_p                                                                                                                | Thermal / ambient or phase-property owner       | Resolve whether lava-only or general temperature                          |
| constitutive | grav, inv_grav, pres, inv_pres, T_ambient, rho_a_amb                                                                                                      | Physical constants / ambient conditions         | Immutable runtime configuration vs derived reciprocals                    |
| constitutive | rho_l, inv_rho_l, kin_visc_l, sp_heat_l, gamma_steam                                                                                                      | Liquid / vaporization                           | Check vaporization and composition clients                                |
| constitutive | sp_heat_a, sp_gas_const_a, kin_visc_a, sp_heat_g, sp_gas_const_g                                                                                          | Gas and phase properties                        | Allocated gas vectors, n_add_gas and restart                              |
| constitutive | rho_c_sub, kin_visc_c, Tref_Suth, muRef_Suth, S_mu                                                                                                               | Carrier gas / liquid / Sutherland               | sp_heat_c removed as unread shared state; retain local mixture outputs    |
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
| parameters   | rheology_flag, rheology_model, liquid_flag, gas_flag, pore_pressure_flag, sutherland_flag                                                                 | Model physics configuration / rheology selector | Rheology selection migration is included in D05                           |
| parameters   | slope_correction_flag, curvature_term_flag, topo_change_flag, water_level                                                                                 | Geometry / spatial physics configuration        | Time-dependent geometry remains updateable                                |
| parameters   | stochastic_flag, stoch_transport_flag, length_spatial_corr                                                                                                | Stochastic model configuration                  | RNG/OU state and restart distinct from input settings                     |
| parameters   | x_source, y_source, r_source, r2_source, azimuth_source, arc_width_source, mfr_source, vel_source, T_source, h_source, alphas_source(:), alphag_source(:) | Source configuration (radial/initial)           | Inspect activation flags, units and source initialization                 |
| parameters   | lateral_source_flag, source_side, x1_source, x2_source, y1_source, y2_source, bottom\_\*\_source_flag, n_fissures, \*\_fissures                           | Source configuration (lateral/fissural)         | Time histories, fixed arrays and source validation                        |
| parameters   | collapsing_volume_flag, x_collapse, y_collapse, r_collapse, T_collapse, h_collapse, alphas_collapse(:)                                                    | Initial/collapse source configuration           | Shared initialization and startup state                                   |
| parameters   | particle_distribution_flag, fractal_dim, diam_min, diam_max, rho_s_min, rho_s_max, sigmoid_k, sigmoid_d0                                                  | Particle-size distribution input                | Separate PSD parameters from derived per-class arrays                     |
| parameters   | bc type, bcW(:), bcE(:), bcS(:), bcN(:)                                                                                                                   | Boundary-condition descriptors                  | Ownership and allocation/restart semantics                                |
| parameters   | dt_output, verbose_level, n\_\*\_levels, thickness_levels(:), dyn_pres_levels(:)                                                                              | Input/configuration and model layout             | Output timers moved to `runtime_state_type`; retain setup/diagnostic ownership review |

### D04 — Eliminazione reologia 8

**Test proporzionati (nuove modifiche a questo step):** T1 rifiuto reologia 8, riferimenti e parsing; T2 un esempio con reologia supportata; T3 regressione altre reologie solo se toccate. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO (gate locale GNU/macOS)

Questa decisione è definitiva.

Il modello 8 era un percorso sperimentale che integrava analiticamente l’ODE
della frizione quadratica, anziché risolvere l’equazione algebrica dello stage
DIRK. La decisione di rimuoverlo è stata implementata.

La verifica del repository non ha trovato input di esempio o generatori che
selezionassero `RHEOLOGY_MODEL = 8`; il parser ora include un controllo
esplicito del suo rifiuto.

Sono stati rimossi:

Il bypass ODE in `nonlinear_solver_2d.f90`.

La procedura `integrate_friction_term` in `equation_terms_2d.f90`.

Il trattamento speciale in `time_integration_2d.f90`.

La selezione e validazione del modello 8 nell'input; il parser lo respinge come
modello non supportato.

Il ramo per `shearVel`/Rouse specifico del modello 8 in `inpout_2d.f90`.
Le variabili NetCDF restano nello schema e conservano il comportamento nullo
degli altri modelli.

I tre fixture ODE legacy dalla caratterizzazione corrente; i risultati storici
D01/D02 restano archiviati come evidenza della baseline precedente.

La futura soluzione analitica della frizione quadratica apparterrà a
`rheology_quadratic_2d.f90` e risolverà correttamente lo stage DIRK.

Non verrà mantenuto un metodo speciale di integrazione dell’ODE esterno alla nuova architettura.

#### Caso \#8 ed identificatori testuali

Completato: rimossi i rami `rheology_model==8` da solver, integratore e output,
oltre alla routine `integrate_friction_term`; gli input con modello 8 sono
respinti. Il modello 6 conserva la propria fisica Bursik--Woods. Il futuro
dispatcher D05 risolverà una volta gli identificatori esterni basati sui
suffissi `rheology_<nome>.f90`; il modello 7 resta da auditare prima di
assegnargli un nome/file. Nessun generatore esistente selezionava il modello 8.

Decisione approvata e implementata: eliminare completamente il modello
reologico 8, mantenere il modello 6 per la legge Bursik--Woods; la soluzione
analitica dello stage è una delle due modalità di gruppo, insieme al solve
Newton. `IM_GROUP_METHOD_ANALYTIC` indica la soluzione chiusa dell’equazione
implicita DIRK del gruppo;
`IM_GROUP_METHOD_NEWTON` indica il solve implicito numerico, con SI dove
previsto.

Validazione D04:
- build production optimized e strict-debug riuscite;
- 29/29 casi parser passati in entrambi i profili, inclusi il rifiuto del
  modello 8 e il controllo positivo per il modello 6;
- 255/255 fixture di caratterizzazione passati in strict-debug e optimized;
- stage quadratic, Coulomb-stop e lava-coupled passati con uno e quattro thread
  in entrambi i profili.

### D05 — Estrazione moduli reologici

**Test proporzionati (nuove modifiche a questo step):** T1 parsing dei nomi, contratti dei moduli reologici e REAL/COMPLEX; T2 un caso per modello attivo e confronto pre/post; T3 nessun cambiamento di traiettoria per il solo spostamento del codice. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Prerequisiti: D04 completato e sotto-audit mirato dei parametri e delle interfacce reologiche da estrarre. La chiusura dell'audit generale D03 non è richiesta. Durante l’estrazione dei moduli reologici, trasferire agli owner di legge mu/xi, parametri lava, lahar, plastic, quadratic, Froude, Lucas e mu(I); non trascinare proprietà fisiche non reologiche in un generico modulo reologico. Registrare read/write e mantenere la numerica REAL/COMPLEX. La selezione testuale descritta in D08 è inclusa in D05: l’input usa i nomi canonici e il dispatcher risolve il nome una sola volta all’inizializzazione. I vecchi ID numerici sono respinti esplicitamente.

Stato ufficiale: COMPLETATO

Regola di struttura: un modello reologico per file e per modulo. Tutti i
contributi di una legge devono risiedere nel suo modulo, anche quando oggi
sono implementati in routine solver diverse: termini semi-impliciti,
impliciti e, se prevista, la soluzione analitica dello stage. I moduli condivisi
sono ammessi per parametri o utility comuni, ma non per accorpare kernel di
leggi distinte. In particolare `mui` e `mui_voellmy`, pur avendo attualmente
formule uguali, mantengono moduli separati. `rheology_parameters_2d` è soltanto
l’owner condiviso dei coefficienti ancora usati da namelist e consumer
esistenti; non contiene leggi reologiche.

| ID interno | Nome input | Modulo previsto |
|------------|------------|-----------------|
| 0          | `none`     | nessuna reologia |
| 1          | `voellmy`  | `rheology_voellmy_2d.f90` |
| 2          | `plastic`  | `rheology_plastic_2d.f90` |
| 3          | `lava`     | `rheology_lava_2d.f90` |
| 4          | `lahar`    | `rheology_lahar_2d.f90` |
| 5          | `empirical`| `rheology_empirical_2d.f90` |
| 6          | `quadratic`| `rheology_quadratic_2d.f90` |
| 9          | `froude`   | `rheology_froude_2d.f90` |
| 10         | `lucas`    | `rheology_lucas_2d.f90` |
| 11         | `mui`      | `rheology_mui_2d.f90` |
| 12         | `mui_voellmy` | `rheology_mui_voellmy_2d.f90` |

Gli ID 7 (non supportato) e 8 (rimosso in D04) non hanno nomi input.

Questa è la parte operativa: cosa deve cambiare rispetto all’implementazione attuale.

| File                           | Modifiche richieste                                                                                                                                     |
|--------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------|
| parameters_2d.f90              | Memorizzare il nome canonico e l’ID interno della reologia risolto una sola volta.                                                                       |
| inpout_2d.f90                  | Leggere i nomi, respingere input numerici/sconosciuti, validare i parametri e stampare il nome selezionato.                                               |
| constitutive_parameters_2d.f90 | Mantenere i parametri globali attuali durante l’estrazione; trasferire solo i coefficienti con ownership reologica non condivisa.                         |
| equation_terms_2d.f90          | Conservare l’assemblaggio generale e delegare ai moduli reologici i contributi specifici delle leggi.                                                    |
| rheology_2d.f90                | Risolvere il nome in un ID al setup e dispatchare ai kernel della legge, senza confronti tra stringhe nei cicli numerici.                                 |
| rheology\_\*.f90               | Un file/modulo per legge, con tutti i contributi SI/IM e l’eventuale soluzione analitica della legge; mantenere invariati i contratti REAL/COMPLEX. |
| stochastic_module.f90          | Usare l’ID interno per le scelte stocastiche dipendenti dalla reologia.                                                                                   |
| EXAMPLES/, TESTS/, generatori  | Aggiornare gli input attivi, gli esempi e gli script; testare nomi, maiuscole, nomi sconosciuti e rifiuto degli ID numerici.                              |
| Makefile.am                    | Registrare i nuovi moduli e le dipendenze Fortran; rigenerare i file Automake derivati.                                                                  |

D05 resta circoscritto all’estrazione delle reologie e alla selezione testuale. Le sorgenti termiche e di pore pressure restano in D06; EX_terms/IM_terms restano in D10; i gruppi e il Jacobiano ridotto restano in D11–D13.

Esito del sotto-audit: i coefficienti reologici sono posseduti dai rispettivi
moduli e letti/validati da `inpout_2d`; `alpha1_coeff` per lahar viene derivato
dall’input durante il setup. L’ID del modello ha inoltre consumer fuori dai
kernel: partizione della temperatura in `init_problem_param`, selezione della
viscosità e del settling lahar, rumore stocastico e diagnostiche di attrito.
Per questo l’input testuale viene risolto una volta in un ID interno. `xi`,
`T_ref` e gli altri parametri condivisi restano in `rheology_parameters_2d`;
le sorgenti termiche e di pore pressure generali restano in D06.

Sono ora presenti moduli separati per tutti i modelli attivi: `voellmy`,
`plastic`, `lava`, `lahar`, `empirical`, `quadratic`, `froude`, `lucas`, `mui`
e `mui_voellmy`. Ogni modulo contiene i contributi del modello che prima erano
dispersi tra le routine implicite e semi-implicite; il modulo lahar contiene
anche la relazione di viscosità usata dal settling e la propria regola di
deposizione. I coefficienti specifici sono stati trasferiti ai rispettivi
moduli; `xi`, `T_ref` e i coefficienti già condivisi restano nel modulo
parametri, che non contiene kernel. I consumer di setup, diagnostica e rumore
stocastico restano nei rispettivi percorsi generali, ma leggono gli owner
espliciti. La caratterizzazione passa 255/255 fixture in strict e optimized,
il parser passa 32/32 casi e la build completa passa in strict-debug e
optimized.

Esito T2/T3 — confronto pre/post tra baseline `379301f1` e candidato
`664c82bc`, con GNU Fortran 15.2 e la fixture dinamica
`TESTS/TEST_DAMBREAK/TEST_STEP_FINH`. È stato eseguito un caso per ciascuno dei
dieci modelli attivi, con output a 0, 0.005 e 0.01 s; per Lahar sono state
impostate la fase liquida e `SED_VOL_PERC=50`. Tutti i run hanno raggiunto la
fine e prodotto stati finiti e dinamici. In strict (`-O0`, controlli runtime,
un thread), i 30 snapshot confrontati sono byte-identici. In optimized
(`-Ofast`), sia a uno sia a quattro thread, la massima norma `Linf` per campo,
normalizzata con `max(1, Linf(baseline))` e presa su snapshot e variabili di
stato, resta sotto la soglia `1e-10` per ogni modello:

| Modello | Strict, 1 thread | Massima norma normalizzata optimized, 1/4 thread |
|---------|------------------|---------------------------------------------------|
| `voellmy` | identico | 5.751e-13 |
| `plastic` | identico | 0 |
| `lava` | identico | 1.724e-20 |
| `lahar` | identico | 0 |
| `empirical` | identico | 0 |
| `quadratic` | identico | 0 |
| `froude` | identico | 0 |
| `lucas` | identico | 0 |
| `mui` | identico | 0 |
| `mui_voellmy` | identico | 3.553e-17 |

Il gate T3 non rileva cambiamenti di traiettoria oltre la tolleranza stabilita;
D05 è completato. D06 resta il prossimo step.

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
├── physics/
│   ├── equation_terms_2d, mass_exchange_laws_2d,
│   │   pore_pressure_laws_2d, thermal_laws_2d
│   └── rheology_laws/ rheology_2d, rheology_voellmy,
│       rheology_plastic, rheology_lava,
│       rheology_lahar, rheology_empirical,
│       rheology_quadratic, rheology_froude,
│       rheology_lucas, rheology_mui,
│       rheology_mui_voellmy \[\*.f90\]
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

thermal_laws_2d

Pore pressure

pore_pressure_laws_2d

Schema delle responsabilità e delle interfacce concettuali, non una sequenza obbligatoria di chiamate per la soluzione analitica.

La caratteristica fondamentale della nuova struttura sarà la separazione fra il modello fisico, il trattamento numerico delle sorgenti e l’integratore IMEX.

In particolare, aggiungere una nuova reologia dovrà richiedere principalmente un nuovo file rheology\_\<nome\>.f90 e la sua registrazione nel dispatcher, senza modificare il ciclo RK o l’algoritmo Newton.

Questo è il piano consolidato che utilizzerei come riferimento per l’implementazione e per il successivo aggiornamento della Development Roadmap.

### D06 — Sorgenti termiche e pore pressure

**Test proporzionati (nuove modifiche a questo step):** T1 sorgenti termiche/pore pressure e limiti sottili; T2 casi lava/raffreddamento; T3 energia e stiff RK se cambiano i termini. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Perimetro dell’estrazione:

- `thermal_laws_2d` possiede le leggi per l’energia delle iniezioni generali, la vaporizzazione e il raffreddamento radiativo/convettivo della lava. Il termine energetico del mass exchange viene assemblato direttamente in `equation_terms_2d::eval_mass_exchange_terms`, insieme agli altri termini dello scambio. `rheology_lava_2d` resta proprietario dei soli termini reologici; parametri termici dichiarati ma non usati non costituiscono nuove leggi di conduzione o riscaldamento.
- `pore_pressure_laws_2d` possiede diffusione implicita (inclusa la sensibilità COMPLEX), perdita Darcy di gas e accoppiamenti conservativi con iniezioni e mass exchange.
- Le due formulazioni della permeabilità dinamica restano provider distinti: diffusione con pesi di sfericità `alpha*d^2`, gas loss con `alpha/d`. L’estrazione conserva e verifica entrambe; non ne presume l’equivalenza né ne uniforma la fisica.

Stato ufficiale: COMPLETATO (estrazione e gate locale GNU/macOS). T1 strict/optimized con gli oracoli indipendenti e i test termici, e T2 lava/raffreddamento, sono superati localmente. Decisione: mantenere distinti i provider di permeabilità (`alpha*d^2` per la diffusione, `alpha/d` per la Darcy gas-loss); un’eventuale armonizzazione resta una modifica scientifica futura, non parte di D06. T3 non era necessario perché i termini osservabili non sono cambiati.

### D07 — Migrazione progressiva dei parametri non reologici

**Test proporzionati (nuove modifiche a questo step):** T1 ownership parametri, inizializzazione e default; T2 input/restart rappresentativi; T3 solo se mutano effetti numerici. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: IN PAUSA (in attesa della ripresa e del gate D03).

#### Ambito e dipendenze

Eseguire per famiglie, DOPO il gate D03 e coordinandosi con D05/D06 (reologie e sorgenti). D07 resta in pausa finché D03 non viene ripreso e il suo gate non è chiuso. Non spostare ripetutamente le variabili da un contenitore globale a un altro. La riallocazione deve portare ownership semantica, API esplicite e minor accoppiamento verificato; E01 restera un audit architetturale indipendente post-refactoring.

#### Sequenza operativa

\(1\) Parametri ambientali, proprieta di gas/liquidi/solidi e loro reciprocali; (2) trasporto/deposizione e mass exchange; (3) pressione interstiziale e coefficienti di inibizione; (4) configurazione sorgenti, condizioni al contorno, output e runtime da parameters_2d; (5) metadati layout, solver e geometria. La migrazione termica va sincronizzata con D06 e non duplicata. Ogni tranche deve identificare gli USE, le inizializzazioni, le allocazioni e la proprieta dello stato.

#### Policy dello stato

Parametri letti da namelist convalidati una volta; grandezze derivate calcolate al momento corretto e aggiornabili se la fisica richiede variabilita; tempo corrente, diagnostica e workspace di stage non devono diventare configurazioni globali. Non assumere geometria/topografia statica: il DEM puo evolvere nel tempo. Gli array allocabili devono avere ownership esplicita, senza copie divergenti.

#### Acceptance per tranche

Compilazione GNU e ifx, OpenMP 1/4 thread, test di input/restart, confronto numerico strict/optimized e casi fisici mirati. Registrare symbol manifest before/after, dependency graph, confronto di default e prestazioni; ogni modifica alla semantica delle fonti/risultati e un task scientifico separato. Aggiornare Doxygen per ogni simbolo e procedura alterati.

### D08 — RHEOLOGY_MODEL da numeri a stringhe (assorbito in D05)

La migrazione dell’input a nomi canonici, il rifiuto degli ID numerici e i
relativi test sono stati trasferiti in D05 per evitare due step duplicati. D08
non è più un’attività indipendente; i test di gruppi ANALYTIC/NEWTON e
OpenMP restano assegnati a D11–D14.

### D09 — Correzione del solo attrito

**Test proporzionati (nuove modifiche a questo step):** T1 isolamento frizione e segno dei termini; T2 casi dissipativi/non dissipativi; T3 conservazione e regressione reologica se cambia forza effettiva. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO (gate locale GNU/macOS).

La procedura:

correct_friction_momentum(...)

è contenuta in:

src/physics/rheology_laws/rheology_2d.f90

Sostituisce i due controlli ad hoc dopo SI e dopo Newton. Riceve esplicitamente
il momento prima e dopo il sottopasso di frizione e proietta il candidato sulla
direzione iniziale con un fattore compreso fra zero e uno. Se il momento iniziale
è nullo, il risultato resta nullo.

La correzione deve preservare la direzione e non aumentare il modulo della quantita di moto SOLO per il sottopasso in cui agisce esclusivamente la frizione. La stessa proprieta non si applica allo stato complessivo di uno stage DIRK o di un gruppo in cui siano presenti gravita, termini termici, pore pressure o altre sorgenti che possono ruotare o accelerare il momento.

Pertanto, per un aggiornamento esclusivamente reologico:

La correzione non deve essere applicata a contributi spaziali, gravitativi o ad altre forze che possono modificare la direzione del moto.

La sostituzione elimina l’incoerenza precedente fra OR e AND nei controlli
componente per componente.

Gate locale superato: i test unitari verificano allineamento, proiezione di un
candidato ruotato, arresto per inversione, momento iniziale nullo e limite
superiore al modulo. Le regressioni di integrazione strict/optimized verificano
quadratica, arresto Coulomb e lava accoppiata a 1/4 thread. Le formule delle
reologie non sono state modificate; la correzione si applica solo alle coppie
di stati che delimitano gli aggiornamenti semi-implicito e implicito di frizione,
non allo stato globale dello stage.

#### Correzione esclusivamente reologica

correct_friction_momentum vive in rheology_2d.f90 ed opera solo sul contributo della frizione isolato: m_after = alpha\*m_before, 0 \<= alpha \<= 1, con arresto ammissibile e nessuna rotazione/inversione. L'interfaccia rende espliciti momentum_before_friction e momentum_after_friction e non proietta lo stato globale dello stage. L'integrazione attuale usa il correttore sugli aggiornamenti SI e NEWTON; ANALYTIC e i futuri gruppi misti richiedono una prova separata dell'isolamento del contributo prima di applicarlo.

### D10 — EX_terms / IM_terms

**Test proporzionati (nuove modifiche a questo step):** T1 contratto EX_terms/IM_terms e consistenza delle somme; T2 uno step IMEX con tutte le componenti; T3 casi stiff se cambia partizione. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO (gate locale GNU/macOS)

Il ciclo RK conosce e combina soltanto:

spatial_terms  
EX_terms  
IM_terms

La separazione NH e SI_NH non è esposta all’integratore, poiché i due contributi utilizzano lo stesso tableau implicito.

Il ciclo costruisce la parte nota dello stage con spatial_terms ed EX_terms tramite il tableau esplicito, e con IM_terms degli stage precedenti tramite il tableau implicito. La diagonale implicita corrente viene risolta nel solver locale.

Per uno stage con $a_{ii} \neq 0$, il solver locale restituisce lo stato aggiornato e il contributo IM effettivo, ricostruito dopo le correzioni applicate.

Il trattamento SI, Newton e le correzioni di frizione restano interni al solver locale. Il caso predefinito a singolo gruppo mantiene la sequenza SI+Newton; i tableau IMEX non sono stati modificati.

Gate locale superato: `make check` strict e optimized; i test degli stage di produzione strict e optimized per quadratica, arresto Coulomb e lava accoppiata, con 1 e 4 thread; suite cell-local con 255 fixture strict e optimized. Nei test stage `known_stage2`, `effective_increment` e `final_projection` hanno errore zero.

#### SI legacy, correzioni e incremento implicito effettivo

Nel default un unico gruppo NEWTON deve riprodurre in sequenza la parte SI, lo stato intermedio q_si, la soluzione Newton e le correzioni attuali, salvo la successiva correzione scientifica esplicitamente testata del criterio di arresto. Non confondere una sorgente fisica valutata dopo clipping con il contributo effettivo RK. Dopo le correzioni ricostruire:

Con tau=0, lo stage non effettua un solve implicito: evitare divisioni per zero e definire il contributo/stage coerentemente con il tableau e la sua necessità nel passo. I tableau IMEX e la distinzione sorgenti esplicite/operatore spaziale non vengono modificati.

### D11 — IM_group e metadati

**Test proporzionati (nuove modifiche a questo step):** T1 metadati e copertura IM_group; T2 selezioni e subset; T3 regressioni RK quando cambia il raggruppamento. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO (gate locale GNU/macOS)

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

#### Implementazione ed evidenza D11

`equation_partition_type%initialize` accetta opzionalmente la mappa dei gruppi
e i relativi metodi. Le chiamate legacy che forniscono solo `implicit_mask`
mantengono un unico gruppo NEWTON; una partizione tutta esplicita non crea
gruppi. L'inizializzazione valida dimensioni, corrispondenza tra maschera e
gruppi, numerazione densa, metodi supportati e mappe compatte.

`init_problem_param` assegna gruppi ANALYTIC separati a energia lava e
sovrapressione quando attive, seguiti dal gruppo NEWTON dei momenti. I gruppi
disattivati vengono omessi e rinumerati senza buchi. Questi campi restano
metadati: il solver, il tableau e il dispatch non li consumano in D11.

Verifica locale: `TESTS/TEST_EQUATION_METADATA/run_test.sh` supera i casi
strict per partizioni legacy, gruppi multipli, subset espliciti/impliciti,
metodi e invarianti. `make check` supera il build ottimizzato e le regressioni
di integrazione locali; il runner dei metadati è incluso in `check-local`.

### D12 — Newton per gruppo e Jacobiano ridotto

**Test proporzionati (nuove modifiche a questo step):** T1 Jacobiani ridotti, complex step e convergenza Newton; T2 problemi stiff; T3 full/split e profili se cambiano traiettorie. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: IN CORSO — implementata la sweep ordinata dei gruppi e il
Newton ridotto alle sole equazioni del gruppo; gate d'integrazione completo
ancora aperto.

#### Jacobiano Newton

La dimensione del Jacobiano Newton dipende dal gruppo attivo, non dal totale
`n_implicit`. Il residuo, le righe/colonne del Jacobiano complex-step e la
line-search includono soltanto le equazioni del gruppo; le altre componenti
restano fisse durante la sua iterazione. I sistemi 2×2/3×3 riutilizzano i
risolutori pivotati esistenti. Le mappe compatte sono precalcolate nei
metadati; il solver usa workspace locali automatici, senza allocazioni
esplicite per cella.

#### Inizializzazione e validazione dei gruppi

equation_partition_type conserva n_equations, n_explicit, n_implicit, implicit(:), implicit_map(:) e explicit_map(:). Aggiunge n_IM_groups, IM_group(:) per equazione, IM_group_method(:) per gruppo e le mappe compatte precomputate. La relazione IM_group(i)=0 se e solo se .NOT. implicit(i) è un’invariante. I gruppi devono essere numerati da 1 a n_IM_groups senza vuoti e ogni metodo deve essere valido. Ogni sorgente IM del gruppo deve essere supportata dal suo trattamento: Newton per il residuo implicito oppure integrazione analitica a coefficienti congelati, con formula chiusa ove disponibile e inversione numerica controllata dell’integrale temporale negli altri casi. Il mapping globale non può essere sovrascritto con indici locali temporanei. implicit(:) deve avere una sola fonte di verità, evitando configurazioni indipendenti incoerenti.

#### Contratto algebrico dei gruppi

Sia tau = dt\*a_diag e q_known lo stato noto dopo l’assemblaggio dei termini RK dagli stage precedenti. Un gruppo g risolve soltanto le proprie componenti, sempre con il termine noto dello stage q_known(i); le sorgenti possono leggere lo stato aggiornato dai gruppi precedenti e il valore ancora noto dai gruppi successivi. Questa strategia è una singola sweep block Gauss–Seidel, non identica alla soluzione simultanea se esiste coupling intergruppo. Vanno testati ordine di convergenza temporale, stabilità e sensibilità all’ordine dei gruppi. Nessuna maschera separata di lagging. Il gruppo unico NEWTON costituisce un gate di equivalenza con il solver monolitico preesistente, con confronto di residuali, stati e convergenza agli stage secondo le tolleranze scientifiche approvate. I gruppi multipli ordinati e la singola sweep Gauss-Seidel costituiscono invece una MODIFICA NUMERICA, non un refactoring equivalente: devono superare un distinto gate di accuratezza, convergenza rispetto al coupling e dipendenza dall'ordine dei gruppi, anche quando l'esito non coincide bitwise col solve simultaneo.

Per NEWTON, costruire il residuo locale R_g(x_g) e il Jacobiano J_g = dR_g/dx_g con dimensione n_g × n_g, n_g = COUNT(IM_group==g). Complex-step perturba solo colonne del gruppo e mantiene le altre variabili fisse. Riutilizzare gli ottimizzati solver lineari 2×2/3×3 quando appropriato, calcolare gradienti e line search nello spazio locale e garantire nessuna modifica a componenti esterne. Preallocare workspace per thread, dimensionati per il massimo gruppo Newton, **non** allocare per cella.

#### Implementazione ed evidenza D12

`solve_implicit_groups` attraversa i gruppi nell’ordine dichiarato;
`solve_rk_step` resta il wrapper compatibile. Le sole modalità di gruppo sono
`IM_GROUP_METHOD_ANALYTIC` e `IM_GROUP_METHOD_NEWTON`: la prima usa
il flusso nel tempo dei termini IM con coefficienti non cinematici congelati
(formule chiuse quando disponibili, integrazione numerica controllata
dell’integrale temporale altrimenti); la seconda risolve il residuo con
Newton e include il contributo semi-implicito per la quantità di moto. L’input assegna ciascuna
famiglia (velocità, energia, pressione dei pori) a un ID di gruppo da 1 a 3:
ID uguali combinano le famiglie nello stesso gruppo e gli ID crescenti
stabiliscono l’ordine dello sweep. `IMEX_GROUP_METHOD(1:3)` sceglie `ANALYTIC`
o `NEWTON` per l’ID corrispondente. Per esempio, assegnare tutte le famiglie
al gruppo 1 e impostare `IMEX_GROUP_METHOD(1)='NEWTON'` risolve tutte le
equazioni implicite in un unico gruppo Newton. Per un gruppo analitico misto,
i termini sono integrati nell’ordine delle equazioni conservative,
ricostruendo lo stato prima del termine successivo.

Le variabili appartengono al namelist `NUMERIC_PARAMETERS`; i default
mantengono l’ordine energia, pressione dei pori, velocità:

```fortran
imex_energy_group = 1
imex_pore_pressure_group = 2
imex_velocity_group = 3
imex_group_method = 'ANALYTIC', 'ANALYTIC', 'NEWTON'
```

Per ogni gruppo Newton, il solver riduce residuo, Jacobiano e sistema lineare
al relativo mapping e mantiene invariate le componenti esterne. Il test
`TEST_RHEOLOGY_CELL_LOCAL/test_grouped_implicit_stage.f90` verifica l’ordine e
la composizione configurati, un gruppo Newton congiunto, i solve DIRK chiusi e
il rate effettivo netto in strict/optimized. La line-search conserva
esattamente le componenti esterne al gruppo attivo.

Il confronto spaziale e gli stage IMEX usano ora l’oracolo NumPy indipendente
KT corrente; i core HP-PCCU storici non sono più usati per residui, velocità
caratteristiche o CFL. `TEST_SPATIAL_OPERATOR` passa in strict e optimized:
otto fixture, quattro configurazioni geometriche, esecuzione a 1/4 thread e
sette controlli negativi. `TEST_IMEX_STAGES` passa i confronti degli stage e
le due policy CFL in entrambi i profili; i controlli analitici degli stage
N_RK=2/3/4 passano anch’essi in strict e optimized. Il riferimento KT applica
ora la soglia dry anche alla profondità usata per le velocità d’onda; il suo
checksum aggiornato è registrato nel contratto Gate-H v4 e i sei test
autonomi del riferimento passano.

Il runner non collega più i moduli spaziali storici al driver attuale per il
test dei confini: i relativi workspace hanno layout incompatibili, quindi la
doppia allocazione osservata nel binario baseline era un problema ABI del test,
non un difetto riprodotto nel candidato. Il driver verifica ora direttamente
le tracce analitiche del candidato. La sezione `stages` e i controlli dei
confini passano in strict e optimized, inclusi assi collassati e thread 1/4.

La verifica precedente si arrestava in `order_suite` su
`manufactured mass remains constant` (errore scalato `0.31518397499667228`,
limite `7.2759576141834259e-12`): la fixture usava flusso uniforme ma
condizioni di default che bloccano l’ingresso a ovest/sud e consentono
l’uscita a est/nord.

Aggiornamento 23 ottobre 2026: il controllo è stato riallineato ai bordi
outflow-only. Il fixture d’ordine usa una griglia 16x16; gli oracoli analitici
confrontano la regione interna, separata dai bordi da sei celle, mentre il
bilancio globale verifica che la massa non aumenti e non impone più una massa
costante. Il controllo del drag verifica l’equazione implicita locale rispetto
allo stato noto di ciascuno stage, senza assumere nullo il contributo esplicito
di bordo. La sezione `order` passa in strict e optimized per thermal, drag e
mixed (`N_RK=2`). D12 resta IN CORSO fino alla riesecuzione del runner
integrato.

### D13 — Integrazione analitica esatta dei termini IM

**Test proporzionati (nuove modifiche a questo step):** T1 confronto quantitativo con le soluzioni esatte disponibili; T2 regimi limite, soglie e monotonia; T3 regressioni del tableau per N_RK 2/3/4. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO (validazione locale GNU/macOS). I gruppi
ANALYTIC integrano nel tempo i termini IM con i coefficienti non integrati
congelati all’ingresso del gruppo. Le formule di flusso chiuse sono usate
quando disponibili; gli altri termini sono trattati tramite inversione
controllata dell’integrale scalare del tempo. Le suite strict e optimized
verificano soluzioni esatte, regimi limite, tableau `N_RK=2/3/4` e traiettorie
complete con gruppi multipli (D14), con esecuzioni a uno e quattro thread.

#### Flusso congelato per la velocità

Le routine di modulo integrano nel tempo il termine sorgente combinato per le
due componenti orizzontali della velocità. Profondità, temperatura, densità,
frazioni solide, gravità, stress normali, geometria/curvatura, pressioni e
rapporti tra velocità totale e orizzontale sono ingressi espliciti congelati
durante l’integrazione. La direzione orizzontale viene preservata; si evolve
il modulo scalare della velocità e lo si arresta a zero quando il flusso si
ferma.

Empirica, plastic, lava e quadratic drag usano le rispettive formule esplicite
di flusso. Per gli altri modelli si inverte l’integrale separato
\f$\int dr/f(r)\f$ con quadratura adattiva e bisezione, senza sotto-passi
temporali. I coefficienti non cinematici restano congelati; l’integrale
numerico ha controlli espliciti di convergenza.

#### Energia e pressione interstiziale

Ogni quantità segue lo stesso contratto di metodo: `ANALYTIC` integra il
termine IM nel tempo con i coefficienti non cinematici congelati; `NEWTON`
risolve il residuo implicito dello stage e include SI per la quantità di moto
dove previsto. Per la velocità, il solve analitico comprende già i termini
semi-implicito e implicito, quindi non si applica SI separatamente.

Per il raffreddamento lava, `solve_lava_analytic_temperature` è l’unico solve
selezionato da `ANALYTIC` e integra l’ODE di raffreddamento a coefficienti
congelati. Non esiste un solve algebrico DIRK alternativo.

Per la diffusione della pressione interstiziale, `solve_analytic_pore_pressure`
è l’unico solve selezionato da `ANALYTIC` e usa la soluzione esponenziale
esatta con coefficienti congelati. Se energia e pressione condividono un
gruppo, la pressione usa la temperatura aggiornata nella sequenza del gruppo.
`advance_lava_cooling_temperature` e `advance_pore_pressure_diffusion`
mantengono i rispettivi nomi di compatibilità.

I test diretti verificano le traiettorie con build strict e optimized. Il test
`test_analytic_velocity` controlla la proprietà di semigruppo per tutti i
modelli scalari, inclusi attraversamenti delle soglie Froude/Lucas, e verifica
gli endpoint contro una soluzione di riferimento RK4 indipendente e gli
estremi di arresto e intervallo nullo. Controlla inoltre i limiti di tempo
nullo, velocità nulla e coefficiente di drag nullo per tutte le leggi scalari,
i rami degeneri della legge polinomiale e i casi di arresto finito. I test
energia e pressione verificano anche equilibrio alla temperatura ambiente,
tempo nullo, limite a lungo tempo, pressione nulla, stato senza volume poroso
e permeabilità nulla, oltre alla composizione temporale e al limite a breve
intervallo.
La suite verifica gruppi ANALYTIC e NEWTON, famiglie condivise, ordine
configurato, aggiornamenti visibili ai gruppi successivi e coerenza del rate
effettivo con l’incremento conservativo netto. La regressione
`TEST_IMEX_STAGES/compare_stages.py --section analytic` verifica inoltre,
per `N_RK=2/3/4`, gli stati noti e risolti degli stage e l’assemblaggio finale
contro flussi quadratici indipendenti; il test controlla anche la soluzione
esatta dell’endpoint nel caso a due stage e la convergenza al raffinamento
temporale per i tableau a tre e quattro stage, con esecuzioni a uno e quattro
thread. La traiettoria completa dell’intero tableau con gruppi multipli è
verificata dalla suite D14.

#### Contratto d’integrazione

Il metodo `ANALYTIC` integra **l’ODE del termine IM su `dt*a_diag` con i
coefficienti non cinematici congelati**, non risolve l’equazione algebrica
DIRK. Quando esiste una formula chiusa dello stato finale, si usa quella;
negli altri casi si inverte l’integrale del tempo con un metodo numerico
controllato. Non utilizzare Newton implicitamente dentro una routine
ANALYTIC. Il test
`../TESTS/TEST_RHEOLOGY_CELL_LOCAL/test_analytic_velocity.f90` verifica i
limiti fisici e le formule esatte disponibili.

### D14 — Gruppi misti e validazione

**Test proporzionati (nuove modifiche a questo step):** T1 gruppi misti e fallback; T2 combinazioni rappresentative; T3 confronto delle traiettorie se cambia la strategia implicita. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: COMPLETATO (validazione locale GNU/macOS). `TEST_RHEOLOGY_CELL_LOCAL/run_test.sh`
passa 255/255 fixture in strict e optimized, inclusi i contratti dei gruppi,
le combinazioni rappresentative e il solve Newton monolitico. Inoltre,
`TEST_IMEX_STAGES/compare_stages.py --section groups` verifica traiettorie
complete `N_RK=2/3/4` con gruppi ANALYTIC e NEWTON ordinati in entrambi i modi.
Gli oracoli controllano il solve Newton della temperatura, il flusso esatto
della velocità con coefficienti congelati, la ricorrenza degli stage e
l’assemblaggio finale; l’ordine dei gruppi produce una differenza misurabile.
Quest’ultima suite passa in strict e optimized con thread 1/4.

### D16 — Leggi fisiche, rate e assemblaggio IMEX

**Stato ufficiale: COMPLETATO nel perimetro del gate locale e della CI ifx documentata; il contratto definitivo delle sorgenti è D17.** D16 mantiene la responsabilità della separazione delle leggi e della consistenza dei rate, mentre D17 definisce catalogo canonico, quattro interfacce uniformi, gruppi, scheduler e transazione mass exchange. I gate D16 non costituiscono da soli l'accettazione scientifica D15 della nuova architettura complessiva.

#### Obiettivi scientifici e numerici

Le formule dei processi fisici e reologici non sono proprietà dell'integratore né dei risolutori lineari. Ogni legge valuta contributi all'evoluzione della miscela multifase nello stato fisico locale, preservando unità, convenzioni dei segni, dipendenze su densità/temperatura/frazioni, limiti di ammissibilità e interpretazione dei regimi costitutivi. La reologia deve poter variare mantenendo invariati il solve Newton per gruppi e il tableau RK. Gli addendi costitutivi di una stessa legge rimangono all'interno della famiglia fisica, senza diventare termini configurabili indipendenti.

La conversione fra conservative e primitive e l'insieme di proprietà condivise rimangono in `state_conversion_2d`. Il percorso di derivazione delle sorgenti conserva i kernel REAL/COMPLEX e la valutazione complex-step dei Jacobiani, compreso il trattamento esplicito delle funzioni non analitiche, dei rami a soglia e dei punti degeneri. Le chiamate LAPACK restano tipizzate e l'ordine dei blocchi deriva dagli indici di equazione, non dal nome di una reologia.

#### Organizzazione vigente della fisica

- Reologie: `src/physics/rheology_laws/` e dispatcher di famiglia `rheology_2d`; una voce canonica per legge completa, incluse componenti interne di attrito/drag. In una nuova legge si modificano esclusivamente modulo fisico, catalogo e collegamento statico della famiglia.
- Termica e pressione interstiziale: `src/physics/thermal_laws_2d.f90` e `pore_pressure_laws_2d.f90`; la sorgente di raffreddamento è attivata da `thermal_cooling_flag` indipendentemente da `rheology_model`, con default retrocompatibile quando il parametro non compare negli input più vecchi.
- Iniezioni di massa e geometria delle sorgenti: moduli in `src/physics/mass_injection_terms/`; supportano sorgenti circolari, fessurali, radiali/laterali, finestre temporali e maschere senza ridefinire il contributo idrostatico dell'operatore HP-PCCU.
- Trasferimenti flusso/letto: processi distinti in `src/physics/mass_exchange/`, coordinati da `mass_exchange_laws_2d`; la perdita indotta dalla deposizione è interna alla deposizione e il limitatore sul carrier loss viene applicato a contributi sommati. Le variazioni conservative e topografiche sono assemblate insieme senza aggiornare direttamente lo stato fisico.
- Assemblaggio: `source_terms_2d` espone rate espliciti, semi-impliciti, impliciti e operazioni analitiche; `equation_terms_2d` conserva i soli contributi spaziali pertinenti. Le regole di cadenza, gruppi e schemi sono riportate in D17.

#### Sequenza numerica e regole di compatibilità

Per ciascuno stadio RK, gli addendi spaziali, espliciti e impliciti vengono assemblati senza doppio conteggio e senza moltiplicare arbitrariamente i rate per coefficienti o intervalli temporali a livelli fisici. Il solver Newton compone tutti i contributi delle equazioni del gruppo e applica le correzioni specifiche dei blocchi nei punti post-SI e post-Newton prescritti; la finalizzazione del passo appartiene a `state_2d`. Dopo il completamento RK, l'esecutore STEP tratta i trasferimenti di massa come una transazione indivisibile, applica una sola volta rate/inventari e aggiorna `B_vertex` mediante `geometry_2d`.

Non cambiare per effetto di questo consolidamento: relazioni costitutive approvate, coefficienti dei tableaux IMEX RK, trattamento REAL/COMPLEX, convenzioni wet/dry, proprietà termodinamiche, maschere di sorgente e bilanci di erosione/deposizione. Ogni nuova proposta scientifica avrà una propria motivazione e un gate proporzionato al rischio.

#### Evidenze, interfacce e dipendenze successive

**Accettazione D16:** il consolidamento delle leggi, dei rate e del solve locale ha superato i gate locali e la CI ifx pertinenti documentati nel repository; le evidenze sono limitate alle revisioni e ai casi effettivamente eseguiti. **Accettazione D17:** il successivo contratto implementato è riassunto in D17, con due Actions verdi sul commit `0e545ac9`. **Accettazione complessiva:** D15 rimane il gate per confronti end-to-end selezionati con la baseline A11, restart/output, configurazioni miste e documentazione della validità scientifica. Non duplicare nel piano D16 nuove campagne di test già considerate nel gate D15.

### D17 — Sorgenti fisiche configurabili e gruppi di integrazione omogenei

**Stato ufficiale: IMPLEMENTATO — refactoring architetturale completato, con gate tecnico GNU/CI superato; accettazione scientifica end-to-end ancora da documentare in D15.** Riferimento: branch `feature/d17-source-refactor`, revisione `0e545ac9766f9633ad115fe240f5b3216f86578a`, due Actions manuali verdi confermate il 25 ottobre 2026. La successiva revisione `bc38fbc16c1c647f84a0654056e5573addfc4114` aggiorna unicamente la specifica. La chiusura tecnica non estende automaticamente la baseline scientifica A11 e non sostituisce il gate D15.

**Obiettivo.** Le sorgenti fisiche sono leggi dotate di identità canonica e capacità numeriche dichiarate; il sistema di integrazione conosce gruppi, equazioni, cadenze e metodi, senza codificare l'identità numerica della singola legge. Una reologia è una sola sorgente configurabile anche quando contiene più contributi costitutivi interni. I termini di trasferimento di massa che condividono vincoli di disponibilità e limitazione costituiscono un batch atomico. Le routine numeriche e le leggi rimangono separabili e verificabili.

**Perimetro implementato.** Il catalogo attuale contiene **24 sorgenti configurabili**, ognuna con **quattro procedure concrete** (`explicit`, `semi_implicit`, `implicit`, `analytic`), per **96 interfacce**. Un metodo non supportato dispone di una procedura di rifiuto esplicito: non è una capacità simulata da un risultato nullo. Le combinazioni abilitate sono dichiarate dai metadati e validate prima dell'avvio; la cadenza `STEP` ammette attualmente il mass exchange `EXPLICIT` come batch unico. Il kernel `rheology_newmodel` documenta l'estensibilità tecnica, non costituisce una nuova legge scientificamente calibrata.

#### D17.0 — Architettura normativa e responsabilità

- `source_catalog_2d`: elenco statico di leggi fisiche configurabili, ciascuna con `name`, `target`, `has_explicit`, `has_newton`, `has_analytic`. Nessuna identità fisica trasportata come ID numerico persistente, né ruoli/bundle numerici di sotto-termine. Gli indici numerici restano ammessi per equazioni, gruppi e metodi numerici, mai come identità fisica della legge.
- `source_registry_2d`: costruzione immutabile, copia controllata e validazione dei gruppi; `source_group_type` espone ID del gruppo, `legacy_id` ove necessario, `cadence`, `method`, `phase`, `execution_order`, `source_names(:)`, `equation_ids(:)`. Nessun `term_ids` e nessun sotto-solver schedulato nel gruppo.
- `source_input_2d` e `source_configuration_2d`: conversione dell'input testuale e della configurazione legacy in un registro validato una sola volta. I nomi delle leggi sono canonici e univoci; nomi inesistenti o incompatibili sono rifiutati prima della simulazione.
- `source_terms_2d`: ingressi comuni per l'assemblaggio delle sorgenti, senza governare i tableaux RK né modificare gli inventari del letto. `equation_terms_2d` conserva solo i termini appartenenti all'operatore spaziale e non gestisce le leggi fisiche configurabili.
- `nonlinear_solver_2d`: risoluzione locale dei blocchi di equazioni mediante indici, Jacobiana REAL/COMPLEX e complex-step, senza dispatch numerico della specifica reologia. `time_integration_2d` gestisce schemi e successione degli stadi, non controlli fisici hard-coded.
- `state_2d`: stato secco/quasi secco, controlli e finalizzazione; le conversioni fisiche rimangono in `state_conversion_2d`. `step_source_executor_2d` applica una sola volta la transazione di massa, gli inventari e l'aggiornamento topografico post-RK.

#### D17.1 — Contratto fisico uniforme e capacità

Ogni legge configurabile possiede quattro procedure pubbliche con nomi e responsabilità uniformi:

| Trattamento | Convenzione | Semantica |
|---|---|---|
| Esplicito | `eval_<nome>_explicit` | Restituisce il rate senza mutare lo stato |
| Semi-implicito | `eval_<nome>_semi_implicit` | Restituisce i termini SI necessari al metodo Newton |
| Implicito | `eval_<nome>_implicit` | Restituisce rate/residui IM, coerenti nel ramo REAL e COMPLEX |
| Analitico | `solve_<nome>_analytic` | Aggiorna le sole equazioni target secondo un solve analitico ammesso |

La capability `NEWTON` richiede congiuntamente le interfacce semi-implicita e implicita, anche se per una legge specifica la parte SI risulta identicamente zero. Le procedure di trattamenti non supportati generano errore in caso di invocazione: il validatore preventivo impedisce di schedularle. I nomi rimangono coerenti in reologie, termica, pore pressure, iniezioni, curvatura e trasferimenti di massa; la struttura delle chiamate non deve essere duplicata negli stadi numerici.

**Modello di estensione ordinaria.** Per aggiungere `rheology_newmodel` o una nuova reologia: creare la legge completa in `src/physics/rheology_laws/`, fornire le quattro procedure e i rami REAL/COMPLEX pertinenti, aggiungere la voce canonica al catalogo e il collegamento statico nella famiglia reologica. Non cambiare `source_registry_2d`, `source_configuration_2d`, `source_terms_2d`, `nonlinear_solver_2d` né `time_integration_2d`. L'assenza di nuove eccezioni in questi moduli è un criterio architetturale verificabile.

#### D17.2 — Regole dei gruppi e cadenze

Un gruppo ha un solo metodo (`EXPLICIT`, `NEWTON` o `ANALYTIC`), una sola cadenza (`STAGE` o `STEP`), un ordine definito e l'unione coerente delle equazioni target. La stessa legge fisica attiva non può essere presente in più gruppi. Sono vietati nomi derivati, duplicati, target incompatibili e frazionamenti che cambiano il significato di un termine fisico indivisibile.

- `EXPLICIT/STAGE`: contribuisce a ciascuno stadio secondo i coefficienti dello schema, includendo iniezioni, curvature ed eventuali termini espliciti ammessi. I termini idrostatici dell'operatore spaziale non vengono reinseriti come sorgenti.
- `NEWTON/STAGE`: un **unico sistema Newton accoppiato per gruppo** e tutte le equazioni incluse; il gruppo può comprendere contemporaneamente velocità e temperatura. Gruppi Newton distinti sono ammessi quando configurati e validati, ma i loro solve sono sequenziali e non vengono dichiarati equivalenti al solve monolitico.
- `ANALYTIC/STAGE`: esattamente un solve analitico per gruppo. Reologia, raffreddamento e pore pressure possono richiedere gruppi distinti.
- `EXPLICIT/STEP`: batch post-RK atomico, oggi per il mass exchange, eseguito come ultimo gruppo nella configurazione legacy. Le altre combinazioni `STEP` non supportate sono respinte esplicitamente, non implementate implicitamente.

Il parser `SOURCE_GROUPS` deve conservare ordine, metodo, target ed esclusività. Senza `SOURCE_GROUPS`, la normalizzazione legacy conserva per i casi supportati la partizione dei blocchi implicitamente accoppiati, l'ordine e la scelta dei metodi IMEX; i gruppi misti analitici sono separati in solve singoli. Le liste di nomi di sorgenti e di equazioni sono immutabili durante il passo.

#### D17.3 — Assemblaggio Newton accoppiato, Jacobiane e correzioni numeriche

Il solve Newton usa `group_equations(:)` per costruire un residuo e una Jacobiana congiunti. Tutte le leggi appartenenti al gruppo devono contribuire alle componenti pertinenti; non è ammesso risolvere indipendentemente più sottoblocchi del medesimo gruppo. Si conserva il solve indicizzato della baseline, incluso il percorso linearizzato, le interfacce LAPACK tipizzate, il Jacobiano complex-step e il riuso dei kernel fisici REAL/COMPLEX. Sono coperti numericamente i blocchi 2×2 e 3×3 e l'accoppiamento reologia/raffreddamento nel perimetro dei test disponibili.

La correzione di direzione/verso della quantità di moto è un'operazione sul **blocco di equazioni**, non un'opzione dipendente dal nome del modello reologico. Si applica in due momenti distinti: **dopo SI**, con riferimento allo stato noto dello stadio, e **dopo Newton**, con riferimento allo stato post-SI/pre-Newton. Rimane in `state_2d` ed è richiamata dal solver nei punti prescritti, senza attivazione automatica nella finalizzazione generica. Conservare dissipazione, limite di arresto, casi SI nullo e assenza di inversioni spurie della quantità di moto.

#### D17.4 — Operatori espliciti, analitici e raffreddamento termico

`source_terms_2d::eval_source_explicit` è l'ingresso pubblico comune: un'interfaccia generica Fortran con procedure concrete distinte per firme STAGE e STEP. L'assemblaggio STAGE deve includere iniezioni di fondo circolari/fessurali, sorgenti radiali/laterali e curvatura soltanto nei casi pertinenti, rispettando finestre, frazioni di cella e maschere, senza doppio conteggio del contributo idrostatico. L'assemblaggio STEP delega al coordinatore mass exchange la preparazione dei rate conservativi e accessori.

Il percorso `ANALYTIC` ammette un unico operatore coerente col metodo e con le equazioni target; gli accoppiamenti che richiedono Newton non vengono scomposti silenziosamente. Le leggi termiche e della pressione interstiziale rimangono indipendenti dalla legge reologica selezionata.

**Attivazione termica:** `thermal_cooling_flag` in `TEMPERATURE_PARAMETERS` controlla il raffreddamento indipendentemente dalla scelta lava/non lava. L'assenza del flag nei vecchi input conserva il comportamento di default compatibile; le combinazioni non supportate vanno rifiutate all'avvio. I contributi radiativo, convettivo, la capacità termica della miscela e la natura depth-averaged del termine devono rimanere coerenti anche con spessori estremamente piccoli; la soluzione implicita/analitica locale non implica di per sé positività di ogni stadio RK. Eventuali sviluppi dei metodi devono preservare i bilanci energetici e le capability dichiarate.

#### D17.5 — Trasferimenti di massa e coordinatore STEP

Il mass exchange è organizzato per processo in `src/physics/mass_exchange/`: erosione, deposizione, perdite di fase vettore, entrainment, vaporizzazione, assemblaggio di bilancio e coordinatore `mass_exchange_laws_2d`. Deposizione e relativa perdita indotta costituiscono **un unico termine fisico configurabile**, `mass_exchange_deposition`; il contributo `deposition_induced_carrier_loss` è calcolato internamente e **non** ha una voce separata nel catalogo. Restano configurabili gli altri processi indipendenti effettivamente supportati.

L'ordine delle operazioni fisiche è vincolante: ottenere i rate per classe e vettore; rispettare l'inventario erodibile e la disponibilità delle fasi; comporre le perdite di fase continua da liquid loss, perdita indotta dalla deposizione e Darcy; applicare **una sola volta** il limitatore post-somma e gli eventuali altri vincoli condivisi; assemblare massa totale, quantità di moto, energia, masse delle specie e tasso di variazione del fondo. Il batch non può essere spezzato in gruppi indipendenti che cambierebbero l'ordine dei rate o il limite dipendente da `dt`.

`step_source_executor_2d` chiama `source_terms_2d::eval_source_explicit` una sola volta per cella risolta, applica il rate conservativo, aggiorna gli inventari `deposit`, `erosion`, `erodible` e delega a `geometry_2d` il gather nodale Q1 e il refresh della geometria. Evitare doppio aggiornamento di `B_cent`, doppia valutazione dei rate e applicazione prima del completamento dello stato RK. Le celle secche e le maschere delle sorgenti non devono consumare inventario senza aggiornare coerentemente flusso e topografia.

La topografia nodale `B_vertex` resta autoritativa e viene aggiornata dopo il gather delle proposte locali limitate; le quote/pendenze/derivate vengono poi rigenerate. Il batch legacy rimane l'ultimo `STEP/EXPLICIT` e deve preservare l'ordine temporale rispetto alle iniezioni, alla finalizzazione del passo e all'avanzamento del tempo simulato.

#### D17.6 — Stato fisico e ordine completo del timestep

La sequenza di riferimento è: configurazione validata all'avvio → ricostruzione e contributi spaziali → per ciascuno stadio RK, contributi espliciti e gruppi impliciti/analitici previsti dal tableau → correzioni post-SI/post-Newton ove applicabili → finalizzazione dello stato RK e conversioni primitive → esecuzione ordinata dei gruppi STEP, incluso il batch atomico mass exchange per ultimo → aggiornamento geometrico nodale/derivate e osservatori pertinenti → avanzamento del tempo e operazioni di output/restart secondo la policy globale.

`time_integration_2d` controlla tableau e scheduling; `state_2d` conserva controlli numerici e riparazioni fisicamente definite su stato secco, solidi, temperatura e roundoff; `state_conversion_2d` conserva `q↔qp`. Inizializzazione di primitive secche, halo e stati intermedi usa la medesima semantica. Le correzioni numeriche nei due momenti del solve implicito restano distinte dalla finalizzazione dello stato.

Questa fase **non** implementa time-step rejection/rollback o nuove policy CFL: tali funzioni appartengono alla fase E. Non cambia i tableaux RK approvati né introduce genericamente una capability STEP per metodi impliciti o analitici.

#### D17.7 — Default compatibili e regole di input

La configurazione legacy supportata, quando `SOURCE_GROUPS` non è presente, si traduce in gruppi omogenei con le sorgenti effettivamente abilitate dai flag, senza nuove opzioni imposte all'utente. La matrice di default da mantenere e documentare è:

| Famiglia | Fonte di attivazione legacy | Target | Metodo/cadenza |
|---|---|---|---|
| Reologia completa (una voce per legge) | `rheology_flag` e `rheology_model` | quantità di moto x/y | `NEWTON/STAGE` o `ANALYTIC/STAGE` secondo gruppo legacy |
| Raffreddamento | `thermal_cooling_flag`, fallback di input legacy | energia | `NEWTON/STAGE` o `ANALYTIC/STAGE` se ammesso |
| Dissipazione pressione interstiziale | `pore_pressure_flag` | equazione pore pressure | `NEWTON/STAGE` o `ANALYTIC/STAGE` |
| Curvatura e iniezioni | flag fisici e finestre locali | equazioni conservative pertinenti | `EXPLICIT/STAGE` |
| Erosione del letto | `erosion_coeff` | massa, specie, energia, letto | `EXPLICIT/STEP`, batch unico |
| Deposizione, inclusa perdita indotta | `settling_flag` e porosità | massa, specie, energia, letto | `EXPLICIT/STEP`, batch unico |
| Perdita liquida, Darcy, entrainment, vaporizzazione | flag/perdite dei processi | variabili conservative e accessori pertinenti | `EXPLICIT/STEP`, batch unico |

Per le reologie senza componente fisica attiva e per i layout disabilitati non si devono inventare rate né introdurre sorgenti fittizie. Un gruppo Newton può contenere velocità e temperatura **in un solo sistema accoppiato**, se l'input lo specifica. Un raggruppamento di metodi/cadenze non supportato, una sorgente ignota, un duplicato, una fonte derivata schedulata o una separazione del batch mass exchange è una configurazione non valida. Il parser deve emettere una diagnostica prima dell'avvio del tempo fisico.

#### D17.8 — Struttura del codice e invarianti dei moduli

| Modulo o directory | Responsabilità stabile | Divieto architetturale |
|---|---|---|
| `src/numerics/source_catalog_2d.f90` | Catalogo statico 24 sorgenti, target, tre capability | ID numerici fisici e sotto-componenti come fonti indipendenti |
| `src/numerics/source_registry_2d.f90` | Gruppi, ordine, indici di equazione e validazione | Sotto-solver Newton per gruppo |
| `src/numerics/source_input_2d.f90` | Input `SOURCE_GROUPS` | Valutazione fisica nei parser |
| `src/physics/source_configuration_2d.f90` | Normalizzazione legacy e attivazione sorgenti | Dipendenze dalla formula interna di ogni reologia |
| `src/physics/source_terms_2d.f90` | Assemblatori e dispatch comune per trattamenti | Aggiornare inventari fisici o duplicare lo STEP |
| `src/physics/rheology_laws/` | Leggi reologiche complete | Scheduling RK o solve di altri modelli |
| `src/physics/mass_exchange/` | Processi, `mass_exchange_laws_2d`, rate e limitatori | Applicazione diretta del passo alla griglia |
| `src/numerics/nonlinear_solver_2d.f90` | Newton indicizzato e correzioni collegate a blocchi | Selezione per ID di legge fisica |
| `src/numerics/time_integration_2d.f90` | Stadi e scheduling | Riparazioni costitutive hard-coded |
| `src/numerics/step_source_executor_2d.f90` | Transazione STEP, inventari, chiamata di geometria | Ricalcolo delle formule costitutive |
| `src/state/state_2d.f90`, `src/state/state_conversion_2d.f90` | Ammissibilità, dry/halo, `q↔qp` | Modelli reologici hard-coded per correzioni di blocco |
| `src/geometry/geometry_2d.f90` | Topografia Q1, gather e refresh | Quote aggiornate per vie indipendenti |
| `src/physics/equation_terms_2d.f90` | Residui/termini propri dell'operatore spaziale | Routing delle sorgenti fisiche configurabili |

Le directory indicano l'organizzazione attuale del codice. Gli step futuri G/H devono operare su layout e domini astratti e non reintrodurre dipendenze numeriche dall'identità della sorgente.

#### D17.9 — Piano di implementazione e stato operativo

Questa tabella è la **sequenza normativa di realizzazione**: descrive i componenti della roadmap e il loro stato, non la cronologia dei tentativi.

| Passo | Attività implementativa e accettazione | Stato |
|---|---|---|
| 1 | Inventario delle sorgenti, catalogo e target; confronto con baseline e riferimenti approvati | **IMPLEMENTATO**; riferimento scientifico A11/D15 separato |
| 2 | Contratto un gruppo/un metodo/una cadenza/un solve e gruppi immutabili | **COMPLETATO — CI** |
| 3 | Registro e parser senza sotto-solver e senza ID fisici | **COMPLETATO — CI** |
| 4 | Normalizzazione legacy e validazione preventiva dei nativi | **COMPLETATO — CI**, per configurazioni supportate |
| 5 | Assemblaggio SI/IM congiunto su tutte le sorgenti del gruppo | **COMPLETATO — CI**, percorsi coperti |
| 6 | Newton indicizzato unico, 2×2/3×3, complex-step e correzioni post-SI/post-Newton | **COMPLETATO — CI**, casi mirati |
| 7 | Solve analitico unico per gruppo, senza miscelare metodi | **IMPLEMENTATO — CI** |
| 8 | Sorgenti `EXPLICIT/STAGE`, boundary e finestre; raffreddamento indipendente | **IMPLEMENTATO — CI**, matrice end-to-end non esaustiva |
| 9 | Quattro trattamenti fisici uniformi, rifiuto dei non supportati, estensibilità `rheology_newmodel` | **COMPLETATO — CI** |
| 10 | Processi mass exchange per modulo, deposizione/perdita interna, limitatore unico, coordinatore e ingresso EXPLICIT comune | **COMPLETATO — CI** |
| 11 | Separazione `state_2d` e finalizzazione post-RK/STEP, revisione hot loop | **COMPLETATO — CI**, regressioni observer mirate |
| 12 | Regression testing risk-based su default/nativo, RK, bilanci, stato e source registry | **IMPLEMENTATO PER GATE TECNICO — CI**; equivalenza end-to-end scientifica da documentare in D15 |
| 13 | Accettazione scientifica della nuova architettura contro reference approvate, output/restart e casi end-to-end | **GATE D15 — DA COMPLETARE** |
| 14 | Manutenzione di compatibilità D17 di tutti i test storici, degli esempi, dei runner, dei generatori e dei template di input; migrazione delle API abolite senza alterare asserzioni o casi scientifici | **IN CORSO — MIGRAZIONE CODICE EFFETTUATA IN PARTE; VERIFICA COMPILAZIONE/ESECUZIONE APERTA** |
| 15 | Inventario esaustivo degli input generati, smoke test degli esempi, regression testing selettivo di mass exchange/raffreddamento, output e restart | **DA COMPLETARE — NESSUNA CERTIFICAZIONE GLOBALE** |

#### D17.9a — Compatibilità delle suite storiche, esempi e generatori

L'architettura D17 è implementata; la compatibilità dell'intero ecosistema di verifica e distribuzione costituisce una **attività di manutenzione necessaria prima dell'accettazione complessiva D15**. Tutti i nuovi test devono usare nomi canonici e `source_names`, gli entry point di `source_terms_2d`, un solo solve per gruppo, il coordinatore `mass_exchange_laws_2d` e i runner allineati all'elenco dei sorgenti di `src/Makefile.am`. Non ripristinare `source_solver_type`, `term_ids` o vecchi `SOURCE_TERM_*` come compatibilità interna dei test.

**Compatibilità termica delle fixture:** in `TEST_RHEOLOGY_CELL_LOCAL/probe_sources.f90` e `test_grouped_implicit_stage.f90` il raffreddamento termico è ora abilitato esplicitamente tramite `thermal_cooling_flag=.TRUE.`. D17 separa questo switch dalla scelta reologica; ometterlo lascerebbe i test di sorgente termica a verificare un contributo nullo invece della legge originariamente oggetto del confronto. Resta necessaria l'esecuzione strict/optimized delle suite per confermare i valori.

**Stato al 25 ottobre 2026:** aggiornati i test `TEST_IMEX_STAGES/test_stage_count.f90`, `test_imex_stages.f90`, `TEST_SPATIAL_OPERATOR/test_spatial_operator.f90`, `TEST_CURVATURE_SOURCE/test_curvature_source.f90`, `TEST_MASS_EXCHANGE/test_mass_exchange.f90`, `TEST_RHEOLOGY_FROUDE/test_rheology_froude.f90`, `TEST_RHEOLOGY_CELL_LOCAL/probe_sources.f90`, `test_grouped_implicit_stage.f90` e il relativo `run_test.sh`. La prova del gruppo implicito verifica adesso un unico Newton accoppiato, conservando filtri per sorgente, residuali e controlli REAL/COMPLEX. **Queste modifiche non sono ancora validate mediante una compilazione ed esecuzione completa**: non attribuire loro l'esito verde delle precedenti Actions D17 sul commit `0e545ac9`.

**Piano operativo residuo (criterio di uscita esplicito):** completare scansione e adeguamento di tutti i test, `run_test.sh`, script Python/Shell, generatori e input statici nelle directory `TESTS/`, `EXAMPLES/`, `tools/`, `UTILS/` e workflow pertinenti; confrontare gli input generati con i template per individuare placeholder o parametri D17 non supportati; conservare il fallback legacy dove previsto e controllare che il raffreddamento lavico e gli scambi di massa non cambino involontariamente; compilare ed eseguire le suite risanate, quindi lanciare smoke test selettivi di ciascuna famiglia di esempi, con casi negativi del parser, bilanci e restart per le superfici effettivamente modificate. Registrare per ogni suite/caso revisione, comando, toolchain, configurazione, risultato **PASS/FAIL/NON ESEGUITO** e scostamenti fisici; risolvere i fallimenti senza sopprimere verifiche scientifiche o alterare silenziosamente i default. I confronti end-to-end richiesti per l'approvazione finale rimangono nel gate D15.

**Report tecnico collegato:** `docs/acceptance/D17_EXAMPLES_TESTS_INPUT_GENERATORS_AUDIT_2026_10_25.md`. Il report contiene gli errori rinvenuti e la loro tracciabilità; **questa roadmap è la fonte autoritativa dello stato e delle attività da completare**. Gli errori rilevati nella baseline del report non devono essere interpretati automaticamente come ancora presenti nei file già migrati.

**Esiti delle sei suite storiche (25 ottobre 2026; run locale sul commit e4080c9b):** `TEST_MASS_EXCHANGE` **PASS** in strict/optimized a 1/4 thread, con 18 record diagnostici per profilo e continuità dei checkpoint; `TEST_IMEX_STAGES` **FAIL compilazione** (import boundary e costruttori di stringhe); `TEST_RHEOLOGY_CELL_LOCAL` **FAIL link** per elenco parziale degli oggetti; `TEST_CURVATURE_SOURCE` e `TEST_RHEOLOGY_FROUDE` **FAIL build** per moduli/percorso precedenti; `TEST_SPATIAL_OPERATOR` **BLOCCATO dall'ambiente Python** (NumPy non disponibile), dopo compilazione/link Fortran. Gli script e i sorgenti sono stati corretti per le cause individuate nel successivo commit di bonifica. **Non dichiarare PASS prima di nuove esecuzioni sui file corretti.** I runner curvature e Froude usano ora la lista dei sorgenti canonici da `src/Makefile.am`; il link reologico include tutte le dipendenze D17; i generatori Python possono selezionare `PYTHON` dall'ambiente. La visione scientifica finale e gli altri gate restano invariati.

**Seconda verifica delle suite storiche — 25 ottobre 2026 (commit `66fd58eb`):** `TEST_MASS_EXCHANGE`, `TEST_CURVATURE_SOURCE` e `TEST_RHEOLOGY_FROUDE` **PASS**; mass exchange comprende le modalità strict/optimized e 1/4 thread e Froude le modalità strict/optimized. `TEST_IMEX_STAGES` ha compilato i driver Fortran ma si è interrotto nell'import del riferimento Python, prima dei confronti IMEX; il compilatore ha inoltre segnalato tre costruttori di stringhe con lunghezza disomogenea. `TEST_SPATIAL_OPERATOR` ha incontrato lo stesso errore Python: il riferimento usa `np.ndarray | None`, incompatibile con l'interprete Python dell'ambiente corrente. `TEST_RHEOLOGY_CELL_LOCAL` è stato interrotto da un import rimasto in `test_source_modules.f90` (`correct_friction_momentum` da `rheology_2d` anziché da `state_2d`). Le tre cause sono state corrette nei file di test/riferimento, mantenendo invariati algoritmi, tolleranze e baseline; le tre suite richiedono **nuova esecuzione**, non sono ancora PASS. La compatibilità della sintassi dei type hint è ottenuta con `from __future__ import annotations`, senza alterare l'implementazione numerica.

**Terzo ciclo di verifica (25 ottobre 2026, commit `4267bcc8`):** `TEST_IMEX_STAGES` e `TEST_SPATIAL_OPERATOR` terminano nel controllo SHA-256 del riferimento HP-PCCU 1D: l'aggiunta dell'import `from __future__ import annotations`, necessaria al Python locale, cambia l'hash del file pur senza cambiare il suo algoritmo. Il controllo d'integrità resta obbligatorio; l'hash approvato è stato aggiornato al valore SHA-256 verificato `c35f8ec13f1b0c8e93eb1e31c02e6308bcf679bf182bd00eb265ca3f54ee7f38` per questa sola modifica di compatibilità delle annotazioni. `TEST_RHEOLOGY_CELL_LOCAL` avanza alla compilazione del test analitico ma trova un ID reologico numerico `HUGE(1)` usato come caso negativo; il test conserva il controllo di rifiuto sostituendolo con il nome non supportato `rheology_unsupported`, coerente col dispatcher D17. Queste tre suite restano **DA RIESEGUIRE**; non si dichiarano PASS. Le prove mass exchange, curvature e Froude erano PASS nei log precedenti, non sono state ripetute in questo ciclo.

**Quarto ciclo di verifica (25 ottobre 2026, commit `b904a943`):** `TEST_SPATIAL_OPERATOR` **PASS** col riferimento HP-PCCU e checksum attivo. `TEST_IMEX_STAGES` supera i controlli RK e il rifiuto di `N_RK` non supportati; fallisce ancora la fixture Python `source-groups-native` perché emette `SOURCE_GROUPS` nel vecchio formato a sotto-solver, abolito da D17. La fixture è stata migrata ai target di gruppo `source_group_target_count` / `source_group_targets` e a nomi di sorgenti canonici mantenendo tre gruppi e tutte le verifiche di fase/bilancio. `TEST_RHEOLOGY_CELL_LOCAL` fallisce nel confronto analitico perché passava nomi reologici legacy al dispatcher D17: aggiornati i nomi dei casi in `rheology_*` preservando le attese numeriche e i test negativi. Le due suite modificate sono **DA RIESEGUIRE**, non PASS.

**Quinto ciclo di verifica (25 ottobre 2026, commit `c521b9e8`):** `TEST_IMEX_STAGES` **PASS** nell'esecuzione caricata (inclusi il formato nativo `SOURCE_GROUPS`, i controlli `N_RK`, l'ordine degli stadi e i bilanci). `TEST_RHEOLOGY_CELL_LOCAL` ha superato i controlli di flusso reologico analitico, raffreddamento lavico e diffusione della pressione interstiziale, ma **NON PASS**: la compilazione di `test_grouped_implicit_stage.f90` si interrompe su due array constructor `CHARACTER` con type-spec all'interno della concatenazione. Le due righe e gli stessi errori latenti in `probe_sources.f90` sono stati corretti con un unico type-spec esterno `[CHARACTER(LEN=64) :: active_names, '...']`. Nessuna soglia o formula fisica è stata modificata. `TEST_RHEOLOGY_CELL_LOCAL` necessita nuova esecuzione integrale strict/optimized; per le altre cinque suite rimangono validi i PASS documentati nei cicli recenti, non sono state ripetute qui.

**Sesto ciclo (25 ottobre 2026, commit `327a45b7`):** `TEST_RHEOLOGY_CELL_LOCAL` supera compilazione, test analitici reologici, raffreddamento lavico e diffusione della pressione interstiziale, ma si ferma nel controllo `joint Newton group residual` del test raggruppato. L'audit ha individuato che l'asserzione ricostruiva l'equazione di stage dalla sola componente implicita, mentre `solve_imex_stage` applica anche il contributo semi-implicito della reologia lavica (prima del Newton), e può correggere la quantità di moto. La ricostruzione del residuo della fixture è stata aggiornata per includere esplicitamente il rate SI del gruppo. Il controllo resta attivo, usa una tolleranza relativa `2e-5` già usata dagli altri test di residuo accoppiato di questa suite, ed emette il valore per equazione in caso di fallimento. **Da confermare con esecuzione strict/optimized; non dichiarare PASS.** I cinque PASS delle altre suite rimangono documentati.

**Settimo ciclo di verifica (25 ottobre 2026, commit `2369371a`):** `TEST_RHEOLOGY_CELL_LOCAL` ora compila e supera i test esatti di flusso reologico, raffreddamento, diffusione di pressione e il controllo Newton/DIRK congiunto. Il caratterizzatore indipendente `characterize.py` segnala però **162/255 fixture PASS** in strict e optimized, con 93 fallimenti concentrati su `implicit_oracle_error` (oltre a verifiche di stage e a riposo). La diagnosi ha isolato un errore sistematico della migrazione della fixture: `probe_sources.f90` forzava `thermal_cooling_flag=.TRUE.` per ogni modello, mentre l'oracolo scientifico indipendente include il raffreddamento solo quando `rheology_model='lava'`. Dopo la lettura della namelist la fixture ora configura esplicitamente il raffreddamento per la sola reologia lavica, preservando il comportamento fisico originariamente caratterizzato senza reinserire dipendenze nel codice di produzione. **Nessuna tolleranza o oracolo modificato.** Nuova esecuzione strict/optimized ancora necessaria; accettazione della suite **APERTA**.

**Ottavo ciclo di verifica (25 ottobre 2026, correzione verificata):** l'audit del test `TEST_RHEOLOGY_CELL_LOCAL` ha individuato nel percorso senza gruppo di `eval_source_implicit` in `src/physics/source_terms_2d.f90` la selezione `selected_thermal=rheology_flag`, incompatibile con la separazione termica prevista da D17.4. Corretto il percorso a `selected_thermal=thermal_cooling_flag`; la configurazione per gruppo resta invariata. La modifica non reintroduce ID reologici numerici né adattatori legacy e non modifica formule fisiche, fixture, oracoli o tolleranze. **Verifica completata sul commit `a54f7596` (log locale `TEST_RHEOLOGY_CELL_LOCAL(7).log`, 25 ottobre 2026):** `D02/D03 strict: 255/255 fixtures pass; acceptance=PASS` e `D02/D03 optimized: 255/255 fixtures pass; acceptance=PASS`. Tutti i 93 fallimenti precedenti risultano risolti in entrambi i profili. Superati anche i controlli dei flussi analitici reologici, raffreddamento esatto della lava, diffusione esatta della pressione interstiziale e contratti grouped DIRK/effective IMEX stage-rate. Restano solo warning di link per `-lnetcdf` duplicato. L'utente ha inoltre riferito che entrambe le Actions sono verdi; non si attribuisce a queste Actions la copertura della suite locale. **Gate `TEST_RHEOLOGY_CELL_LOCAL`: PASS strict/optimized; accettazione complessiva D17 ancora APERTA.**

**Nono ciclo — gate CI riproducibile delle suite storiche (25 ottobre 2026, PENDING CI):** aggiunto `.github/workflows/d17-historical-regressions.yml` per eseguire come job indipendenti le sei suite `TEST_IMEX_STAGES`, `TEST_SPATIAL_OPERATOR`, `TEST_MASS_EXCHANGE`, `TEST_CURVATURE_SOURCE`, `TEST_RHEOLOGY_FROUDE` e `TEST_RHEOLOGY_CELL_LOCAL`. Ciascun job usa lo `run_test.sh` originale, senza disattivare assert, cambiare soglie o introdurre percorsi legacy; installa GNU Fortran, netCDF, LAPACK, Python/NumPy, registra SHA e log come artifact anche in caso di fallimento. `fail-fast: false` consente di raccogliere tutti gli esiti. Il workflow è avviabile **esclusivamente manualmente** tramite `workflow_dispatch`; nessun trigger `push` o `pull_request` è configurato. **Gate CI VERIFICATO:** run [40175887639](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/40175887639), commit `cb129420d94b8ffd913a5723be3349cac735ca14`, conclusione `success`: ciascuno dei sei job (`TEST_IMEX_STAGES`, `TEST_SPATIAL_OPERATOR`, `TEST_MASS_EXCHANGE`, `TEST_CURVATURE_SOURCE`, `TEST_RHEOLOGY_FROUDE`, `TEST_RHEOLOGY_CELL_LOCAL`) ha stato `completed/success`, con step `Run existing suite (all profiles)` e caricamento artifact `success`. **Stato del gate delle sei suite: PASS (6/6 job CI).** Il run è stato avviato da `push` prima della successiva modifica del workflow che ha rimosso tale trigger; il workflow aggiornato è ora esclusivamente `workflow_dispatch`. Il risultato rimane riferito al commit testato, non costituisce prova di nuovi run su commit successivi. I PASS locali precedentemente documentati restano evidenze storiche distinte. Il gate non sostituisce simulazioni end-to-end e criteri D15.\n\n**Matrice CI delle sei suite storiche (25 ottobre 2026; run 40175887639):** `TEST_IMEX_STAGES` **PASS**; `TEST_SPATIAL_OPERATOR` **PASS**; `TEST_MASS_EXCHANGE` **PASS**; `TEST_CURVATURE_SOURCE` **PASS**; `TEST_RHEOLOGY_FROUDE` **PASS**; `TEST_RHEOLOGY_CELL_LOCAL` **PASS**. In tutti e sei i job è superato lo step `Run existing suite (all profiles)`. Evidenza specifica del commit `cb129420` su Ubuntu 24.04 con GNU Fortran; le sei suite usano i runner esistenti. Questo chiude il gate di compatibilità delle sei suite storiche, **non** il piano completo di audit di tutti i test/esempi/generatori o l'accettazione scientifica D15/D17.

**Decimo ciclo — audit trasversale esempi e generatori (25 ottobre 2026, STATIC/PARTIAL):** inventario dell'albero Git (618 elementi, 122 file sotto `EXAMPLES/`) e controllo mirato delle dipendenze D17 di generatori, template, lancio ensemble Etna e precedenti evidenze C06. Report: `docs/acceptance/D17_EXAMPLES_GENERATORS_FULL_AUDIT_2026_10_25.md`. Individuati rischi concreti nel generatore ensemble (`create_inputfiles.py` riscrive ogni file regolare anche senza placeholder, senza filtro testuale), propagazione degli errori in `templatedir/run.sh`, provisioning dell'eseguibile e deriva documentale del comando DAMBREAK. Nessun run di esempio sul D17 HEAD è stato eseguito: **audit statico parziale, gate dinamico esempi/generatori APERTO**. I 20 esempi PASS del candidato C06 precedente restano storici e non certificano D17. Prossimo passo: bonifica coordinata dell'ensemble e matrice smoke test manuale con tracciabilità input, output e bilanci.

**Undicesimo ciclo — bonifica ensemble Etna (25 ottobre 2026, PENDING TEST):** `create_inputfiles.py` sostituisce i placeholder nei soli file testuali previsti e solo se cambiati, preservando symlink, file binari e permessi; usa sostituzione atomica, blocca token sconosciuti e mantiene `ENSEMBLE_vx/vy` per il preprocessing successivo. Le directory incomplete sono rimosse quando la generazione fallisce. `launch_jobs.py` controlla prima i permessi/esistenza di script ed eseguibile e conserva l'attesa dei processi e i codici di uscita; `templatedir/run.sh` propaga i fallimenti di preprocessing con `set -euo pipefail`. **Verifica funzionale comunicata dall'utente sul commit `1d60cc8f`: ESITO POSITIVO** per il workflow ensemble dopo la bonifica. Non sono stati forniti log, elenco dei membri eseguiti o risultati separati dei casi negativi: l'esito non certifica automaticamente i test dei guasti intenzionali, la copertura di tutti i membri o una nuova campagna numerica end-to-end. Tali gate rimangono APERTI.

**Dodicesimo ciclo — gate comune dei 19 esempi rimanenti (25 ottobre 2026, PENDING CI):** creati `tools/d17/audit_example_inputs.py` e `.github/workflows/d17-example-input-audit.yml`, una Action **solo manuale** con 19 job indipendenti (compresi i quattro BUMP, escluso l'ensemble Etna già trattato). Ogni job controlla README, presenza del template/input, struttura delle namelist, placeholder residui negli input statici, simboli interni aboliti e sintassi Python dei generatori; salva JSON e log per caso. **Questo gate è dichiaratamente STATIC_ONLY**: non esegue generatori o solver, non controlla l'equivalenza numerica, la semantica dei fallback né i bilanci; i casi da eseguire restano NOT_RUN finché non saranno disponibili risultati dinamici sul commit D17. Non cambiare configurazioni scientifiche, default o tolleranze per ottenere PASS.

#### D17.10 — Verifiche e criteri di accettazione

**Evidenza disponibile sul refactoring:** Actions di compilazione GNU e test source registry verdi sul commit `0e545ac9`; 24 sorgenti/96 interfacce; parser e registro; prova reologica estensibile; Newton accoppiato locale, SI/IM REAL/COMPLEX e complex-step; correzioni di direzione; bilanci conservativi locali; test di stato e dry/halo. Le verifiche mirate dei coefficienti RK non costituiscono un confronto generale dell'intero tableau su griglie reali.

**Regressione reologica locale D02/D03 (25 ottobre 2026):** `TEST_RHEOLOGY_CELL_LOCAL` supera 255/255 fixture in `strict` e 255/255 in `optimized` dopo la correzione di `selected_thermal` in `source_terms_2d.f90` (`a54f7596`). Evidenza: log locale `TEST_RHEOLOGY_CELL_LOCAL(7).log`; il risultato non sostituisce i gate scientifici/end-to-end ancora aperti.

**Gate scientifico ancora pertinente (D15), senza campagne indefinite:** selezionare secondo T2/T3 i casi realmente rappresentativi della modifica, precisando revisioni, input, tolleranze, toolchain e baseline A11/precedente codice approvato. Valutare almeno, se richiesto dalle superfici effettivamente toccate: evoluzione di stati `Q/QP` e dry/wet, `n_RK=2,3,4`, finestre circolari/fessurali/radiali e BC, OU con restart, termica lava/non lava, mass exchange con limitatore saturo/non saturo, bilanci del letto, output full/split/restart e determinismo con thread pertinenti. Non attribuire a D17 i gate scientifici di nuove leggi o i futuri confronti verticali/subgrid.

**Criteri D17 tecnici:** (i) 24 nomi canonici, nessun ID fisico; (ii) quattro entry point per legge e capability validate; (iii) un solve Newton o analitico secondo metodo per gruppo; (iv) residuale/Jacobiana su tutte le equazioni target con complex-step; (v) ingresso EXPLICIT comune STAGE/STEP; (vi) deposizione e perdita indotta inseparabili, un limiter aggregato; (vii) inventari e topografia aggiornati una sola volta; (viii) input legacy/nativo supportato e rigetto preventivo dei casi illegali; (ix) le due correzioni della quantità di moto e la gestione di dry state sono preservate; (x) codice compilato e due Actions verdi. **Questi criteri sono soddisfatti nel perimetro tecnico CI verificato.**

**Limite accettativo:** le regressioni complete end-to-end e l'approvazione delle baseline numeriche del nuovo HEAD non risultano documentate come gate scientifico finale. Fino alla decisione di D15, utilizzare lo stato **IMPLEMENTATO (GATE SCIENTIFICO APERTO)** e non dichiarare equivalenza numerica universale né chiusura senza riserve di ogni caso legacy.

**Dipendenza successiva:** D15 valuta la nuova architettura nel suo insieme; D03/D07 ed E01–E03 hanno scopi diversi (ownership globale e readiness) e non riaprono la struttura delle sorgenti se le API D17 restano valide.

### D18 — Configurazione per sorgente di metodo, cadenza e gruppo; metodi impliciti e analitici a STEP

**Stato: IN CORSO — PARZIALMENTE IMPLEMENTATO; GATE D18 APERTO (audit 26 ottobre 2026, HEAD `558fd94`).** Le due GitHub Actions di compilazione e regressioni storiche sono state comunicate come superate dall'utente dopo `558fd94`. Questo verifica i percorsi coperti dalle Actions, **non** tutte le nuove combinazioni D18: il caso esplicito STAGE introdotto è presente nel codice, ma non è dimostrato da un test dedicato alla nuova configurazione. D18 rimane il riferimento **normativo** per le modifiche successive: la tabella di avanzamento riportata sotto separa ciò che esiste da ciò che è ancora obbligatorio. Nessuna formula costitutiva o tableau RK va modificato per completare D18.

**Contratto gruppi concordato il 28 ottobre 2026 (prevale sulle tabelle storiche D18.2, D18.4b e D18.6):** a STAGE, EXPLICIT e NEWTON possono condividere lo stesso ID numerico GROUP, ma restano due sottogruppi numerici distinti sotto il tableau IMEX. I gruppi impliciti/analitici sono eseguiti in ordine crescente di GROUP: ANALYTIC/STAGE seguito da NEWTON/STAGE fornisce al Newton lo stato gia aggiornato dalla mappa analitica. Una sorgente ANALYTIC deve essere l'unico membro del suo GROUP, a STAGE e STEP, e non puo condividere quel numero con alcun'altra sorgente della stessa CADENCE. A STEP un dato GROUP ammette un solo METHOD; metodi differenti con lo stesso GROUP vengono rifiutati, e i gruppi sono applicati in ordine numerico crescente senza priorita predefinite per METHOD o mass exchange. In un gruppo EXPLICIT tutti i rate usano lo stesso stato prima dell'aggiornamento; gruppi STEP diversi vedono gli stati aggiornati dai precedenti. Le regole storiche che raggruppano per (CADENCE,METHOD,GROUP) senza restrizioni o forzano mass exchange per ultimo sono superate. **Attenzione:** il default GROUP=0 del mass exchange comporta che, sotto questo contratto, esso preceda eventuali gruppi STEP con ID maggiore; i casi che richiedono ordine diverso devono configurarlo esplicitamente. La conformita dei default storici e degli inventari resta da verificare con CI.

**Regola di semplicità (vincolante per ogni commit D18):** riutilizzare `source_configuration_2d`, `source_registry_2d`, `equation_partition_type`, `source_terms_2d`, `step_source_executor_2d` e i moduli fisici già presenti; evolverli verso **un solo** percorso di configurazione normalizzato. Non introdurre registry paralleli, scheduler secondari, sottosolver per singola sorgente, nuovi identificativi fisici numerici o adattatori permanenti. Ogni nuovo ramo temporaneo deve avere una rimozione prevista. Eliminare campi e percorsi ridondanti **dopo** aver trasferito i loro utilizzatori. L'obiettivo non è aggiungere capacità sopra D17, ma semplificare l'architettura sostituendo progressivamente la duplicazione legacy.

#### D18.0 — Idea, obiettivo e motivazione

La configurazione deve partire dalla **sorgente fisica**, non da una definizione separata di gruppi. Chi usa il modello seleziona la reologia, il raffreddamento, la dissipazione della pressione interstiziale, le iniezioni o gli scambi di massa nelle **namelist fisiche che già ne contengono i parametri**. Per ogni sorgente fisica attiva, in quella stessa namelist, deve poter specificare direttamente tre proprietà opzionali. **Eccezione strutturale obbligatoria:** i sei processi di mass exchange sono configurati attraverso una sola nuova namelist `MASS_EXCHANGE_PARAMETERS`, con una sola tripletta numerica comune a tutti i processi attivi del batch:

1. **METHOD**: come risolvere la sorgente: EXPLICIT, NEWTON oppure ANALYTIC. La componente SEMI_IMPLICIT delle quattro interfacce fisiche è attualmente parte del trattamento NEWTON e non deve essere presentata come un quarto metodo di scheduling indipendente finché non venga definito e implementato un integratore autonomo.
2. **CADENCE**: STAGE per partecipare a ciascuno stadio RK secondo il tableau IMEX, oppure STEP per un'integrazione **una sola volta dopo il completamento degli stadi RK**, sull'intero intervallo temporale del passo.
3. **GROUP**: intero positivo che stabilisce **sia l'accoppiamento delle sorgenti sia l'ordine di esecuzione del gruppo**. A parità di cadenza, i gruppi sono eseguiti in ordine crescente di GROUP.

**Non sono richiesti PHASE, execution_order o una seconda dichiarazione esplicita dei gruppi.** PHASE è ridondante ed è già stato rimosso dal registro nella revisione precedente; execution_order va eliminato come dato indipendente perché GROUP ne determina l'ordine. Il gruppo è ricostruito automaticamente dalle proprietà delle singole sorgenti. Questo è il modello canonico di D18; eventuali input storici SOURCE_GROUPS o IMEX_* hanno soltanto un percorso di migrazione chiaramente delimitato.

**Invariante globale:** ogni processo fisico attivo compare **una sola volta in tutto lo schema temporale**. La stessa sorgente non può essere sia STAGE sia STEP, né appartenere a più gruppi. Le interfacce explicit, semi_implicit, implicit e analytic di una legge rappresentano metodi alternativi per **una stessa identità fisica**, non sorgenti diverse. Non è dunque lecito raffreddare la lava sia negli stage sia a fine passo.

Lo scopo è permettere accoppiamenti Newton quando necessari, metodi analitici locali quando disponibili, splitting espliciti e combinazioni personalizzate, senza duplicare sorgenti, senza dispatcher specifici per le singole reologie e senza configurazioni ridondanti.

#### D18.1 — Caso scientifico guida: lava con reologia NEWTON/STAGE e cooling ANALYTIC/STEP

La configurazione richiesta deve permettere:

- rheology_lava: METHOD=NEWTON, CADENCE=STAGE, GROUP=1;
- thermal_cooling: METHOD=ANALYTIC, CADENCE=STEP, GROUP=1.

Il numero GROUP=1 può essere riutilizzato nelle due cadenze: **(STAGE,1) e (STEP,1) identificano due gruppi distinti**. La reologia viene risolta implicitamente ad ogni stage IMEX; il raffreddamento viene applicato analiticamente **una sola volta** allo stato assemblato alla fine del passo, per il valore completo di dt. Il cooling non deve comparire contemporaneamente nei rate di stage.

Esempio della **sintassi già riconosciuta dal parser** nelle namelist fisiche. **Attenzione:** le combinazioni NEWTON/ANALYTIC/STEP mostrate sotto non sono ancora tutte eseguibili; la lettura respinge quelle non abilitate:

~~~fortran
&RHEOLOGY_PARAMETERS
  RHEOLOGY_MODEL = 'lava',
  RHEOLOGY_METHOD = 'NEWTON',
  RHEOLOGY_CADENCE = 'STAGE',
  RHEOLOGY_GROUP = 1,
  ! ulteriori parametri fisici della legge ...
/

&TEMPERATURE_PARAMETERS
  THERMAL_COOLING_FLAG = T,
  THERMAL_METHOD = 'ANALYTIC',
  THERMAL_CADENCE = 'STEP',
  THERMAL_GROUP = 1,
  ! ulteriori parametri termici ...
/

&PORE_PRESSURE_PARAMETERS
  PORE_PRESSURE_METHOD = 'ANALYTIC',
  PORE_PRESSURE_CADENCE = 'STAGE',
  PORE_PRESSURE_GROUP = 2,
  ! ulteriori parametri fisici ...
/
~~~

L'ultimo esempio vale soltanto quando la relativa fisica è effettivamente abilitata. I nomi dei parametri nelle altre famiglie saranno definiti dopo l'audit delle namelist. In Fortran i nomi delle variabili di un modulo non sono isolati dal nome della namelist: occorrono prefissi univoci per processo, e più processi della medesima famiglia fisica possono richiedere più triplette nello stesso blocco di input.

**Nota scientifica vincolante:** un solve analitico esatto del solo raffreddamento non rende esatto l'accoppiamento temperatura–viscosità sull'intero timestep. Lo splitting fa sì che il contributo del raffreddamento applicato a fine passo influenzi gli stage del passo successivo; accuratezza, ordine globale, dissipazione e stabilità devono essere verificati. D18 non assume equivalenza numerica fra ANALYTIC/STAGE, NEWTON/STAGE e ANALYTIC/STEP.

#### D18.2 — Definizione dei gruppi derivati e vincoli di coerenza

Sia S_active l'insieme delle sorgenti fisiche abilitate. Ciascuna sorgente s possiede esattamente una tripla efficace (M_s,C_s,G_s), ottenuta combinando i parametri facoltativi con i default compatibili. Un gruppo è l'insieme delle sorgenti con uguale tripla **(CADENCE,METHOD,GROUP)**. Questa definizione permette a gruppi di metodi differenti di condividere il valore numerico GROUP (incluso zero), restando gruppi separati. Tutte le sorgenti appartenenti a uno stesso gruppo hanno necessariamente lo stesso METHOD.

Formalmente, per due sorgenti attive i,j, esse appartengono allo stesso gruppo se e solo se \((C_i,M_i,G_i)=(C_j,M_j,G_j)\).

La condizione è necessaria ma **non sufficiente**: un gruppo deve anche essere risolvibile fisicamente e numericamente. Validare all'inizializzazione tutti i seguenti invarianti:

1. Ogni sorgente fisica attiva ha un'unica identità canonica e una sola tripla; sorgenti inattive e sotto-componenti derivate non possono essere schedulate.
2. Una stessa sorgente non può essere ripetuta in più gruppi, **nemmeno fra STAGE e STEP**. Duplicati e assenze sono errori di configurazione prima della simulazione.
3. Gruppi con identica tripla (CADENCE,METHOD,GROUP) raccolgono le sorgenti interessate; la medesima cifra GROUP è indipendente sia fra cadenze sia fra metodi differenti.
4. Tutte le sorgenti devono supportare realmente il metodo: NEWTON richiede il contributo SI/IM valido, ANALYTIC richiede la mappa analitica valida; uno stub di rifiuto non costituisce una capability.
5. Un gruppo NEWTON con più termini risolve **un solo sistema accoppiato** su tutte le equazioni target, con residuo e Jacobiana congiunti, non vari sottosolver indipendenti.
6. Un gruppo ANALYTIC con più leggi è ammesso soltanto se esiste **un operatore analitico congiunto**; la solvibilità analitica dei singoli addendi non basta.
7. La deposizione comprende la perdita della fase vettore da essa indotta. **Tutti i processi mass exchange attivi appartengono obbligatoriamente allo stesso gruppo, con identici METHOD e CADENCE, assegnati esclusivamente dalla singola tripletta `MASS_EXCHANGE_PARAMETERS`**. Il batch con limiter condiviso, inventari e topografia è indivisibile anche quando soltanto un sottoinsieme delle sei leggi è attivo; non è configurabile separare erosione, deposizione, perdite, entrainment e vaporizzazione in gruppi differenti.
8. L'ordine crescente di GROUP opera all'interno della stessa cadenza, ma non può scavalcare dipendenze obbligatorie del tableau IMEX o introdurre valutazioni esplicite in punti numericamente errati: dichiarare e validare queste restrizioni senza reintrodurre PHASE.

I gruppi devono avere identificativi non negativi (GROUP=0 valido) e un'unione coerente dei target conservativi. Il registro normalizzato è costruito **una sola volta**; nel calcolo per cella non si devono confrontare stringhe, ricostruire liste o schedulare duplicati. Le diagnostiche di errore devono indicare sorgente, namelist, METHOD richiesto e coppia (CADENCE,GROUP) coinvolta.

#### D18.3 — Default opzionali, input legacy e migrazione

I tre attributi sono **facoltativi**. Gli input già validi, senza alcun campo nuovo, devono continuare a riprodurre **la stessa fisica, gli stessi accoppiamenti Newton, le stesse cadenze, lo stesso ordine e gli stessi risultati numerici della baseline approvata**, entro le tolleranze esistenti.

Non basta assegnare indiscriminatamente GROUP=1 o un gruppo separato a ogni sorgente: i vecchi input ammettono gruppi impliciti accoppiati, gruppi analitici separati e un batch STEP indivisibile. La matrice dei default deve essere derivata dall'inventario di tutte le sorgenti canoniche e dalle configurazioni IMEX storiche, non limitata ai tre casi reologia/termica/pore pressure.

Matrice di riferimento da precisare tramite audit esaustivo:

| Sorgente o famiglia | CADENCE legacy | METHOD legacy | Regola per GROUP omesso |
|---|---|---|---|
| Reologia attiva | STAGE | NEWTON oppure ANALYTIC, secondo configurazione IMEX | Ricostruire il blocco legacy, anche accoppiato |
| thermal_cooling attivo | STAGE | NEWTON oppure ANALYTIC, secondo configurazione IMEX | Preservare gruppo e precedenze legacy |
| pore_pressure_dissipation attiva | STAGE | NEWTON oppure ANALYTIC | Preservare gruppo e precedenze legacy |
| Iniezioni e curvatura esplicite | STAGE | EXPLICIT | Preservare trattamento esplicito del tableau |
| Mass exchange, inclusi i componenti del batch | STEP | EXPLICIT | Un unico batch coordinato e coerente |

Per specifiche **parziali** (solo METHOD, solo CADENCE, solo GROUP, oppure due campi), definire una risoluzione esplicita campo-per-campo e validare il gruppo risultante. Se i valori dichiarati confliggono con un accoppiamento legacy implicito non si deve cambiare silenciosamente l'accoppiamento: generare errore oppure applicare una regola di disaccoppiamento esplicitamente documentata e testata, da decidere prima dell'implementazione.

Le namelist realmente presenti includono RHEOLOGY_PARAMETERS, TEMPERATURE_PARAMETERS e PORE_PRESSURE_PARAMETERS; iniezioni e curvatura devono essere associate alle loro namelist fisiche pertinenti. **Per gli scambi di massa si introduce invece la nuova `MASS_EXCHANGE_PARAMETERS`, che è la proprietaria unica dei parametri di attivazione e della tripletta comune dei sei processi.** Questa scelta non introduce una lista globale di sorgenti o gruppi: è una namelist fisica di un batch indivisibile. Se una namelist diversa governa più processi indipendenti, utilizzare nomi con prefisso per processo canonico.

#### D18.3a — Nuova MASS_EXCHANGE_PARAMETERS: contenuti, ownership e contratto

Creare una sola namelist `MASS_EXCHANGE_PARAMETERS` che contenga **una** tripletta opzionale `METHOD`, `CADENCE`, `GROUP` riferita all'intero batch, con default **EXPLICIT/STEP** e con GROUP scelto in modo da riprodurre esattamente la posizione attuale del mass exchange fra i gruppi STEP. La tripla è letta e validata una sola volta, poi assegnata a tutte e sole le sorgenti mass exchange effettivamente attive. Non introdurre `EROSION_METHOD`, `DEPOSITION_GROUP` o altre triplette per componente; non permettere di overrideare singole sorgenti mass exchange con `SOURCE_GROUPS`.

Il catalogo **mantiene sei identità canoniche distinte**: `mass_exchange_erosion`, `mass_exchange_deposition` (inclusa la perdita della fase vettore indotta), `mass_exchange_carrier_liquid_loss`, `mass_exchange_carrier_darcy_loss`, `mass_exchange_entrainment` e `mass_exchange_vaporization`. La nuova namelist raccoglie i **parametri di processo** prima distribuiti fra le sezioni di trasporto e pressione interstiziale; restano nei rispettivi moduli fisici le proprietà del materiale che hanno significato autonomo fuori dal mass exchange (per esempio densità, viscosità e parametri della dissipazione pore pressure). Per i parametri usati da più leggi — inclusi permeabilità, proprietà di settling ed eventuali coefficienti termodinamici — identificare un **unico owner** dopo audit delle dipendenze e documentare la loro provenienza; niente duplicazione di parametri con due valori concorrenti.

Inventario iniziale da verificare nelle routine di lettura e nelle leggi prima di spostare le dichiarazioni:

| Processo canonico | Parametri/processi da ricollocare dalla namelist attuale | Provenienza attuale | Vincolo |
|---|---|---|---|
| `mass_exchange_erosion` | `EROSION_COEFF`, parametri direttamente connessi a erodibilità/erosione | `SOLID_TRANSPORT_PARAMETERS` | Non spostare automaticamente le proprietà di composizione del materiale erodibile se usate altrove |
| `mass_exchange_deposition` | `SETTLING_FLAG`, opzioni di settling e deposizione direttamente legate allo scambio | `SOLID_TRANSPORT_PARAMETERS` | Deposizione e perdita indotta restano un solo processo canonico |
| `mass_exchange_carrier_liquid_loss` | `LOSS_RATE` | `LIQUID_TRANSPORT_PARAMETERS` | La presenza della fase liquida è un prerequisito fisico, non una seconda configurazione numerica |
| `mass_exchange_carrier_darcy_loss` | `GAS_LOSS_FLAG` e parametri specifici della perdita Darcy, dopo audit | `PORE_PRESSURE_PARAMETERS` | Non confondere la perdita di massa Darcy con la sorgente `pore_pressure_dissipation` |
| `mass_exchange_entrainment` | `ENTRAINMENT_FLAG` e relativi parametri di processo | `GAS_TRANSPORT_PARAMETERS` | Preservare i controlli di compatibilità di fase già esistenti |
| `mass_exchange_vaporization` | `LIQUID_VAPORIZATION_FLAG` e relativi parametri di processo | `GAS_TRANSPORT_PARAMETERS` | Preservare dipendenze da fase gas e trasferimenti termici |

Esempio **del contratto finale, non ancora accettato integralmente dal parser**: al momento MASS_EXCHANGE_PARAMETERS riconosce soltanto i campi prefissati `mass_exchange_method`, `mass_exchange_cadence`, `mass_exchange_group`, non i sei parametri fisici elencati di seguito né i nomi non prefissati METHOD/CADENCE/GROUP. Questo esempio è **progettuale**, non un input da eseguire:

~~~fortran
&MASS_EXCHANGE_PARAMETERS
  EROSION_COEFF = 0.01,
  SETTLING_FLAG = T,
  LOSS_RATE = 0.0,
  GAS_LOSS_FLAG = F,
  ENTRAINMENT_FLAG = F,
  LIQUID_VAPORIZATION_FLAG = F,
  METHOD = 'EXPLICIT',
  CADENCE = 'STEP',
  GROUP = 1,
/
~~~

La dichiarazione della namelist non crea sorgenti inesistenti: ogni processo rimane abilitato soltanto secondo i propri flag e prerequisiti fisici. Se **nessun processo mass exchange è attivo**, non schedulare un gruppo vuoto. Se ne è attivo anche uno solo, esso entra nel batch comune. La nuova collocazione dei parametri richiede un piano di migrazione degli input legacy: lettura compatibile temporanea o conversione esplicita, precedenza non ambigua per parametri trovati in vecchie e nuove namelist, aggiornamento di template, generatori, esempi, documentazione e test. Rifiutare valori doppi discordanti invece di selezionare silenziosamente una sezione.

**Metodi supportati:** nel default attuale il batch usa `EXPLICIT/STEP`. `EXPLICIT/STAGE` è una futura opzione da supportare **solo dopo** aver garantito assemblaggio atomico dei rate, limiter e aggiornamenti conservativi coerenti con i coefficienti RK. `NEWTON` e `ANALYTIC`, sia a STAGE sia a STEP, **non sono automaticamente capacità del mass exchange**: restano vietati dal catalogo finché non saranno implementati e verificati un residuo/Jacobiana congiunti oppure una mappa analitica congiunta valida. La configurazione uniforme non giustifica risultati nulli o separazioni implicite del batch.


La vecchia SOURCE_GROUPS, il relativo SOURCE_GROUP_ORDER e le opzioni legacy IMEX_* richiedono una migrazione pianificata: audit degli input presenti, parser di compatibilità di sola lettura ove necessario, diagnosi delle configurazioni miste ambigue, conversione di esempi/template/generatori e successiva rimozione degli adattatori obsoleti. **Non mantenere permanentemente due modelli di configurazione equivalenti né due fonti dell'ordine.** Anche i restart devono validare che schema di gruppi e metodi resti coerente.

#### D18.4 — Semantica dei metodi STAGE e STEP: requisiti imprescindibili

**EXPLICIT/STAGE** significa che i rate di **tutte** le sorgenti del medesimo gruppo sono valutati sullo **stesso stato di stage congelato** e poi sommati, prima di contribuire alla parte esplicita del sistema IMEX con i coefficienti del tableau Runge–Kutta esplicito: non è un aggiornamento di Eulero separato. È vietato aggiornare Q o QP dopo la valutazione di una sorgente e usare tale stato aggiornato per valutare la sorgente successiva dello stesso gruppo. NEWTON/STAGE e ANALYTIC/STAGE seguono invece il trattamento implicito/analitico prescritto per gli stage. Non cambiare tableau, numero di stadi, incremento IM effettivo, correzioni post-SI/post-Newton o ordine matematico.

**STEP significa una sola applicazione post-RK sul passo completo dt**: prima si assembla lo stato conservativo dopo tutti gli stage, poi si eseguono i gruppi STEP in ordine GROUP, aggiornando stato conservativo e primitivo in modo coerente. La collocazione post-RK è una conseguenza di CADENCE=STEP, non un attributo PHASE separato.

**Tutte e sei le combinazioni EXPLICIT/NEWTON/ANALYTIC × STAGE/STEP costituiscono le modalità richieste di scheduling**, compatibilmente con le capability reali di ogni legge e dei gruppi accoppiati. EXPLICIT/STAGE deve alimentare il rate esplicito del tableau; **EXPLICIT/STEP deve avanzare esplicitamente per dt dopo RK**; NEWTON/STEP e ANALYTIC/STEP richiedono veri solutori post-RK e non soltanto l'accettazione dei parametri nell'input:

- **EXPLICIT/STEP:** eseguire realmente un **avanzamento esplicito una sola volta dopo RK per l'intero dt**. Per un gruppo ordinario, dato lo stato Q_in dopo RK o dopo il precedente gruppo STEP, **valutare indipendentemente tutti i rate delle sorgenti del gruppo sul medesimo Q_in e sulle primitive QP_in coerenti, senza alcun aggiornamento intermedio dello stato**, sommarli in S_G(Q_in) e soltanto dopo aggiornare Q_out = Q_in + dt*S_G(Q_in), sulle equazioni conservative target, con conversione Q/QP coerente. Non inserire il rate STEP nei coefficienti RK e non dividerlo automaticamente fra gli stage. Quando il processo possiede vincoli condivisi, limiter, aggiornamenti di inventario o geometria, l'operazione resta un singolo incremento conservativo e atomico vincolato: il mass exchange e il refresh topografico rimangono indivisibili e nello stesso ordine relativo del default.
- **NEWTON/STEP:** implementare un reale solve implicito sull'intero dt, con residuo e Jacobiana **del gruppo completo**, gestione SI pertinente, controlli di convergenza, fallimento e ammissibilità. Definire la discretizzazione temporale di STEP (per esempio backward Euler) senza riutilizzare coefficienti RK di stage come se fossero validi per il passo completo. Lo schema scelto va documentato, verificato e incluso nel confronto scientifico.
- **ANALYTIC/STEP:** implementare la mappa di flusso analitica del gruppo su dt, **una sola volta**, con coefficienti congelati e target aggiornati secondo la legge. Per la lava riusare la routine di raffreddamento esatto già presente, attraverso il dispatcher fisico comune e non tramite chiamate specifiche nell'integratore. Estendere il percorso alla dissipazione della pressione interstiziale e ad altre leggi analiticamente risolvibili, ove compatibili.

La presenza dell'interfaccia analitica di una legge non garantisce da sola che possa essere composta analiticamente con altre leggi. Tutti i gruppi STEP devono specificare lo stato iniziale (post-RK o aggiornato dal gruppo STEP precedente), l'intervallo dt, i coefficienti congelati, il trattamento dry/wet, l'aggiornamento dei target e la propagazione degli errori. Due gruppi STEP che toccano lo stesso target sono ordinati e in generale **non commutano**; l'ordine GROUP è quindi scientificamente significativo.

D18 non introduce automaticamente subcycling, rejection/rollback né nuove politiche CFL: l'interfaccia con il futuro controllo di robustezza appartiene alla fase E. La composizione di stage e step è uno splitting da validare rispetto all'accuratezza temporale, non una riformulazione matematicamente equivalente per definizione.

#### D18.4a — Implementare EXPLICIT usando le componenti NEWTON già esistenti

**Attività nuova da realizzare in D18:** per le sorgenti che già dispongono di un trattamento NEWTON ma non hanno una implementazione EXPLICIT utilizzabile, completare la normale interfaccia \`eval_<source>_explicit\` riutilizzando le componenti fisiche esistenti SI e IM; nessuna nuova categoria di sorgenti o nuova interfaccia.

Quando la decomposizione effettiva della legge rappresenta il rate completo:

\[
S_{\mathrm{explicit}}(Q) =
S_{\mathrm{semi\_implicit}}(Q) +
\operatorname{Re}\left[S_{\mathrm{implicit}}(Q)\right].
\]

Le due componenti **devono ricevere lo stesso stato conservativo Q e le stesse primitive coerenti**, senza mutazioni o incrementi intermedi. In questa valutazione non si esegue Newton, non si applicano correzioni post-SI/post-Newton e non si modifica lo stato. L'estrazione della parte reale segue le interfacce REAL/COMPLEX esistenti e **non modifica** il percorso complex-step del Jacobiano.

Prima di abilitare una sorgente, verificare che SI+REAL(IM) sia davvero **l'intero rate fisico**, senza termini mancanti, sovrapposti o con segno errato; una procedura SI o IM presente ma incompleta non basta. In caso contrario mantenere EXPLICIT non disponibile finché il contratto non sia corretto. Le sorgenti che già implementano correttamente EXPLICIT non cambiano. Aggiornare \`has_explicit\` in \`source_catalog_2d\` **solo dopo** che l'interfaccia è operativa e testata.

La medesima interfaccia EXPLICIT è impiegata indifferentemente per \`EXPLICIT/STAGE\` (rate nel tableau esplicito RK) ed \`EXPLICIT/STEP\` (un solo incremento di ampiezza dt post-RK). Quando più sorgenti appartengono allo stesso gruppo EXPLICIT, **anche quelle costruite da SI+IM sono valutate tutte sullo stato comune immutato del gruppo**, prima di sommare i rate.

**Verifica di correttezza:** controllare che EXPLICIT rappresenti SI+REAL(IM) sullo stesso stato e che Newton, ANALYTIC e Jacobiane REAL/COMPLEX restino invariati. **Decisione operativa (26 ottobre 2026):** l'utente ha chiesto di non introdurre test numerici dedicati alle reologie esplicite; per ora mantenere i gate di compilazione e regressione esistenti e concentrarsi sulla semplificazione architetturale. Tale scelta non equivale a certificare una nuova accuratezza scientifica, che resta un gate separato.

#### D18.4b — Fotografia storica dell'implementazione (26 ottobre 2026; SUPERATA, non stato corrente)

**Evidenza:** codice su `master` dopo `558fd94`; l'utente ha confermato entrambe le Actions verdi per tale commit. Le Actions non sostituiscono una dimostrazione delle modalità nuove non esercitate dai fixture. Stati: **IMPLEMENTATO** = presente nel codice; **PARZIALE** = solo sottoinsieme o transizione; **DA FARE** = requisito non implementato; **VALIDATO CI** = compilazione/regressioni storiche superate sui percorsi coperti.

| Ambito / riferimento | Stato audit | Evidenza presente e lacuna concreta |
|---|---|---|
| D18.0 modello sorgente-centrico | **PARZIALE** | Triplette opzionali per reologia, cooling, pore pressure e batch mass exchange; rimangono `SOURCE_GROUPS`, opzioni `IMEX_*`, `execution_order` e due percorsi di inizializzazione. Le namelist di iniezioni e curvatura non espongono ancora triplette proprie. |
| D18.1 esempio lava NEWTON/STAGE + cooling ANALYTIC/STEP | **DA FARE** | Parser e scheduler respingono il mixed-mode e non eseguono cooling analitico post-RK. Non confondere l'esempio normativo con un input funzionante. |
| D18.2 gruppi derivati | **PARZIALE** | `initialize_explicit_physical_sources` raggruppa il solo caso interamente EXPLICIT/STAGE per ID. Nessuna costruzione unificata dei gruppi misti; `execution_order` e `SOURCE_GROUP_ORDER` restano attivi per compatibilità. |
| D18.3 default e parsing delle triplette | **PARZIALE** | Input legacy e `SOURCE_GROUPS` restano utilizzabili; campi fisici letti in `inpout_2d` e conservati in `source_input_2d`. Solo EXPLICIT/STAGE è accettato per le tre sorgenti fisiche; selezioni miste respinte. Manca risoluzione campo-per-campo degli input parziali e migrazione di esempi/generatori. |
| D18.3a mass exchange unico | **PARZIALE** | Namelist opzionale `MASS_EXCHANGE_PARAMETERS` legge una tripletta prefissata; EXPLICIT/STEP con GROUP positivo modifica l'ID del batch legacy. I parametri fisici di erosione/deposizione/perdite ecc. **non** sono ancora trasferiti; non c'è ancora un unico owner per tali parametri. Metodi e cadenze alternativi non supportati. |
| D18.4a interfacce EXPLICIT locali | **IMPLEMENTATO; VALIDATO CI per build/regressioni legacy** | Tutte le 11 reologie implementano `eval_rheology_<model>_explicit` nel proprio file riutilizzando SI+REAL(IM); `rheology_2d` fa dispatch; thermal e pore pressure hanno interfacce EXPLICIT nei moduli proprietari; `source_terms_2d` le assembla. Il valutatore generale `eval_newton_source_explicit` è stato **rimosso**. Mancano verifiche funzionali sistematiche delle nuove combinazioni; su richiesta dell'utente **non aggiungere test numerici dedicati alle reologie esplicite** in questa fase. |
| D18.4 STAGE EXPLICIT | **PARZIALE; build/regressioni legacy OK** | `initialize_explicit_physical_sources` crea gruppi EXPLICIT/STAGE e una partizione senza target IM. Non ancora verificata con un gate funzionale nuovo; passaggio del valore stocastico locale EXPLICIT incompleto (attualmente zero). |
| D18.4 gruppi misti NEWTON/ANALYTIC/EXPLICIT STAGE | **DA FARE** | Serve una sola normalizzazione che costruisca **insieme** registry e partizione implicita, preservando solve congiunti e tableaux. |
| D18.4 esecutore STEP generale | **DA FARE** | `step_source_executor_2d` è ancora centrato sul batch EXPLICIT/STEP; NEWTON/STEP, ANALYTIC/STEP e gruppi STEP ordinari multimodali non sono eseguiti. |
| D18.5 migrazione e rimozione ridondanze | **DA FARE** | Eliminare progressivamente `SOURCE_GROUP_ORDER`, `execution_order` e i doppioni legacy solo dopo parità degli input e verifica dei loro utilizzatori. |
| D18.6 accettazione scientifica finale | **APERTO** | Actions storiche verdi non dimostrano splitting STEP, metodi misti, restart con schema nuovo né equivalenza numerica delle configurazioni nuove. |

**Audit dei commit D18 già confluiti:** `d1929d4` (interfacce iniziali); `93ac675`, `8bc7297` (riparazioni regressioni storiche); `b8e6c13`, `09c49c9` (**tentativo superato e poi rimosso** del valutatore SI+IM generale e relativo test); `ee97ee2`, `66f0f3f`, `07a820e` (ownership EXPLICIT negli undici moduli reologici e rimozione valutatore generale); `0b2c9d3` (assemblatore EXPLICIT); `10fbcfc`, `94bc6a8` (input fisici e mass exchange); `0750e4d`, `659d41f` (ID del batch e reset input); `558fd94` (primo percorso completamente EXPLICIT/STAGE). **Non trattare i commit storici superati come funzioni tuttora presenti.**

**Priorità per il prossimo sviluppo, minimizzando il codice:** (1) sostituire la biforcazione legacy/tutto-EXPLICIT con **una** normalizzazione che calcola triplette efficaci per ogni sorgente attiva e costruisce contestualmente `equation_partition_type` e `source_registry_type`; (2) riutilizzare gli esecutori STAGE esistenti per gruppi misti, senza nuove strutture di scheduling, con controllo della coerenza METHOD per `(CADENCE,GROUP)`; (3) estendere l'esecutore STEP esistente ai metodi realmente disponibili, senza decomporre mass exchange né duplicare il motore Newton; (4) migrare parametri mass exchange e le altre namelist fisiche; (5) eliminare gradualmente i campi ridondanti dopo il passaggio. **Non introdurre un terzo percorso di inizializzazione.**

**Semplificazione di D18 — 26 ottobre 2026 (CI da verificare):** eliminata l'interfaccia pubblica `initialize_explicit_physical_sources` e rimossa la biforcazione corrispondente in `simulation_2d`. Gli input fisici non nativi passano da `initialize_configured_sources`, che applica internamente la scelta esplicita quando richiesta, utilizzando le strutture già presenti. Non sono stati aggiunti scheduler, registry o moduli. Resta aperta la normalizzazione unificata per configurazioni miste EXPLICIT/NEWTON/ANALYTIC e la rimozione delle duplicazioni interne legacy. Stato: PARZIALE, non D18 completato.

**Gruppi a metodi differenti — intervento parziale (26 ottobre 2026, CI da verificare):** la procedura privata `apply_explicit_physical_selection` nel modulo esistente ora preserva le sorgenti non promosse nei gruppi impliciti legacy, rimuove le sole sorgenti optate EXPLICIT/STAGE e ricostruisce insieme partizione IMEX e registry. Conserva i gruppi NEWTON residui accoppiati, senza introdurre strutture o sottosolver. I gruppi espliciti con stesso `(STAGE,GROUP)` vengono aggregati; collisioni METHOD sono respinte. **Non è ancora la normalizzazione D18 completa**: NEWTON/ANALYTIC espliciti nelle triplette non riscrivono il metodo legacy; nuovi GROUP su gruppi impliciti e metodi STEP restano aperti. Verifica CI necessaria.

#### D18.4c — Audit dei contributi di curvatura e identificazione dei gruppi di default (26 ottobre 2026)

**Riscontro sul codice:** `geometry_curvature` non è un alias della correzione reologica. In `source_terms_2d::eval_geometry_curvature_explicit`, la derivata seconda del fondo forma `a_c = B_xx u² + 2 B_xy uv + B_yy v²`; la sorgente autonoma di quantità di moto vale `S_x = -G_c (rho h) a_c B_x` e `S_y = -G_c (rho h) a_c B_y`. La chiamata avviene una volta nel gruppo EXPLICIT/STAGE quando `curvature_term_flag` è attivo. Nelle leggi reologiche (ad esempio Voellmy e mu(I)) la **stessa accelerazione geometrica** entra invece nel calcolo dello sforzo normale efficace e quindi dell'attrito: è una dipendenza costitutiva della sorgente reologica, non la stessa forza tangenziale addizionata una seconda volta. Le due espressioni sono diverse e hanno bersagli fisici distinti; **non rimuovere automaticamente** `geometry_curvature` per il solo fatto che una reologia usa `a_c`. L'audit statico non costituisce da solo una dimostrazione completa dell'assenza di doppi conteggi rispetto a ogni variante delle equazioni HP: conservarne la verifica scientifica nel gate D15.

**Regola di configurazione approvata per D18, non ancora implementata nel registry:** identificare ogni gruppo con `(CADENCE, METHOD, GROUP)`, non con `(CADENCE, GROUP)`. `GROUP=0` deve essere un identificativo **valido** e il default universale per un campo GROUP omesso. Questo consente contemporaneamente `(STAGE,NEWTON,0)`, `(STAGE,EXPLICIT,0)` e `(STEP,EXPLICIT,0)`; le sorgenti di una stessa tripla sono raccolte nello stesso gruppo. All'interno di una stessa classe `(CADENCE,METHOD)`, l'ordine è crescente per GROUP. In caso di METHOD diversi con lo stesso numero GROUP non c'è conflitto: sono **gruppi distinti**, non gruppi eterogenei. L'ordine fra metodi diversi durante STAGE rimane imposto dal tableau IMEX, non dal numero GROUP.

**Default fisici da rendere espliciti nel catalogo esistente `source_catalog_2d`:** per le reologie, thermal cooling e pore-pressure dissipation `NEWTON/STAGE/0` quando attive, formando **un unico solve Newton accoppiato** in assenza di override; per iniezioni e `geometry_curvature`, `EXPLICIT/STAGE/0`; per tutto il batch mass exchange, `EXPLICIT/STEP/0`. Le interfacce supportate rimangono le attuali capability della sorgente; `semi_implicit` è una componente di NEWTON e non un METHOD autonomo. Non aggiungere un secondo catalogo, né triplette indipendenti per ciascuna delle sei leggi mass exchange.

**Compatibilità da verificare prima del cambio dei default:** il percorso IMEX legacy ammette opzioni esplicite `IMEX_*`, e il registry corrente considera validi solo identificativi GROUP positivi. Il comportamento richiesto di default è quello dichiarato per le simulazioni legacy con unico Newton accoppiato; preservare comunque i vecchi input che impostano esplicitamente gruppi/metodi diversi. Normalizzare tali input *una volta* all'avvio, poi eseguire lo stesso scheduler delle namelist fisiche. Non conservare due interpretazioni operative permanenti e non introdurre cambi di METHOD, CADENCE o GROUP durante la simulazione.

**Esito audit curvatura:** `geometry_curvature` è una sorgente geometrica EXPLICIT/STAGE distinta dall'impiego di curvatura nel calcolo costitutivo reologico. Il catalogo deve riconoscerla separatamente, ma la formula ausiliaria `curvature_acceleration` va riutilizzata, non duplicata. Nessuna modifica alle formule o alle interfacce fisiche è giustificata da questo solo audit.

#### D18.4d — Audit strutturale di semplificazione (26 ottobre 2026)

**Scopo e natura:** audit del codice al commit `3861548` (sulla base dei file letti da `master`), senza modifiche alle formule o agli operatori numerici. Questo paragrafo ha **valore progettuale per D18**, non prescrive di mantenere i moduli D17. Non considerare un modulo necessario solo perché è già presente. Il criterio di scelta è il **numero minimo di dati indipendenti, percorsi di inizializzazione e trasformazioni** compatibile con le capacità fisiche e numeriche richieste. Il commit `d7da3e0` è stato precedentemente confermato solo a livello di pubblicazione, non di Actions; non attribuirgli un gate CI che non risulta comunicato.

**Inventario reale (linee indicative al momento dell'audit):**

| Modulo | Dimensione | Responsabilità effettiva | Decisione D18 |
|---|---:|---|---|
| `source_catalog_2d` | 123 righe | Identità canoniche, target e tre capability; lookup e associazione nomi reologici | **MANTENERE LA FUNZIONE, NON NECESSARIAMENTE IL MODULO:** un'unica definizione delle capability e dei default è necessaria per validare l'input. Aggiungere METHOD/CADENCE default, GROUP=0 comune; evitare un secondo elenco di sorgenti. Accorpare soltanto se riduce davvero dipendenze e righe. |
| `source_registry_2d` | 525 righe | `source_group_type`, copie dei gruppi, validazione e grande traduzione legacy `initialize_legacy_registry` | **RIDURRE RADICALMENTE:** mantenere soltanto i gruppi normalizzati indispensabili al runtime e una validazione, spostando/rimuovendo la traduzione legacy dopo la migrazione. Valutare unione con la configurazione **solo dopo** la rimozione dei percorsi doppi. |
| `source_configuration_2d` | 574 righe | Accensione sorgenti, produzione gruppi legacy, adattamento dei metodi EXPLICIT, parsing/traduzione `SOURCE_GROUPS`, ricostruzione della partizione | **UN SOLO ALGORITMO:** risolvere le triplette efficaci per sorgente attiva e costruire direttamente gruppi e partizione; eliminare il passaggio `initialize_legacy` → `apply_explicit_physical_selection` → nuova validazione. La lettura legacy deve produrre soltanto gli stessi input normalizzati. |
| `source_input_2d` | 187 righe | Storage fisico, array `SOURCE_GROUPS`, scanner e reader di `MASS_EXCHANGE_PARAMETERS` | **RIDURRE:** nessun doppio storage permanente per configurazioni equivalenti. Concentrare gli input e convertire quelli storici all'avvio; non mantenere `SOURCE_GROUP_ORDER` come proprietà canonica. |
| `equation_metadata_2d` | 323 righe | `equation_partition_type`, mappe e blocchi impliciti usati dal solver | **MANTENERE LA FUNZIONE:** la partizione e la mappa dei target sono dati numerici necessari; non ricrearle nel registry o nei gruppi. Evitare memorizzazioni ridondanti derivabili dalla partizione. |
| `time_integration_2d` | 674 righe | Tableau IMEX e valutazione dei gruppi STAGE espliciti, accesso a partizione e gruppi | **MANTENERE IL MOTORE:** non introdurre un altro scheduler RK; consumare gruppi già validati. Ottimizzare in seguito i confronti sui nomi e le copie di gruppi, non le formule. |
| `step_source_executor_2d` | 222 righe | Transazione mass exchange, limiter condiviso, inventari e topografia dopo RK | **MANTENERE IL KERNEL ATOMICO:** il loop attuale esegue il batch STEP; è scorretto sostituirlo con incrementi indipendenti per ogni sorgente. Estendere in questo percorso gli altri operatori STEP, senza duplicare il batch. |
| `simulation_2d` | 193 righe | Inizializzazione contesto e scelta tra `SOURCE_GROUPS` e percorso fisico | **SEMPLIFICARE:** una chiamata alla configurazione normalizzata indipendente dalla sintassi di input; nessuna logica dei metodi nel contesto. |

**Ridondanze dimostrate dall'ispezione:** (1) attivazione delle stesse sorgenti elencata sia in `initialize_configured_sources` sia in `initialize_native_source_configuration`; (2) traduzione `initialize_legacy_registry` seguita, per scelta EXPLICIT, da una seconda costruzione `apply_explicit_physical_selection`; (3) ordine espresso da `source_group_order`, `execution_order` e ID GROUP, nonostante D18 richieda solo quest'ultimo; (4) `legacy_id` persistente nei descrittori di esecuzione; (5) `source_groups` e `source_input` mantengono in parallelo impostazioni legacy e fisiche; (6) `stage_groups()` e `step_groups()` producono copie allocabili in più call site; (7) `time_integration_2d` ricostruisce un sottoinsieme di gruppi espliciti, e i rate per cella interrogano stringhe di sorgenti; (8) `step_source_executor_2d` conta i gruppi ma applica una sola transazione globale, il che conferma che **la semplice esistenza del registry non costituisce un esecutore STEP multimodale**.

**Informazioni minime necessarie dopo l'inizializzazione:**
- Per ogni **sorgente attiva**, identità canonica (o un riferimento stabile a essa), metodi realmente supportati e DEFAULT_METHOD/DEFAULT_CADENCE; non salvare sei copie di proprietà identiche del batch mass exchange.
- Per ciascun **gruppo effettivo**, una sola chiave `(CADENCE,METHOD,GROUP)` con GROUP non negativo, membri fisici univoci e target strettamente necessari alla valutazione. I target impliciti e le mappe di accoppiamento sono posseduti da `equation_partition_type`, non duplicati in rappresentazioni concorrenti. L'ordine è derivato da METHOD/tableau e GROUP; nessun `execution_order` indipendente.
- Per il runtime, una rappresentazione costruita **una sola volta** e immutabile: nessuna scelta dinamica di metodo, nessuna ricostruzione di gruppi e nessuna riconciliazione legacy nel ciclo di tempo. Il batch mass exchange è sempre una sola transazione STEP (salvo nuove capability validate).
- Non serve necessariamente un oggetto chiamato `source_registry_type`: **serve** la rappresentazione dei gruppi e la sua validazione. Decidere se eliminare/accorpare il modulo soltanto dopo aver enumerato le dipendenze di solver, test e restart.

**Architettura obiettivo minima (logica, non nuovi moduli):**
`namelist fisiche / input storici → risoluzione unica dei default e delle sorgenti attive → validazione capability e gruppi (CADENCE,METHOD,GROUP) → costruzione simultanea di gruppi e partizione IMEX → esecuzione con operatori STAGE/STEP esistenti`.
Non produrre prima gruppi legacy e poi «spostare» sorgenti; ogni sorgente riceve la propria configurazione **fissa all'inizializzazione**. Identità di default `GROUP=0` per ciascuna coppia METHOD/CADENCE. Reologia, cooling e pore pressure di default formano `NEWTON/STAGE/0` **come unico Newton accoppiato**; iniezioni/curvatura `EXPLICIT/STAGE/0`; mass exchange `EXPLICIT/STEP/0`. Compatibilità con input IMEX espliciti va risolta *prima* di costruire i gruppi, non con un secondo scheduler.

**Decisione sui moduli:** mantenere **per ora** il piccolo catalogo delle capability, la partizione delle equazioni e i kernel di esecuzione: hanno compiti distinti. Non decidere anticipatamente di conservare tutto `source_registry_2d` o tutto `source_configuration_2d`; la loro attuale separazione è riesaminabile. L'eliminazione di un modulo è un risultato possibile, non un requisito estetico. La misura di successo è meno stato mutabile, meno traduzioni e minore codice complessivo, non soltanto meno file.

**Roadmap di riduzione, in modifiche coordinate:**
1. **Contratto e default:** estendere il metadato fisico *esistente* con default METHOD/CADENCE, usare GROUP=0 senza valore sentinella e aggiornare la validazione della tripla; separare l'assenza del campo dai valori espliciti durante il solo parsing. Risolvere esattamente l'input legacy a singolo Newton e quelli con IMEX_* espliciti.
2. **Unificazione in un'unica costruzione:** ottenere tutte le sorgenti attive con una sola routine, risolvere le triplette, creare direttamente gli stessi gruppi e la partizione, rimuovendo `apply_explicit_physical_selection` e il rebuild da `initialize_legacy_registry` nel percorso di produzione. Conservare temporaneamente una funzione di lettura legacy che **traduce i valori**, non i gruppi.
3. **Sfoltimento delle strutture:** eliminare `legacy_id`, `execution_order`, `SOURCE_GROUP_ORDER` dal contratto normalizzato dopo migrazione dei loro utilizzatori; rimuovere copie, ordinamenti e validazioni duplicate. Valutare se accorpare/eliminare `source_registry_2d` una volta che la traduzione legacy non vi risiede più.
4. **STEP multimodale senza duplicazioni:** estendere il passo post-RK nel medesimo esecutore; riusare i solver fisici, il Newton congiunto e l'operatore analitico esistenti; conservare atomico il mass exchange. Non dichiarare supportati NEWTON/STEP o ANALYTIC/STEP finché non eseguiti realmente.
5. **Chiusura:** migrare il contenuto fisico di `MASS_EXCHANGE_PARAMETERS`, gli altri input e la documentazione; aggiornare startup, output e restart; eseguire regressioni legacy e gate scientifici. Eliminare gli adattatori residui dopo la parità verificata.

**Criteri obbligatori per ogni commit di semplificazione:** descrivere *che cosa viene eliminato* (campi, routine, duplicazioni e call site), mantenere costanti formule, tableaux, batch topografico e accoppiamento Newton, non cambiare oracoli per far passare CI, aggiornare gli stati della roadmap, verificare build e regressioni. **Non aggiungere test numerici dedicati alle reologie EXPLICIT**, come richiesto; questo non autorizza a dichiarare scientificamente verificati metodi che non sono ancora stati esercitati.

**Avanzamento del contratto a triplette (26 ottobre 2026; CI non ancora eseguita):** in `source_catalog_2d` ogni descrittore fisico dispone ora anche di `default_method` e `default_cadence`; NEWTON/STAGE per le sorgenti implicite, EXPLICIT/STAGE per gli apporti geometrici, EXPLICIT/STEP per le leggi mass exchange. La validazione `source_registry_2d` accetta GROUP=0 e identifica i duplicati con la tripla `(CADENCE,METHOD,GROUP)`. La ricostruzione EXPLICIT/STAGE usa per default il gruppo 0 e non confonde gruppi con metodi differenti. **Transizione incompleta:** i percorsi legacy e SOURCE_GROUPS usano ancora `legacy_id`, `execution_order`, TARGETS e ID positivi; il catalogo nuovo non è ancora l'unico punto da cui la normalizzazione determina i default. Non dichiarare conclusa la riduzione del registry né la migrazione completa prima della rimozione di tali residui. Il presente commit non modifica formule o metodi di integrazione. Compilazione e regressioni richieste prima di ulteriori eliminazioni.

**Riduzione effettiva di codice (26 ottobre 2026; compilazione/regressioni da verificare):** l'enumerazione delle sorgenti fisiche attive, duplicata nei percorsi `initialize_configured_sources` e `initialize_native_source_configuration`, è stata accentrata in `collect_active_source_names`. La selezione termica resta parametrizzata per non cambiare l'attivazione legacy legata a `partition%implicit(4)` rispetto all'attivazione fisica usata dall'input `SOURCE_GROUPS`. Questo **non** elimina ancora il percorso intermedio `initialize_legacy_registry` → `apply_explicit_physical_selection`, né i campi `legacy_id` ed `execution_order`; il gate di semplificazione D18 resta aperto. Non dichiarare completata la normalizzazione diretta delle triplette.

**Refactoring della costruzione fisica diretta (26 ottobre 2026; CI da verificare):** il percorso `initialize_configured_sources` ora enumera le sorgenti attive, risolve METHOD/CADENCE/GROUP da metadati e input fisici, aggrega direttamente i gruppi finali e costruisce la partizione IMEX. Rimossa integralmente la procedura intermedia `apply_explicit_physical_selection`; il percorso di produzione non chiama più `registry%initialize_legacy`. Il metodo e i target della partizione legacy sono ancora usati **solo come default di compatibilità** in questa transizione. `initialize_legacy_registry` rimane temporaneamente nel registry per i test D17; `SOURCE_GROUPS` ha ancora il proprio lettore. Sono ancora residui `legacy_id` ed `execution_order` nel tipo e restano da eliminare. I metodi STEP non realmente eseguibili sono rifiutati. Le regression e la corretta conservazione dei target-only legacy devono essere accertate prima di considerare questa migrazione conclusa.

**Semplificazione runtime dei descrittori (26 ottobre 2026; Actions da verificare):** rimossi `legacy_id` ed `execution_order` da `source_group_type`, la loro inizializzazione, le verifiche ridondanti sul registry e le assegnazioni nei test esistenti. L'ordine di esecuzione è ormai l'ordine dell'array dei gruppi, prodotto una sola volta dall'inizializzazione; il metodo e il numero GROUP restano nella tripla di identità. Il lettore `SOURCE_GROUPS` conserva temporaneamente `source_group_order` **solo come compatibilità dell'input storico** per produrre tale ordine, non come campo runtime. Ancora da eliminare dopo la migrazione dell'input storico; i casi legacy con ordine di dipendenza non devono cambiare. Questo commit non altera i kernel fisici né la sequenza risolta dal solver.

**Rimozione SOURCE_GROUP_ORDER (26 ottobre 2026, verifica CI richiesta):** eliminati `source_group_order` dall'input `SOURCE_GROUPS`, dalla lettura namelist e dall'ordinamento nel percorso nativo. Le configurazioni native sono ora ordinate per classe di metodo e `GROUP`; le identità sono convalidate come triple `(CADENCE,METHOD,GROUP)`, quindi `GROUP=0` è ammesso e ID uguali con metodi diversi non collidono. Aggiornati i test che costruiscono input nativi. La dipendenza legacy dal vecchio ordine IMEX viene mantenuta nel percorso di compatibilità distinto. **Incompatibilità intenzionale**: input nativi che contengono ancora `SOURCE_GROUP_ORDER` devono rimuovere quel campo; un ordine numerico passato non viene silenziosamente rispettato. L'eliminazione completa dei percorsi doppi e dei default legacy resta aperta.

**Aggiornamento implementativo (27 ottobre 2026; PR #3, NON ancora validata):** dopo la rimozione di `SOURCE_GROUPS`, il branch `d18-remove-legacy-imex-config` elimina la chiamata di produzione a `init_problem_param`: `initialize_configured_sources` costruisce i gruppi e la partizione IMEX direttamente dalle triplette fisiche e dai default del catalogo. La namelist `NUMERIC_PARAMETERS` non accetta piu i quattro campi `IMEX_VELOCITY_GROUP`, `IMEX_ENERGY_GROUP`, `IMEX_PORE_PRESSURE_GROUP`, `IMEX_GROUP_METHOD`. Sono iniziati gli adeguamenti dei fixture `TEST_IMEX_STAGES` e `TEST_SOURCE_REGISTRY`. **Non dichiarare completata D18:** `init_problem_param` e gli array `IMEX_*` sono ancora presenti in altri moduli e usati nei test di accoppiamento; vanno eliminati solo insieme alla migrazione dei casi numerici esistenti. I nuovi default possono raggruppare le sorgenti in modo diverso dai vecchi input; la convalida dei gruppi ANALYTIC e dei gruppi NEWTON accoppiati deve essere esplicita, senza alterare le formule. Non e stata eseguita la compilazione locale o la CI di questa PR. Il branch e una modifica in lavorazione, non pronta al merge.

**Pulizia addizionale del branch (27 ottobre 2026; non validata):** rimossi i quattro globali `IMEX_*` da `parameters_2d`, eliminata `init_problem_param` da `equation_terms_2d`, aggiornati i test D17 di Newton, `TEST_RHEOLOGY_CELL_LOCAL`, `TEST_MASS_EXCHANGE`, `TEST_SOURCE_REGISTRY` e `TEST_IMEX_STAGES` in direzione del solo percorso fisico. Nei test reologici eliminato il rebuild `refresh_source_registry` basato su `initialize_legacy`. La compatibilità con i vecchi input `IMEX_*` e `SOURCE_GROUPS` non e un requisito: gli esempi e fixture che li contengono vanno migrati. **Gate aperto:** i nuovi default/casi analitici con gruppi condivisi, la compilazione di tutti i target e le regressioni devono ancora essere verificati. Questa nota documenta le modifiche effettive del branch, non dichiara la chiusura della PR o di D18.

**Audit di rimozione adattatori (27 ottobre 2026; branch PR #3, CI da rieseguire):** eliminata l'interfaccia `source_registry_type%initialize_legacy` e la routine `initialize_legacy_registry` da `src/numerics/source_registry_2d.f90`. Eliminati anche gli helper esclusivamente legacy `target_equation_ids_first`, `names_for_targets` e `make_group`; mantenuti `target_equation_ids` e `target_index` per la validazione corrente dei gruppi. I tre utilizzatori di `initialize_legacy` in `TEST_SOURCE_REGISTRY`, `TEST_IMEX_STAGES` e `TEST_RHEOLOGY_CELL_LOCAL/probe_sources.f90` sono stati convertiti alla generazione tramite `initialize_configured_sources` e triplette fisiche. Rimosso il traduttore `resolve_imex_group_method` e i suoi test. Rimangono da verificare con CI i nuovi fixture e da migrare eventuali input, descrizioni o generatori obsoleti; nessuna compatibilita con input storici viene richiesta. I metadati interni della partizione IMEX rimangono necessari e non sono un'interfaccia di configurazione alternativa.

#### D18.4e — Audit delle 24 identita canoniche e delle namelist (28 ottobre 2026; HEAD verificato prima delle modifiche)

Questo audit confronta il codice `source_catalog_2d`, `source_configuration_2d`, `source_input_2d` e le dichiarazioni NAMELIST in `inpout_2d` su `master`. **CI verde dei test esistenti non equivale a completezza del parser.** Tutti i 24 descrittori del catalogo hanno `default_method` e `default_cadence`; la normalizzazione assegna `GROUP=0` quando non sovrascritto. Le capability dichiarate richiedono ancora verifica fisica caso per caso.

| Identita fisica (catalogo) | Default METHOD/CADENCE/GROUP | Input di scheduling esposto | Esito |
|---|---|---|---|
| `rheology_voellmy`, `rheology_plastic`, `rheology_lava`, `rheology_lahar`, `rheology_empirical`, `rheology_quadratic`, `rheology_froude`, `rheology_lucas`, `rheology_mui`, `rheology_mui_voellmy` (10 alternative) | NEWTON/STAGE/0 | `RHEOLOGY_PARAMETERS`: `rheology_method/cadence/group` comune alla reologia selezionata | Disponibile; singola reologia attiva per volta |
| `rheology_newmodel` (fixture estensibile) | NEWTON/STAGE/0 | La stessa tripletta reologica | Solo NEWTON dichiarato; non trattare come legge con explicit/analytic disponibili |
| `thermal_cooling` | NEWTON/STAGE/0 | `TEMPERATURE_PARAMETERS`: `thermal_method/cadence/group` | Disponibile |
| `pore_pressure_dissipation` | NEWTON/STAGE/0 | `PORE_PRESSURE_PARAMETERS`: `pore_pressure_method/cadence/group` | Disponibile |
| `bottom_circular_injection` | EXPLICIT/STAGE/0 | Nessuna tripletta in `RADIAL_SOURCE_PARAMETERS` | **Manca configurabilita per sorgente** |
| `bottom_fissural_injection` | EXPLICIT/STAGE/0 | Nessuna tripletta in `FISSURAL_SOURCE_PARAMETERS` | **Manca configurabilita per sorgente** |
| `radial_injection` | EXPLICIT/STAGE/0 | Nessuna tripletta in `RADIAL_SOURCE_PARAMETERS` | **Manca configurabilita per sorgente** |
| `geometry_curvature` | EXPLICIT/STAGE/0 | Nessuna tripletta in `EXPL_TERMS_PARAMETERS` o altra namelist fisica | **Manca configurabilita per sorgente** |
| `boundary_injection` | EXPLICIT/STAGE/0 | Nessuna tripletta in `LATERAL_SOURCE_PARAMETERS` | **Manca configurabilita per sorgente** |
| `mass_exchange_erosion`, `mass_exchange_deposition`, `mass_exchange_carrier_liquid_loss`, `mass_exchange_carrier_darcy_loss`, `mass_exchange_entrainment`, `mass_exchange_vaporization` (6) | EXPLICIT/STEP/0 | Una tripletta `mass_exchange_method/cadence/group` in `MASS_EXCHANGE_PARAMETERS` | Scheduling unico disponibile; **migrazione parametri fisici incompleta** |

**Risultato numerico dell'inventario:** 11 reologie alternative + 2 sorgenti termo/pore + 5 sorgenti geometriche/di iniezione + 6 leggi mass exchange = **24**. La selezione della tripletta per target in `source_configuration_2d` tratta esplicitamente `velocity`, `temperature`, `pore_pressure` e `mass_exchange`; tutti i target `explicit` ricadono nel `CASE DEFAULT` e non leggono un override METHOD/CADENCE/GROUP da input. I default nel catalogo non risolvono questa lacuna.

**Ownership mass exchange ancora da migrare:** `MASS_EXCHANGE_PARAMETERS` in `source_input_2d` dichiara solo `mass_exchange_method`, `mass_exchange_cadence`, `mass_exchange_group`. In `inpout_2d`, `erosion_coeff` e `settling_flag` sono ancora in `SOLID_TRANSPORT_PARAMETERS`, `loss_rate` in `LIQUID_TRANSPORT_PARAMETERS`, `entrainment_flag` e `liquid_vaporization_flag` in `GAS_TRANSPORT_PARAMETERS`, `gas_loss_flag` in `PORE_PRESSURE_PARAMETERS`. Separare in un audit di dipendenza i veri parametri di processo dalle proprieta condivise prima di spostare le dichiarazioni; aggiornare parser, test, esempi, generatori e documentazione in una modifica coerente. L'attuale parser NON accetta ancora la sintassi progettuale non prefissata `METHOD/CADENCE/GROUP` nel blocco mass exchange.

**Contratto di gruppo attuale approvato:** ANALYTIC ha sempre un gruppo esclusivo a STAGE e STEP; STAGE permette la stessa cifra fra EXPLICIT e NEWTON, applicati secondo IMEX; STEP vieta metodi diversi nello stesso GROUP e ordina gli operatori per GROUP crescente. Questo criterio prevale sulle tabelle storiche che affermano il contrario.

**Piano correttivo prioritario (prima di ulteriori test di splitting):** (1) introdurre triplette opzionali con prefissi distinti per le cinque sorgenti di geometria/iniezione, risolverle nell'unico `initialize_configured_sources` e coprire parser e capability; (2) migrare la proprieta dei parametri del batch in `MASS_EXCHANGE_PARAMETERS` senza duplicazioni, conservando atomicita e limiter; (3) aggiornare tutti gli input e i generatori distribuiti; (4) verificare la compatibilita dei default con ordine crescente dei gruppi STEP e proseguire con i gate scientifici realmente residui. **Nessuna formula, tableau RK o trattamento fisico deve cambiare per questi lavori.**

**Aggiornamento implementativo 28 ottobre 2026 (commit `70b6b6f`, `38e2723`; CI da verificare):** le cinque sorgenti `bottom_circular_injection`, `bottom_fissural_injection`, `radial_injection`, `geometry_curvature` e `boundary_injection` dispongono ora di triplette opzionali prefissate nelle rispettive namelist `RADIAL_SOURCE_PARAMETERS`, `FISSURAL_SOURCE_PARAMETERS`, `RADIAL_SOURCE_PARAMETERS`, `EXPL_TERMS_PARAMETERS`, `LATERAL_SOURCE_PARAMETERS`. La selezione passa attraverso il ramo `explicit` di `initialize_configured_sources`, senza un secondo scheduler, e il test `TEST_SOURCE_REGISTRY` verifica cinque ID individuali e il rifiuto di metodi non disponibili. Le selezioni `EXPLICIT/STAGE` restano il default. **Il batch MASS_EXCHANGE_PARAMETERS non e stato ancora migrato**: per evitare un cambiamento silenzioso degli input esistenti, il trasferimento dei parametri di processo dalle vecchie namelist deve essere accompagnato da un aggiornamento coordinato degli esempi, delle letture e delle regole di duplicazione. Prima di considerare operativa una selezione alternativa `EXPLICIT/STEP` per una sorgente geometrica, occorre verificarne il contratto con il percorso di iniezione e il dominio; non basta che la capability EXPLICIT sia dichiarata nel catalogo.

**Verifica CI del 28 ottobre 2026:** sul commit `1bb373ad51d2313ed146527f552313e6f08fc959`, `D17 historical physics regressions` risulta verde dopo la correzione del generatore condiviso `TEST_PCCU_LAKE_AT_REST/generate_case.py`: `ENTRAINMENT_FLAG` e `LOSS_RATE` sono ora emessi solo in `MASS_EXCHANGE_PARAMETERS`. Il fallimento precedente di `TEST_IMEX_STAGES` avveniva durante il parsing di `GAS_TRANSPORT_PARAMETERS`, prima della verifica di `N_RK`, e non dimostrava un errore del tableau IMEX. Questo risultato chiude la regressione del generatore, **non** i gate integrati D18 su splitting, restart di produzione e mass exchange simultaneo.\n\n**Migrazione mass exchange implementata il 28 ottobre 2026 (commit `cfe95fa`–`e64bd76`; CI ancora da eseguire).** `MASS_EXCHANGE_PARAMETERS` possiede ora `erosion_coeff`, `settling_flag`, `loss_rate`, `entrainment_flag`, `gas_loss_flag`, `liquid_vaporization_flag` oltre alla tripletta indivisibile `mass_exchange_method/cadence/group`. I sei parametri sono stati rimossi dai NAMELIST `SOLID_TRANSPORT_PARAMETERS`, `LIQUID_TRANSPORT_PARAMETERS`, `GAS_TRANSPORT_PARAMETERS`, `PORE_PRESSURE_PARAMETERS`; i moduli che memorizzano il dato restano invariati per non duplicare il dato fisico. La lettura avviene prima dei controlli sui trasporti; il backup emette una sola namelist mass exchange. Sono stati migrati i 18 file `.inp` distribuiti mantenendone i valori. **Gli input esterni che usano i vecchi blocchi richiedono migrazione**: non e prevista alcuna compatibilita con la vecchia posizione dei parametri. La CI e la verifica dei generatori di input non ancora coperti restano gate aperti. Il mantenimento dell'ordine e della conservativita del batch fisico non modifica il relativo kernel.

#### D18.5 — Piano di implementazione coordinato

1. **Inventario e matrice dei default.** Identificare ogni processo fisico canonico, attivazione, namelist esistente, target e dipendenze conservative; per mass exchange censire tutti i parametri oggi distribuiti nelle namelist solid/liquid/gas/pore pressure, separando parametri di processo da proprietà fisiche condivise e definendo un solo owner; registrare default effettivi legacy per METHOD/CADENCE/GROUP, accoppiamenti, batch, gruppi vuoti ed eccezioni. Registrare le configurazioni già presenti in test, esempi e generatori.
2. **Completamento delle interfacce EXPLICIT — IMPLEMENTATO per reologie, cooling e pore pressure.** La costruzione risiede nei singoli moduli fisici; mantenere questa proprietà. Non reinserire il valutatore generale SI+IM. Nessun test numerico reologico EXPLICIT aggiuntivo richiesto in questa fase; regressioni esistenti preservate.
3. **Parametri opzionali nelle namelist fisiche.** Estendere la lettura in inpout_2d e nei moduli input pertinenti, usando nomi univoci prefissati per processo **salvo la sola tripletta METHOD/CADENCE/GROUP della nuova MASS_EXCHANGE_PARAMETERS**; inizializzare default prima del parsing. Controllare valori consentiti, interi positivi, input parziali, flag disabilitati e presenza di vecchie opzioni IMEX/SOURCE_GROUPS.
4. **Normalizzazione e registry.** In source_configuration_2d costruire una tripla efficace per ogni sorgente attiva; in source_registry_2d raggruppare per (CADENCE,GROUP), verificare METHOD identico, identità globalmente unica, target e capability, ordinare tramite GROUP. Eliminare execution_order come campo indipendente, SOURCE_GROUP_ORDER come impostazione canonica e ogni riferimento residuale a PHASE nel contratto corrente.
5. **Validazione numerica dei gruppi.** Forzare tutte le sorgenti mass exchange attive nello stesso gruppo indivisibile e rifiutare configurazioni che separino il batch o ne duplichino il metodo/cadenza; preservare il solve NEWTON congiunto e il sistema REAL/COMPLEX/complex-step per l'unione dei target; non suddividere il gruppo in sottosolver. Vietare gruppi ANALYTIC senza mappa congiunta e metodi non realmente supportati. Proteggere il batch mass exchange e il relativo aggiornamento del letto.
6. **Percorso STAGE senza regressioni.** Adeguare lo scheduler e i solver a consumare gruppi normalizzati senza alterare i tableaux RK, la semantica delle sorgenti esplicite, le correzioni numeriche o i passaggi SI/IM. Separare regole del tableau da sequenze di gruppi controllate dall'utente.
7. **Esecutore STEP multimodale.** Estendere step_source_executor_2d oltre l'attuale percorso esplicito: dispatch per EXPLICIT/NEWTON/ANALYTIC e per GROUP; nei gruppi EXPLICIT/STAGE ed EXPLICIT/STEP separare rigorosamente fase di valutazione di **tutti** i rate sullo stesso stato, assemblaggio del rate di gruppo e applicazione dell'aggiornamento, vietando aggiornamenti sequenziali tra sorgenti dello stesso gruppo; integrazione implicita completa su dt e mappa analitica su dt, aggiornamento coerente di Q e QP, gestione degli stati secco/quasi secco e propagazione errori. Riusare i kernel fisici e il motore Newton, evitando formule duplicate.
8. **Integrazione con runtime e output.** Collegare la nuova esecuzione STEP fra assemblaggio RK, applicazione delle transazioni finali, controlli di stato, aggiornamento topografico e avanzamento/output/restart; verificare finestre temporali, topografia evolutiva, stato stocastico e impossibilità di doppia applicazione.
9. **Migrazione documentale e codice morto.** Aggiornare tutorial, input, template, generatori, test, manifest e documentazione del formato effettivo; rimuovere letture/adattatori superati una volta approvato il passaggio. Nessun dispatcher parallelo nascosto per sorgenti legacy.

Ogni passo deve costituire una modifica **sostanziale e coordinata** con i test pertinenti; non disperdere la migrazione in piccole modifiche senza chiusura dell'intero contratto. Non cambiare formule fisiche, tolleranze o oracoli per ottenere PASS.

#### D18.MX - PRIORITA' 1: architettura unificata del mass exchange (28 ottobre 2026)

**Stato: SPECIFICA DEFINITA / IMPLEMENTAZIONE DA COMPLETARE.** Questa sezione descrive direttamente l'architettura obiettivo per tutte le sorgenti di scambio di massa. E' il contratto normativo per i lavori futuri; D17.5, D17.8-D17.10 e D18.0-D18.6 documentano anche implementazioni storiche, i cui stati CI non certificano ancora questa architettura. Priorita' immediata prima dell'estensione delle leggi fisiche.

**Principi strutturali.**
- **Un solo coordinatore/assemblatore:** `src/physics/mass_exchange/mass_exchange_2d.f90` seleziona i processi attivi, richiama le rispettive interfacce, risolve i vincoli condivisi e accumula direttamente le sorgenti conservative **processo per processo**. Non esiste un modulo separato `mass_exchange_balance_2d` nell'architettura finale; assorbire le sue responsabilita' e quelle di `mass_exchange_laws_2d` nel coordinatore unico, eliminando duplicazioni e gli adattatori non necessari.
- **Una responsabilita' fisica per processo:** ogni modulo dedicato (erosione, deposizione granulare, liquid loss, Darcy gas loss, air entrainment, vaporizzazione da acqua del substrato, deposizione lavica) calcola l'intera fisica del proprio meccanismo, inclusi i contributi di massa e di energia. La legge non scrive nel vettore `eqns_term`, non aggiorna direttamente letto/inventari e non implementa le formule di accumulo conservative generiche.
- **Un solo applicatore temporale:** l'executor numerico riceve le sorgenti assemblate e i dati per inventari/geometria; applica STEP e gli aggiornamenti previsti dal contratto temporale una sola volta, senza selezionare leggi o contenere formule specifiche dei singoli processi.
- **Estensibilita' e minimo numero di file:** aggiungere un nuovo *modello* a un processo gia' esistente deve idealmente modificare solo il modulo di quel processo (oltre a input e test). Un processo fisico nuovo richiede la registrazione D17 e la connessione nel coordinatore. Nessun secondo registry, scheduler, contenitore di ID fisici numerici o gerarchia di tipi obbligatoria.

**Modalita' fisiche.** `lava_flag=T` identifica una fase continua lavica, con `gas_flag=F` e `liquid_flag=F`; la lava e' distinta dalla modalita' liquida. Con `lava_flag=F` rimangono valide le configurazioni gas, liquido e gas+liquido (gas continuo e liquido disperso). La selezione della reologia, del raffreddamento e delle sorgenti rimane indipendente dai flag di fase: `rheology_model='lava'` non determina `lava_flag`. Validare le combinazioni incompatibili prima del timestep.

**Matrice dei processi e delle compatibilita'** (G=gas; L=liquido; GL=gas continuo e liquido disperso; V=lava; 'ammesso' non significa 'attivo'; 'da verificare' non autorizza implicitamente la combinazione):

| Processo fisico | G | L | GL | V | `mixture_mass_rate` del processo | Altri rate/output fisici |
|---|---|---|---|---|---|---|
| Erosione granulare | ammesso | ammesso | ammesso | escluso inizialmente | **positivo**: massa erosa, solidi e materiale interstiziale | `solid_erosion_rate(:)`, `continuous_phase_erosion_rate`, `energy_rate` dell'ingresso dal substrato |
| Deposizione granulare | ammesso | ammesso | ammesso | escluso inizialmente | **negativo**: solidi e fase interstiziale depositati | `solid_deposition_rate(:)`, `continuous_phase_deposition_rate`, `energy_rate` sottratto |
| Perdita prescritta di liquido | escluso | ammesso | ammesso | escluso | **negativo**: liquido sottratto alla miscela | rate della componente liquida se dispersa; `energy_rate`; nessun deposito |
| Perdita gas/Darcy | con pressione interstiziale | escluso | da verificare | escluso | **negativo** se il gas abbandona il flusso | rate delle componenti pertinenti e `energy_rate`; nessun deposito |
| Entrainment aria | ammesso | **escluso** | da verificare | escluso inizialmente | **positivo**: aria ambiente entrante | `energy_rate` alla temperatura ambiente; nessun deposito |
| Vaporizzazione da sedimentazione in acqua del substrato | ammesso | escluso | ammesso | escluso inizialmente | **positivo**: massa del vapore entrante | **`added_gas_rate(1)` positivo** e `energy_rate`; stessa massa in equazione totale e prima specie gassosa |
| Deposizione/immobilizzazione lava | escluso | escluso | escluso | ammesso | **negativo**: massa lavica immobilizzata | `continuous_phase_deposition_rate`, `energy_rate` sottratto; aggiornamento di `lava_deposit` e topografia a valle |

**Contratto uniforme delle interfacce fisiche.** Ogni legge restituisce *soltanto* gli output pertinenti fra:
1. **`mixture_mass_rate` [kg m^-2 s^-1]**, con segno: positivo per ingresso nel flusso, negativo per uscita dal flusso (inclusa la deposizione). E' il contributo di massa totale del singolo processo, non un rate volumetrico e non un accumulato globale. La routine converte le quantita' volumetriche usando le densita' appropriate dei materiali coinvolti; non si sommano direttamente rate volumetrici di fasi con densita' diverse.
2. **`energy_rate`**, con segno e nelle unita' della corrispondente equazione conservativa: ogni processo lo calcola dalle proprieta' termiche e dalle temperature fisicamente pertinenti, comprese quelle del substrato o dell'aria entrante. Non usare automaticamente temperatura/densita' della miscela per tutte le sorgenti. L'energia effettiva deve essere coerente con la massa effettivamente trasferita.
3. **Rate delle componenti** interessate (classi solide, liquido disperso, gas aggiuntivi), per aggiornare le relative equazioni senza sommare nuovamente la loro massa all'equazione totale. Le componenti seguono le convenzioni conservative e i segni della loro equazione.
4. **Rate di erosione/deposizione** ove il processo modifica substrato o deposito: `solid_erosion_rate(:)`, `solid_deposition_rate(:)`, `continuous_phase_erosion_rate`, `continuous_phase_deposition_rate`, con unita' volumetriche [m s^-1] dove definite in questo modo. Questi output rappresentano il *destino* del materiale e servono per inventari e topografia; sono gia' contabilizzati in `mixture_mass_rate` e non devono generare una seconda sorgente massica.
5. Altri output soltanto dove effettivamente necessari a variabili ausiliarie della legge o a una cinematica diversa da quella standard. **Non** restituire per ogni processo due sorgenti di momento, un vettore completo `eqns_term` o campi fittizi privi di significato.

**Accumulo conservativo nel solo `mass_exchange_2d`.** Per ciascun processo con output *definitivi dopo i limiter*, il coordinatore aggiorna subito:
```fortran
eqns_term(1) = eqns_term(1) + mixture_mass_rate
eqns_term(2) = eqns_term(2) + r_u * MIN(0.0_wp, mixture_mass_rate)
eqns_term(3) = eqns_term(3) + r_v * MIN(0.0_wp, mixture_mass_rate)
eqns_term(4) = eqns_term(4) + energy_rate
! Aggiungere anche i rate delle sole componenti coinvolte.
! Accumulare separatamente i dati per substrato, deposito e topografia.
```
Questa convenzione di momento descrive la **cinematica attuale**: i materiali entranti da erosione, entrainment e vaporizzazione non portano momento orizzontale iniziale; le masse uscenti lo sottraggono alla velocita' locale del flusso. L'accumulo va effettuato **per singolo processo prima di eventuali cancellazioni fra rate positivi e negativi**. Per esempio, erosione positiva e deposizione negativa possono avere massa netta zero ma produrre un termine di momento negativo. Se si introduce una legge con momento non nullo all'ingresso, la cinematica supplementare deve essere resa esplicita nel contratto e nei test: la formula con `MIN` non e' un'identita' generale della fisica.

**Sequenza e limiter.** Una legge valuta il proprio meccanismo e i propri vincoli locali. Il coordinatore risolve i limiter condivisi *prima* di rendere definitivi i contributi; se piu' processi consumano lo stesso inventario, non possono impegnarlo due volte ne' dipendere arbitrariamente dall'ordine delle chiamate. Non accumulare nel bilancio output provvisori che saranno successivamente modificati: in presenza di limiter accoppiati, valutare prima i rate candidati, applicare il vincolo congiunto e solo allora accumulare i contributi finali dei processi. **Massa, energia, rate delle componenti e rate del deposito devono corrispondere alla medesima quantita' trasferita dopo i limiter**; ricalcolare i contributi dipendenti se cambia il rate limitato. Il coordinatore non replica le leggi fisiche.

**Regole per i singoli processi e gli inventari.**
- **Erosione granulare:** la routine comprende solidi e fase continua del substrato e calcola anche l'energia introdotta alle proprieta'/temperatura del substrato. La frazione prelevata va conteggiata negli inventari erodibili e non costituisce una seconda sorgente di massa.
- **Deposizione granulare:** la routine produce sia i rate per classe solida sia `continuous_phase_deposition_rate` per il materiale interstiziale, con trattamento fisicamente corretto di `coeff_porosity` e `maximum_solid_packing`. L'incremento geometrico considera la porosita' granulare e l'inventario appropriato.
- **Perdita liquida e Darcy:** le perdite che non formano deposito hanno `mixture_mass_rate<0` e non incrementano gli inventari. Anche la perdita di liquido *disperso* in GL riduce la massa totale della miscela e modifica la relativa equazione di componente; non presupporre che esista un'equazione separata per la fase continua.
- **Entrainment:** la sorgente introduce massa ed energia dell'aria ambiente. In solo liquido (`gas_flag=F,liquid_flag=T,lava_flag=F`) e' **vietato**, con rifiuto della configurazione quando `entrainment_flag=T`; GL resta da verificare.
- **Vaporizzazione:** l'acqua inizialmente si trova nel substrato (`B_zone`) e viene vaporizzata dal deposito di particelle calde; la massa del vapore contribuisce una volta a `mixture_mass_rate>0` e contemporaneamente a `added_gas_rate(1)`, senza doppio conteggio della massa totale. Richiede `gas_flag=T` e `n_add_gas>=1`, non `liquid_flag=T`. `added_gas(1)` e' **sempre** il vapore acqueo in questa configurazione; mantenere `eqns_term(5+n_solid)`, `sp_heat_g(1)`, `sp_gas_const_g(1)`, senza introdurre un indice configurabile della specie.
- **Deposizione lavica:** la legge fisica restituisce `mixture_mass_rate<0`, `energy_rate` e `continuous_phase_deposition_rate`; `lava_deposition_rate` puo' restare un **parametro del modello** ma non e' il nome dell'output. `lava_deposit` e' un **inventario** aggiornato a valle, non un output della legge. L'incremento topografico lavico non usa la porosita' granulare. La disponibilita' lavica va verificata sui consumi fisicamente concorrenti, senza sottrarre arbitrariamente un carrier_loss generico.

**Separazione fra bilancio e stato del fondo.** `mass_exchange_2d` restituisce i termini conservativi assemblati, i rate di erosione/deposizione e l'informazione minima per gli aggiornamenti del substrato, degli inventari `deposit(:)`/`lava_deposit` e della geometria. L'executor applica gli aggiornamenti una sola volta nella cadence prescritta, rispettando le maschere delle sorgenti e la politica dry/wet; non contiene formule fisiche laviche, granulometriche o di perdita. L'incremento di topografia **non** puo' essere dedotto dal solo `mixture_mass_rate`: occorrono rate di deposito, densita' e porosita' pertinenti.

**D17 e compatibilita' numerica.** Mantenere il registro canonico delle sorgenti, le quattro interfacce EXPLICIT/SEMI_IMPLICIT/IMPLICIT/ANALYTIC, capability e diagnostiche dei metodi non supportati, METHOD/CADENCE/GROUP, gruppi atomici e trattamento STEP/STAGE senza duplicare scheduler o registry. I cambiamenti dei moduli fisici e della struttura del coordinatore non devono alterare di nascosto il significato delle sorgenti conservative o l'ordine temporale. Le routine di processo non introducono dipendenze specifiche dello STEP executor.

**Realizzazione coordinata.**
1. Audit dei processi attuali: flag, input, densita', unita' dei rate, equazioni delle componenti, energia, momento, limiter, deposito/topografia, inventari, restart e gate CI; preservare quantitativamente la fisica legacy nel solo refactoring.
2. Introdurre `lava_flag` con validazioni dei modi fisici e delle combinazioni di sorgenti; aggiornare matrici di test delle modalita' G/L/GL/V.
3. Portare in ciascun modulo fisico i suoi `mixture_mass_rate`, `energy_rate`, rate delle componenti e di erosione/deposizione. Separare liquid loss, Darcy e deposizione granulare prima aggregati in carrier_loss; eliminare le formule specifiche dei processi dagli altri moduli.
4. Consolidare in `mass_exchange_2d.f90` il coordinamento e l'accumulo conservativo; eliminare `mass_exchange_balance_2d.f90`, `mass_exchange_laws_2d.f90` e `add_lava_deposition_balance` una volta migrati tutti i chiamanti, senza mantenere wrapper permanenti superflui. Alleggerire lo STEP executor.
5. Eseguire test unitari delle leggi e dei limiter, test con processi simultanei e massa netta nulla ma momento non nullo, oracoli indipendenti su energia e specie, bilanci massa-energia-momento-inventari-topografia, vaporizzazione in `added_gas(1)`, lava, casi dry/wet, restart e regressioni storiche; verificare GitHub Actions. Qualsiasi modifica della fisica va motivata e validata separatamente.

**Criterio di completamento:** D18.MX e' COMPLETATO soltanto quando il codice implementa l'unico coordinatore/assemblatore, ogni processo possiede i propri output fisici coerenti dopo i limiter, sono eliminati gli assemblatori duplicati, e passano tutti i test indicati. I risultati CI storici rappresentano la baseline, non l'accettazione di questa specifica.

#### D18.6a — Audit di chiusura del 27 ottobre 2026 (master; nuovi test in attesa della CI)

**Implementazione:** triplette METHOD/CADENCE/GROUP operative in produzione; percorsi STAGE esplicito, Newton e analitico e STEP esplicito, Newton e analitico collegati ai kernel esistenti; batch mass exchange esplicito/STEP mantenuto atomico. Rimossi i principali adattatori legacy. La configurazione guida lava NEWTON/STAGE + cooling ANALYTIC/STEP e stata verificata a livello di partizione, Newton di stage e mappa analitica cell-local (CI segnalata verde dopo `eb3b460`).

**Audit condivisione (commit `1d6ed4e`, `fa242ca`):** adattatore geometrico unico `eval_cell_source_explicit` per STAGE e STEP; contesto fisico/stocastico analitico STEP popolato con `cell_analytic_source_context`; i due coordinatori temporali restano distinti perche una sorgente STAGE contribuisce al tableau mentre STEP aggiorna lo stato completato. Questi commit richiedono ancora conferma dalle Actions.

**Nuovo test di produzione (commit `493dc7e`, `375cf82`; CI da ripetere):** `TEST_MASS_EXCHANGE/test_mass_exchange.f90::test_analytic_cooling_step_executor` passa attraverso `execute_step_source_groups` e confronta l'energia con l'oracolo del cooling esatto su dt. Controlla aggiornamento del primitivo, conservazione delle altre variabili conservative, celle asciutte e celle non attive. La fixture e isoterma solo nel senso della legge applicata: resta una fase liquida per una conversione fisica coerente. Questo test copre la fase post-RK, ma **non** attraversa in una singola simulazione il precedente RK Newton della lava.

**Avanzamento 27 ottobre 2026 (commit `2e875d9`; CI da verificare):** il test `test_analytic_cooling_step_executor` include ora un solve reale con `solve_imex_stage` per la reologia della lava `NEWTON/STAGE` su una cella attiva, verifica che l'energia resti invariata nello stage e usa lo stato risultante nella chiamata di produzione `execute_step_source_groups` con cooling `ANALYTIC/STEP`, confrontando l'energia di tutte le celle attive con la mappa esatta. Si tratta di un test di **composizione di uno stage implicito piu executor post-RK**, non ancora di un avanzamento completo attraverso tutti gli stage del tableau RK e il ciclo `simulation_2d`; non va presentato come validazione completa del time integrator.\n\n**Gate ancora aperti (non dichiarare COMPLETATO prima della CI e delle prove pertinenti):**
1. Integrazione completa di produzione lava NEWTON/STAGE + cooling ANALYTIC/STEP nello stesso avanzamento e verifica che il cooling sia applicato una sola volta.
2. Oracoli indipendenti per gruppi EXPLICIT multi-sorgente STAGE/STEP sullo stesso stato iniziale, ordine invariante all'interno del gruppo; NEWTON/STEP multi-sorgente con residuo di sistema accoppiato.
3. Regressioni di default del mass exchange atomico, conservazione/topografia e restart; test di splitting e raffinamento temporale per lava termodipendente.
4. Audit finale di input, rifiuti capability, namespace e ownership dei gruppi; aggiornamento delle tabelle storiche D18.4b, da mantenere etichettate come fotografia del 26 ottobre e non come stato corrente; registrazione della CI verde su master.

**Avanzamento di validazione (27 ottobre 2026, commit `4251b3d` e `e7d2367`; CI da eseguire):** aggiunto in `TEST_RHEOLOGY_CELL_LOCAL` un test locale di un unico gruppo NEWTON/STEP che comprende reologia della lava, cooling e pressione interstiziale, con controllo del residuo conservativo sulle quattro equazioni implicite e dell'invarianza delle equazioni non target. In `TEST_MASS_EXCHANGE` il test production STEP cooling controlla ora anche la proprieta di composizione dell'operatore analitico (due intervalli dt contro soluzione esatta su 2 dt). Il test storico `run_restart` continua a verificare bit-per-bit stato, geometria e inventari del mass exchange dopo checkpoint, ma non e ancora un restart del caso lava misto.

**Gate non completati:** (i) avanzamento vero attraverso tutti gli stage RK di `time_integration_workspace_type%advance` e successivo STEP cooling in un caso fisico della lava; (ii) gruppi EXPLICIT/STAGE ed EXPLICIT/STEP con piu sorgenti, usando uno stato comune e un oracolo di somma indipendente; (iii) restart della simulazione lava con split e regressione quantitativa contro configurazioni legacy; (iv) raffinare dt e documentare accuratezza e splitting. Non interpretare i test locali sopra come copertura di questi gate.

**Gate full RK aggiunto (27 ottobre 2026; commit `6401132`, `bb18ed0`; CI non ancora verificata):** il driver gia usato da `TEST_IMEX_STAGES` esegue ora, su richiesta della fixture, tutta la sequenza `time_integration_workspace_type%advance` con lava `NEWTON/STAGE`, seguita dalla chiamata reale `execute_step_source_groups` per cooling `ANALYTIC/STEP`. Per isolare il problema termico il registry post-RK del test contiene soltanto il gruppo cooling, senza mass exchange: questo non sostituisce un test con tutti i processi attivi. Si controllano mappa analitica esatta per cella, invarianza delle variabili non target e presenza di stage Newton. `compare_stages.py` esegue questa fixture per `n_RK=2,3,4`, due passi e 1/4 thread, confrontando byte-per-byte i risultati paralleli. Resta da raccogliere l'esito della CI.

**Ancora aperti per chiusura D18:** test EXPLICIT STAGE/STEP multi-sorgente con oracoli di somma indipendenti, restart della lava mixed-mode, raffinatezza temporale e accuratezza dello splitting, regressione con mass exchange simultaneo, audit del parser e dell'ownership delle configurazioni.

**Avanzamento EXPLICIT multi-sorgente (28 ottobre 2026; commit `7bb6154`; CI da verificare):** `TEST_MASS_EXCHANGE` esercita un gruppo con cooling termico esplicito e curvatura geometrica, ciascuno con contributo indipendente alla sorgente. Per entrambe le cadenze STAGE e STEP confronta la valutazione congiunta con la somma delle valutazioni separate sullo stesso stato conservativo, includendo il controllo di invarianza rispetto all'ordine dei nomi; per STEP confronta inoltre l'incremento attraverso `execute_step_source_groups` con `dt` per il tasso calcolato dallo stato iniziale. Il test STAGE verifica direttamente l'adattatore di produzione comune all'integratore RK, ma non e un test integrale con piu sorgenti esplicite attraverso l'intero tableau. Non dichiarare il gate chiuso prima della CI. Restano restart della lava, splitting/raffinamento temporale e regressioni combinate mass exchange.

**Ripresa del caso lava (28 ottobre 2026; commit `26c173f`, `f11e24e`; CI da verificare):** il driver `TEST_IMEX_STAGES` conserva e ripristina lo stato conservativo e primitivo e il tempo fisico usando un file binario della fixture. `compare_stages.py` confronta bit-per-bit gli endpoint di due passi continui con un primo passo seguito da un secondo processo che ricrea il registry e l'integratore e prosegue dal checkpoint, per `n_RK=2,3,4` nella configurazione lava NEWTON/STAGE + cooling ANALYTIC/STEP. **Limite importante:** il formato di checkpoint e specifico del test; non copre ancora `write_restart_file/read_restart_file` della produzione, ne topografia/inventari. Il gate di restart ufficiale rimane aperto; il test e una regressione di continuita del metodo e di ricostruzione del registry.

**Avanzamento 28 ottobre 2026 (commit `f512f4e`, `7b8be5b`, CI da verificare):** esteso `TEST_MASS_EXCHANGE` con un test dell'effettivo formato di restart di produzione (`write_restart_file/read_restart_file`) mentre il cooling e configurato `ANALYTIC/STEP`; dopo ripristino si ricostruisce il registry, si rigenera il primitivo e si confronta bit-per-bit l'endpoint con la continuazione non interrotta. La fixture non contiene ancora l'intero avanzamento RK della lava con mass exchange. In `TEST_IMEX_STAGES` aggiunta regressione di raffinamento a tempo fisso (2/4/8/16 passi) dell'intera sequenza lava NEWTON/STAGE + cooling ANALYTIC/STEP per 2/3/4 stage, con evidenza delle distanze fra soluzioni successive. Questa e una **verifica di Cauchy**, non una dimostrazione dell'ordine teorico dello splitting. Restano aperti il restart ufficiale con l'intero solver RK della lava e un benchmark indipendente dell'errore di splitting, oltre alla regressione simultanea con mass exchange attivo.

**Avanzamento restart mass exchange + cooling (28 ottobre 2026, commit `94ce2351`; CI da eseguire):** esteso il test ufficiale `TEST_MASS_EXCHANGE::test_official_cooling_restart` per attivare l'erosione durante il cooling `ANALYTIC/STEP` e confrontare bit-per-bit stato conservativo, topografia nodale, depositi, erosione e inventario erodibile al checkpoint e dopo la continuazione. La fixture richiede erosione non nulla; la sequenza utilizza `write_restart_file/read_restart_file`, non un checkpoint specifico del test. **Limite:** il caso non percorre tutti gli stage RK della lava e non dimostra ancora il restart integrato lava NEWTON/STAGE + mass exchange + cooling. Nessun gate viene chiuso prima delle Actions verdi.\n\n#### D18.6b — Audit aggiornato al 28 ottobre 2026 (HEAD `66de2acd`; CI verificata)

**Stato complessivo D18: PARZIALE AVANZATO — validazione integrata lava completata per i gate implementati, accettazione scientifica generale D18 ancora APERTA.** Questo audit prevale sugli stati storici D18.4b e sulle note con `CI da verificare` relative ai commit successivamente validati. Non reinterpretare la fotografia del 26 ottobre come elenco di difetti ancora presenti.

| Area / gate | Stato verificato | Evidenza e limite |
|---|---|---|
| Sorgenti canoniche, capability e quattro interfacce | **IMPLEMENTATO; REGRESSIONE VERDE** | 25 sorgenti nel catalogo; `mass_exchange_lava_deposition` con quattro entry point fisici (solo EXPLICIT implementato; gli altri rifiutano). `D17 source registry tests` verde a `7dd1d2ca`, run `40741563524`. |
| Build e dipendenze del modulo fisico lava | **IMPLEMENTATO; BUILD VERDE** | `D17 compile diagnostics` verde a `95fb529b`, run `40740965220`; successivo aggiornamento test non modifica le dipendenze Fortran. |
| STEP GROUP ordinati e batch mass exchange atomico | **IMPLEMENTATO; TESTATO** | Esecutore post-RK con `mass_exchange` EXPLICIT/STEP GROUP=0 e cooling ANALYTIC/STEP GROUP=1; regressioni legacy e nuove verdi. |
| Deposizione di lava immobile | **IMPLEMENTATO; REGRESSIONE VERDE** | `LAVA_DEPOSITION_FLAG` off per default; modelli CONSTANT/TEMPERATURE/VISCOSITY, viscosita condivisa con reologia; sottrazione conservativa dal flusso, inventario indipendente `lava_deposit`, incremento del letto senza porosita granulare; limite sullo spessore disponibile. Il test attualmente copre soprattutto configurazioni pure-liquid. |
| Formato di restart e inventari | **IMPLEMENTATO; REGRESSIONE VERDE CON LIMITI** | Formato binario aggiornato da v5 a v6, con inventario `lava_deposit`; `TEST_MASS_EXCHANGE` verifica round-trip ufficiale e continuazione bitwise con cooling/erosione e separatamente con deposito di lava. I file v5 non sono compatibili con v6. |
| Composizione lava NEWTON/STAGE + deposizione VISCOSITY/STEP + cooling ANALYTIC/STEP | **VALIDATA** | `TEST_MASS_EXCHANGE::test_integrated_lava_deposition_cooling_restart` usa Newton reale di stage e due passi STEP, includendo feedback termico, bilanci e restart ufficiale; **non** e un doppio avanzamento RK completo. |
| Avanzamento full RK/IMEX n_RK=2/3/4 | **VALIDATO PER LA FIXTURE** | `TEST_IMEX_STAGES`: `integration%advance` seguito dall'executor STEP reale; confronto continuo versus checkpoint della fixture, stato Q/QP, topografia e `lava_deposit`. `D17 historical physics regressions` verde su `66de2acd`, run `40746443339` (include test full-RK aggiunto). |
| Due passi full RK + **restart binario ufficiale** nello stesso esperimento | **APERTO (rafforzamento)** | La fixture full-RK utilizza un checkpoint stream specifico del test; i test con `write_restart_file/read_restart_file` percorrono invece la transazione STEP con stage Newton limitato. Non sommare i due esiti dichiarando che questo particolare caso e stato eseguito. |
| Accuratezza dello splitting e confronto scientifico quantitativo | **APERTO** | Le verifiche di raffinamento/Cauchy non costituiscono una dimostrazione dell'ordine asintotico dello splitting o un benchmark indipendente dell'errore/costo. |
| Acceptance D18 generale: parser, migrazione completa, default e tutte le combinazioni | **APERTO** | Continuano a valere i requisiti della matrice D18.6 su ownership dei parametri, input legacy, gruppi multi-sorgente, compatibilita e oracoli indipendenti; le tre Actions verdi non dimostrano da sole l'intero contratto scientifico D18. |

**Evidence Actions (GitHub):** compile run `40740965220` (`95fb529b`, SUCCESS); source registry run `40741563524` (`7dd1d2ca`, SUCCESS); historical physics regressions run `40746443339` (`66de2acd`, SUCCESS). Lo stato dell'HEAD documentale resta da validare se sono introdotte ulteriori modifiche al codice.

**Decisione di chiusura:** considerare completato il *gate integrato di deposizione lava attualmente implementato*, senza marcare **D18 COMPLETATO**. Prossimi interventi prioritari: (1) test end-to-end full RK con restart ufficiale; (2) audit puntuale input/default/ownership e compatibilita restart; (3) studio quantitativo di splitting, accuratezza e costo prima dell'accettazione scientifica finale.

#### D18.6 — Matrice di regressione e criteri di accettazione

| Caso | Risultato obbligatorio |
|---|---|
| lava NEWTON/STAGE GROUP=1; cooling ANALYTIC/STEP GROUP=1 | Due gruppi distinti; Newton agli stage e un solo cooling esatto post-RK |
| lava e cooling NEWTON/STAGE con stesso GROUP | Un **solo** sistema Newton accoppiato quantità di moto–energia |
| cooling ANALYTIC/STAGE oppure NEWTON/STAGE | Modalità stage originale preservata |
| gruppo multi-sorgente **EXPLICIT/STAGE** | Tutti i rate valutati sullo stesso stato di stage, sommati e inseriti una sola volta nel contributo EX del tableau; nessun aggiornamento intermedio |
| gruppo multi-sorgente **EXPLICIT/STEP** | Tutti i rate calcolati sul medesimo Q_in post-RK/post-gruppo precedente, somma prima dell'aggiornamento e un solo incremento conservativo sull'intero dt |
| cooling **NEWTON/STEP** | Solve implicito reale sul passo intero e residuo verificato |
| cooling **ANALYTIC/STEP** | Mappa esatta su dt, aggiornamento conservativo e nessuna chiamata negli stage |
| pore pressure **ANALYTIC/STEP** | Dissipazione post-RK realmente eseguita, quando attiva |
| gruppo NEWTON/STEP multi-sorgente valido | Un sistema accoppiato, non solve sequenziali |
| stessa coppia CADENCE/GROUP ma metodi diversi | Errore preventivo con nomi delle sorgenti |
| stessa sorgente assegnata a STAGE e STEP | Errore: processo duplicato |
| stesso GROUP con due cadenze diverse | Valido: namespace distinti |
| gruppo ANALYTIC senza operatore congiunto | Errore preventivo |
| sei processi mass exchange attivi con unica tripla in MASS_EXCHANGE_PARAMETERS | Unico gruppo, stesso METHOD/CADENCE/GROUP, rate sullo stesso stato e limiter post-somma |
| uno o più processi mass exchange disattivati | Un solo gruppo per il sottoinsieme attivo, nessuna sorgente fittizia |
| EXPLICIT/STAGE richiesto per il batch prima della sua validazione fisica | Rifiuto preventivo; nessun falso supporto |
| NEWTON o ANALYTIC richiesto per il batch senza un operatore congiunto implementato | Rifiuto preventivo |
| parametri mass exchange duplicati fra namelist legacy e nuova con valori discordanti | Errore esplicito; nessuna precedenza silenziosa |
| mass exchange spezzato violando limiter/bed | Errore preventivo |
| input legacy senza i nuovi parametri | Accoppiamenti, ordine e soluzioni della baseline riprodotti |

**Vincolo di valutazione sincrona per gruppi EXPLICIT:** per ogni gruppo G, sia STAGE sia STEP, definire lo stato di ingresso Q_G e ricavarne una sola vista primitiva coerente QP_G. Calcolare S_G(Q_G) = somma delle sorgenti S_s(Q_G), con ogni legge valutata su Q_G/ QP_G e senza che una legge modifichi ciò che leggono le altre. Assemblare e limitare il contributo aggregato secondo il contratto fisico del gruppo (per il mass exchange, un solo limiter dopo la composizione dei rate). Solo dopo l'assemblaggio completo applicare il contributo nel tableau EX di stage oppure un unico incremento esplicito su dt a STEP. Gruppi diversi possono invece vedere stati diversi se la semantica dello scheduler prevede l'aggiornamento fra gruppi: **non confondere il sequenziamento dei gruppi con quello, vietato, dei rate interni a un gruppo**.

**T1 — Contratto:** per ciascuna sorgente NEWTON inizialmente priva di EXPLICIT, test dell'identità EXPLICIT = SI + REAL(IM), stato comune immutato, segni e target, e regressione delle Jacobiane REAL/COMPLEX;  test del parser, dei default e dell'unicità globale; rifiuti espliciti, metodi/cadenze incompatibili, indici/target, nomi e gruppi. **T2 — Operatori:** includere test multi-sorgente EXPLICIT/STAGE ed EXPLICIT/STEP che confrontino la somma dei rate sullo stato comune contro un oracolo indipendente e che falliscano se l'ordine delle sorgenti altera il risultato; oracoli analitici per cooling e pore pressure; residui Newton STEP; bilanci conservativi, Q/QP, secco/quasi secco, limiter del mass exchange, contatori di quante volte viene eseguita ciascuna sorgente. **T3 — Integrazione:** test IMEX per n_RK=2/3/4 ove supportati, input legacy, lava con viscosità termodipendente, full/split binary restart, output/manifest, casi rappresentativi OpenMP e toolchain GNU/ifx secondo rischio. **T4 — Scientifico:** studio di raffinamento in dt e confronto fra NEWTON/STAGE, ANALYTIC/STAGE, NEWTON/STEP e ANALYTIC/STEP, quantificando accuratezza di splitting, stabilità e costo: nessuna equivalenza globale automatica.

**Gate D18 COMPLETATO solo quando:** le sorgenti NEWTON per cui SI+IM ricostruisce il rate completo possiedono EXPLICIT concretamente implementato e testato, con capability aggiornate senza regressioni NEWTON/complex-step;  tutte le sorgenti fisiche configurabili hanno una tripla opzionale nella namelist fisica pertinente **oppure, per i sei processi mass exchange, una singola tripla comune nella nuova MASS_EXCHANGE_PARAMETERS con i parametri di processo migrati e ownership fisica non duplicata**; ogni sorgente viene assegnata una volta sola; ogni gruppo ha METHOD omogeneo e GROUP come unico ordine configurabile; **EXPLICIT/STAGE** contribuisce correttamente al tableau EX, **EXPLICIT/STEP** applica un unico avanzamento esplicito su dt, e **NEWTON/STEP e ANALYTIC/STEP** sono implementati, realmente eseguiti e testati; il default conserva l'accoppiamento e il comportamento numerico legacy; nessun PHASE/execution_order indipendente sopravvive al contratto finale; mass exchange, topografia, conservazione, input e restart sono coerenti; test mirati e gate scientifici richiesti risultano documentati e superati.

---

### D15 — Acceptance finale della nuova architettura

**Test proporzionati (nuove modifiche a questo step):** gate scientifico finale sul codice consolidato da D16: matrice reologie e IMEX, conservazione, REAL/COMPLEX, restart, baseline A11 e variazioni numeriche autorizzate; decisione scientifica autonoma. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

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

Relocation sorgenti (completata) → refactoring x/y → riorganizzazione I/O e OUTPUT/\<RUN_NAME\>/ → politica Autotools B e packaging con due archivi → separazione reologie e riorganizzazione solver IMEX → D16 (consolidamento) → D15 (accettazione scientifica finale) → E03 (accettazione architetturale) → nuovi sviluppi scientifici pianificati (G01–H05 e altre estensioni).

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

File interessati: src/state_2d.f90; src/domain_2d.f90; src/simulation_2d.f90; src/spatial_operator_2d.f90; src/time_integration_2d.f90; src/physics/mass_exchange/; src/numerics/step_source_executor_2d.f90; TESTS/TEST_MODEL_LAYOUT/\*.

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

File interessati: src/domain_2d.f90; src/hp_reconstruction_2d.f90; src/spatial_operator_2d.f90; src/physics/mass_exchange/; src/numerics/step_source_executor_2d.f90; nuovi test tile.

Criteri di accettazione: Risultati identici per tile sizes multipli, gestione stencil e rispetto bounds, OpenMP stabile.

Build rectangular tile descriptors over the owned local domain while retaining compact active lists as the correctness reference.

Tile halo/stencil reach must cover HP reconstruction, two-cell LS derivative geometry and curvature/source needs.

Treat tiles only as scheduling/locality units; they must not change wet/dry activity or vertex mass-exchange semantics.

Add tile-vs-compact-list regressions before removing any old traversal path.

### H03 — Decomposizione MPI

**Test proporzionati (nuove modifiche a questo step):** T1 MPI exchange e partizioni; T2 1/2/4 rank; T3 conservazione e regressioni rank-invariant selezionate. La selezione concreta segue la regola generale; non cambia lo stato o le evidenze storiche.

Stato ufficiale: DA FARE

Attività e verifiche operative: Decomposizione cartesiana MPI e halo exchange generico; B_vertex ownership condivisa e gather nodale across rank; continuità correlazioni stocastiche ai rank boundaries.

File interessati: Nuovo communication/MPI module; src/domain_2d.f90; src/geometry_2d.f90; src/physics/mass_exchange/; src/numerics/step_source_executor_2d.f90; src/stochastic_module.f90; nuovi test MPI.

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
| D08      | P4                        | RHEOLOGY_MODEL da numeri a stringhe (assorbito in D05) |
| D09      | P5                        | Correzione del solo attrito                          |
| D10      | P6                        | EX_terms / IM_terms                                  |
| D11      | P7                        | IM_group e metadati                                  |
| D12      | P7                        | Newton per gruppo e Jacobiano ridotto                |
| D13      | P8                        | Integrazione analitica esatta dei termini IM         |
| D14      | P8                        | Gruppi misti e validazione                           |
| D16      | —                         | Consolidamento leggi fisiche e assemblaggio IMEX     |
| D17      | —                         | Struttura unificata dei termini sorgente e separazione tra fisica e numerica |
| D15      | P9                        | Acceptance finale della nuova architettura           |
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