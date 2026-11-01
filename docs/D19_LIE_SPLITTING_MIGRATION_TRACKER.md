# D19 — Architettura Lie definitiva e stato della migrazione

**Branch di lavoro:** `master`  
**Stato:** revisione architetturale dei percorsi Lie principali conclusa; verifica esaustiva di firme/rami Fortran e validazione eseguibile ancora aperte.  
**Ultima revisione:** 2026-11-01  
**Scopo:** questo è il riferimento operativo unico per l'architettura desiderata, lo stato verificabile dell'implementazione e i lavori ancora necessari. Le sezioni descrivono esclusivamente il progetto finale: non contengono versioni storiche o regole alternative.

## 1. Architettura definitiva

### 1.1 Integrazione temporale Lie

L'unico integratore temporale globale è il **Lie splitting di primo ordine**. Ogni passo conserva lo stato originale `q_n` e lo stato fisico corrispondente.

- **Gruppo 0 (sempre esplicito):** valutare l'operatore spaziale su `q_n`, applicare il suo incremento per ottenere `q_sp` e verificarne l'ammissibilità. **Non applicare limiter ai rate spaziali.** Valutare tutti i rate sorgente del gruppo 0 sullo stato originale `q_n` e applicare i rispettivi incrementi allo stato destinatario `q_sp`. I limiter delle sorgenti utilizzano le disponibilità effettive dello stato destinatario e degli inventari.
- **Gruppi successivi:** elaborare nell'ordine configurato; ogni gruppo riceve lo stato completo e la geometria aggiornati al termine del precedente e avanza per l'intero `dt`. I metodi disponibili secondo le capacità fisiche del modello sono EXPLICIT, ANALYTIC e NEWTON.
- Ogni gruppo produce un risultato completo e ammissibile prima che gli aggiornamenti siano pubblicati. Il tempo globale avanza **una volta sola** per passo Lie.
- Se un gruppo fallisce, ripristinare l'intero passo: stato fluido, inventari, topografia e grandezze geometriche derivate.

### 1.2 Responsabilità dei moduli

| Modulo | Responsabilità finale |
| --- | --- |
| `spatial_operator_2d` e moduli dei flussi | Calcolare **solo** operatori spaziali e flussi numerici; nessun limiter fisico sui loro rate. |
| `physics/source_terms_2d` e modelli specifici | Assemblare **solo** i termini sorgente e i trasferimenti fisici; calcolare e applicare i limiter ai **rate sorgente**. Restituire trasferimenti coerenti e già limitati. |
| `state_2d` | Definire e verificare ammissibilità e validità degli stati conservativi e primitivi. Le conversioni rimangono nelle interfacce dedicate allo stato. |
| `numerics/time_integration_2d` | Coordinare il passo Lie, l'avanzamento spaziale e la transazione dell'intero passo. |
| `numerics/lie_source_executor_2d` | Coordinare l'ordine dei gruppi, applicare nel tempo i loro incrementi, predisporre e pubblicare in modo consistente gli stati evoluti **alla fine del gruppo**. Non implementare modelli fisici né limiter. |
| `numerics/nonlinear_solver_2d` | Risolvere il gruppo NEWTON senza effetti collaterali su depositi, erodibile o topografia durante residui, Jacobiani e iterazioni. |
| `geometry_2d` | Rappresentare geometria e topografia; calcolare o aggiornare le **quantità derivate** (quote geometriche, pendenze, curvature, coefficienti geometrici) a partire da una topografia già assegnata. **Nessuna integrazione nel tempo e nessuna dipendenza da `dt` o dal metodo Lie.** |

Depositi, erodibile e topografia sono grandezze che possono **evolvere** nel tempo. Il loro incremento temporale e la relativa accettazione competono al livello di integrazione, non a `geometry_2d`. La proprietà dei dati e la loro organizzazione saranno definite coerentemente con questa responsabilità; non si deve usare `geometry_2d` come sede delle formule `inventario + dt * rate` o `topografia + dt * rate`.

### 1.3 Contratto dei limiter e dei trasferimenti

Per un gruppo esplicito `G`, `source_terms_2d` assembla **tutti** i rate fisici accoppiati, inclusi i trasferimenti con il substrato, e calcola un **solo fattore di limitazione condiviso** `lambda_G` nell'intervallo `[0,1]`, riferito alle disponibilità dello stato e degli inventari destinatari.

Lo stesso fattore scala **una sola volta** il vettore conservativo, le variazioni di depositi, erosione, erodibile, componente carrier, gas, energia, lava e topografia. Nessun trasferimento fisicamente accoppiato può essere lasciato senza scala o scalato due volte. L'interfaccia restituisce una **transazione di rate già limitati**, non un fattore separato che l'esecutore deve applicare a propria discrezione.

I limiter si applicano esclusivamente ai **rate sorgente**, non all'operatore spaziale. Non aggiungere un limiter a posteriori agli incrementi NEWTON o ANALYTIC alterando il metodo; verificarne invece la corretta soluzione e l'ammissibilità.

### 1.4 Confine di gruppo, geometria e Newton

La geometria e gli inventari utilizzati per valutare un gruppo sono quelli **all'inizio del gruppo** e restano fissi per l'intera valutazione del gruppo, inclusi tutti i tentativi e le iterazioni NEWTON.

Per ciascun gruppo:

1. Assemblare i rate fisici usando la geometria e gli inventari iniziali del gruppo.
2. Eseguire il metodo scelto e calcolare, mediante l'integratore, i candidati di stato fluido, inventari e topografia. Verificare l'ammissibilità tramite `state` e i vincoli dei rispettivi dati.
3. Solo dopo l'accettazione di **tutte** le celle, pubblicare insieme gli incrementi temporali di depositi, erosione, erodibile, lava e topografia. Nessun aggiornamento temporale diventa visibile durante la valutazione delle altre celle del medesimo gruppo.
4. **Se la topografia è cambiata**, chiamare `geometry_2d` per ricalcolare le grandezze derivate, incluse pendenze e curvature, **prima del gruppo successivo**. Saltare il ricalcolo se non necessario.
5. Conservare la capacità di rollback **dell'intero passo Lie**, non soltanto dell'ultimo gruppo.

`geometry_2d` **non** calcola variazioni del tipo `dt * dB_dt`, non applica operatori di avanzamento temporale e non deve conoscere gruppi Lie, stadi o solver Newton. La proiezione spaziale dei valori di topografia da celle a vertici, essendo una **operazione geometrica**, può risiedere in `geometry_2d`, ma l'integrazione di tali incrementi nel tempo avviene nel livello numerico.

### 1.5 Newton: equazioni interessate

Per ogni gruppo NEWTON, le incognite sono **tutte e sole** le equazioni conservative i cui rate sorgente di gruppo sono strutturalmente non nulli, come unione dei processi del gruppo; non dipendono da quali rate siano accidentalmente nulli in una cella specifica. Risolvere il sistema backward Euler sul sottosistema completo con Jacobiano complesso-step. Non ripristinare ricostruzioni triangolari, partizioni RK o condizioni fisiche hard-coded nel solver. Durante la soluzione nessun inventario o dato geometrico viene aggiornato.

## 2. Stato di implementazione verificato

Le voci seguenti descrivono il codice **già registrato su `master`**, non attestano compilazione o correttezza numerica.

| Componente | Stato corrente |
| --- | --- |
| Integratore Lie globale, gruppo 0 e gruppi successivi | Implementazione di produzione presente; la conversione dello stato originale prima del flusso spaziale ora controlla esplicitamente l'esito. Completare audit e accettazione finale. |
| Registro gruppi e selezione delle incognite NEWTON | Validazione preventiva delle configurazioni mass-exchange solo EXPLICIT; correggere le capacità del catalogo per `rheology_froude`, `rheology_lucas`, `rheology_mui`: NEWTON non è implementato dai rispettivi dispatcher, mentre EXPLICIT e ANALYTIC restano disponibili. Le equazioni del blocco Newton sono costruite come unione strutturale dei processi ammessi. |
| `source_terms_2d`: assemblaggio esplicito e limiter comune | Presente `assemble_lie_explicit_transaction`; i trasferimenti conservativi e al substrato vengono limitati con un fattore coerente. Tutti i campi del trasferimento mass-exchange sono verificati per finitezza prima e dopo la scalatura fisica. |
| `state_2d`: controlli di ammissibilità | Le funzioni `physical_cell_state_finite` e `conservative_mass_admissible` sono procedure del modulo (non type-bound) e vengono utilizzate dai percorsi Lie. Ripristinata la corretta sezione `CONTAINS` di `state_type` con le sole procedure type-bound. |
| `lie_source_executor_2d`: staging per gruppo | Presenti candidati separati per stato fluido (`q_group`, `qp_group`) e per depositi, erosione, erodibile, lava e rate del fondo. I candidati sono pubblicati soltanto dopo la verifica dell'intero gruppo; la geometria derivata è aggiornata al suo confine. |
| `state_2d`: proprietà degli inventari | `state_type` possiede `deposit`, `deposit_tot`, `erosion`, `erosion_tot`, `erodible`, `lava_deposit` con allocazione/deallocazione nello stato; `geometry_2d` conserva solo topografia e geometria. Integratore e Lie executor ricevono il medesimo `state_type`. I campi aggregati sono aggiornati tramite `state%refresh_substrate_totals()` dopo l'inizializzazione della soluzione, al termine del passo Lie e dopo il caricamento restart; l'output ESRI legge i totali senza modificarli. |
| `geometry_2d`: calcolo dei candidati temporali | Calcolato dall'esecutore Lie a partire dai rate già limitati; nessuna routine temporale di candidatura nel modulo geometrico. |
| `geometry_2d`: aggiornamento del fondo | `project_cell_field_to_vertices` proietta campi geometrici senza conoscere il tempo; `assign_topography_vertices` assegna e valida fondo e cache; `begin/commit/rollback_geometry_transaction` descrivono snapshot geometrici indipendenti da Lie. |
| Transazione di rollback | Snapshot dell'intero passo separati per responsabilità: `time_integration_2d` conserva e ripristina depositi, erosione, erodibile e lava; `geometry_2d` conserva e ripristina soltanto fondo e campi geometrici derivati. |
| Pulizia delle API storiche RK/STEP/STAGE | In corso: eliminate la finalizzazione IMEX e le routine non operative di curvatura; rimosso il wrapper sorgenti semimplicito inutilizzato. Le interfacce fisiche dei modelli ancora necessarie restano disponibili. I test storici saranno aggiornati nella fase finale. |
| Compilazione e suite finale | **Non eseguite**: per decisione progettuale si effettuano soltanto dopo tutte le modifiche al codice di produzione. |

## 3. Passi necessari per completare l'implementazione

- [x] **Separare l'evoluzione temporale dalla geometria**: estrarre da `geometry_2d` il calcolo `B + dt * rate`. Conservare lì soltanto rappresentazione/proiezione geometrica e refresh di pendenze, curvature e coefficienti; l'integratore effettua il calcolo temporale sul fondo.
- [x] **Spostare la costruzione temporale degli inventari nel livello di integrazione**: preparazione di depositi, erosione, erodibile e lava tramite `dt` e i rate già limitati; nessuna di queste formule dentro `geometry_2d`.
- [x] **Rimuovere la gestione evolutiva da `geometry_2d`**: sostituire `commit_lie_group_inventories` e le altre interfacce che conoscono il progresso temporale con un'assegnazione coerente dei dati accettati e un'interfaccia puramente geometrica per aggiornare i campi derivati. Definire l'ownership finale delle variabili evolutive senza mischiare fisica e geometria.
- [x] **Differire la pubblicazione degli stati fino a fine gruppo**: candidati fluidi e inventari rimangono separati durante l'elaborazione delle celle. L'aggiornamento del fondo e della geometria derivata avviene dopo la validazione dei candidati. Il refresh geometrico usa ora una copia locale dei campi geometrici: se i campi derivati risultano non finiti, ripristina il fondo e tutte le cache senza esporre un aggiornamento parziale. La nuova topografia viene assegnata e validata prima di pubblicare gli inventari del gruppo; un refresh rifiutato lascia anche gli inventari e i candidati fluidi invariati. Rimangono da verificare nell'accettazione finale i percorsi di errore numerici e l'integrazione dei moduli.
- [ ] **Mantenere la geometria congelata durante Newton** e verificare che nessuna routine invocata nei residui/Jacobiani modifichi inventari o geometria. Audit delle interfacce: `eval_source_implicit` costruisce uno stato fisico complesso locale e inoltra esclusivamente rate alle routine `eval_rheology_implicit` e `eval_thermal_implicit`; il dispatch reologico esaminato riceve lo stato con `INTENT(IN)` e restituisce un rate con `INTENT(OUT)`. Rimane da verificare in dettaglio l'assenza di effetti collaterali all'interno di ogni implementazione specifica dei modelli. Audit statico completato per gli otto modelli reologici chiamati da `eval_rheology_implicit` (Voellmy, mui-Voellmy, plastic, lava, lahar, empirical, newmodel, quadratic) e `eval_thermal_cooling_implicit`: nei corpi esaminati non risultano scritture dirette su inventari/topografia, `INTENT(INOUT)` o dichiarazioni `SAVE`. Rimane necessaria la verifica finale dell'intera catena di chiamate e dei casi di esecuzione.
- [x] **Separare la proprietà degli snapshot di rollback**: l'integratore conserva e ripristina gli inventari temporali; `geometry_2d` conserva e ripristina fondo e cache geometriche. Rimane da verificare la corretta esecuzione dei percorsi di rifiuto nella fase finale, senza compilation/test intermedi.
- [ ] **Completare l'audit del codice di produzione**: eliminare le ulteriori API, variabili, commenti e dipendenze RK/STAGE/STEP residue senza introdurre compatibilità legacy. I due entrypoint espliciti interni a `source_terms_2d` sono ora denominati in base alla funzione fisica (`eval_geometry_cell_source_explicit`, `eval_mass_exchange_source_explicit`), mantenendo il dispatch generico pubblico. Verificata inoltre la corrispondenza dei nomi di apertura/chiusura delle 14 subroutine di `source_terms_2d` e corretto `END SUBROUTINE eval_source_analytic_state`. Eliminata la finalizzazione legacy `finalize_imex_cell_state` che eseguiva clipping post-IMEX; ricerca repository: riferimenti residui soltanto nei test storici, da aggiornare al termine della migrazione. Verificate le due chiamate di produzione a `eval_source_implicit` nel solver Newton: il residuo usa input/output reali, il Jacobiano usa input/output complessi. Entrambi i percorsi sono necessari; eliminato invece lo `STOP` non diagnostico per argomenti incoerenti e sostituito con `fatal_error`. Audit del percorso implicito: `eval_source_implicit` accetta esclusivamente coppie coerenti di input/output reali oppure complessi; rimosse le variabili locali inutilizzate che duplicavano campi del `complex_physical_state_type`. La routine valuta soltanto il sorgente, senza mutare inventari o geometria.
- [x] **Verificare la separazione di responsabilità rimasta**: revisione delle API di proiezione del fondo, delle variabili evolutive ancora esposte da `geometry_2d` e delle eventuali dipendenze nascoste tra gruppi. Eliminati `begin_topography_rate_update`, `record_topography_rate_cell`, `project_lie_bed_rate`, il workspace globale dei rate topografici e il relativo rilascio da `simulation_2d`: i rate temporanei restano nell'esecutore Lie e la geometria fornisce la proiezione generica cella-vertice. Rimane l'audit dello storage evolutivo residuo e dei percorsi di stato. Rinominate le API geometriche senza mantenere alias legacy. I test Fortran rilevanti sono stati aggiornati. Il substrato e i campi aggregati appartengono a `state_type`, non alla geometria. Trasferiti tutti e sei i campi nel tipo `state_type`, compresi i due aggregati; rimossa la proprietà globale dalla geometria. Migrazione delle routine di produzione effettuata senza compatibilità delle variabili globali obsolete; restano le verifiche statiche dei chiamanti e l'adeguamento dei test.
- [ ] **Completare il controllo statico end-to-end** di flussi, sorgenti, registro gruppi, Newton, conversioni e rollback prima di tentare la compilazione. Audit del gruppo 0: `qp_original` alimenta i rate, mentre `q_group` post-spaziale è usato come destinatario dai limiter. Il risultato esplicito è validato in `source_terms_2d` anche dopo la scalatura comune (rate conservativi e tutti i campi dei trasferimenti). Da approfondire ancora gli errori di conversione e l'insieme dei campi di stato condivisi. La ricostruzione pre-spaziale ora rifiuta esplicitamente conversioni non valide, senza ritoccare lo stato o limitare i flussi. Il registry rifiuta già all'inizializzazione un mass-exchange non esplicito. Allineate le capacità NEWTON del catalogo ai soli modelli con valutazione implicita effettiva: Froude, Lucas e μ(I) non possono usare NEWTON fino a quando non sarà implementata una corrispondente legge implicita; il metodo di default è EXPLICIT. Corretto l'aggiornamento dei totali deposito/erosione: operazione dello stato, non effetto collaterale di `output_esri`; sincronizzazione anche dopo restart. Corretto il posizionamento illegale delle due funzioni di validazione dentro il tipo derivato `state_type`: ora appartengono alla sezione `CONTAINS` del modulo. Completata la sincronizzazione degli aggregati anche dopo `read_solution`, e verificata esplicitamente la porosità del substrato prima del calcolo dei totali (porosità in [0,1)).
- [ ] **Soltanto alla fine**, compilare e correggere i problemi di integrazione, poi aggiornare/eseguire i test Lie e la suite di regressione completa.

## 4. Criteri di accettazione finale

- Il gruppo 0 calcola i rate spaziali e sorgente sullo stato originale ma limita **solo** i rate sorgente rispetto allo stato destinatario.
- Ogni transazione sorgente limitata mantiene identica la scala di tutti i trasferimenti fisici accoppiati.
- Depositi, erodibile e topografia cambiano **solo a fine gruppo**; il gruppo successivo usa slope e curvature aggiornate.
- L'integrazione temporale di tali grandezze non compare in `geometry_2d`.
- Newton usa geometria/inventari congelati durante le iterazioni ed esattamente l'insieme strutturale delle equazioni del gruppo.
- Tutti i gruppi condividono lo stesso `dt`, il tempo globale avanza una sola volta, e il fallimento di qualunque gruppo ripristina l'intero passo.
- Codice di produzione completo prima della compilazione e dei test.

## 5. Manutenzione del documento

Mantenere **una sola versione** del contratto in §1. Quando una decisione viene precisata, aggiornare direttamente la descrizione dell'architettura e la checklist, **senza aggiungere una nuova regola contraddittoria o una cronologia di decisioni**. Registrare i progressi soltanto nello stato di implementazione (§2) e nella checklist (§3). Un commit non implica compilazione o test superati.

### Audit inventari Lie (2026-11-01)

Il gruppo mass-exchange verifica finitezza e nonnegativita di deposito, erosione, erodibile e lava immobilizzata sia allo stato di ingresso sia per i candidati; nessun clipping.

### Migrazione organica dei test strutturali (2026-11-01)

- `TEST_MASS_EXCHANGE/test_projection_stress.f90`: usa `state_type` per gli inventari e la firma corrente di `execute_lie_source_groups(state,...)`; mantiene il confronto dei bilanci fisici e della geometria, senza facciate globali di deposito/erosione/lava.
- Transazioni di geometria: test aggiornati a `begin_geometry_transaction` / `rollback_geometry_transaction`.
- `TEST_WORKSPACE_LIFECYCLE`: eliminate le chiamate al workspace topografico temporale rimosso da `geometry_2d`.
- Non eseguiti compilazione o test; gli altri test legacy e i percorsi di produzione rimangono oggetto della verifica statica finale.

### Transazione degli inventari nello stato (2026-11-01)

- Introdotto `substrate_snapshot_type` in `state_2d`: snapshot con deposito, erosione, erodibile e lava immobilizzata.
- `state%snapshot_substrate` e `state%restore_substrate` possiedono la cattura e il ripristino degli inventari; il ripristino sincronizza anche i campi aggregati.
- `time_integration_2d` non conserva piu quattro copie indipendenti: orchestra la transazione dello stato e quella della geometria, mantenendo separata la topografia.
- Adeguate le diagnostiche test del mass exchange. Restano la verifica integrale dei test legacy e la compilazione finale concordata.

### Correzione critica di atomicita dei gruppi (2026-11-01)

- Audit del dispatcher Lie: il ramo ANALYTIC scriveva direttamente nel vettore `q` autorevole, anziche nel candidato `q_group`. Corretto per aggiornare solo `q_group(1:n_eqns,j,k)` prima dei controlli di ammissibilita e della conversione; il commit dei campi `q` e `qp` rimane alla fine del gruppo, come per EXPLICIT e NEWTON.
- Verificato staticamente che `solve_backward_euler_group` opera sul `q_group` passatogli per riferimento.
- Nessuna compilazione o test ancora eseguiti, secondo il piano.

### Chiusura del percorso analitico del gruppo 0 (2026-11-01)

- Il ramo ANALYTIC usa ora `q_original` per valutare il contributo del gruppo 0, come gli altri termini fisici del gruppo iniziale, e applica il contributo al candidato post-spaziale `q_group`.
- Per i gruppi successivi ANALYTIC usa il proprio stato di ingresso `q_group`. Il salto di aggiornamento non viene pubblicato prima della validazione dell'intero gruppo.
- L'audit statico delle vecchie API `begin_lie_geometry_transaction`, `rollback_lie_geometry_transaction`, `release_topography_workspace` e `project_lie_bed_rate` non trova chiamanti Fortran di produzione residui; rimane una menzione documentale storica.
- Non attribuire lo stato «validato» all'architettura finche non siano stati eseguiti compilazione e test di integrazione nella fase conclusiva.

### Eliminazione dell'equation partition legacy (2026-11-01)

- Rimosso da `simulation_context_type` il campo `equation_partition_type`, che non esiste piu nella sorgente.
- La chiamata `initialize_configured_sources` riceve ora soltanto registro e stringa d'errore, coerentemente con la sua interfaccia in `source_configuration_2d`.
- Le partizioni dei gradi di liberta impliciti sono definite dai metadati dei gruppi sorgente e non da una partizione globale parallela.

### Coerenza transazionale degli inventari per gruppo (2026-11-01)

- `substrate_snapshot_type` include ora anche `deposit_tot` e `erosion_tot`; il rollback ripristina esattamente i sei campi originali, senza ricalcolare durante il ripristino.
- Ogni gruppo mass-exchange pubblicato sincronizza i totali nello stesso confine di gruppo in cui pubblica deposito, erosione, erodibile e lava. I gruppi successivi osservano quindi uno stato di substrato internamente consistente.
- Lo stato completo rimane pubblicato dall'integratore solo alla conclusione del passo; in caso di errore, geometria e inventari tornano allo snapshot pre-step.
- Ancora non compilato o testato: la verifica dinamica resta fase successiva alla chiusura dell'audit statico.

### Audit statico delle dipendenze build (2026-11-01)

- `time_integration_2d` usa `lie_source_executor_2d`: corretta in `src/Makefile.am` la dipendenza oggetto `time_integration_2d.$(OBJEXT)` -> `lie_source_executor_2d.$(OBJEXT)` e riordinati i sorgenti.
- Verificati i riferimenti legacy principali in produzione: la chiamata `initialize_configured_sources(registry,error)` coincide con la firma, e gli inventari fisici non sono piu importati dalla geometria nei moduli centrali controllati.
- Audit statico non dichiarato esaustivo finche non vengono esaminati sistematicamente tutti i moduli e gli altri file di build. Compilazione e test non eseguiti come concordato.

### Esito ricognizione statica trasversale, 1 novembre 2026

**Ambito effettivamente verificato.** Lettura dei 64 file Fortran sotto `src/` e dei 29 programmi/test Fortran direttamente appartenenti alle directory di test (escluso il riferimento storico `TEST_HYDROSTATIC_ROUNDOFF/reference/hyperbolic_2d.f90`). Ricerca nei sorgenti delle interfacce rimosse (`begin_lie_geometry_transaction`, `rollback_lie_geometry_transaction`, `release_topography_workspace`, `project_lie_bed_rate`, `equation_partition_type`) e degli import di inventari da `geometry_2d`: nessuna occorrenza nei file controllati. Verificate le due chiamate di produzione del dispatcher e la firma `initialize_configured_sources(registry,error)`. In `src/Makefile.am` la dipendenza esecutore Lie -> integratore e l'ordine dei due sorgenti sono stati corretti.

**Confine dell'audit:** questa e una verifica statica delle interfacce obsolete e delle responsabilita dei moduli, non una dimostrazione automatica della correttezza di tutti i nomi importati con `USE ONLY`, di ogni ramo del codice o della compatibilita delle configurazioni. Il controllo completo della compilazione, la risoluzione degli errori di integrazione e la regressione restano necessariamente aperti; non devono essere segnati come verificati prima dell'esecuzione. L'audit semantico di ogni catena di chiamate non e stato effettuato esaustivamente e rimane aperto nel controllo statico end-to-end.


### Verifica incrociata dei percorsi centrali - 1 novembre 2026

**Revisione esaminata:** `2070ce698250adf8395f7cbd0ac38e84bd7e8d25`, prima dell'aggiornamento di questo documento. **Stato:** audit statico parziale, non ancora chiuso; nessuna compilazione o suite eseguita.

- Controllati in lettura i moduli di stato, geometria, esecutore Lie, integratore, Newton, configurazione dei sorgenti e `source_terms_2d`, con ricognizione delle clausole `USE` e delle procedure principali. Questa lettura non equivale a risoluzione esaustiva dei simboli esportati, delle firme dei chiamanti e dei test.
- In `lie_source_executor_2d` il ramo ANALYTIC usa il vettore originale per calcolare il contributo del gruppo 0 e aggiorna `q_group`; EXPLICIT valuta il candidato di gruppo e la conversione; NEWTON lavora sul candidato. Prima della pubblicazione sono controllate ammissibilita conservativa e conversione a primitive. L'assegnazione autorevole `q=q_group`, `qp=qp_group` avviene alla fine del gruppo.
- Per gruppi con mass exchange gli inventari vengono preparati come candidati, il nuovo fondo e le cache geometriche vengono accettati prima della pubblicazione degli inventari e viene invocato `state%refresh_substrate_totals()` al confine del gruppo.
- `time_integration_2d` controlla la conversione della soluzione originaria prima del flusso spaziale, rifiuta stati post-spaziali inammissibili senza limiter spaziale, crea snapshot degli inventari e della geometria e ripristina entrambi se l'esecutore Lie segnala fallimento. Questo controllo e' statico: il rollback numerico e la corretta copertura dei rami di errore restano da testare.
- In `nonlinear_solver_2d` la struttura Newton e' guidata dagli `equation_ids` del gruppo e usa `eval_source_implicit` in residuo reale e Jacobiano complesso, da verificare ancora in tutti i modelli.
- **Punto emerso:** `domain_2d` importa e usa `n_RK` come numero di passate per l'espansione delle maschere di celle attive. Non eliminarlo meccanicamente: verificare se rappresenta la larghezza dell'alone spaziale, eventualmente rinominandolo in base alla semantica, e controllare tutti gli altri usi del parametro. Non e' evidenza che vengano ancora eseguiti stage temporali RK.

**Verifiche aperte, da completare prima della compilazione finale:** (1) censimento automatico *integrale* dei moduli e dei loro simboli pubblici, incluse generic interfaces, procedure type-bound, rinomine e `USE ONLY`; (2) confronto tra firme dei chiamanti e delle definizioni, intent, rank e keyword arguments; (3) esame di tutti i percorsi e dei possibili effetti collaterali per gruppo 0, EXPLICIT, ANALYTIC, NEWTON, substrato, topografia e rollback; (4) aggiornamento degli eventuali test legacy incompatibili. Non segnare queste voci come completate in assenza di evidenza di scansione integrale. La compilazione resta l'ultimo passaggio di integrazione dopo la verifica statica.

### Eliminazione della dipendenza temporale dalle maschere spaziali (2026-11-01)

- Audit di `domain_2d`: il modulo importava `n_RK` da `parameters_2d` benché il simbolo non sia più definito; la compilazione sarebbe fallita. La sua unica funzione era determinare il numero di espansioni della maschera attiva, non l'integrazione temporale.
- Sostituito con la costante locale `active_halo_layers=2`: due strati geometrici, indipendenti dai gruppi Lie, per il calcolo delle interfacce e la ricostruzione spaziale a tre punti. Il significato è ora spaziale, non un retaggio RK. La scelta di due strati va confermata tramite i test del bordo asciutto e delle worklist di ricostruzione.
- **Non è ancora una prova dell'assenza di altri simboli mancanti**. La risoluzione esaustiva di ogni `USE ONLY`, delle firme, dei percorsi semantici e la compilazione finale rimangono aperte. Nessuna compilazione o test eseguiti in questa verifica: il clone diretto GitHub dall'ambiente di esecuzione non ha accesso DNS.


### Controllo aggiuntivo delle interfacce critiche (2026-11-01)

- Confermato nel sorgente aggiornato che `domain_2d` non importa piu `n_RK` e usa la costante locale `active_halo_layers` per l'espansione delle maschere; gli altri moduli centrali esaminati non usano `n_RK`. La ricerca testuale di GitHub può restituire frammenti indicizzati della versione precedente: per questo l'esito e' fondato sul contenuto recuperato direttamente dal branch.
- Verificato che il generic `eval_source_explicit` in `source_terms_2d` e' dichiarato con due `MODULE PROCEDURE` effettivamente presenti (`eval_geometry_cell_source_explicit`, `eval_mass_exchange_source_explicit`): le due chiamate controllate si risolvono a interfacce differenti, non sono riferimenti a una routine inesistente.
- Confermate in `time_integration_2d` le chiamate ai metodi `state%snapshot_substrate` e `state%restore_substrate`, oltre alle transazioni geometriche; confermato in `lie_source_executor_2d` che gli aggiornamenti `q` e `qp` autorevoli seguono la validazione dei candidati.
- **Blocco operativo constatato**: `git ls-remote https://github.com/demichie/IMEX_SfloW2D_dev.git HEAD` nell'ambiente di esecuzione fallisce per DNS (`Could not resolve host: github.com`). Non esiste quindi una copia locale verificata dalla quale eseguire la compilazione; il connettore GitHub permette di consultare e registrare il codice, ma non equivale a compilazione.
- Restano aperti il confronto esaustivo di tutti i simboli `USE ONLY` e delle firme, la verifica di tutti i rami semantici e i test. Nessuna chiusura fittizia dell'audit statico o della compilazione.

### Verifica integratore e indice canonico (2026-11-01)

- Esaminata integralmente la routine `advance_lie`: lo stato originale viene copiato nei workspace, il flusso spaziale aggiorna solo `q_trial`, le conversioni e l'ammissibilita' sono controllate prima dei gruppi; il rollback di substrato e topografia avviene in caso di errore sorgente, mentre `state%q` e `state%qp` vengono pubblicati solo dopo il completamento di tutti i gruppi.
- Corretta una dipendenza implicita dall'ordinamento delle variabili: la ricostruzione pre-spaziale controllava `q_before(1,j,k)` per distinguere celle bagnate/asciutte, ora usa l'indice canonico `idx_totMassEqn`, gia' importato dal modulo. Non cambiano ne' il metodo temporale ne' la fisica del flusso.
- Restano da verificare esaustivamente le altre firme/import e i singoli percorsi numerici; compilazione non eseguita, poiche' il clone Git diretto continua a fallire per risoluzione DNS. Non dichiarare chiuso l'audit.


### Controllo statico mirato del registro e dei gruppi Lie (2026-11-01)

- Esaminata l'implementazione completa di `initialize_lie_registry` in `source_registry_2d`: richiede il gruppo iniziale ID 0 con METHOD EXPLICIT; impone ID strettamente crescenti, nomi sorgente unici, metodi ammessi dal catalogo, insieme NEWTON uguale all'unione strutturale delle equazioni coinvolte e mass exchange soltanto EXPLICIT e in un unico gruppo. L'inizializzazione risulta transazionale: il registro viene inizialmente finalizzato e `initialized=.TRUE.` soltanto dopo tutte le validazioni.
- Verificata la valutazione in `lie_source_executor_2d`: nel gruppo ID 0 `qpentry` e' tratto da `qp_original`; nei gruppi successivi da `qp_group`. I rate EXPLICIT sono applicati al destinatario post-spaziale `q_group`; l'aggiornamento di ciascuna cella passa per i controlli di ammissibilita' e conversione. I rate di mass exchange e gli inventari sono preparati prima del commit di gruppo. Il ramo ANALYTIC legge `q_original` solo per il gruppo zero, sebbene la configurazione corrente richieda METHOD EXPLICIT per tale gruppo.
- Le procedure type-bound `source_registry%groups()`, `%is_initialized()`, `%finalize()` e `state%refresh_substrate_totals()` risultano dichiarate nei rispettivi tipi. Confermata la presenza del generic `eval_source_explicit` con due implementazioni specifiche.
- **Ambito di questa evidenza:** lettura e confronto statico di questi percorsi; non costituisce ancora una verifica esaustiva delle 64 unita' Fortran, di tutti gli `USE ONLY`, delle firme e dei rami fisici. Nessuna compilazione e nessuna GitHub Action avviate, per indicazione esplicita. D19 resta aperto finche' la verifica statica integrale non e' completata.

### Audit statico di reinizializzazione del registro (2026-11-01)

- Controllata la catena `initialize_configured_sources` -> `collect_active_source_names` -> validazione metodo/fase/gruppo -> `registry%initialize`. Prima della correzione i rami di errore anteriori a `registry%initialize` potevano lasciare utilizzabile un registro già inizializzato in una precedente configurazione.
- Inserito `CALL registry%finalize()` all'ingresso di `initialize_configured_sources`: anche gli errori precoci invalidano il registro e impediscono l'uso di metadati stale. `initialize_lie_registry` resta responsabile della successiva validazione transazionale.
- Esaminati i metadati del catalogo relativi ai metodi: Froude, Lucas e mui non dichiarano NEWTON; i mass exchange dichiarano il solo EXPLICIT. La configurazione e il registro controllano le capacità richieste.
- Nessuna compilazione o Action eseguita per accordo. Il controllo integrale di ogni firma e percorso, inclusi i test, resta aperto; questa nota non dichiara concluso D19.

### Chiusura del vincolo sulle masse di componente (2026-11-01)

- Controllo end-to-end di `conservative_mass_admissible` e dei campi di composizione in `state_conversion_2d`: le masse conservative dei solidi e dei gas addizionali sono entrambe parti della massa totale della miscela. Prima il guard verificava la nonnegativita' di ciascuna componente ma confrontava solo la somma dei solidi con la massa totale. Uno stato con solidi + gas addizionali superiori alla massa totale poteva superare questo guard e fallire soltanto durante la conversione o la chiusura delle frazioni.
- Corretto il controllo per imporre `SUM(solid masses) + SUM(additional gas masses) <= total mixture mass` senza clipping e senza limiter. Per gli insiemi di componenti vuoti la somma Fortran e' zero. Gli indici canonici sono gia' importati dal modulo.
- La correzione vale per stato post-spaziale e candidati dei gruppi in quanto condividono lo stesso guard. Non attesta superamento dei test: compilazione e GitHub Actions restano riservate all'utente dopo D19.
- Il censimento integrale delle firme/import e degli ulteriori percorsi resta aperto; audit non ancora dichiarato concluso.

### Limiter EXPLICIT e chiusura delle masse di componente (2026-11-01)

- Riscontrata incoerenza fra il controllo di ammissibilita' conservativa aggiornato e `limit_explicit_source_rates`: quest'ultimo limitava le diminuzioni di massa dei solidi, ma non quelle dei gas addizionali; inoltre definiva il carrier residuo sottraendo soltanto i solidi alla massa totale.
- Corretto il limiter in `source_terms_2d` per controllare separatamente i rate negativi di tutte le masse solide e dei gas addizionali con i rispettivi indici canonici. La disponibilita' e il consumo del carrier derivano ora da massa totale meno solidi meno gas addizionali.
- La scalatura resta unica per tutti i rate sorgente e i campi trasferiti tramite `scale_mass_exchange_transfer`; nessun limiter viene applicato al flusso spaziale. I nuovi controlli si allineano alla stessa condizione di ammissibilita' `m_s + m_g <= m_tot`.
- Verificato staticamente il percorso `assemble_lie_explicit_transaction` -> `limit_explicit_source_rates` -> controllo di finitezza dei rate/trasferimenti, ma resta da controllare la totalita' dei sottoprocessi e dei test. Nessuna compilazione/Action eseguita.

### Rettifica fisica del limiter della fase gassosa (2026-11-01)

- Su revisione della semantica fisica, la modifica al limiter del commit `f466782` e' stata corretta: i gas aggiuntivi sono componenti della fase gassosa complessiva, non un inventario carrier indipendente da sottrarre alla sua massa disponibile. La disponibilita' del carrier e' `m_tot - SUM(m_solidi)`.
- In `mass_exchange_deposition_2d` il rate negativo dei gas aggiuntivi e' proporzionale alle loro frazioni rispetto alla massa della fase gassosa; `mass_exchange_vaporization_2d` aggiunge invece massa al gas. Si torna quindi a limitare esplicitamente massa totale, singole masse solide e disponibilita' della fase carrier come insieme, **senza limiter aggiuntivi per ciascun gas**.
- Rimane il distinto controllo di ammissibilita' in `state_2d`: ogni massa trasportata e' nonnegativa e la somma di solidi e gas aggiuntivi non supera la massa totale. Tale controllo accerta la consistenza delle frazioni massiche e non scala i rate ne' corregge lo stato.
- Nessuna compilazione o Action eseguita. La verifica completa delle firme/import e dei restanti percorsi resta aperta.

### Passaggio trasversale sull'intero inventario di produzione (2026-11-01)

- Recuperati direttamente dal branch `master` e letti tutti i **64 file `.f90`** sotto `src/`, in quattro lotti da 16 file. Il censimento copre i moduli core/IO, geometria e stato, integratore Lie, registro/catalogo, operatori spaziali, mass exchange (coordinatore e sei processi), injection (coordinatore e varianti), reologie, termica, pore pressure e stocastico.
- Estratte le clausole `USE ... ONLY` mediante ricongiungimento delle continuazioni e censite le dichiarazioni `PUBLIC ::` esplicite, inclusi i riferimenti tra coordinatori e implementazioni fisiche. Esaminati i percorsi `time_integration_2d` -> `lie_source_executor_2d` -> `source_terms_2d`, la validazione dei metodi del catalogo e il confine `nonlinear_solver_2d` -> `eval_source_implicit`. Non e' emersa in questa scansione una ulteriore correzione verificabile che giustificasse un intervento cieco.
- Verificata la differenza fra `USE ONLY` con semplice simbolo e rinomina, come gli alias `darcy_*` in `pore_pressure_2d`: un confronto basato solo su nomi identici darebbe falsi errori. Analogamente un modulo senza lista `PUBLIC ::` puo' esportare simboli per default; le procedure type-bound e le generic interfaces devono essere risolte separatamente.
- **Limite esplicito della verifica**: questa e' una scansione integrale delle dichiarazioni `USE ONLY`, non ancora una dimostrazione meccanica che ciascun import sia accessibile e che ogni chiamata concordi per tipo, rank, optional e keyword. Non sono stati confrontati automaticamente gli alberi sintattici di tutte le definizioni/callsite; ne' eseguiti i test. D19 resta aperto per tale confronto e per il completamento dei rami semantici. Non eseguire compilazione o GitHub Actions: sara' l'utente a lanciarle manualmente dopo la chiusura statica.

### Verifica incrociata delle firme nei confini Lie/NEWTON/ANALYTIC (2026-11-01)

- Confrontati direttamente gli header delle procedure nei moduli `nonlinear_solver_2d`, `source_terms_2d`, `rheology_2d`, `thermal_laws_2d`, `state_conversion_2d` e le rispettive chiamate dai coordinatori. `solve_backward_euler_group(source_group,q_cell,dt,Bprime_x,Bprime_y,stochastic_value,iterations_used,converged,linear_info,line_search_failed)` corrisponde alla chiamata di `lie_source_executor_2d`; `eval_source_analytic(group,q_before,context,interval,rate,status)` corrisponde alla chiamata del ramo analitico.
- L'interfaccia `qc_to_qp(qc,qp,p_dyn,context,accepted)` corrisponde alle chiamate per candidati EXPLICIT e fine-gruppo. Le due forme di `eval_source_implicit(Bprimej_x,Bprimej_y,Zij,c_qj,c_nh_term_impl,r_qj,r_nh_term_impl,group)` sono invocate con keyword `r_qj/r_nh_term_impl` nel residuo e `c_qj/c_nh_term_impl` nello Jacobiano, senza mischiare varianti reali e complesse. I dispatch delle routine specifiche `eval_rheology_implicit(models,state,use_stochastic,stochastic_value,rate)` ed `eval_thermal_implicit(models,state,radiative,convective,ambient,rate)` usano lo stesso ordine dei parametri.
- Dall'esame del solver Newton, il gruppo risolve soltanto le componenti strutturalmente attive, congelando le altre al valore di ingresso; in caso di non convergenza la routine ripristina `q_cell=q_in` prima del ritorno. Il solver delega i residui e lo Jacobiano alle due varianti senza pubblicare lo stato globale.
- **Nessun difetto nuovo confermato** in questi confini, quindi nessuna modifica artificiale al codice. Il tentativo di far passare tutti e 64 i file in un solo controllo automatico ha superato il limite di chiamate del connettore: tale esecuzione non conta come verifica riuscita. L'inventario dei 64 `USE ONLY` dei passaggi precedenti resta valido, ma la risoluzione automatica globale e l'analisi esaustiva di rank/optional per ciascuna callsite sono ancora aperte. D19 non e' chiuso.
- Per accordo con l'utente: nessuna compilazione, nessun workflow GitHub avviato; la prima compilazione sara' una Action manuale successiva alla chiusura statica.

### Audit trasversale sul commit/rollback e sui dati del gruppo Lie (2026-11-01)

- Confrontata la sequenza `simulation_2d%advance_timestep` -> `time_integration_2d%advance_lie` -> `lie_source_executor_2d%execute_lie_source_groups` usando sorgenti recuperati dal branch. La simulazione prepara i coefficienti temporali e, se abilitato, il campo stocastico, passa lo stesso `step_start` e `step_length` all'integratore, quindi avanza `runtime%t` una sola volta dopo il completamento dell'intera sweep Lie.
- Corretto in `time_integration_2d` il nome del dummy `Z`, che denotava in realta' `stochastic%effective_Z`, in `stochastic_values`. La modifica non altera l'interfaccia posizionale, i valori, la fisica o la proprieta' della topografia: esplicita il confine fisico fra campo stocastico e geometria ed evita ambiguita' nei successivi audit.
- Esaminati i rami principali di `lie_source_executor_2d`: ogni gruppo opera su copie `q_group/qp_group`; il fallimento della cella interrompe la sweep senza pubblicare il gruppo corrente; per mass exchange i candidati dell'inventario vengono accumulati per tutte le celle, e solo dopo i controlli vengono pubblicati; la topografia candidata viene proiettata e validata prima del commit degli inventari. Il chiamante ha effettuato `snapshot_substrate` e `begin_geometry_transaction`, e su fallimento richiama `restore_substrate` e `rollback_geometry_transaction`; `state%q/qp` resta invariato fino al successo globale.
- Confermata la corrispondenza degli argomenti di `assemble_lie_explicit_transaction` e del suo utilizzo nel gruppo EXPLICIT, inclusi `qp_eval` originale per gruppo 0 e `q_recipient` post-spaziale. Il metodo ANALYTIC aggiorna il candidato locale; NEWTON applica solo le equazioni strutturalmente attive e ripristina la cella all'input se non converge.
- Questo e' un controllo statico dei confini e dei rami principali, **non** la prova automatica di ogni firma, ogni path fisico e tutti gli import. Rimangono ancora aperti i controlli statici esaustivi indicati nella sezione dei criteri di chiusura. Nessuna compilazione o Action avviata, come richiesto.

### Correzione fisica trasversale: stocasticita' nel gruppo EXPLICIT (2026-11-01)

- Rilevato un difetto nel percorso `simulation_2d -> time_integration_2d -> lie_source_executor_2d -> assemble_lie_explicit_transaction -> eval_cell_source_explicit -> eval_geometry_cell_source_explicit -> eval_rheology_explicit`: il valore `stochastic%effective_Z(j,k)` veniva preparato e letto come `noise` dall'esecutore ma il dispatcher reologico EXPLICIT riceveva sempre `0.0_wp`. La reologia EXPLICIT ignorava quindi il campo stocastico anche quando abilitato, mentre ANALYTIC e NEWTON usavano il valore della cella.
- Aggiunto argomento esplicito obbligatorio `stochastic_value` alla catena dei quattro adapter sopra elencati, con inoltro dal valore `noise` della cella in esecuzione. Il dispatcher reologico ora riceve lo stesso valore stocastico previsto dai metodi ANALYTIC/NEWTON. L'assemblaggio di mass exchange rimane immutato.
- Questa modifica riguarda i sorgenti fisici, non gli operatori spaziali o i limiter. Da verificare successivamente che eventuali callsite diretti dei test usino le firme aggiornate; nessuna compilazione o Action avviata per accordo.

### Allineamento test dopo propagazione stocastica EXPLICIT (2026-11-01)

- Aggiornati i 10 callsite diretti di `eval_source_explicit` in `TESTS/TEST_CURVATURE_SOURCE/test_curvature_source.f90` e il callsite in `TESTS/TEST_SPATIAL_OPERATOR/test_spatial_operator.f90` per fornire esplicitamente `stochastic_value=0.0_wp`, evitando chiamate con la vecchia firma. Nei due test non viene esercitata la componente stocastica, perciò il valore nullo mantiene la fisica attesa.
- L'esecutore Lie passa invece il campo stocastico effettivo della cella fino al modello reologico EXPLICIT, risolvendo il bug identificato nella scansione trasversale. Non sono state effettuate compilazioni o Action.

### Sostituzione del gate legacy D17 con un gate statico Lie (2026-11-01)

- Individuato un test di accettazione non piu' utilizzabile: `TESTS/ACCEPTANCE/test_d17_physical_interfaces.py` tentava di leggere `src/numerics/step_source_executor_2d.f90`, rimosso dalla migrazione; cercava catalog entries con `SOURCE_CADENCE_STAGE/STEP`, non piu' previste, e procedure `eval_stage_source_explicit`, anch'esse eliminate. Questo test falliva indipendentemente dalla correttezza del codice e impediva una chiusura credibile della migrazione.
- Sostituito il controllo legacy con un gate statico per i 24 nomi di sorgente del catalogo corrente, verifica dei trattamenti dichiarati nei moduli reologia/mass-exchange/termica/pore pressure, condizioni di default METHOD, generic EXPLICIT, contratti di registro, ingresso della stocasticita' nel percorso EXPLICIT, snapshot/rollback e assenza dei principali identificatori STAGE/STEP di produzione.
- Nessun compilatore, processo di test o Action avviato. Il test e' stato aggiornato come sorgente ma non eseguito; eventuali problemi del gate e i restanti controlli di firma/rank vanno verificati prima della chiusura definitiva di D19.

- Verifica testuale post-commit del nuovo gate: corretto l'escaping delle espressioni regolari Python, che durante la generazione della prima versione del file aveva duplicato i backslash nelle raw string. Ora i pattern riconoscono whitespace, parentesi e nomi di procedure come previsto. Nessuna esecuzione del test.

### Audit dei gate legacy di test e manifest distribuzione (2026-11-01)

- Sostituito `TESTS/ACCEPTANCE/test_d17_state_finalization.py`: il vecchio test leggeva `step_source_executor_2d.f90`, `finalize_imex_cell_state` e `imex_RK_solver`, tutti rimossi. Il nuovo gate statico verifica la sequenza Lie di avanzamento spaziale, ammissibilita', snapshot/rollback e commit di substrato/geometria, oltre all'inizializzazione delle celle asciutte. Nessun riferimento a STEP/RK.
- Rimossi dal manifest `Makefile.am` i riferimenti alle fixture di sorgente non piu' presenti (`TEST_IMEX_STAGES`, `TEST_EQUATION_METADATA`, `probe_sources.f90`, `characterize.py`) e il relativo target `check-local` per `TEST_EQUATION_METADATA`. Non sono state avviate compilazioni ne' Action.
- I controlli aggiornati non sono ancora stati eseguiti: la verifica del loro esito verra' fatta quando l'utente avviera' la pipeline manuale. D19 non e' dichiarato chiuso, poiche' restano da riconciliare gli altri test legacy e le firme in tutti i callsite.

### Separazione definitiva della condizione al bordo dalla sweep Lie (2026-11-01)

- Verificato il percorso effettivo: `reconstruction_2d` applica `eval_boundary_injection_explicit` come condizione di Dirichlet sulle tracce delle facce; `spatial_operator_2d` passa l'override `boundary_source_enabled` a tale ricostruzione. Non si tratta di un rate volumetrico di cella. Tuttavia `source_configuration_2d` inseriva `boundary_injection` fra i nomi di sorgenti Lie, e la sweep EXPLICIT non lo applicava come sorgente volumetrico: gruppo spurio, possibile configurazione METHOD/GROUP fuorviante.
- Eliminata l'attivazione di `boundary_injection` dal registro Lie e dal dispatch per gruppi; l'iniezione laterale continua nel percorso spaziale. In `reconstruction_2d` il modulo fisico del bordo viene ora importato direttamente da `mass_boundary_injection_2d` invece di passare attraverso il re-export improprio di `source_terms_2d`; questi quattro re-export sono rimossi.
- La voce `boundary_injection` del catalogo resta descrittiva della famiglia fisica, ma non e' schedulata come gruppo Lie; i vecchi parametri `boundary_injection_method/group` non devono pilotare il termine di bordo spaziale. Eventuale rimozione dell'input legacy va coordinata con il parser di input. Nessuna compilazione o Action eseguita.

### Chiusura della separazione boundary/spatial (2026-11-01)

- Rimossi `boundary_injection_method` e `boundary_injection_group` dalla dichiarazione di configurazione in `source_input_2d`, dagli import, default e namelist `/lateral_source_parameters/` in `inpout_2d`. Non sono parametri dei sorgenti volumetrici Lie: il bordo resta una condizione spaziale.
- Rimosso l'import inutilizzato `lateral_source_flag` da `source_configuration_2d`. Si completa l'intervento iniziato nel commit `d9d9da5`: l'iniezione laterale non e' piu' schedulata tra i sorgenti Lie.
- Non modificati i test, nessuna compilazione o Action eseguita.

### Catalogo Lie limitato a sorgenti realmente schedulabili (2026-11-01)

- Consolidata la separazione tra condizione al contorno laterale e sorgenti volumetrici: rimosso `boundary_injection` da `source_catalog_2d%SOURCE_CATALOG`, ora di 23 sorgenti. Il modulo fisico `mass_boundary_injection_2d` resta utilizzato direttamente da `reconstruction_2d` per le tracce alle interfacce.
- Il registro Lie non puo' ora piu' accettare un gruppo fittizio `boundary_injection` tramite il catalogo. Eliminato anche un commento sul vecchio concetto di cadence STAGE/STEP. La configurazione del bordo rimane esclusivamente spaziale.
- Nessuna modifica ai test, nessuna compilazione o Action. I controlli di firma e i restanti percorsi fisici non sono ancora dichiarati esaustivamente conclusi.



### Rettifica: zone topografiche invarianti (2026-11-01)

- Su indicazione esplicita del responsabile scientifico, le zone topografiche `B_zone` **non devono essere modificate** dalla migrazione Lie. Ripristinato integralmente `src/geometry/geometry_2d.f90` alla versione precedente ai commit `1ad9018` e `c684437`, annullando sia il ricalcolo nelle transazioni dinamiche sia gli interventi sulle dimensioni/inizializzazione degli array della classificazione. L'architettura dei sorgenti Lie rimane separata dalla gestione originale delle zone.
- Non eseguiti test, compilazione o Actions.

### Verifica finale delle interfacce di produzione - 2026-11-01

- Controllati i 64 sorgenti Fortran presenti in `src/`: gli import `USE` sono associati a moduli del repository o a dipendenze esterne previste (`netcdf`, `omp_lib`, moduli intrinseci). Questo controllo NON equivale a una compilazione e non prova tutte le firme di procedura.
- Verificata la catena `simulation_2d` -> `time_integration_2d` -> `lie_source_executor_2d` -> `source_terms_2d` e il percorso `source_terms_2d` -> `mass_exchange_2d`, ora con `selected_names` propagato esplicitamente.
- Le sorgenti Lie sono raggruppate con un solo METHOD per gruppo; ANALYTIC ammette un solo sorgente; GROUP 0 e' EXPLICIT. `source_registry_2d` verifica anche il footprint delle equazioni (NEWTON: insieme esatto; ANALYTIC: vuoto; EXPLICIT: tutte le equazioni conservative).
- I commit `72fdfb7` e `a735293` hanno reso il coordinatore di mass exchange dipendente dalle sorgenti selezionate nel gruppo e allineato le compatibilita' di fase per perdita del carrier liquido e vaporizzazione.
- Sono invarianti `B_zone` e il calcolo delle zone topografiche. Non sono stati eseguiti test, compilazioni o workflow.
- **Stato: audit statico delle principali dipendenze e dei dispatcher svolto; D19 NON certificato concluso sul piano della correttezza Fortran/numerica.** La verifica completa delle firme in ogni callsite e la compilazione finale sono ancora non eseguite, per scelta concordata.

### Dipendenza deposizione-vaporizzazione verificata (2026-11-01)

- `mass_exchange_vaporization` usa esclusivamente la deposizione calcolata e limitata nella medesima transazione Lie. In precedenza la vaporizzazione configurata senza deposizione riceveva un rate nullo, divenendo silenziosamente inefficace.
- Validazione introdotta in `source_configuration_2d` e `source_registry_2d`: vaporizzazione richiede la presenza di deposizione; il vincolo preesistente di un solo gruppo mass exchange garantisce che i processi siano nello stesso gruppo EXPLICIT. `mass_exchange_2d` verifica anche il contratto direttamente al proprio ingresso.
- Non modificati gruppi analitici, `B_zone`, test o workflow; nessuna compilazione eseguita. La chiusura statica completa di D19 resta subordinata alla revisione finale delle interfacce non ancora accertate.

### Unicita' globale del modello reologico nella registry Lie (2026-11-01)

- Rilevata una incongruenza tra il registro, che poteva accettare due sorgenti `rheology_*`, e `group_rheology_name`, che restituisce intenzionalmente un solo modello. L'input di produzione attuale seleziona gia' una sola reologia, ma un registro costruito direttamente poteva essere accettato e poi eseguito solo parzialmente.
- `source_registry_2d%initialize` ora rifiuta piu' di una sorgente reologica complessivamente, anche se assegnate a gruppi distinti. La regola vale per i metodi EXPLICIT, NEWTON e ANALYTIC; restano obbligatori gruppi omogenei e un solo sorgente ANALYTIC per gruppo.
- Nessun cambiamento ai modelli, alle zone `B_zone`, a test, workflow o compilazione. La verifica statica delle chiamate di tutte le leggi fisiche non e' certificata esaustiva.

### Chiusura del fallback implicito fuori registro (2026-11-01)

- Rimosso in `source_terms_2d%eval_source_implicit` il percorso alternativo che, in assenza del descrittore Lie, selezionava reologia e raffreddamento da flag globali. L'argomento `group` e' ora obbligatorio e il metodo deve essere `SOURCE_METHOD_NEWTON`; la selezione delle leggi segue esclusivamente i nomi validati del gruppo.
- Verificati i due callsite di produzione in `nonlinear_solver_2d` (residuo reale e Jacobiano complex-step): entrambi gia' passano `group=source_group`. Nessuna modifica necessaria al solver.
- La ricerca globale mostra anche chiamate dirette nei test storici: il loro eventuale adeguamento e' rinviato alla fase concordata di compilazione e test, senza condizionare l'architettura finale.
- Non modificati `B_zone`, test o workflow; non eseguita compilazione. D19 resta aperto per le verifiche statiche ancora non esaustive.

## Esito consolidato dell'audit architetturale D19 (2026-11-01)

**Concluso per i percorsi di produzione verificati:**

1. `simulation_2d` -> `time_integration_2d` -> `lie_source_executor_2d`: avanzamento spaziale senza limiter, controllo di ammissibilita', gruppo zero valutato sullo stato pre-spaziale e applicato allo stato post-spaziale, successivi gruppi in sequenza.
2. Registro `source_registry_2d`: gruppo zero EXPLICIT, metodi omogenei, massimo una sorgente ANALYTIC per gruppo, footprint coerente con il metodo, una sola reologia, mass exchange nello stesso gruppo EXPLICIT.
3. `source_terms_2d` e `mass_exchange_2d`: valutazione dei soli processi selezionati, limiter accoppiati sul gruppo, dipendenza deposizione-vaporizzazione, nessun ricalcolo indipendente della deposizione nel modello della vaporizzazione.
4. Newton: residuo reale e Jacobiano complex-step ricevono il gruppo Lie validato, nessun fallback ai flag globali nel dispatcher implicito.
5. Stato e geometria: pubblicazione conservativa finale, snapshot/rollback degli inventari e dei campi geometrici derivati; zone topografiche `B_zone` non modificate dalla migrazione.
6. Ricognizione delle dipendenze `USE` dei 64 sorgenti di `src/` svolta in piu' blocchi; non individuati import verso moduli di progetto rimossi. Il confronto esaustivo di ogni `USE ONLY` e firma non e' dimostrato da questa ricognizione.

**Non attestato e quindi non segnato come completato:** confronto automatico/esaustivo di tutte le firme, rank, optional, keyword e rami fisici; compilazione; suite di test e regressione. Questi aspetti restano una fase di verifica distinta. Il completamento dell'audit architetturale sopra elencato NON certifica che ogni unita' Fortran compili, ne' che la correttezza numerica sia dimostrata.

**Decisione operativa:** nessun ulteriore cambiamento strutturale individuato nei percorsi critici verificati. La prossima fase e' la validazione di integrazione, che non e' stata avviata su richiesta del responsabile scientifico. Non usare il termine 'D19 completamente chiuso' finche' le verifiche residue richieste dal tracker non sono effettivamente eseguite.

## Verifica mirata firme dei modelli specifici (2026-11-01)

### Svolto

- Censite le dichiarazioni degli entrypoint degli 11 modelli reologici (`voellmy`, `plastic`, `lava`, `lahar`, `empirical`, `quadratic`, `newmodel`, `froude`, `lucas`, `mui`, `mui_voellmy`) e le chiamate nel dispatcher `rheology_2d`, per le quattro famiglie di interfaccia esposte quando disponibili. La ricognizione ha identificato le routine attive e la loro struttura di argomenti; non prova ancora per tutte le chiamate il rank e il tipo effettivo.
- Censite e confrontate nominalmente le sei interfacce mass exchange EXPLICIT (`erosion`, `deposition`, `carrier_liquid_loss`, `entrainment`, `vaporization`, `lava_deposition`) con il coordinatore, incluso l'argomento `deposition` passato alla vaporizzazione.
- Censite e confrontate nominalmente le dodici interfacce dei tre modelli di mass injection (bottom circular, bottom fissural e radial) con il relativo dispatcher per EXPLICIT, SEMI_IMPLICIT, IMPLICIT e ANALYTIC.
- Censite le quattro interfacce del raffreddamento (`thermal_cooling`) e le quattro del modello di pressione di poro Darcy, insieme ai corrispondenti dispatcher.
- Riletti i due callsite di produzione di `eval_source_implicit` per Newton reale e complex-step: entrambi inoltrano il gruppo validato.

### Rimane da verificare prima di certificare tutte le firme

- Confronto completo delle liste di argomenti multilinea dei dispatcher reologici con le implementazioni, inclusi tipi, `INTENT`, rank, optional e chiamate per ogni ramo raggiungibile (non limitato ai nomi dei primi argomenti).
- Confronto completo di tipi/rank/optional delle firme termiche, Darcy, mass exchange e injection gia' censite nominalmente.
- Conversioni e chiamate di geometria/stato non comprese nel censimento delle famiglie fisiche, incluse chiamate indirette e funzioni type-bound.
- Verifica automatica esaustiva di tutte le chiamate e degli `USE ONLY` nelle 64 unita' di produzione; compilazione e test restano una fase successiva concordata.

**Esito:** nessuna incongruenza certa rilevata nel solo censimento nominale di questo passaggio. Non trasformare il censimento nominale in una certificazione dei tipi/rank o dell'audit statico integrale. Nessuna compilazione, test, workflow o modifica alle zone `B_zone`.

### Verifica effettiva tipo/rank: sei processi mass exchange EXPLICIT (2026-11-01)

Confrontate le dichiarazioni degli argomenti e le chiamate in `mass_exchange_2d` con le implementazioni `mass_exchange_erosion`, `deposition`, `carrier_liquid_loss`, `entrainment`, `vaporization` e `lava_deposition`. Nessuna incongruenza tipo/rank rilevata in questi sei callsite: `physical_state_type` e `mass_exchange_rates_type` sono i tipi trasmessi; `dt` e' scalare `REAL(wp)`; `erodible` e `deposition` sono vettori `REAL(wp)` di `n_solid`; `B_zone` e' `INTEGER` scalare. Il coordinatore riceve `selected_names(:)` carattere e dispatcha solo i nomi selezionati. L'interfaccia sorgente `eval_mass_exchange_source_explicit` inoltra `B_zone`, `erodible`, `dt` e `selected_names` con i corrispondenti tipi e rank. Gli argomenti `r` dei modelli sono `INTENT(INOUT)` e vengono inizializzati dal coordinatore mediante le routine di reset.

Limite preciso: questa verifica certifica **solo i sei callsite EXPLICIT e la catena del coordinatore mass exchange**, non le interfacce storiche semi-implicite/implicite inutilizzate dal percorso Lie. Rimane da confrontare allo stesso livello ogni callsite reologico, termico, Darcy, injection e stato/geometria; nessuna compilazione o test e nessun intervento su `B_zone`.

### Controllo tipi/rank: mass injection, termica e Darcy (2026-11-01)

- Verificati i tipi e i rank delle dodici chiamate del dispatcher `mass_injection_2d` verso le implementazioni bottom circular, bottom fissural e radial nei quattro metodi. Le interfacce reali impiegano `physical_state_type` e rate `REAL(wp)(n_eqns)`; le interfacce implicite impiegano `complex_physical_state_type` e rate `COMPLEX(wp)(n_eqns)`. Tempo, fattori geometrici e aree sono scalari reali; il vettore delle frazioni di fessura ha rank uno. Non emergono incongruenze tra le dichiarazioni e le chiamate esaminate.
- Verificate le quattro firme del dispatcher `thermal_laws_2d` e delle implementazioni `thermal_cooling_2d`: stato reale per EXPLICIT, SEMI_IMPLICIT e ANALYTIC, stato complesso per IMPLICIT, rate/scalari termici reali salvo il rate IMPLICIT complesso. Non emergono incongruenze tipo/rank.
- Confrontate le quattro firme di `pore_pressure_2d` e `pore_pressure_darcy_2d`, inclusi gli array `rate(n_eqns)`, `q_before(n_eqns)`, il risultato logico `admissible` e l'opzionale `include_gas_loss` passato dal dispatcher. I rami EXPLICIT e ANALYTIC hanno stato reale, IMPLICIT complesso. Non emergono incongruenze tipo/rank nei confini esaminati.
- Ancora aperta la verifica completa delle singole chiamate reologiche e delle catene residue di stato/geometria; la verifica statica globale di tutte le unita' Fortran non e' stata certificata. Non eseguiti test, compilazione o workflow.

### Verifica sistematica delle chiamate reologiche (2026-11-01)

- Eseguito confronto strutturale delle **36 chiamate** presenti in `rheology_2d` alle dichiarazioni delle routine nei **11 modelli specifici**. Lo scanner rispetta le continuazioni di riga e le parentesi annidate degli elementi di array `process_rate(1)` e `process_rate(2)`. Risultato: **36/36 con numero di argomenti coerente**, nessuna routine di destinazione mancante nell'insieme esaminato. Il controllo iniziale con semplice split sulle virgole aveva prodotto 15 falsi positivi per gli indici `process_rate`; eliminati con il parsing delle parentesi.
- Ripartizione: 11 EXPLICIT, 7 SEMI_IMPLICIT, 8 IMPLICIT, 10 ANALYTIC. NEWTON Lie utilizza le chiamate IMPLICIT, mentre le SEMI_IMPLICIT rimangono presenti nel dispatcher ma non costituiscono un metodo di gruppo Lie.
- Resta da certificare il confronto integrale di tipi, KIND, rank, OPTIONAL e INTENT per ciascuna delle 36 chiamate, oltre alle interfacce di conversione/stato/geometria e alla verifica automatizzata globale dell'intero repository. Nessuna compilazione, test o workflow eseguito; `B_zone` invariato.

### Controllo statico tipi/rank delle 36 chiamate reologiche (2026-11-01)

- Confrontate le dichiarazioni dei dummy argument negli 11 file `rheology_*_2d` con le chiamate e i tipi locali del dispatcher `rheology_2d` per 11 EXPLICIT, 7 SEMI_IMPLICIT, 8 IMPLICIT e 10 ANALYTIC. Per EXPLICIT e ANALYTIC viene trasmesso `physical_state_type`, con rate/velocita' `REAL(wp)` di dimensione 2 e parametri reali/logici scalari. Per IMPLICIT viene trasmesso `complex_physical_state_type`; i due `rate_x/rate_y` complessi scalari ricevono `process_rate(1:2)` come elementi, con parametri stocastici reali/logici. SEMI_IMPLICIT usa stato reale e due rate reali scalari. Non sono emerse incompatibilita' di tipo e rank nei 36 collegamenti del dispatcher esaminati.
- Il modello `newmodel` e' correttamente escluso dai rami ANALYTIC e NEWTON non supportati dal catalogo; Froude, Lucas e mui non sono registrabili in NEWTON. L'interfaccia SEMI_IMPLICIT rimane nel dispatcher ma non e' un metodo Lie.
- Questo esito certifica il confronto dei confini dispatcher-modello reologico esaminati; non sostituisce la compilazione ne' la verifica automatica esaustiva dell'intero progetto, incluse funzioni type-bound, chiamate indirette e conversioni stato/geometria.

### Verifica contratti integratore, stato e geometria (2026-11-01)

- Confrontati i callsite di `time_integration_2d` e `lie_source_executor_2d` con le firme di `execute_lie_source_groups` e `qc_to_qp`: campi conservativi `REAL(wp)(n_vars,comp_cells_x,comp_cells_y)`, primitivi `REAL(wp)(n_vars+2,...)`, stato precedente, passo scalare, registry, dominio, stocastico, esito `LOGICAL` e stringa diagnostica. `qc_to_qp` riceve in ogni callsite esaminato `qc(n_vars)`, `qp(n_vars+2)`, pressione scalare e gli opzionali `context` carattere e `accepted` logico.
- Esaminati `begin_geometry_transaction`, `commit_geometry_transaction`, `rollback_geometry_transaction`: nessun argomento; snapshot e restore di `state%snapshot_substrate` / `state%restore_substrate` utilizzano `substrate_snapshot_type`; metodo type-bound `refresh_substrate_totals` senza argomenti espliciti.
- Verificate le firme di `project_cell_field_to_vertices`: input cella `REAL(wp)(comp_cells_x,comp_cells_y)` e output ai vertici `REAL(wp)(comp_interfaces_x,comp_interfaces_y)`; `assign_topography_vertices` riceve matrice reale candidata e `accepted` logico; `initialize_dry_physical_cell` riceve array reale di variabili primitive.
- Non trovate incompatibilita' di tipo/rank nei confini diretti verificati. Questa e' la chiusura della verifica mirata delle interfacce stato/geometria/Lie, **non** una certificazione automatizzata di ogni chiamata Fortran del repository.
- Restano come limiti non attestati: risoluzione globale ed esaustiva di tutti i `USE ONLY`, generic, procedure type-bound e invocazioni indirette in 64 sorgenti, e compilazione/regressione ancora non autorizzate in questa fase. Nessun test o workflow eseguito, nessuna modifica a `B_zone`.

### Censimento globale dei moduli di produzione, tutti i 64 sorgenti (2026-11-01)

- Letti tutti i 64 file `.f90` presenti sotto `src/` tramite 8 blocchi separati per rispettare i limiti del connettore GitHub. In ciascun file censite le dichiarazioni `MODULE` e le istruzioni `USE`. I moduli di progetto importati risultano rappresentati nei sorgenti della stessa revisione, inclusi i moduli state, geometry, spatial, physics, numerics, io, stochastic e i singoli modelli. Distinte dipendenze esterne/intrinseche (`iso_fortran_env`, `ieee_arithmetic`, `ieee_exceptions`, `iso_c_binding`, `omp_lib`, `netcdf`).
- **Chiuso: controllo globale dell'esistenza dei moduli importati.** Il controllo e' stato effettuato su tutti i 64 file, non su un campione. Il tentativo iniziale in una singola chiamata ha superato il limite del connettore ed e' stato sostituito dalla lettura in blocchi.
- **Non ancora chiuso: confronto esaustivo di ogni simbolo importato con `USE ONLY`, accessibilita' pubblica/privata, renaming, generic, type-bound e di ogni chiamata e funzione indiretta.** Questa distinzione non puo' essere convertita in un esito positivo senza un controllo specifico; i confini Lie/modelli principali hanno gia' verifiche mirate nelle sezioni precedenti.
- Non eseguite compilazioni, test o workflow, e nessuna modifica a `B_zone`.

### Risoluzione mirata USE ONLY nei moduli Lie centrali (2026-11-01)

- Analizzati `lie_source_executor_2d`, `time_integration_2d`, `source_terms_2d`, `source_registry_2d`, `geometry_2d`, `state_2d`, `state_conversion_2d`, `nonlinear_solver_2d`. Confrontati **67 riferimenti di import `USE ONLY` interni a questo insieme** con i testi dei moduli fornitori della stessa revisione: nessun nome non trovato. La ricerca per presenza testuale del simbolo NON verifica `PRIVATE`/`PUBLIC`, rinomine complesse e generic interface: risultato mirato, non verifica semantica del compilatore.
- Rimane aperta la certificazione esaustiva dei simboli `USE ONLY` dell'intero repository, delle chiamate indirette e dell'accessibilita' pubblica, in particolare nei moduli non compresi nell'insieme sopra. La compilazione e i test restano differiti secondo accordo.
- Nessuna modifica a `B_zone`, test, workflow o al codice di produzione.

### Verifica accessibilita' degli import centrali aggiuntivi (2026-11-01)

- Censite **80 voci `USE ONLY`** nei quattro chiamanti `simulation_2d`, `source_configuration_2d`, `source_input_2d` e `time_integration_2d` verso 18 moduli fornitori distinti. Questo numero riguarda le voci importate nei quattro file, non una verifica esaustiva dell'intero repository.
- Controllate le dichiarazioni di accessibilita' dei fornitori `parameters_2d`, `geometry_2d`, `state_2d`, `source_registry_2d`, `mass_exchange_lava_deposition_2d`, `source_input_2d`, `model_layout_2d` e `runtime_2d`: moduli con `PRIVATE` generale esportano esplicitamente tipi e procedure richiesti; i sette parametri pubblici della deposizione lava importati dall'input sono dichiarati con attributo `PUBLIC`, incluso `lava_deposition_flag`. Nessun errore confermato in questi contratti.
- Da non confondere con una risoluzione semantica esaustiva: restano i simboli dei fornitori non inclusi nel controllo di accessibilita', le generic/type-bound indirette e tutte le chiamate non coperte dagli audit mirati. Nessuna compilazione, test o modifica alle zone `B_zone`.

### Visibilita' effettiva dei simboli centrali (2026-11-01)

- Controllate le liste esplicite `PUBLIC ::` dei moduli `lie_source_executor_2d`, `source_terms_2d`, `rheology_2d`, `nonlinear_solver_2d`, `source_registry_2d`, `state_2d`, `state_conversion_2d`, `thermal_laws_2d`: questi moduli usano `PRIVATE` generale e pubblicano i dispatcher e gli entrypoint Lie richiesti.
- Incrociati gli import `USE ONLY` di produzione tra `source_terms`, `rheology`, `nonlinear_solver` ed `lie_source_executor`: `eval_source_implicit`, i dispatcher reologici, `assemble_lie_explicit_transaction`, `eval_source_analytic`, `cell_analytic_source_context`, `analytic_term_context_type`, `solve_backward_euler_group`, `qc_to_qp`, `initialize_dry_physical_cell`, `conservative_mass_admissible` e `physical_cell_state_finite` risultano espressamente pubblici.
- Controllato anche che `physical_state_type`, `complex_physical_state_type` in `state_conversion_2d` e `state_type` in `state_2d` siano dichiarati con attributo `PUBLIC` direttamente sulla definizione del tipo, evitando falsi errori dovuti alla loro assenza dalle liste `PUBLIC ::`.
- Le type-bound procedure `initialize`, `finalize`, `refresh_substrate_totals`, `snapshot_substrate` e `restore_substrate` sono dichiarate su `state_type`; `groups`, `initialize`, `finalize`, `is_initialized` e `has_implicit_sources` sono dichiarate su `source_registry_type`.
- Nessuna incoerenza di accessibilita' confermata in questi collegamenti. Restano fuori da questo controllo i simboli e le chiamate degli altri moduli del repository e la validazione semantica con il compilatore, non ancora autorizzata. Nessuna modifica ai sorgenti, test, workflow o `B_zone`.

### Audit import dei moduli spaziali, stocastici e dispatcher fisici (2026-11-01)

- Censite le istruzioni `USE ONLY` in `mass_exchange_2d`, `mass_injection_2d`, `pore_pressure_2d`, `thermal_laws_2d`, `spatial_operator_2d`, `reconstruction_2d`, `hyperbolic_2d` e `stochastic_module`.
- Verificata la visibilita' dichiarata dei fornitori `mass_boundary_injection_2d`, `equation_terms_2d`, `nonconservative_2d`, `pccu_2d`, `reconstruction_2d`, `hyperbolic_2d` e `stochastic_random_2d` (moduli con `PRIVATE` generale): gli entrypoint richiesti sono dichiarati `PUBLIC` e i tipi workspace sono `TYPE, PUBLIC`; le costanti `PATH_DIR_X`, `PATH_DIR_Y` e `pccu_speed_tolerance` sono dichiarate `PUBLIC` direttamente.
- Le rinomine `darcy_explicit`, `darcy_semi_implicit`, `darcy_implicit`, `darcy_analytic` in `pore_pressure_2d` si riferiscono ai quattro simboli pubblici corrispondenti di `pore_pressure_darcy_2d` (gia' censiti nel controllo delle interfacce Darcy).
- Nessuna incompatibilita' di visibilita' confermata nei collegamenti esaminati. Non dichiarare completato l'audit semantico di tutti gli `USE ONLY` o di tutte le chiamate indirette: restano le verifiche esaustive e, nella fase successiva autorizzata dall'utente, la compilazione. Nessun intervento su sorgenti di produzione, `B_zone`, test o workflow.

### Correzione concreta: dispatcher sorgenti Lie indipendente da `rheology_flag` (2026-11-01)

Durante l'audit residuo sono state identificate tre verifiche ridondanti di `rheology_flag` in `source_terms_2d`: nei percorsi ANALYTIC, NEWTON e EXPLICIT impedivano l'esecuzione di una sorgente reologica gia' selezionata dal gruppo Lie. Particolarmente critico il ramo ANALYTIC, che restituiva `status=0` anche se il flag globale aveva impedito ogni aggiornamento. Sono state eliminate queste tre condizioni: l'abilitazione della reologia e' responsabilita' di `source_configuration_2d`, la selezione per ciascun passo di `group%source_names` e dei metodi validati dal registry. Restano i controlli fisici (per es. cella asciutta) e le limitazioni numeriche pertinenti. La modifica non tocca `B_zone`.

**Stato:** bug di selezione corretto nel codice di produzione; audit statico dei collegamenti principali eseguito, ma non certificato come equivalente a compilazione o controllo semantico esaustivo. Nessuna compilazione, test o workflow eseguito.

### Allineamento ulteriore dei dispatcher ANALYTIC ed EXPLICIT (2026-11-01)

Dopo la correzione dei gate `rheology_flag`, eliminati anche i doppi gate `thermal_cooling_flag` nei rami ANALYTIC ed EXPLICIT e il controllo `pore_pressure_flag` nella selezione ANALYTIC. La scelta dei processi viene effettuata in `source_configuration_2d` e validata per gruppo nel registry; i dispatcher invocati con il gruppo eseguono il processo selezionato, senza annullarlo silenziosamente a seguito di modifiche a flag globali. I flag che parametrizzano la fisica rimangono invariati. Nessuna compilazione o test eseguito.
