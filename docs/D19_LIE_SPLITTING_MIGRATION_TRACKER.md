# D19 — Architettura Lie definitiva e stato della migrazione

**Branch di lavoro:** `master`  
**Stato:** migrazione del codice di produzione in corso; architettura finale non ancora completamente implementata.  
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
- [ ] **Verificare la separazione di responsabilità rimasta**: revisione delle API di proiezione del fondo, delle variabili evolutive ancora esposte da `geometry_2d` e delle eventuali dipendenze nascoste tra gruppi. Eliminati `begin_topography_rate_update`, `record_topography_rate_cell`, `project_lie_bed_rate`, il workspace globale dei rate topografici e il relativo rilascio da `simulation_2d`: i rate temporanei restano nell'esecutore Lie e la geometria fornisce la proiezione generica cella-vertice. Rimane l'audit dello storage evolutivo residuo e dei percorsi di stato. Rinominate le API di assegnazione e rollback geometrico per rimuovere il legame nominale con Lie; i test storici che ancora chiamano i vecchi nomi vanno adeguati nella fase finale. Restano nello storage geometrico le variabili fisiche di substrato, il cui passaggio a un contenitore di stato richiede audit dei chiamanti. Trasferiti tutti e sei i campi nel tipo `state_type`, compresi i due aggregati; rimossa la proprietà globale dalla geometria. Migrazione delle routine di produzione effettuata senza compatibilità delle variabili globali obsolete; restano le verifiche statiche dei chiamanti e l'adeguamento dei test.
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
