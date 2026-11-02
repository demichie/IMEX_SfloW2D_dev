# D20 — Ripristino e verifica della parallelizzazione OpenMP del passo Lie

**Repository:** `demichie/IMEX_SfloW2D_dev` · **branch:** `master`  
**Data dell'analisi:** 2026-11-01  
**Stato:** implementazione iniziale su master; i nuovi cicli paralleli richiedono compilazione, verifica delle chiamate transitive e regressioni prima dell'accettazione.  
**Dipendenza architetturale:** [D19 — Architettura Lie definitiva](D19_LIE_SPLITTING_MIGRATION_TRACKER.md).

## 1. Obiettivo e vincoli non negoziabili

Ripristinare il parallelismo **per cella** nel passo Lie, mantenendo l'ordine sequenziale dei gruppi e la semantica delle transazioni. La versione finale deve fornire risultati indipendenti dal numero di thread, sia in caso di accettazione sia in caso di rifiuto dello stato.

- L'ordine temporale rimane: snapshot iniziale, valutazione spaziale su `q^n`, avanzamento conservativo spaziale, verifica di ammissibilità **senza limiter spaziali**, gruppi sorgente Lie ordinati, commit completo o rollback dell'intero passo.
- Il **gruppo 0** calcola i rate a partire dallo snapshot pre-spaziale, ma ne limita/applica gli incrementi al candidato post-spaziale. I gruppi successivi utilizzano lo stato e la geometria accettati dal gruppo precedente.
- L'ordine **fra gruppi** è seriale; solo le celle indipendenti **dentro un gruppo** sono parallelizzabili.
- Ogni modello fisico conserva le **quattro interfacce** EXPLICIT, SEMI_IMPLICIT, IMPLICIT, ANALYTIC; le interfacce non supportate continuano a generare un errore esplicito. Non aggiungere wrapper legacy o routine che inoltrino semplicemente le chiamate.
- Le funzioni in `physics/` ricevono **valori fisici locali**, non indici `j,k`, né riferimenti a celle o dominio. L'estrazione dei dati di cella resta in `numerics/`.
- Lo stato possiede depositi, erosione, erodibile, lava e aggregati; la geometria possiede topografia e grandezze derivate. Nessun modulo fisico deve scrivere in questi array.
- Le modifiche riguardano parallelizzazione e correttezza concorrenziale, **non** formule, limiter, scelta dei modelli, disponibilità dei metodi o ordine degli operatori.

## 2. Stato del codice: riscontri dalla lettura di `master`

| Area | Riscontro statico | Conseguenza |
| --- | --- | --- |
| `numerics/time_integration_2d.f90`, `advance_lie` | I due cicli su `domain%reconstruction_cells` e `domain%solve_cells` non hanno direttive OpenMP. Le variabili `j,k,pd,conversion_ok` sono locali alla routine ma attualmente condivise se si aggiunge un semplice `PARALLEL DO`. | Parallelizzare con variabili private e gestione appropriata degli errori di conversione. |
| `spatial/spatial_operator_2d.f90` | Già presenti `!$OMP PARALLEL DO` e riduzione per le velocità massime. | Conservare questo parallelismo e controllare la catena di chiamate spaziali. |
| `spatial/hyperbolic_2d.f90` | Già presenti diverse regioni `PARALLEL DO` e una regione `PARALLEL` con due `DO`. | Evitare regioni annidate accidentali, verificare `PRIVATE`, `NOWAIT` e scritture sugli scratch condivisi. |
| `numerics/lie_source_executor_2d.f90`, `execute_lie_source_groups` | Il ciclo `DO group_index` è seriale (corretto); `cell_loop: DO cell=1,domain%solve_cells` è seriale (da cambiare). Nel corpo sono presenti `EXIT cell_loop`, `CYCLE`, `reject_cell` e scritture dirette a `step_failed/failure_reason`. | Serve separare valutazione locale, staging e aggregazione degli esiti prima di applicare OpenMP. |
| `lie_source_executor_2d` | `q_group=q` e `qp_group=qp` creano buffer temporanei per il gruppo. I risultati per cella vengono scritti in `q_group(:,j,k)` e `qp_group(:,j,k)`. | Indici univoci nell'active list sono un prerequisito di race-freedom; verificare che nessuna cella compaia due volte. |
| `lie_source_executor_2d` | I trasferimenti mass-exchange vengono depositati in `staged_deposit/erosion/erodible(:,cell)`, `staged_lava(cell)`, `staged_bed_rate(cell)`. `group_j/group_k` conservano la mappatura. | Gli slot distinti per `cell` sono già adatti allo staging parallelo, purché tutte le letture dei dati iniziali restino congelate. |
| `lie_source_executor_2d` | La proiezione `project_cell_field_to_vertices`, `assign_topography_vertices`, il refresh dei totali e il commit inventari sono **dopo** il ciclo sulle celle. | Mantenere questa barriera; la proiezione/assegnazione geometrica non deve entrare nella regione parallela delle sorgenti. |
| `physics/source_terms_2d.f90` | `assemble_lie_explicit_transaction` riceve stato fisico, inventario, geometria e descrittore gruppo locali; produce rate e trasferimenti. Assembla anche un descrittore `non_exchange` temporaneo con `PACK`. | Contratto concettualmente compatibile con le celle parallele, ma audit di tutte le chiamate e delle allocazioni automatiche obbligatorio. |
| `numerics/nonlinear_solver_2d.f90` | `solve_backward_euler_group` lavora su `q_cell` e array locali, ma il modulo conserva `h` e `one_by_h` come variabili globali inizializzate da `initialize_nonlinear_solver`. | Devono essere solo lette dopo una inizializzazione seriale garantita; preferibile eliminare lo stato mutabile se non necessario. |
| `state/state_conversion_2d.f90` | Conversioni e strutture fisiche locali allocano array; alcune procedure possono terminare tramite `fatal_error`. | Verificare che non vi siano scratch `SAVE`/globali mutabili e definire quali fallimenti sono recuperabili. |
| `stochastic/stochastic_module.f90` | Possiede già regioni OpenMP e chiama il RNG in un percorso distinto. | Non produrre nuovi campioni dentro i cicli Lie: leggere il campo `stochastic_values` già preparato. |

Questi riscontri sono **statici**, non attestano ancora reentrancy di tutta la catena di chiamate, né equivalenza numerica. Non dedurre che una routine sia thread-safe dalla sola assenza di `j,k`.

## 3. Rischi tecnici da chiudere prima di aggiungere OpenMP

### 3.1 Workspace e variabili private

All'interno del ciclo delle sorgenti i seguenti dati sono **per cella** e non possono restare `SHARED`:

- coordinate `j,k`, indice `cell`, `qpentry`, `q_candidate`, `qp_candidate`, `total_rate`, `analytic_rate`, `analytic_input`, `pressure`, `noise`, `analytic_context`, `transfers`;
- `deposit_candidate`, `erosion_candidate`, `erodible_candidate`, `lava_candidate`, `ok`, `conversion_ok`, `failed`, `status`, `iterations`, `linear_info`, `context`;
- ogni nuovo buffer di diagnostica/rifiuto e gli scratch interni di Newton e delle conversioni.

Le strutture `groups`, lo stato all'ingresso del gruppo, la geometria in ingresso, la configurazione immutabile e il campo stocastico sono **SHARED in sola lettura**. `q_group`, `qp_group` e gli array staged sono SHARED soltanto per **scritture disgiunte**, da dimostrare con unicità degli indici.

Si raccomanda `!$OMP PARALLEL DO DEFAULT(NONE)` con dichiarazione esplicita di `SHARED`/`PRIVATE` (o uno schema a subroutine locale con stato automatico + `DEFAULT(NONE)`). Se si preferisce usare un tipo `cell_workspace_type`, deve essere una struttura numerica **per thread/per invocazione**, non un nuovo stato fisico trasportato in `physics/`.

### 3.2 Errori, rifiuto e determinismo

**Non** utilizzare `EXIT cell_loop`, `RETURN` dalla regione o scritture condivise su `failure_reason`. Per ogni cella raccogliere:

- `cell_failed(cell)` o un codice `cell_status(cell)`;
- identificatore del gruppo, indice della cella e dettagli diagnostici locali (es. `linear_info`, numero di iterazioni);
- eventuale testo diagnosticato, senza operazioni I/O in thread concorrenti.

Al termine del `PARALLEL DO` (barriera implicita) scegliere il **primo fallimento nell'ordine stabile dell'active list**, non il primo thread a concludere. Produrre il messaggio una volta, in seriale. Nessun rate o staging viene pubblicato se anche una sola cella fallisce.

Distinguere **errore recuperabile della cella** (candidate non finito, conversione non valida, fallimento Newton) ed **errore di configurazione/programmazione** (`fatal_error`, trattamento fisico non supportato). Le routine che possono eseguire `ERROR STOP` o `fatal_error` necessitano di un audit: dove possibile validare prima della regione parallela oppure propagarne lo stato senza interrompere un worker. Non mascherare errori fatali con valori nulli.

### 3.3 Immutabilità dello stato di ingresso e transazioni

Per il gruppo in elaborazione:

1. snapshot `q_group`/`qp_group` o buffer distinti di ingresso e uscita; leggere solo il gruppo precedente o i campi originali congelati;
2. valutare i candidati delle celle indipendentemente senza rendere visibili rate, inventari e geometria accettati;
3. controllare l'intero gruppo e raccogliere gli esiti;
4. in seriale proiettare le variazioni del fondo e verificare la geometria candidata;
5. in seriale pubblicare inventari e aggregati, quindi i campi fluidi candidati;
6. iniziare il gruppo successivo soltanto dopo il commit completo.

**Punto delicato:** il codice attuale aggiorna `q_group(:,j,k)` nel ciclo, mentre altri thread potrebbero leggere `q_group` per il calcolo. Serve verificare che ciascuna routine riceva esclusivamente la *slice locale* destinataria, e che l'active list non contenga duplicati. Se non si può dimostrare, usare buffer separati `q_group_input`/`q_group_candidate` (e analoghi per `qp`). Preferire questo secondo schema per rendere l'invariante verificabile e ridurre accoppiamenti accidentali.

La proiezione cella-vertice **non** è indipendente per cella: è esclusa dalla prima fase di parallelizzazione. L'ordine delle somme o la modalità di proiezione non deve cambiare involontariamente.

### 3.4 Newton e componenti condivisi

Verificare `nonlinear_solver_2d` in modo transitivo: `solve_backward_euler_group` → `solve_newton_group`, residui, `eval_jacobian`, LAPACK, `eval_source_semi_implicit` e `eval_source_implicit` → modelli reologici/termici → `state_conversion_2d`.

- Accertare che nessuna routine usi `SAVE`, allocabili globali temporanei, cache aggiornate, `RANDOM_NUMBER`, diagnostica interattiva o output condiviso.
- Inizializzare `h`/`one_by_h` prima dell'avvio dei thread; valutarne la sostituzione con valori locali/parametri per rendere il solver intrinsecamente reentrant.
- LAPACK richiamato da thread diversi deve ricevere matrici/vettori distinti; controllare anche il threading BLAS per evitare oversubscription, senza modificare la precisione o l'algoritmo.
- Un fallimento Newton deve tornare al chiamante come stato locale; non modificare le altre celle o il gruppo prima dell'esito aggregato.

### 3.5 Domini e celle attive

Verificare che `domain%j_cent/k_cent` identificano celle **uniche** e che `domain%j_reconstruction/k_reconstruction` non contengono duplicati prima di rendere parallele scritture su `q_trial`/`qp_trial`. Dove sono presenti halo, garantire che la ricostruzione non scriva in locazioni condivise con altre celle.

Il controllo di unicità può essere eseguito durante la creazione/refresh del dominio (non dentro il ciclo fisico), con un test dedicato che costruisca liste errate e verifichi il rifiuto.

## 4. Modifiche pianificate per modulo

| Priorità | File | Modifica |
| --- | --- | --- |
| P0 | `src/numerics/lie_source_executor_2d.f90` | Estrarre il calcolo di **una** cella in un worker numerico reale (non un wrapper di inoltro): input congelati, candidati, esito e diagnostica esplicita; nessuna scrittura globale. Parallelizzare il ciclo sulle celle con `DEFAULT(NONE)`. |
| P0 | `src/numerics/lie_source_executor_2d.f90` | Introdurre codici di esito per cella e aggregazione seriale deterministica; separare buffer di input/candidato se necessario; mantenere seriali proiezione, commit inventari, refresh geometrici e gruppi. |
| P0 | `src/numerics/time_integration_2d.f90` | Parallelizzare conversione pre-spaziale e aggiornamento+ammissibilità post-spaziale. Privatizzare tutti gli scratch; raccogliere errori per cella e chiamare `fatal_error` **fuori** dal `PARALLEL DO`. |
| P0 | `src/geometry/domain_2d.f90` | Verificare/garantire unicità delle liste active/reconstruction e corretta separazione delle zone di scrittura. |
| P0 | `src/numerics/nonlinear_solver_2d.f90` | Eliminare o rendere immutabili i parametri condivisi complesso-step; verificare reentrancy di Newton, residui, Jacobiani e LAPACK. |
| P1 | `src/physics/source_terms_2d.f90` | Audit transitivo dei dispatcher e degli stati allocabili locali; il rate resta solo un risultato, non viene pubblicato dal modulo fisico. |
| P1 | `src/physics/{rheology_laws,thermal_laws,mass_exchange,mass_injection_terms,pore_pressure}/` | Audit delle quattro interfacce e delle loro chiusure: niente stato temporaneo mutabile globale durante la valutazione; preservare espliciti errori di metodi non supportati. |
| P1 | `src/state/state_conversion_2d.f90` | Verificare tutte le conversioni reali/complex, incluso comportamento d'errore in OpenMP; nessun array globale scratch. |
| P1 | `src/spatial/`, `src/stochastic/` | Audit delle regioni OpenMP già presenti; evitare generazione stocastica dentro Lie e nested regions involontarie. |
| P2 | `TESTS/` e workflow CI | Regressioni seriale/parallelo, race checks e casi di rollback; test con combinazioni di metodi/gruppi. |
| P2 | `docs/D19_LIE_SPLITTING_MIGRATION_TRACKER.md` | Collegare a D20 e aggiornare lo stato dell'accettazione OpenMP al termine dei gate. |

## 5. Sequenza di implementazione: quattro fasi sostanziali

### Fase A — Invarianti e audit di thread safety

1. Mappare l'intero call graph del percorso Lie, incluse routine `physics` e `state_conversion`, e classificare ogni dato condiviso: read-only, output disgiunto, accumulazione o mutabile non sicuro.
2. Rimuovere o confinare eventuali scratch globali; rendere verificabile l'unicità delle celle nel dominio.
3. Preparare la strategia comune di diagnostica per-cell e i test del primo fallimento deterministico.
4. Fissare baseline seriale e casi di test **prima** di cambiare l'ordine di esecuzione.

**Gate A:** nessuna scrittura condivisa nascosta e nessuna condizione di errore non gestita nei worker, con elenco puntuale delle eccezioni.

### Fase B — Integratore e preparazione dello stato

1. Introdurre `PARALLEL DO DEFAULT(NONE)` sui due cicli principali di `advance_lie`.
2. Rendere privati `j,k,pd,conversion_ok` e gli altri scratch; mantenere read-only snapshot e operatore spaziale all'interno dei rispettivi confini.
3. Convertire i rifiuti a status per cella e aggregare dopo la barriera; non applicare mai limiter al rate spaziale.

**Gate B:** confronto `OMP_NUM_THREADS=1,2,4` sul passo spaziale, per tutti i campi e su celle asciutte/bagnate, senza regressioni.

### Fase C — Parallelismo dei gruppi Lie

1. Riorganizzare **in modo sostanziale** `execute_lie_source_groups`: per-cell evaluation, staging, errori locali, verifica, commit seriale. Non aggiungere `CRITICAL` attorno alla fisica: annullerebbe il parallelismo e nasconderebbe dati condivisi.
2. Abilitare EXPLICIT (inclusi gruppo 0, iniezione e mass exchange) con convalida dei trasferimenti e della topografia candidata.
3. Abilitare ANALYTIC (reologia, termico, pressione di poro) mantenendo input originali del gruppo 0 e verifica di ammissibilità.
4. Abilitare NEWTON soltanto dopo la conferma dell'assenza di scratch globali e della gestione locale delle iterazioni.
5. Mantenere la pubblicazione seriale di `state%deposit/erosion/erodible/lava`, aggregati, topografia e geometria derivata alla barriera del gruppo.

**Gate C:** identità seriale/parallelo sui candidati e sugli inventari; nessuna race; stessa scelta deterministica dell'errore.

### Fase D — Accettazione finale e prestazioni

1. Eseguire checker statici, build OpenMP, regressioni Lie/HP-PCCU/stochastic/restart e varianti della fisica.
2. Ripetere con `OMP_NUM_THREADS=1,2,4` e, ove disponibile, 8; eseguire almeno una configurazione con `OMP_DYNAMIC=FALSE`. Provare `OMP_SCHEDULE=static` e `dynamic` senza influenzare i risultati numerici.
3. Confrontare stato conservativo, primitivo, massa, energia, depositi, erosione, erodibile, lava, quota del fondo, grandezze geometriche, tempi/output e diagnostica dei rifiuti; usare confronto bitwise dove l'ordine aritmetico è invariato, altrimenti tolleranze esplicite e giustificate.
4. Testare rifiuti artificiali di più celle simultanee, rollback dell'intero passo, gruppo 0 e gruppi successivi, casi dry/wet, topografia attiva/inattiva e solver Newton che non converge.
5. Usare strumenti per race detection ove disponibili (es. ThreadSanitizer/Intel Inspector), senza considerare l'assenza di segnalazioni prova sufficiente.
6. Misurare speedup separato di spaziale, conversioni e gruppi EXPLICIT/ANALYTIC/NEWTON, evitando falsi confronti dovuti a BLAS multi-thread e allocazioni per cella.

**Gate finale:** correttezza scientifica preservata, race freedom verificata e parallelizzazione significativa delle celle. La prestazione non giustifica mai deroghe all'atomicità del passo.

## 6. Matrice minima di test

| Scenario | Proprietà da verificare |
| --- | --- |
| Solo operatore spaziale + gruppo 0 esplicito | Legge del gruppo 0 e nessun limiter spaziale. |
| EXPLICIT con scambio di massa, topografia attiva | Budget di fase, staging disgiunto, proiezione e geometria seriali. |
| EXPLICIT con iniezioni e celle asciutte | Disponibilità e convertibilità, assenza di race nelle sorgenti. |
| ANALYTIC reologia | Valutazione su stato e geometria congelati del gruppo. |
| ANALYTIC raffreddamento | Identico avanzamento termico per numero di thread diverso. |
| ANALYTIC pressione di poro/Darcy | Conservazione e ammissibilità degli stati con/senza perdita gas. |
| NEWTON misto | Equazioni attive corrette, stesso esito/iterazioni, nessun scratch condiviso. |
| Rifiuto simultaneo in più celle | Diagnosticare la prima cella *nell'ordine dell'active list*; nessun commit parziale. |
| Geometria candidata non finita | Rollback completo: stato, inventari, quota e cache geometriche. |
| Riavvio binario + evoluzione | Stessa soluzione e continuità dei campi persistenti. |
| Lista active con duplicati (test negativo) | Rifiuto prima dell'esecuzione parallela. |

## 7. Condizioni di completamento e registro delle verifiche

- [x] Lettura statica dell'integratore Lie e dell'esecutore: identificati i cicli seriali e i punti di commit.
- [x] Individuati i rischi: variabili locali condivise, errori con `EXIT`, parametri Newton globali, output disgiunti, topografia.
- [ ] Audit transitivo e matrice completa delle chiamate thread-safe.
- [ ] Implementazione dello staging e diagnostica deterministica per cella.
- [x] OpenMP attivato nei due cicli di `advance_lie` con raccolta seriale dell'esito. **Da compilare e validare**.
- [x] Prima implementazione OpenMP nel ciclo celle di `execute_lie_source_groups` per EXPLICIT, ANALYTIC e NEWTON, con risultati per cella e scelta seriale deterministica dell'errore. **Da compilare e verificare per thread safety transitiva**.
- [x] Verificata staticamente la costruzione ordinaria delle liste di celle mediante scansione delle maschere in `domain_2d` (una visita per coppia di indici). **Resta da testare** l'invariante e l'eventuale modifica esterna delle liste.
- [ ] Confronti di correttezza `1/2/4/8` thread e test di rollback.
- [ ] Analisi race e validazione prestazionale.
- [ ] Integrazione finale con D19 e chiusura del piano.

## 8. Registro implementazione (2026-11-02)

- `time_integration_2d.f90`: `PARALLEL DO` su ricostruzione pre-spaziale e avanzamento spaziale, scratch `PRIVATE`, codici di errore per cella, controllo seriale dopo la barriera. Commit `10e224f`.
- `lie_source_executor_2d.f90`: `PARALLEL DO` nel gruppo Lie, temporanei `PRIVATE`, diagnostiche per slot `cell_failure(cell)` e selezione deterministica della prima cella fallita; rimangono seriali gruppi, proiezione, aggiornamento degli inventari e refresh. Commit `fad351b`, correzione `f6d47e8`.
- **Non certificato**: reentrancy dell'intera catena di chiamate (conversioni, fisica, Newton), univocità effettiva delle scritture in ogni configurazione, compilazione con `-fopenmp`, test 1/2/4/8 thread, casi di rollback e race detection. Questi restano gate obbligatori prima di considerare attiva e sicura la parallelizzazione.

- Audit Newton: le scritture verbose e la pausa interattiva della line-search sono state disabilitate quando la chiamata avviene dentro una regione OpenMP, mantenendo il comportamento diagnostico nell'esecuzione seriale. Le routine non devono usare I/O interattivo concorrente. **Resta da verificare** la sicurezza di ogni chiusura transitiva e delle interfacce del convertitore. 

- Audit statico aggiuntivo: in `mass_exchange_2d`, `rheology_2d`, `thermal_laws_2d` e `pore_pressure_2d` i dispatcher possono raggiungere `fatal_error` per configurazioni non supportate. L'esecutore valida ora la presenza della lista e il metodo del gruppo **prima** del ciclo OpenMP, ma le altre combinazioni devono essere garantite dal registry/catalogo e coperte dai test; non si assume che tutte le chiusure siano automaticamente thread-safe. `qc_to_qp` restituisce uno status opzionale, ma percorre anche `enforce_mass_fraction_closure` con diagnostica fatale: questa dipendenza richiede ancora un controllo dei limiti di admissibilità in tutti i casi.

- Conversione stato, commit `2df999f`: `enforce_mass_fraction_closure` accetta ora lo stato opzionale `accepted` per segnalare frazioni materialmente negative senza invocare un errore fatale. `qc_to_qp`, quando richiesto l'esito `accepted`, utilizza questo percorso recuperabile; i chiamanti tradizionali senza stato mantengono la diagnostica fatale. La verifica dei dati prima di `r_phys_var` e la reentrancy dei percorsi transitivi restano soggette a compilazione e test.

- **Newton reentrant (2026-11-02):** eliminati i due scalari globali mutabili di perturbazione complex-step; `eval_jacobian` determina localmente `complex_step_perturbation` e il reciproco dalla configurazione `n_vars`, preservando la formula originale. Rimossi i corrispondenti entrypoint di inizializzazione/finalizzazione ormai superflui e le chiamate in `simulation_2d`. Nessuna modifica alle quattro interfacce dei termini fisici. **Validazione compilatore e confronto numerico ancora da eseguire.**

**Nota di metodo:** questo piano deriva dalla lettura del codice di produzione alla data sopra indicata. Le proprietà di thread safety transitive sono *ipotesi da validare*, non garanzie; non è stata effettuata compilazione né esecuzione dei test come parte della redazione.
