# IMEX_SfloW2D — Roadmap di riprogettazione del sistema degli stati

**Specifica tecnica operativa** · 5 novembre 2026 · Repository `demichie/IMEX_SfloW2D_dev` · Branch `master`.

> **Ambito:** progettazione clean-slate della gestione e implementazione degli stati e del loro utilizzo in numerica e fisica. Non è una riprogettazione dell’intero codice: leggi fisiche, discretizzazioni, integratori e consumatori vengono modificati solo nella misura necessaria per usare i nuovi contratti di stato. La destinazione deve essere ricavata dai flussi di dati, dai costi computazionali e dagli invarianti numerici, non dalla suddivisione attuale dei file o dalla compatibilità con organizzazioni storiche. Il codice esistente è una fonte di formule, invarianti e casi da preservare, non il vincolo architetturale di partenza. Questo documento è una roadmap tecnica specialistica subordinata alla [roadmap generale](IMEX_SfloW2D_Development_Roadmap.md): non introduce un piano generale concorrente. Aggiornare checklist, prove e rischi qui; riportare nella roadmap generale solo stato aggregato e priorità.

> **Stato corrente al 5 novembre 2026 (baseline architetturale `c5d6eb2`):** il cutover del codice di produzione e completato: unica directory `src/state/`, moduli `state_types`, `state_fields`, `state_composition`, `state_mixture`, `state_conversions`, `state_physical` e facciata `state`; rimossi `state_2d`, `state_conversion_2d` e `src/state_new/`. I consumatori di produzione e il Makefile usano i nuovi moduli. L'utente ha confermato una Action **Build** verde dopo il cutover del 5 novembre; verificare lo SHA del relativo run prima di attribuirgli anche la successiva rinomina. **Traguardo architetturale raggiunto; traguardo computazionale NON ancora certificato:** esistono riusi mirati e kernel minimi, ma non e stato completato l'audit di tutte le conversioni, dei calcoli derivati e delle allocazioni nei percorsi caldi. Il test scientifico, il profiling e le regressioni restano gate distinti. Tutte le Actions partono solo manualmente.

**F6/F7 primo blocco di migrazione completa dei chiamanti (5 novembre 2026, CODICE PRONTO / CI DA VALIDARE):** migrati da `state_conversion_2d` a `stateNEW` `time_integration_2d` (conversioni pre-spatial e post-spatial), `lie_source_executor_2d` (conversioni dei candidati espliciti/Newton e proprietà fisiche per buoyancy), `init_2d` e tutti i punti di `inpout_2d` ancora dipendenti dalle vecchie conversioni. Le conversioni inverse frazioni volumetriche -> frazioni massiche e primitive -> frazioni volumetriche sono ora kernel canonici nel modulo `state_new_mixture`, mantenendo le formule storiche e la normalizzazione dei tracciati. Gli stati fisici di input e output sono predisposti una volta per routine, non allocati per ogni cella. Nei due operatori Lie si usa `conservative_to_primitive(numerical_reconstruction=.TRUE.)` per la semantica near-dry originale. Nessuna modifica agli inventari di settling, rate, limiter o adozione di legacy wrapper. Aggiornato `src/Makefile.am` con le dipendenze canoniche; esteso il test di equivalenza sui percorsi I/O e collasso. **Build, CI e confronto output restano gate aperti**; la vecchia implementazione `state_conversion_2d` resta temporaneamente per test e altri chiamanti non inclusi in questo blocco. La riorganizzazione dell'ownership `state_2d` e il cleanup finale del backend sono esclusi da questo change set e appartengono al secondo blocco.

**F5 correzione del contratto di validazione (5 novembre 2026):** dopo il fallimento della regressione `TEST_MASS_EXCHANGE/test_projection_stress` (fatal error `Invalid mass exchange physical state`), annullati i controlli di ammissibilita introdotti in `source_terms_2d` nei percorsi explicit, mass exchange, analytic, semi-implicit e Newton, e annullata la conseguente propagazione di `accepted` nel risolutore Newton. Ripristinati inoltre i comportamenti originari dei costruttori `state_new_physical` e della chiusura composizionale `state_new_composition`, inclusi i test associati al cambiamento poi annullato. La verifica dell'ammissibilita del candidato globale deve essere governata dall'esecutore Lie / driver al confine del time-step, senza fatal error aggiunti nei singoli evaluator fisici. Restano i controlli di dominio intrinseci alle singole leggi (es. argomenti fisicamente impossibili) e i controlli esistenti del limitatore e del rollback, che non sono stati modificati in questo intervento. I commit precedenti su queste verifiche sono storici e **sostituiti** da questa decisione. La verifica CI di questa correzione e ancora da eseguire; F5 resta IN CORSO.

**F5 - contratto transazionale globale Lie (5 novembre 2026; chiarimento definitivo).** Il rollback appartiene esclusivamente al driver del time-step ed e **globale**: nessuna cella e nessun gruppo Lie viene ripristinato o accettato separatamente. Lo stato autorevole del dominio prima del time-step rimane immutato durante la costruzione del candidato (flusso, substrato, geometria); i campi candidati sono pubblicati solo dopo la verifica finale di tutto il dominio. Se il time-step fallisce, si scarta l'intero candidato e si ripristinano insieme tutti i campi eventualmente mutati. Gli scratch locali per worker e le copie locali di cella sono dettagli di implementazione, non transazioni con commit o rollback autonomi. Eliminare i controlli di ammissibilita intermedi aggiunti a spatial/source-group e ricondurre la verifica dell'ammissibilita fisica al confine finale del time-step; non confondere gli errori intrinseci di una legge o l'esito di un solver con un'accettazione anticipata. **Audit corrente:** `time_integration_2d` controlla ancora la massa dopo il solo spazio (righe 170-175) e `lie_source_executor_2d` controlla candidati tra gruppi (righe 460 e 570); queste verifiche anticipate devono essere riorganizzate senza perdere la possibilita di valutare in sicurezza i gruppi successivi. Questa e una prescrizione per il prossimo change set, non una modifica al comportamento gia effettuata. Prove richieste: nessuna pubblicazione parziale quando fallisce una qualsiasi cella; integrita del dominio, degli inventari e della geometria; singola validazione finale con OpenMP 1/2/4 thread. F5 IN CORSO.

**F5 analytic dispatcher deduplication (5 novembre 2026):** `source_terms_2d` costruisce lo stato fisico una sola volta per l'unico modello analitico selezionato, invece di ripetere il codice in tre rami. Le reologie e la termica mantengono gli stessi dati geometrici; Darcy mantiene il no-op se la massa e zero. Non sono stati aggiunti controlli `accepted`, validazioni intermedie o cambiamenti al rollback Lie; eliminata la qualificazione `RECURSIVE` non necessaria. Compilazione e test regressivi del nuovo commit ancora da verificare.

**F5 regressione descrittori granulometrici (5 novembre 2026):** `TEST_STATE_NEW/test_state_new.f90` esercita Sauter e densita media con frazioni diverse da zero, esattamente nulle e inferiori alla soglia `EPSILON`, controllando i risultati attesi delle formule storiche dopo la migrazione. La nuova Action e ancora necessaria per convalidare questo test. Non sono state alterate formulazioni o logiche di accettazione Lie.

**Gate aggiuntivo F6/F7 (5 novembre 2026):** la workflow `State architecture validation` include ora un job `strict-debug-driver`, con `autoreconf`, `configure --enable-strict-debug` e `make -j2` dell'intero eseguibile. E necessario per verificare i cambiamenti di `IMEX_SfloW2D.f90` e delle dipendenze di build non coperti dai test unitari; il nuovo job non e ancora stato eseguito sul commit finale.

**F5 — ripristino semantica numerica dei budget di deposizione Lie (5 novembre 2026, DA VALIDARE):** i dati `D19_SOURCE_AUDIT` dell'utente mostrano una sorgente radiale correttamente calcolata su 31 celle, ma il settling consuma praticamente tutte le particelle appena iniettate e il sink finale `red_grav<=0` elimina la corrente. Il passaggio `ea213c5` aveva sostituito il budget conservativo del gruppo 4 (stato pre-spaziale `q`) con lo stato destinatario `q_local` dopo l'iniezione; questo ha autorizzato per la prima volta la deposizione *nello stesso passo* dei solidi nuovi, modificando lo splitting storico. La soluzione conserva entrambe le condizioni: i rate negativi possono prelevare soltanto inventari presenti nello stato pre-spaziale **e** in quello destinatario aggiornato; la validazione del candidato resta effettuata **esclusivamente sul destinatario** `q_local`. Il contratto fisico riceve un secondo vettore locale opzionale `q_preexisting` e non indici di cella; il vecchio `state_conversion_2d` non viene ripristinato. In una cella inizialmente asciutta, la deposizione nello stesso passo viene annullata come nella semantica precedente, pur mantenendo l'ammissibilita esatta dello stato destinatario. Aggiunto `check_preexisting_settling_budget` in `TEST_MASS_EXCHANGE`. Il nuovo limiter contiene anche la correzione indipendente del margine floating-point per il D19 sul bordo del solido. **Non dichiarare risolta la regressione originale post-F5 senza CI e confronto TAAL con settling**.

**Gate riproducibile TAAL settling (5 novembre 2026):** le modifiche al kernel di composizione non hanno risolto il rollback, che resta riproducibile a `t=0.0226` esclusivamente con `SETTLING_FLAG=T`. Per evitare ulteriori correzioni speculative la CI `strict-debug-driver` esegue ora uno smoke test reale di `EXAMPLE_TAAL` su DEM originale e input originale con solo `T_END=0.05`, mantenendo settling e gli altri gruppi fisici. Il test raccoglie un log diagnostico e fallisce in presenza di D19. L'esecutore Lie, solo quando `IMEX_LIE_DIAGNOSTICS=1`, registra per la prima cella fallita gruppo/cella, stato conservativo destinatario, rate sorgente gia scalato, candidato e verifica della ricostruzione `q+dt*rate`, **dopo** il barrier OpenMP e senza stampe per cella nel ciclo. Questi dati consentono di discriminare mismatch limiter/esecutore, mutazioni dello stato e violazioni della chiusura. Non modificati fisica, soglie, numerica Lie o rollback. **Risultato smoke test ancora da vedere; F5 non si considera validata su TAAL**.

**Audit differenziale F5 del settling, vecchio state vs stateNEW (5 novembre 2026, DA VALIDARE):** analisi diretta del vecchio `physical_state_from_qp` alla revisione `cbcdc174ad` e del nuovo `build_physical_state_from_qp_real`. Confermata una regressione nelle tracce con componenti quasi esauriti: la vecchia `enforce_mass_fraction_closure` conserva tutte le frazioni **positive** (anche inferiori a `4096*EPSILON`) e corregge solo i valori negativi di roundoff; `close_composition_real` in `state_new_composition` azzerava anche le frazioni positive con `ABS(xs)<=tol` o `ABS(xg)<=tol`. Questo modifica `state%xs`, le frazioni volumetriche `state%alphas` e i rate di sedimentazione; non e una semplice differenza formale. Corretto il kernel REAL canonico usando `MAX(xs,0)` e `MAX(xg,0)` dopo il controllo di negativita materiale, come il costruttore storico, senza reintrodurre le vecchie API di produzione. Il kernel COMPLEX e la soglia di rifiuto negativa restano invariati. La regressione `TEST_STATE_EQUIVALENCE` verifica componenti solide/gas positive sotto la tolleranza, uguaglianza esatta di `xs/xg/alphas` vecchio/nuovo e del rate derivato dal settling, oltre alla perturbazione COMPLEX. Formule di settling, Lie, rollback, limiter e bilanci conservativi invariati in questo change set. **La root cause differenziale e dimostrata nel costruttore; la sua responsabilita esclusiva per i precedenti rollback TAAL/CF richiede verifica sul run completo con settling attivo.**

**Correzione strutturale al budget dei gruppi Lie successivi (5 novembre 2026, DA VALIDARE):** EXAMPLE_TAAL avanza fino a t=0.0226 e quindi rifiuta il candidato sorgente; disattivando `SETTLING_FLAG` avanza. Ispezionando `lie_source_executor_2d` e emerso un errore distinto dal solo arrotondamento: `assemble_lie_explicit_transaction` riceveva `q(1:n_eqns,j,k)` quale riferimento di disponibilita, sebbene il rate del gruppo successivo venga applicato a `q_local`, gia mutato dai gruppi precedenti (es. iniezione del gruppo 0). Il limiter valutava pertanto un budget *stale*. Corretto il chiamante per passare `q_local(1:n_eqns)`, preservando la valutazione del rate a `qpentry` (per il gruppo 0 e `qp_original`), tutte le formule e i trasferimenti di energia/massa. Aggiunto test deterministico `check_late_group_recipient_budget`: dimostra un caso in cui un limite riferito allo stato precedente e inammissibile sul vero destinatario e verifica che il limite sul destinatario rispetti la conservazione rigorosa. La salvaguardia ULP introdotta nella correzione precedente resta una difesa floating-point secondaria. Non cambiano l'ammissibilita rigorosa, le regole di accettazione Lie o il rollback. **Richiesti nuova Action e run reali di CF/TAAL con settling attivo per confermare la risoluzione.**

**Correzione D19 deposizione in EXAMPLE_CF e EXAMPLE_TAAL (5 novembre 2026, DA VALIDARE):** al primo time-step CF rifiutava il candidato del gruppo EXPLICIT 4 alla cella (228,101): massa totale positiva `1.4739E-04`, massa solida `-3.4660E-20`. E un rifiuto rigoroso dovuto ad arrotondamento di un bilancio limitato esattamente a zero, non un fallimento Newton. In `limit_explicit_source_rates` il rapporto analitico rimane il limite fisico, ma viene ora controllata la **stessa operazione floating-point del candidato Lie** `q_reference+dt*(limiter*rate)`. Solo se quel candidato non e ammissibile si abbassa il fattore del rate di un ULP alla volta (massimo 256 tentativi); aggiunta la protezione delle singole specie gas, oltre ai solidi e carrier. Se non esiste un fattore verificato la transazione viene rifiutata normalmente. Non si applica clipping delle masse, non si cambia la funzione di ammissibilita ne il rollback. Il predicato conservativo unico appartiene a `state_new_composition` ed e re-esportato da `state_2d` senza duplicazioni. Esteso `TEST_MASS_EXCHANGE/test_projection_stress.f90` con 1000 casi deterministici near-depletion e controllo esatto del candidato. Mancano compilazione/Action, riesecuzione di EXAMPLE_CF e EXAMPLE_TAAL e confronto bilanci.

**F6/F7 blocco di migrazione coordinata delle conversioni (5 novembre 2026):** i consumatori CFL (`spatial_operator_2d`) e diagnostici del driver (`IMEX_SfloW2D`) usano ora il kernel REAL `conservative_to_primitive` di `state_new_conversions` con `numerical_reconstruction=.TRUE.`. Il ramo standard rigoroso resta invariato; quello numerico mantiene il contratto legacy di stato secco per massa <= EPSILON, temperatura ambiente per massa <= eps_sing o temperatura non positiva, e normalizzazione delle frazioni di uscita **dopo** il calcolo delle proprieta sullo stato originale. I fallimenti di conversione vengono aggregati dai thread OpenMP e segnalati dopo la barriera. In `init_2d` lo stato vuoto inizializza direttamente `q=0`, mentre il collasso usa l'assemblaggio numerico condiviso, rinominato `numerical_primitive_to_conservative`; il percorso advettivo PCCU e i test usano la stessa interfaccia, senza alias legacy. `TEST_STATE_EQUIVALENCE` confronta la conversione numerica old/new e la pressione dinamica su casi dry, wet, near-dry, con test aggiuntivi per temperatura di fallback e soglia sub-epsilon. `time_integration_2d` conserva **deliberatamente** le attuali chiamate legacy e i controlli intermedi fino al change set separato sul rollback globale dopo F5. La nuova Action CI e il confronto scientifico multi-esempio sono ancora da eseguire su questo commit; non dichiarare F6/F7 complete.

**Gate confermato prima del nuovo change set:** la workflow `State architecture validation` e terminata `success` su `f6d6d3b2e5387038c7266899b4cf527092996559`, run `42045188411`. Include i primi change set F6/F7, l'equivalenza vecchio/nuovo e il ripristino del permesso eseguibile dello script.

**Correzione eseguibilita CI F6/F7 (5 novembre 2026):** la workflow `State architecture validation` falliva con exit code 126 prima di eseguire `TEST_STATE_EQUIVALENCE`: il relativo `run_test.sh` era registrato nel tree Git come `100644` ma chiamato direttamente. Ripristinato il modo eseguibile `100755` senza modificare il contenuto dello script; la workflow usa inoltre esplicitamente `sh`, come gli altri job, per evitare dipendenze dal modo del file. Questo errore non e una regressione fisica; la compilazione e l'equivalenza dei nuovi kernel F6 restano da verificare al prossimo run.

**F7 primo audit del grafo Make (5 novembre 2026):** rimosse in un unico blocco 26 dipendenze obsolete di compilazione verso `state_conversion_2d` dalle regole `physics/` di reologie, termica, Darcy, mass exchange, iniezioni e dispatcher, incluse due regole costituite solo dal vecchio prerequisito. Le dipendenze `state_new` dei moduli fisici restano esplicite. Eliminati i prerequisiti ridondanti senza modificare alcuna routine fisica; build parallel/CI ancora da validare. F7 e ancora parziale: restano le API legacy realmente usate nei driver, in `spatial_operator_2d`, `time_integration_2d`, `init_2d`, I/O e gli helper storici da rimuovere solo dopo migrazione e regressione dei chiamanti.

**F6 secondo change set coordinato (5 novembre 2026):** migrata anche la conversione delle tracce advettive in `hyperbolic_2d`, rimuovendo l'uso di `qp_to_qc` da quella routine. La nuova API `face_primitive_to_conservative` in `state_new_conversions` applica la soglia dry storica `h<=EPSILON` e la normalizzazione numerica di faccia, mentre la conversione di stato ordinaria mantiene la sua chiusura rigorosa. Le tre valutazioni REAL/COMPLEX/face condividono ora il *medesimo* `conservative_from_primitive.inc`: niente triplicazione della formula massa/momento/energia, nessun wrapper di inoltro. Estesa la regressione `TEST_STATE_EQUIVALENCE` per confrontare direttamente `qp_to_qc` legacy e face conversion nei casi wet gas/solid, wet liquid/solid, dry, near-dry, composizione sovrariempita e negativo entro roundoff, oltre ai test idrostatici/velocita. Aggiornate dipendenze e distribuzione di `.inc`. **CI/strict-debug e confronti HP-PCCU del nuovo blocco da verificare prima di proseguire alle conversioni del driver**. Rollback Lie non modificato.

**F6 primo change set strutturale HP-PCCU (5 novembre 2026):** `equation_terms_2d` calcola densita e gravita ridotta per il coefficiente idrostatico e per le velocita caratteristiche direttamente con `evaluate_mixture_density` da `state_new_mixture`, senza `mixt_var` e senza valutare cp, Richardson e frazioni volumetriche non richieste. La chiusura numerica delle composizioni di faccia ha ora un'unica routine `normalize_face_composition` in `state_new_composition`, utilizzata direttamente anche da `reconstruction_2d`; per preservare HP-PCCU mantiene il vecchio contratto di clipping dei negativi entro roundoff e normalizzazione di ogni overfill, distinto dalla chiusura rigorosa `close_composition` delle conversioni. Nessuna modifica a flux splitting, path, wet/dry, limiter o rollback. `TEST_STATE_EQUIVALENCE` verifica gamma, gravita e wave bounds x/y contro `mixt_var` per gas/solid, liquid/solid, dry, near-dry, overfill di faccia e negativi entro roundoff; build e CI richieste sul commit corrente. Dopo questo blocco F6 resta da migrare `qp_to_qc` nelle tracce advettive, le conversioni del driver e le ulteriori costruzioni spaziali, con test dedicati. F6 IN CORSO.

**F5 Action e baseline verificata (5 novembre 2026):** workflow `State architecture validation` verde su `9618d3294282dda5c367e0ce709311b18c138dde` (run `42038943703`). Questo certifica i nuovi kernel REAL Newton e la migrazione dei descrittori delle particelle, non il successivo change set F6 ancora da validare.

**Stato dei gate CI (5 novembre 2026):** ultima baseline verde confermata `9618d32`; i due commit F6 `99dec9a` e `46dffed` e la pulizia F7 dei prerequisiti Make richiedono un nuovo run completo prima di dichiarare chiuse le relative migrazioni. La validazione globale del rollback resta esplicitamente successiva a F5.

**F5 migrazione descrittori delle particelle (5 novembre 2026):** le formule `sauter_diameter` e `average_density_solids` sono state trasferite senza modifiche al modulo `state_new_mixture`, che ne diventa l'unico proprietario. Aggiornati i chiamanti `rheology_mui_2d`, `rheology_mui_voellmy_2d` e `inpout_2d`; rimosse le definizioni e le esportazioni dal vecchio `state_conversion_2d` senza lasciare wrapper pass-through. Aggiornate le dipendenze Fortran in `src/Makefile.am`. Nessun intervento su Lie, rollback o validazioni intermedie. Necessario verificare CI e la parita numerica dei casi granulometrici, anche con frazioni solide quasi nulle. F5 IN CORSO.

**F5 CI 344160b verificata (5 novembre 2026):** workflow `State architecture validation` run `42036539220` completata con `success` su `344160b3b65a8039304acf242312218e6763dc72`. La nuova migrazione dei descrittori delle particelle NON e coperta da questa Action precedente e richiede un nuovo run.

**F5 audit dipendenze legacy (5 novembre 2026):** rimosso l'import inutilizzato `real_part_physical_state` dal modello `rheology_plastic_2d`, completando la rimozione di questo residuo storico nel modello. Restano le funzioni fisiche `sauter_diameter` e `average_density_solids`, ancora collocate in `state_conversion_2d` e usate da mu(I), mu(I)-Voellmy e I/O; il loro trasferimento va progettato senza duplicare formule ne introdurre wrapper pass-through. Gli operatori spaziali e i controlli Lie restano fuori ambito F5 e sono stati lasciati invariati.

**F5 completamento separazione implicito REAL/COMPLEX (5 novembre 2026):** rimosso il ramo misto legacy da `eval_source_implicit`: ora la sua firma accetta soltanto `c_qj` COMPLEX e restituisce `c_nh_term_impl` COMPLEX per il Jacobiano. Il residuo usa la nuova routine `eval_source_implicit_real` senza promuovere ogni cella REAL a COMPLEX e senza estrarre `REAL(source_term)`. Mantenuto lo stesso metodo Newton, gli stessi gruppi risolti e gli stessi rate fisici per ogni legge, con formula condivisa tramite `.inc`. Non sono cambiati i controlli dello stato o il rollback. CI del nuovo commit ancora da eseguire.

**F5 Newton REAL residual / COMPLEX Jacobian (5 novembre 2026):** introdotto il dispatcher `eval_rheology_implicit_real_id` che usa i kernel impliciti REAL delle otto reologie Newton, tutti condivisi con i corrispondenti COMPLEX tramite gli include gia esistenti (incluso nuovo `rheology_newmodel_linear.inc`). Anche la termica espone `eval_thermal_implicit_real`. Il nuovo `eval_source_implicit_real` costruisce direttamente uno stato REAL e valuta il rate del residuo; `nonlinear_solver_2d` riceve lo scratch REAL gia pronto dal worker e lo passa a `eval_f` attraverso Newton/line search. Il Jacobiano mantiene il percorso COMPLEX con complex-step e scratch dedicato. I test dei kernel condivisi confrontano direttamente rate REAL e COMPLEX (incluso newmodel). Non cambia il contratto di validazione degli stati, non si aggiungono controlli di ammissibilita e non si tocca Lie/rollback. **Commit di integrazione ancora da validare in CI e regressioni dinamiche**; F5 IN CORSO.

**F5 CI verde sul change set del 5 novembre 2026:** workflow `State architecture validation`, run `42034227956`, **SUCCESS** sul commit `09c0259a9a6984543163c00783620cad9cbdbe1e`. Convalida i nuovi dispatcher analitici, il riuso dello stato fisico nei gruppi misti e i test di regressione inclusi nella matrice. Rimangono aperti l'audit delle conversioni legacy e dei chiamanti Newton, i test di equivalenza numerica estesi e la verifica del costo computazionale. Nessuna modifica al rollback Lie prima del completamento F5.

**F5 combined-group regression coverage (5 novembre 2026):** `TEST_MASS_EXCHANGE/test_projection_stress.f90` usa ora un gruppo 0 misto `mass_exchange_erosion` + `geometry_curvature`, che esercita il riuso dello scratch fra il dispatcher fisico e il mass exchange mantenendo il controllo dei budget di massa, inventari e rollback globale esistente. I limiti e la semantica Lie non sono modificati. Eseguire la Action per verificare il caso combinato.

**F5 combined explicit source scratch reuse (5 novembre 2026):** quando un gruppo Lie esplicito combina sorgenti non-exchange e mass exchange valutati sulla stessa `qp_eval`, il dispatcher riusa lo stato fisico gia costruito nello scratch del worker invece di ricostruirlo una seconda volta. Il flag `physical_ready` viene restituito solo dopo la costruzione effettiva, quindi i casi dry con ritorno anticipato non riusano stato stale; i chiamanti standalone e quelli senza scratch conservano la costruzione autonoma. Nessuna modifica ai limiter, ai controlli di ammissibilita, alla sequenza Lie o al rollback. Test/CI del change set da eseguire.

**Attivita differita - DOPO F5 (non bloccante per F5).** Riorganizzare il controllo di ammissibilita e il rollback di Lie come una sola transazione globale a fine time-step sull'intero dominio: nessun rollback di cella/gruppo, nessuna validazione anticipata dopo lo spazio o tra gruppi. Intervenire in `time_integration_2d.f90` e `lie_source_executor_2d.f90` dopo aver completato e validato la migrazione F5; aggiornare allora i test di rollback e OpenMP. **Non modificare ora questa logica come parte di F5.**


**F6/F7 — secondo blocco, confine backend storico (5 novembre 2026; VERIFICA CI MANUALE PENDENTE):** confermato che i sorgenti Fortran di produzione non importano piu direttamente `state_conversion_2d` (il modulo resta presente per i test differenziali vecchio/nuovo). Rimossa dal grafo Make la dipendenza obsoleta `spatial/hyperbolic_2d -> state_conversion_2d`, rimasta dopo la migrazione di `hyperbolic_2d` a `state_new_conversions`. Introdotto `TEST_STATE_NEW/check_legacy_boundary.py`, eseguito dalla workflow manuale `State architecture validation`, che blocca la reintroduzione di import legacy nella produzione e del prerequisito rimosso. **Non si elimina ancora `state_conversion_2d` dal target principale**: i test di caratterizzazione estraggono la lista sorgenti da `src/Makefile.am` e richiedono un percorso di build autonomo prima della sua rimozione. Ownership finale di `state_2d`, eliminazione del backend legacy dalla produzione, validazione dell'equivalenza scientifica e controllo del rollback restano aperti. Nessuna formula o semantica fisica modificata; nessuna Action avviata automaticamente.

**F6/F7 — cutover strutturale del 5 novembre 2026 (IMPLEMENTATO NEL SORGENTE, DA VALIDARE):** unica architettura di produzione sotto `src/state_new/`. Il modulo proprietario di `q`, `qp`, inventari di deposito/erosione/lava, diagnostiche e snapshot globali del substrato e stato trasferito da `src/state/state_2d.f90` a `src/state_new/state_new_fields.f90` con modulo `state_new_fields` e import dei chiamanti aggiornati. `state_new.f90` re-esporta il nuovo modulo. Rimossi dal tree `src/state/`, `state_conversion_2d` e i relativi include obsoleti; eliminate le dipendenze legacy dal target in `src/Makefile.am`. Nessuna nuova struttura alternativa, nessuna modifica intenzionale alle formule o alla transazione Lie. **Per richiesta dell'utente, test, regressioni, benchmark e Actions NON eseguiti in questo passaggio**; gli harness storici che importavano `state_2d`/`state_conversion_2d` richiedono successiva migrazione, aggiornamento o rimozione, e non costituiscono ora criteri di blocco. Prima di certificare il cutover occorre almeno una build completa manuale.

**Rinomina finale delle API di stato (5 novembre 2026; validazione manuale pendente):** la sola architettura di produzione passa da `src/state_new/` a `src/state/`. I moduli e i file `state_new*` sono rinominati `state*` (inclusi `state.f90`, `state_fields.f90`, `state_types.f90`, `state_composition.f90`, `state_mixture.f90`, `state_conversions.f90` e `state_physical.f90`); allineati i `USE`, le dipendenze in `src/Makefile.am` e i sorgenti dei test. Nessuna nuova struttura, wrapper legacy o modifica intenzionale alle formule numeriche. I riferimenti storici ai nomi precedenti in questa roadmap descrivono versioni passate; il nome canonico corrente e `state`. Actions esclusivamente manuali, non lanciate in questo intervento.

## Stato operativo corrente e piano di chiusura computazionale (sezione prevalente)

Questa sezione e la **fonte autorevole per le prossime modifiche**. I registri F0-F9, P0-P3 e i resoconti datati sottostanti descrivono l'evoluzione storica, anche quando citano moduli `state_new` o verifiche allora mancanti; non rappresentano una seconda architettura o un backlog concorrente. I nomi di riferimento attuali sono esclusivamente quelli in `src/state/`.

### A. Obiettivo e confine della chiusura

**Obiettivo computazionale:** per ogni *versione* dello stato di una cella o delle tracce ricostruite di faccia, calcolare ogni prodotto necessario una sola volta per la sua durata di validita, passando gli intermedi ai consumatori successivi senza conversioni o ricostruzioni superflue. Una nuova versione di `q` o della composizione richiede legittimamente nuovi prodotti; non si deve sopprimere una conversione necessaria a mantenere `q` e `qp` coerenti.

- Il conservativo `q` e autorevole; `qp` e una rappresentazione numerica sincronizzata. Lo stato fisico REAL/COMPLEX e un *risultato locale* da costruire quando richiesto, non una cache persistente su tutto il dominio.
- Calcolare soltanto i prodotti richiesti: conversione `q -> qp`, conversione `qp -> q`, densita/capacita termica minime di faccia o stato fisico completo per le leggi; non costruire Richardson, curvature o composizione volumetrica completa se il consumatore non li usa.
- Lo scratch vive nel worker o nella chiamata, e viene riutilizzato solo con la **stessa** versione degli input. Invalidare esplicitamente dopo modifiche di massa, specie, energia, momento, variabili trasportate, profondita o geometria da cui il prodotto dipende. La topografia puo cambiare a ogni stage: non precomputarla staticamente.
- In `physics/` restano vietati indici `j,k`, accessi a campi dell'intera griglia e stato globale mutabile. Nessun wrapper di mera compatibilita e nessuna duplicazione delle formule REAL/COMPLEX; preservare complex-step.
- Il lavoro di riuso **non** deve cambiare ordine Lie, regole wet/dry, rate o limiter fisici, bilanci conservativi, semantica di settling, transazione globale/rollback o discretizzazione HP-PCCU. Le eventuali correzioni di questi contratti sono change set separati.
- Non dedurre un guadagno prestazionale dal solo numero di moduli o dal successo della build. Distinguere audit statico, correttezza numerica e miglioramento misurato.

### B. Cosa e gia implementato — e cosa non prova ancora

| Area | Evidenza nel codice corrente | Stato computazionale |
| --- | --- | --- |
| Ownership e backend | `src/state/state_fields.f90` possiede `q`, `qp`, substrato, aggregati e snapshot; `src/state/state.f90` esporta i moduli canonici; nessuna vecchia implementazione in produzione | **ARCHITETTURA COMPLETATA**; non e una misura dei calcoli ripetuti |
| Kernel mirati | `state_conversions` distingue conversioni numeriche; `state_mixture` espone densita minima, capacita termica e valutazione completa | **IMPLEMENTATO**; verificare duplicazioni tra chiamate annidate |
| Stato locale | Tipi REAL/COMPLEX con composizioni preallocate e workspace di worker; `state_physical` espone costruttori da `q`/`qp`, refresh e aggiornamento cinematico | **IMPLEMENTATO**; durata/invalidazione degli scratch da auditare |
| Pre-spatial e CFL | `time_integration_2d:advance_lie` supporta `pre_spatial_primitives_ready` e il riuso di `state%qp` | **RIUSO IMPLEMENTATO** nel percorso che passa il flag; verificare tutti i chiamanti |
| Sorgenti esplicite | `source_terms_2d` passa `physical_ready` da un gruppo combinato all'exchange sullo stesso `qp`; singolo stato per modello analitico | **RIUSO PARZIALE**; altri gruppi e Newton da auditare |
| Conversioni Lie | `time_integration_2d` converte dopo lo spatial; `lie_source_executor_2d` converte candidati espliciti e non espliciti | **NON TUTTE RIDONDANTI**; classificazione per versione ancora aperta |
| Facce HP-PCCU | `equation_terms_2d` usa `evaluate_mixture_density`; `hyperbolic_2d` converte le tracce advettive con kernel canonico | **KERNEL MINIMI IMPLEMENTATI**; riuso tra endpoint/flux/path non dimostrato |
| Termodinamica nei costruttori | `state_physical:build_physical_state_from_q` REAL e COMPLEX chiamano prima `evaluate_mixture_heat_capacity` per ricavare T e poi `evaluate_mixture`, che include la valutazione di cp | **DUPLICAZIONE CONCRETA DA RISOLVERE**, preservando temperatura ed esatta chiusura |

### C. Regola di audit: un prodotto per versione di stato

Prima di eliminare una chiamata classificare **produttore, consumatore, identita dell'input, dipendenze, durata e invalidazione**. Per ciascun cammino annotare i prodotti `qp`, `cp_m`, `rho_m`, `alphas/alphag`, `red_grav`, `Ri`, cinematica e le conversioni/allocazioni associate.

- **Conversione necessaria:** una nuova `q_trial` dopo lo spatial o un candidato Lie accettato deve generare/aggiornare `qp_trial` prima dei consumatori che lo richiedono.
- **Conversione evitabile:** due chiamate ricostruiscono `qp` dalla *stessa* `q` immutata, senza che cambino composizione, energia o momento.
- **Riutilizzo sicuro:** stesso `q`/ `qp`, composizione/temperatura e medesima geometria rilevante; scratch gia inizializzato, input coerente, nessuna copia profonda non necessaria.
- **Riutilizzo vietato:** candidato Newton/line search differente, rate che ha mutato un componente conservativo, nuova traccia ricostruita, topografia modificata o passaggio REAL/COMPLEX che perde la perturbazione.

Una riduzione dei conteggi non e sufficiente se introduce una proprieta *stale*: ogni ottimizzazione deve rendere visibili i punti di invalidazione.

### D. Piano esecutivo — blocchi C0-C5

L'ordine e vincolante per le dipendenze; il lavoro puo essere committato in blocchi sostanziali. **Non introdurre una nuova gerarchia di moduli:** usare quelli gia presenti. Per ogni blocco registrare nel documento file modificati, SHA, risultati effettivi e lavoro residuo. Non eseguire Actions automaticamente.

| Blocco | Priorita / stato | Moduli e intervento concreto | Deliverable e criterio tecnico di chiusura |
| --- | --- | --- | --- |
| **C0 — Inventario dei percorsi caldi** | **P0 / AUDIT STATICO COMPLETATO** | `time_integration_2d`, `spatial_operator_2d`, `hyperbolic_2d`, `equation_terms_2d`, `lie_source_executor_2d`, `source_terms_2d`, `nonlinear_solver_2d`, `state_*`. Tracciare per cella/stato/faccia le conversioni, `evaluate_mixture*`, `build_physical_state*`, workspace, allocazioni e copie. Seguire i cammini CFL -> spatial -> gruppo 0 -> gruppi successivi -> Newton/analitico -> publish. | Tabella **produttore -> consumatore -> versione -> invalidazione -> numero di valutazioni**. Classificare ogni occorrenza come necessaria, riutilizzata o candidata all'eliminazione; distinguere chiamate per step, gruppo, faccia e iterazione Newton. |
| **C1 — Eliminare i ricalcoli nei kernel dello stato** | **P0 / IN CORSO** | `state_physical.f90`, `state_mixture.f90`, `state_conversions.f90`. Risolvere il doppio calcolo cp nella costruzione da `q` REAL/COMPLEX: calcolare cp per derivare T, poi passare cp gia calcolati al kernel materiale completo senza riduplicarne la formula. Condividere la singola formula via gli `.inc` esistenti. Audit di copie e array temporanei nei costruttori, conversioni e refresh; mantenere stati fisici completi e coerenti. | Nessuna rivalutazione cp per la stessa composizione nel percorso da `q`; nessun nuovo derived type o wrapper pass-through; derivazioni REAL/COMPLEX e semantica dry/near-dry/temperatura invariati. |
| **C2 — Conversioni nell'integratore Lie** | **P0 / IN CORSO** | `time_integration_2d.f90`, `lie_source_executor_2d.f90`, chiamanti CFL/driver. Verificare il flag `pre_spatial_primitives_ready` per il vero percorso produttivo e i chiamanti diretti. Mappare `q^n/qp^n`, `q_trial/qp_trial`, i candidati e il publish; riutilizzare viste gia coerenti, senza convertire due volte lo stesso candidato e senza conversione anticipata non consumata. Non rinviare una conversione se serve al gruppo successivo. | Una sola produzione di `qp` per ogni versione di `q` effettivamente consumata; nessuna lettura di primitivi obsoleti; gruppo 0 pre-spatial, altri gruppi e rollback globale preservati. |
| **C3 — Riuso del physical state nei sorgenti e in Newton** | **P0 / IN CORSO** | `source_terms_2d.f90`, `lie_source_executor_2d.f90`, `nonlinear_solver_2d.f90` e dispatcher reologici/termici. Mappare dipendenze tra explicit, mass exchange, semi-implicit, residuo REAL, Jacobiano COMPLEX, analitico. Preparare una sola volta per *identico* candidato e riusare lo scratch con identita/validita esplicita, mai tra perturbazioni complex-step o iterazioni Newton diverse. Evitare ricostruzioni, promozioni REAL->COMPLEX e allocazioni non richieste. | Numero minimo giustificato di costruzioni fisiche per gruppo e candidato; REAL residuo e COMPLEX Jacobiano restano distinti; quattro interfacce di ogni sorgente preservate e nessuna cache globale. |
| **C4 — Endpoint HP-PCCU e proprieta minime** | **P1 / IN CORSO** | `spatial/hyperbolic_2d.f90`, `spatial/reconstruction_2d.f90`, `spatial/pccu_2d.f90`, `physics/equation_terms_2d.f90`, `state_mixture.f90`. Inventariare proprieta delle tracce sinistra/destra definitive, wave bounds, idrostatica, flussi e path. Riutilizzare rho/gravita/gamma quando gli input sono identici; ricomputare dopo qualunque ricostruzione/normalizzazione che cambi la composizione o profondita. Evitare l'intero physical state sulle facce. | Ogni proprieta di endpoint valutata il numero minimo necessario per la traccia finale; nessuna ipotesi di topografia statica o modifica a well-balancing/wet-dry/path. |
| **C5 — Audit memoria, evidenze e chiusura** | **P1 / IN CORSO** | Tutti i percorsi C1-C4. Verificare setup per-worker, array automatici, allocazioni nascoste, intrinsic assignment con `ALLOCATABLE` e conformita OpenMP. Raccogliere contatori per tipo di chiamata, poi profili/timing separati per spaziale, sorgenti, Newton e conversioni in `EXAMPLE_CF` e `EXAMPLE_TAAL` (e almeno un caso multi-fisica). | Tabella prima/dopo di conteggi e costi con SHA e input; regressioni fisiche e strict-debug manuali quando richieste; nessuna allocazione/copia evitabile nei loop caldi e nessun peggioramento scientifico. Solo allora dichiarare **TRAGUARDO COMPUTAZIONALE COMPLETATO**. |

**Confine C3/C4:** C3 governa gli stati fisici completi destinati alle leggi, C4 i prodotti numerici minimi di faccia. Evitare un sistema generale di cache con version counter su tutta la griglia se bastano passaggio esplicito degli intermedi e scratch locali. Uno stesso valore scalare condivisibile non autorizza a riutilizzare uno stato fisico completo con dipendenze geometriche differenti.

### D.1 Registro operativo C0/C1 — prima implementazione (5 novembre 2026)

**C0, audit statico preliminare (IN CORSO):** inventario verificato sui sorgenti correnti; non confondere il numero dei call site con il numero di chiamate effettive per timestep.

| Prodotto / identita | Produttore e posizione | Consumatori / invalidazione | Esito |
| --- | --- | --- | --- |
| `qp^n` per CFL e Lie | `time_integration_2d:advance_lie`, ingresso `pre_spatial_primitives_ready` | Riusa `state%qp` quando CFL ha gia convertito `q^n`; altrimenti converte in `qp_before`. Si invalida dopo update conservativo | Riuso gia presente; ancora da tracciare tutti i chiamanti |
| `qp_trial` post-spatial | `time_integration_2d`, dopo `q_trial=q^n-dt*spatial_rate` | Gruppi Lie; una nuova versione di `q_trial` necessita nuova conversione | Conversione necessaria; non eliminarla arbitrariamente |
| `qp_candidate` dei gruppi Lie | `lie_source_executor_2d`, rami esplicito e non esplicito | Si aggiorna per i candidati accettati e prima dei gruppi successivi; identita diversa per gruppo/candidato | Audit del riuso tra gruppi aperto |
| Stato fisico del gruppo esplicito | `source_terms_2d`, `physical_ready` passato al mass exchange | Riuso soltanto sulla stessa `qp` e stesso contesto geometrico | Riuso parziale gia implementato |
| Stato REAL/COMPLEX di Newton | `source_terms_2d`, costruzione da `q_cell` o `c_qj` | Residuo REAL e Jacobiano COMPLEX, candidati/perturbazioni distinti | Serve inventario per iterazione e per valutazione; non condividere tra versioni diverse |
| Proprietà faccia | `hyperbolic_2d`, `equation_terms_2d` | Tracce left/right ricostruite, onde e flussi; invalidate dopo ricostruzione/normalizzazione | Kernel minimo densita gia presente, riuso per endpoint ancora da verificare |
| Capacita termica nel costruttore da `q` | `state_physical` REAL/COMPLEX: cp per ricavare T, poi valutazione completa precedente | Identica composizione, cp indipendente da T: secondo calcolo superfluo | **Corretto nel codice, verifica dinamica pendente** |
| Volume fractions da primitivi | `state_mixture:primitive_to_volume_fractions` | Prima valutava anche cp pur restituendo solo `alphas/alphag` | **Corretto nel codice, verifica dinamica pendente** |

**C1, intervento nel codice (CODICE MODIFICATO / VALIDAZIONE PENDENTE):** aggiunto in `state_mixture.f90` il prodotto distinto `evaluate_mixture_phase_properties` con implementazione REAL e COMPLEX tramite **l'unico `mixture_from_mass.inc` esistente**. Calcola densita di miscela/carrier e frazioni volumetriche, senza valutare cp. I costruttori `build_physical_state_from_q` REAL/COMPLEX calcolano una volta `evaluate_mixture_heat_capacity` (necessaria per T), quindi invocano il nuovo kernel che non ripete cp. Anche `primitive_to_volume_fractions` evita la capacita termica non richiesta. Il percorso `evaluate_mixture` completo rimane invariato per chi necessita anche cp; i costruttori da `qp` e il `refresh` non cambiano. Non e una cache, non introduce globali mutabili e non cambia lo schema REAL/COMPLEX. Restano da controllare le composizioni/temperature near-dry, il comportamento floating-point in test, le allocazioni e le altre rivalutazioni lungo i consumatori.

**Nessun test, benchmark o Action avviato in questo change set**; i conteggi e gli effetti sui tempi vanno ancora misurati. Prossimo passaggio: completare C0 per gruppi/endpoint e proseguire con C1 (audit interno conversioni, alias e copie) prima di C2-C4.

### D.2 Registro operativo C2 — riuso di qp nei candidati Lie (5 novembre 2026)

**CODICE MODIFICATO / BUILD E REGRESSIONI DEL COMMIT PENDENTI.** In `numerics/lie_source_executor_2d.f90`, dopo il controllo di ammissibilita del candidato, i rami EXPLICIT e ANALYTIC/NEWTON confrontano il vettore conservativo `q_candidate` con `q_local`. Se il candidato e **esattamente invariato**, `qp_local` e gia la vista coerente dello stesso stato: vengono saltate sia `conservative_to_primitive` sia la copia del vettore primitivo. Se cambia anche un solo componente, si mantiene la conversione numerica precedente (o la costruzione dry) e l'esito `conversion_ok` come prima. Il confronto e sul vettore completo `n_vars` e usa uguaglianza esatta: nessuna soglia o classificazione approssimata dei cambiamenti, nessuna decisione basata solo sul rate, nessuna cache persistente. Il ramo exchange continua a validare e aggiornare depositi, erosione, erodibile e lava anche quando il conservativo e invariato; il controllo di ammissibilita resta precedente al riuso. La valutazione della sorgente e il gruppo Lie non sono saltati.

**Evidenza di C0:** la CFL esegue `q -> qp` una volta sul reconstruction halo e il percorso ordinario passa `pre_spatial_primitives_ready`, mentre `advance_lie` mantiene la conversione per chiamanti diretti. La conversione post-spatial serve quando cambia `q_trial`; il prossimo controllo C2 deve documentare tutti i chiamanti del flag, il costo del confronto per gruppi sempre attivi, i casi dry/zero-rate, e la durata delle viste dopo cambiamenti di stato. Le chiamate `conservative_to_primitive` rimaste nell'esecutore sono necessarie per candidati effettivamente diversi. **Nessuna Action avviata automaticamente.**

**C2, secondo intervento (stesso blocco; compilazione da validare):** nelle conversioni `conservative_to_primitive` invocate dall'integratore Lie e dall'esecutore dei gruppi Lie e stato omesso l'argomento opzionale `p_dyn`, poiche i temporanei `pd` e `pressure` venivano calcolati e mai letti dai due chiamanti. Eliminati i temporanei e i relativi elementi delle clausole `PRIVATE` OpenMP, senza alterare l'API opzionale `p_dyn` per gli altri consumatori che la richiedono, i controlli `accepted` o la ricostruzione numerica. Questa e un'eliminazione verificata staticamente di una proprieta non consumata; non una misura delle prestazioni. Nessuna Action automatica.

### D.3 Registro operativo C3 — eliminazione di valutazioni Newton duplicate (5 novembre 2026)

**Primo intervento, CODICE MODIFICATO / VALIDAZIONE PENDENTE.** Audit di `nonlinear_solver_2d:solve_newton_group`: il residuo `right_term` e il merito `scal_f` vengono valutati alla fine di ogni iterazione, sul candidato del passo completo (prima iterazione) oppure sull'ultimo candidato accettato dalla line search (iterazioni successive). L'iterazione seguente in precedenza invocava di nuovo `eval_f` sul medesimo `qj`, ricostruendo inutilmente il physical state REAL e i rate impliciti.

Aggiunto un flag locale `residual_ready` che trasferisce **solo il residuo e il merito gia calcolati**, non un physical state persistente, all'inizio della prossima iterazione. Il primo residuo viene sempre calcolato normalmente; in caso di fallimento/early exit non si riusa un risultato di un candidato diverso. Gli update di `qj` ricostruiscono esclusivamente le stesse equazioni attive usate nella valutazione del trial. Il Jacobiano continua a chiamare l'evaluator COMPLEX per ogni perturbazione, preservando l'intera derivata complex-step; nessuna cache del Jacobiano o del suo stato. Il gruppo e il dominio numerico restano locali, senza indici di cella nella fisica. I rami semi-impliciti congelati e gli aggiornamenti conservativi non sono cambiati.

**Risparmio atteso per struttura del loop:** fino a una chiamata `eval_f` (e costruzione dello stato fisico REAL implicito) evitata per ogni iterazione Newton successiva alla prima che raggiunge il punto di riuso. Non attribuire tempi o speedup senza profiling. Restano audit del riuso per stati identici tra i dispatcher, del costo delle copie e dei fallback locali e verifica scientifica dei percorsi line search/zero-rate/complex-step.

**Situazione verificata dall'utente prima di C3:** build verde ed `EXAMPLE_TAAL` in esecuzione sulle modifiche C2. Non assumere che questo confermi il nuovo commit C3. Nessuna Action avviata automaticamente.

### D.4 Registro C3 — raffreddamento analitico: solo stato termico necessario (5 novembre 2026)

**Codice modificato, compilazione e regressioni sul nuovo commit da verificare.** Audit di `source_terms_2d:eval_source_analytic_state`: nel ramo `thermal_cooling` l'integratore esatto usa esclusivamente la temperatura iniziale e la capacita termica totale congelata; prima veniva costruito un `physical_state_real` completo (miscela, volumi, gravita ridotta, Richardson e cinematica), anche quando nessuna di queste altre proprieta veniva letta.

Ora per massa conservativa **positiva** il dispatcher calcola `xs` e `xg` da `q_before`, invoca il kernel canonico `evaluate_mixture_heat_capacity`, e ricava la temperatura con la stessa espressione `q_before(4)/(q_before(1)*cp_m)` del costruttore da `q`. La capacita termica congelata dell'integrazione analitica mantiene l'espressione precedente `q_before(4)/temperature_before`; il calcolo della temperatura e il rapporto conservativo mantengono quindi l'ordine delle operazioni. Il dispatcher termico `eval_thermal_analytic` riceve direttamente la temperatura scalare invece di uno stato fisico (nessun wrapper o modulo nuovo); l'integratore `eval_thermal_cooling_analytic` era gia scalare e non cambia. Per massa non positiva rimane la costruzione fisica storica, cosi da conservare il comportamento ambientale e di rifiuto; le altre leggi analitiche conservano lo stato fisico completo. Aggiornato `src/Makefile.am` per la nuova dipendenza diretta da `state_mixture`.

**Evitato:** un `build_physical_state_from_q` completo per ogni valutazione wet di raffreddamento analitico, mantenendo il solo calcolo di `cp_m` necessario per ricostruire T. Nessuna modifica intenzionale ai rate, alle formule, agli inventari o alle quattro interfacce dei processi; test scientifici, comparazioni floating-point e timing ancora pendenti. C3 resta **IN CORSO**: completare audit degli altri consumatori REAL/COMPLEX e degli scratch, senza condividere dati fra candidati Newton perturbati distinti.

### D.5 Registro operativo C4 — riuso della pressione idrostatica delle facce (5 novembre 2026)

**Primo intervento C4, CODICE MODIFICATO / BUILD E VERIFICA SCIENTIFICA PENDENTI.** Audit del percorso `spatial_operator_2d:compute_spatial_timestep` -> `hyperbolic_2d:eval_speeds` -> `hyperbolic_2d:eval_flux_PCCU`. Il primo passaggio costruisce le medesime tracce di estremita, calcola `rho_m` e `g'` tramite `eval_spatial_buoyancy` per ogni endpoint e determina le velocita caratteristiche. Il secondo passaggio ricalcolava `rho_m` e `g'` da temperatura e composizione per ottenere `gamma=rho_m*g'` nei due endpoint di ogni interfaccia, con due chiamate aggiuntive al kernel di miscela/densita. La cache delle velocita di faccia era gia proprieta di `hyperbolic_workspace_type`, e il percorso CFL/spatial corrente puo riusarla sullo stato immutato.

**Implementazione:** `equation_terms_2d:eval_local_speeds_direction` espone un output `gamma` opzionale, scritto usando gli `r_rho_m` e `r_red_grav` **gia valutati**. `hyperbolic_2d:eval_interface_speed_bounds` riceve due output gamma e li produce nella stessa visita dei due endpoint; `hyperbolic_workspace_type` conserva quattro array REAL `gamma_interface_xL/xR/yB/yT` affiancati ai quattro array di velocita (allocati e azzerati in initialize, deallocati in finalize). `eval_pccu_interface` legge i gamma memorizzati e non richiama `eval_hydrostatic_coefficient` su ciascuna faccia. **Invalidazione esatta:** i gamma vengono sovrascritti ogni volta che `eval_speeds` rigenera le tracce; nel percorso `reuse_cfl_traces` non vengono modificati, esattamente come le velocita delle stesse facce. Il percorso dry-donor puo azzerare la profondita tra la stima delle onde e il calcolo del path: per i soli endpoint con profondita finale `<=EPSILON` il gamma del path viene impostato a zero prima della stessa `regularize_dry_gamma_pair` preesistente. Temperatura e frazioni massiche restano quelle delle tracce selezionate; non si riusa un coefficiente dopo una loro mutazione.

**Costo e tradeoff da misurare:** eliminate fino a **due rivalutazioni densita/miscela per faccia attiva** nel calcolo del path di interfaccia; introdotti quattro array scalar-per-face aggiuntivi (due endpoint per orientazione), senza nuove allocazioni nel loop. Gli altri path *interni alla cella* restano calcolati indipendentemente: possono differire per composizione/temperatura/ricostruzione e non e sicuro associarli alla cache di interfaccia. Il modello `gamma` non diventa un campo dello stato fisico generale: e scratch numerico localizzato nel workspace HP-PCCU. Le formule per onde, gamma, wet/dry, regularizzazione e path non sono cambiate intenzionalmente. Nessuna topografia statica; flussi e transazioni Lie invariati.

**Prima di dichiarare C4 completato:** eseguire build manuale e confronti HP-PCCU/well-balanced, shoreline wet-dry, termini di sorgente al bordo, multidimensionalita e riuso della ricostruzione CFL con e senza `reuse_cfl_traces`. Misurare eventuale beneficio rispetto al costo di memoria/banda delle nuove matrici. Nessuna Action avviata in questo intervento.

### D.6 Registro C4/C5 — riuso dei coefficienti nel controllo well-balanced (5 novembre 2026)

**Codice modificato, verifica del nuovo commit ancora pendente.** La seconda passata C4/C5 sul percorso `hyperbolic_2d:eval_hyperbolic_terms -> cancel_resting_roundoff` ha individuato fino a quattro ulteriori `eval_hydrostatic_coefficient` per ogni direzione della cella esattamente ferma. La stima della scala di roundoff utilizza le stesse tracce di faccia finali, e `eval_speeds` ha gia salvato `gamma=rho_m*g'` delle tracce selezionate in `hyperbolic_workspace_type` (primo intervento C4).

**Implementazione:** `cancel_resting_roundoff` legge i quattro coefficienti presenti nel medesimo workspace, senza nuovi array o ownership. `hydrostatic_arithmetic_scale` utilizza il gamma salvato **solo quando sia la profondita di traccia effettiva, sia la profondita selezionata per le onde sono wet secondo la soglia interna `EPSILON(1.0_wp)`**. La differenza di profondita non modifica `gamma` se temperatura e frazioni massiche sono invariate. Quando la traccia effettiva e dry, usa esplicitamente il gamma zero del kernel originale; se la traccia effettiva e wet ma quella selezionata e dry, rivaluta il coefficiente con il vecchio `eval_hydrostatic_coefficient`. Un residuo gia esattamente zero restituisce direttamente zero, senza valutazioni termodinamiche. La scala aritmetica, la soglia `64*epsilon`, le tracce, l'ordine della riduzione `MAX` e la decisione sul residuo non sono alterati nei rami generali.

**Risparmio strutturale potenziale:** evitare fino a otto chiamate di densita per una cella ferma 2D (quattro in x, quattro in y) quando i gamma degli endpoint sono riutilizzabili; meno nei casi dry/near-dry con fallback. Non e ancora un risultato misurato; il controllo viene eseguito soltanto su stencil esattamente fermi. I coefficienti restano locali al workspace dell'operatore spaziale, che viene aggiornato a ogni ricostruzione: topografia non resa statica.

**Verifiche ancora aperte:** build, well-balanced e lake-at-rest su topografia complessa; interfacce near-dry, selezione di profondita differente, CFL e spatial in sequenza e chiamanti senza riuso CFL; 1/2/4/8 thread; confronto numerico e wall time sugli stessi input e SHA. Se l'ottimizzazione peggiora prestazioni o robustezza, rivedere il riuso senza cambiare il contratto scientifico. C4 e C5 restano IN CORSO; nessuna Action avviata automaticamente.

### D.7 Audit statico trasversale completo — disponibilita e riuso degli stati (5 novembre 2026)

**Baseline esaminata:** `a7a3ff872836b1366447955e2d9283a44d7376ec` (`master`), confermata dall'utente come compilabile e funzionante su esempi. **Ambito:** driver principale, inizializzazione/I-O, CFL e ricostruzione, HP-PCCU, workspace Lie, gruppi espliciti/mass exchange/analitici, Newton REAL/COMPLEX, stato, proprieta di miscela e leggi fisiche. Audit di sorgenti e call graph, **non** un profilo runtime: le frequenze sotto sono limiti per cella o chiamata, non contatori misurati.

#### D.7.1 Risultati con evidenze e classificazione

| ID | Priorita | Percorso e prova nel codice | Identita/invalidazione | Giudizio |
| --- | --- | --- | --- | --- |
| **A01** | **P0** | `src/IMEX_SfloW2D.f90`, ciclo subito dopo `simulation%advance_timestep` (circa righe 485-506): `advance_lie` ha gia pubblicato `state%q` **insieme** a `state%qp` (in `time_integration_2d.f90`, circa 299-305), ma le diagnostiche effettuano ancora `conservative_to_primitive` per ciascuna cella wet. | Stessa versione post-Lie di `q`; nessuna mutazione tra publish e diagnostica. `p_dyn` e richiesto per `pdynmax`/vulnerabilita, ma **non** giustifica rifare tutta la conversione `q->qp`. | **RICALCOLO CERTO, HOT PATH a ogni timestep.** Separare `p_dyn` da conversione; preservare le frazioni originali usate dal ramo `numerical_reconstruction`, che possono differire dalle frazioni normalizzate memorizzate in `qp`. Non sostituire alla cieca con `mass/qp_h` o densita da `qp` senza prova di equivalenza numerica. |
| **A02** | **P0** | `src/numerics/lie_source_executor_2d.f90` circa 641-651: dopo i gruppi, per ogni cella con `q_local(1)>0`, `build_physical_state_from_qp(qp_local,...)` viene richiamato per leggere soltanto `red_grav` e decidere l'eventuale sink di buoyancy. | `qp_local` appena convertito dopo il gruppo o gia disponibile post-spatial; geometria assente nella chiamata, tutte le altre proprieta dello stato costruito vengono scartate. | **OVERCOMPUTATION CERTA** (miscela completa, cp, volumi, kinematica, Ri non consumati). Calcolare soltanto `red_grav` con la **stessa** chiusura composizionale e semantica dry del costruttore; riusare prodotti esistenti quando validi. Non saltare il sink ne modificare il rollback. |
| **A03** | **P0** | `src/physics/mass_injection_terms/mass_radial_injection_2d.f90`, `mass_bottom_circular_injection_2d.f90`, `mass_bottom_fissural_injection_2d.f90`: `close_composition`, `evaluate_mixture_density`, `evaluate_mixture_heat_capacity` del materiale **di sorgente** sono dentro le routine cell-local, ripetute a parita di `T_source`, `xs_source`, `xg_source` o parametri della fissura. | Proprieta *della sorgente*, non dello stato destinatario; uguali per tutte le celle che condividono lo stesso materiale e configurazione. Distinguere parametri materiali (validita) da finestre/coefficienti temporali `set_source_time_coefficients` (variano a ogni step). Le fissure possono avere temperature diverse. | **RICALCOLO INTER-CELLA CERTO** nelle sorgenti attive. Prevalutare una volta per materiale/fissura **nella fase corretta**, passare i valori read-only ai kernel cell-local; non precomputare geometria/topografia, non introdurre indici in `physics/`. |
| **A04** | **P1** | `src/numerics/time_integration_2d.f90` circa 163-191: dopo lo spatial `q_trial=q^n-dt*spatial_rate`, viene convertita ogni cella wet risolta, anche quando la differenza conservativa e esattamente zero. | `qp_pre` e gia pronto dalla CFL o dalla preflight, ma **solo** nel medesimo stato e con worklist/halo coerenti; se `q_trial/=q^n`, conversione indispensabile. | **OPPORTUNITA CONDIZIONALE:** evitare la conversione sulle sole celle con `q_trial==q^n` e `qp_pre` effettivamente valido, mantenendo admissibility, dry e diagnostica del fallimento. Valutare costo confronto e guadagno sulle zone inattive. |
| **A05** | **P1** | `src/IMEX_SfloW2D.f90` circa 425-487 e `spatial_operator_2d.f90`: `compute_spatial_timestep` converte l'halo di ricostruzione e prepara tracce/onde; `advance_lie` riusa tutto **solo** quando `cfl/=-1`, stocastico e tutte le iniezioni disattivate. | Le iniezioni sono normalmente applicate come gruppi **dopo** lo spatial, ma possono cambiare coefficienti temporali e condizioni di bordo della ricostruzione; `stochastic%prepare_timestep` avviene tra CFL e spatial. | **OPPORTUNITA DA VERIFICARE, NON UN BUG CONFERMATO:** analizzare separatamente `q`, halo, boundary forcing, geometria e forcing stocastico. Non allargare il flag con una condizione booleana senza dimostrare identita delle tracce. Con `cfl==-1` non c'e preflight CFL da riusare. |
| **A06** | **P1** | `src/physics/source_terms_2d.f90` circa 414-465: `eval_source_implicit_real` e `eval_source_implicit` costruiscono uno stato fisico COMPLETO anche quando il gruppo Newton e soltanto `thermal_active`; la legge di cooling consuma essenzialmente T. | Stato REAL del residuo e COMPLEX del Jacobiano **differiscono** per perturbazione: non riusare uno stato fisico REAL per COMPLEX ne tra colonne. Pero ciascuna valutazione puo costruire solo cp/T con le formule tipate condivise. | **OVERCOMPUTATION CERTA per gruppi thermal-only**, da correggere con kernel minimo REAL/COMPLEX senza duplicazione delle formule, mantenendo sensitivity complex-step e contratti dry. |
| **A07** | **P1** | `src/io/inpout_2d.f90` circa 4446-4560: nel setup di sorgenti, dopo `evaluate_mixture` c'e un `build_physical_state_from_qp` completo sulla stessa miscela; altri costruttori completi sono ripetuti, incluso un ramo in cui `qp_source` non cambia. | Materiale della sorgente costante, ma modifiche a `h_source` o `vel_source` invalidano cinematiche e Ri. Densita/cp e gravita ridotta non dipendono da h e u se T/composizione restano uguali. | **DUPLICAZIONE CERTA MA COSTO DI SOLO SETUP.** Separare proprieta materiale da Ri/cinematica; eliminare chiamate su stato identico. Priorita inferiore agli hot path. |
| **A08** | **P1** | `src/spatial/hyperbolic_2d.f90`, circa 990-1083: i path idrostatici **interni** alla cella chiamano ancora `eval_hydrostatic_coefficient` per entrambi gli endpoint. C4 gia riusa gamma per interfacce e scale well-balanced. | I due endpoint di path interno possono avere profondita diversa dalle tracce dei limiti delle onde; T e composizione vanno confrontate con la corrispondente traccia. La gravita dell'endpoint e valutata ogni stage. | **POTENZIALE RIUSO CONDIZIONALE**, non dichiarare duplicazione senza verificare identita delle tracce. I gamma di faccia possono essere usati solo se preservano wet/dry, closure e valori dei coefficienti. |
| **A09** | **P2** | `src/physics/rheology_laws/rheology_2d.f90` circa 314-336: ramo analitico multi-legge copia `current_state=state` (tipo con array allocabili), ricostruisce cinematica fra leggi. | La chiamata produttiva `source_terms_2d` passa un **solo** modello tramite `[selected_model_name]`, ed entra nel ramo gia ottimizzato senza copia. Per molte leggi consecutive lo stato cinematico evolve, quindi l'aggiornamento e legittimo. | **NESSUNA DUPLICAZIONE PRODUTTIVA ACCERTATA.** Possibile ottimizzazione del ramo multi-legge solo dopo averne verificato l'uso, senza perdere semantica sequenziale. |
| **A10** | **P2** | `src/state/state_types.f90` espone `promote_physical_state` con assegnazioni di array allocabili, ma la ricerca dei chiamanti di produzione non ne individua alcuno: attualmente e usato nei test. `refresh_physical_state` e definito in `state_physical.f90`, senza chiamanti di produzione individuati. | La conversione COMPLEX dei Jacobiani e costruita da candidati perturbati distinti; non richiede promozione del REAL. | **NON CONTARE come costo hot-path.** Prima di rimuovere eventuali API verificare i consumatori esterni; nessuna cancellazione in questo audit. |
| **A11** | **P2** | `src/numerics/lie_source_executor_2d.f90`: workspace REAL/COMPLEX ed exchange allocati al bisogno **per worker** fuori dal loop; `src/state/state_physical.f90` esegue reset degli array di composizione a ogni costruzione. | Ogni nuovo stato deve sovrascrivere le proprieta; alcuni zeroing possono essere evitabili solo dimostrando che nessun percorso di rifiuto legge intermedi vecchi. | **WORKSPACE RIUSATI, non allocazioni per cella** nel percorso produttivo; possibile riduzione dei reset, ma serve audit di tutte le vie di early return. Non ottimizzare alla cieca. |
| **A12** | **P2** | `src/IMEX_SfloW2D.f90`, inizializzazione circa 305-345: conversione iniziale completa delle celle wet anche se alcune inizializzazioni producono gia `qp`. | Restart binario, restart legacy, empty-init e collapsing-volume hanno contratti diversi; non assumere cache valida all'avvio. | **DA CLASSIFICARE PER TIPO DI AVVIO**, impatto una tantum, nessuna prova di ridondanza generale. |

#### D.7.2 Lavoro gia evitato e valutazioni NECESSARIE

- **Gia riusati:** `qp` da CFL allo spatial nel percorso opt-in; facce e velocita da CFL nello spatial; `physical_ready` fra processi espliciti ed exchange dello stesso gruppo; residuo del trial Newton accettato alla prossima iterazione; `qp_local` quando un gruppo Lie non cambia `q`; cp nel costruttore da `q`; gamma delle facce per path d'interfaccia e well-balancing; stato termico minimo del raffreddamento analitico wet.
- **Necessari:** conversione del candidato dopo modifiche conservative spatial/Lie; ricostruzione completa nei processi che usano miscela, cinematica, Ri e geometria; Jacobiano COMPLEX per ciascuna perturbazione; conversione speciale dei primitivi alle facce per il flusso conservativo; ricostruzione di velocita caratteristica quando cambiano tracce, topografia o condizioni al bordo.
- **Non sono stati rilevati** duplicazioni di formula nei dispatcher REAL/COMPLEX gia unificati in `.inc`; questo non esclude altri calcoli identici su punti differenti del call graph.
- **C0 — inventario statico CONCLUSO** per il percorso produttivo e i consumatori di stato principali. Il conteggio dinamico e la certificazione di C5 rimangono APERTI: una chiamata statica in un loop non equivale al numero di esecuzioni o al costo.

#### D.7.3 Piano correttivo vincolato agli input

**Sequenza suggerita per rendere computazionale il traguardo senza nuove ristrutturazioni:**

1. **P0/A01, diagnostica post-Lie:** eliminare la conversione completa del driver, mantenendo il `qp` pubblicato. Ricavare `p_dyn` con una **funzione minima scientificamente equivalente** al ramo `conservative_to_primitive(...,numerical_reconstruction=.TRUE.)` per le stesse `q`; confrontare wet, near-dry, overflow/overfill, floating-point e massima pressione. Se l'output `p_dyn` va conservato nel workspace, motivarne il costo memoria rispetto al ricalcolo locale.
2. **P0/A02, buoyancy:** introdurre la valutazione minima di `red_grav` a partire dal `qp_local` coerente e confrontarla con `build_physical_state_from_qp` per ogni caso dry/near-dry e composizione. Non modificare la rimozione di massa.
3. **P0/A03, materiali delle iniezioni:** produrre `rho_m`, `cp_m`, `red_grav` e frazioni chiuse una sola volta per materiale/fissura configurata; passare i dati locali read-only alle leggi, senza `j,k` dentro `physics`. Separare invalidazione dei parametri materiali e coefficienti temporali; tenere topografia aggiornabile a ogni stage.
4. **P1/A06, Newton thermal-only:** percorso minimo da `q` REAL/COMPLEX per cp/T con `.inc` comune, senza stato completo; testare derivata complex-step, energy-only e gruppi combinati. Non condividere stati fra colonne perturbate.
5. **P1/A04/A05/A08:** riuso selettivo nel post-spatial, esame rigoroso del flag CFL/spatial con sorgenti e verifica identita degli endpoint interni HP-PCCU. Misurare overhead di confronti e cache prima di aggiungere altri array.
6. **P1/P2/A07/A09/A11/A12:** semplificazioni setup, reset e API realmente inutilizzate **solo** se evidenza di costo o manutenzione giustifica la modifica.
7. **C5 — gate finale:** contatori `q->qp`, costruttori REAL/COMPLEX, `evaluate_mixture*`, funzioni densita/cp, allocazioni, facce/well-balanced, tempi per operatori e gruppi; confronti `EXAMPLE_CF`, `EXAMPLE_TAAL` e altri casi con fisiche differenti; same SHA/input/build/OpenMP 1/2/4/8; risultati scientifici prima/dopo, senza eseguire Actions automaticamente.

**Regola:** ogni modifica deve separare *esistenza dello stato gia pronto*, *identita dell'input*, *prodotti davvero consumati*, *punto di invalidazione*, *costo di conservarli* e *rischio numerico*. Non introdurre nuovi wrapper di passaggio, un secondo backend, cache globali su tutta la griglia o precomputi geometrici statici.

### D.8 Implementazione dei risultati A01-A03 dell'audit (5 novembre 2026)

**Stato: CODICE MODIFICATO, BUILD/REGRESSIONI/MISURE ANCORA DA VERIFICARE.** La baseline confermata dall'utente prima di questi interventi e `a7a3ff8`. I risultati statici documentati in D.7 restano validi come motivazione; questa sezione prevale per lo stato degli interventi.

| ID | Modifiche effettive | Invarianti e validazione residua |
| --- | --- | --- |
| **A01 / IMPLEMENTATO** | `src/state/state_conversions.f90`: `dynamic_pressure_from_synchronized_state(q_mass,qp)` ricava `rho_m=q_mass/qp_h` e la pressione dinamica direttamente dai due vettori gia sincronizzati. `src/IMEX_SfloW2D.f90`: il loop diagnostico **post-Lie** non richiama piu `conservative_to_primitive` sulle celle wet; sfrutta `state%qp` gia pubblicato e calcola solo `p_dyn`. Conversione e controllo iniziali all'avvio/restarter rimangono invariati. | Almeno una conversione completa in meno **per cella wet e per timestep** nel driver. La densita ricostruita con `mass/h` puo differire per gli ultimi bit dalla densita originale; confrontare `p_dyn`, `pdynmax` e classificazioni vulnerabilita esattamente sulle soglie e sui casi near-dry. Il controllo di ammissibilita non e stato rimosso dall'integratore Lie; nessun nuovo stato persistente. |
| **A02 / IMPLEMENTATO, CRITERIO FISICO AGGIORNATO SU RICHIESTA** | `src/numerics/lie_source_executor_2d.f90`: rimosso il `build_physical_state_from_qp` finale usato solo per buoyancy. Ora si chiudono le frazioni composizionali con `close_composition`, si calcola **solo `rho_mix`** via `evaluate_mixture_density` (o il valore ambientale nel ramo dry) e si produce `red_grav` per le diagnostiche con la stessa espressione scalare. Il sink si applica se e soltanto se `rho_mix<rho_a_amb`. | La nuova condizione e volutamente **stretta**: miscele esattamente neutrali non vengono rimosse; la versione precedente applicava `red_grav<=0`, che includeva anche l'uguaglianza. Il conteggio diagnostico usa ora `red_grav<0`. Preservati controllo di composizione, errore 15, workset per-worker e pubblicazione/rollback. Verificare fronti wet/dry, celle inizialmente asciutte con iniezione e TAAL. |
| **A03 / IMPLEMENTATO** | `src/physics/mass_injection_terms/mass_injection_coefficients_2d.f90`: nuovo **contenitore tipato** `injection_material_set_type` e preparazione `prepare_injection_materials` fuori dal loop OpenMP: chiusura composizionale, densita, cp e gravita ridotta **una volta per materiale comune e una per fissura**, per ogni sweep Lie. `lie_source_executor_2d` possiede il contenitore e lo passa in sola lettura fino a `source_terms_2d`, `mass_injection_2d` e ai tre modelli (radiale, circolare, fissurale). `src/Makefile.am` dichiara le dipendenze Fortran aggiornate. | Nessuna cache `SAVE` globale, nessun indice spaziale in `physics`, nessun precomputo di topografia. Le proprieta vengono ricalcolate **una sola volta prima del sweep** per accettare configurazioni cambiate; i coefficienti temporali e i rate per cella non vengono memorizzati. Le fissure conservano composizione e temperatura proprie. Verificare percorsi senza sorgenti, input per-volume, fissure eterogenee, finestre temporali, OpenMP e disponibilita degli array al setup. |

**Ordine delle verifiche manuali raccomandate:** build su HEAD; `EXAMPLE_TAAL` con diagnostiche A02, casi di iniezione radiale/circolare/fissurale (inclusi dry recipient e valori composizionali limite); confronti di `p_dyn` e soglie vulnerabilita; esempio senza iniezione e restart; infine profili di `conservative_to_primitive`, `build_physical_state_from_qp` e `evaluate_mixture_*` e tempi con 1/2/4/8 thread. Un run verde non dimostra ancora l'assenza di differenze scientifiche al limite delle soglie o un beneficio di wall time.

**Nota sui gate:** A01-A03 sono modifiche al codice, **non ancora chiusure verificate**. I blocchi C4/C5 e le altre opportunita A04-A12 rimangono in roadmap. Non avviare automaticamente Actions o benchmark.

### D.9 Correzione ownership A02 — buoyancy esclusivamente in physics (5 novembre 2026)

**Risolto il rilievo architetturale dell'utente, codice modificato e build da verificare.** Il primo intervento A02 aveva impropriamente inserito in `numerics/lie_source_executor_2d.f90` la chiusura della composizione, il calcolo `rho_mix` e il criterio di galleggiabilita. Questi elementi non appartengono alla numerica.

La valutazione e stata spostata integralmente in `physics/source_terms_2d.f90`, nella routine **cell-local** `evaluate_buoyancy_removal(qp,remove_cell,red_grav,accepted)`. La routine chiude la composizione usando `close_composition`, calcola soltanto la densita con `evaluate_mixture_density`, conserva per le diagnostiche la formula `red_grav=(rho_mix-rho_a_amb)*inv_rho_mix*grav` e decide l'eliminazione tramite **`rho_mix<rho_a_amb`**, come richiesto. Se la composizione non e ammissibile `accepted=.FALSE.`. Nessun indice spaziale, campo di griglia, mutazione di stato globale o dipendenza dalla numerica e introdotto nella nuova interfaccia.

L'esecutore Lie invoca soltanto il criterio fisico, interpreta `accepted` (errore 15), registra le diagnostiche, e applica/committa la rimozione del conservativo e la primitiva dry quando `remove_cell` e vero. Sono stati rimossi da Lie gli import `state_mixture`, `state_composition`, `constitutive_parameters_2d` usati per A02 e i relativi intermedi REAL e array di composizione. L'ottimizzazione A02 e preservata: **nessuno stato fisico completo ricostruito per la buoyancy**. Neutrali conservati per scelta esplicita.

**Gate:** verifica manuale del compilatore, EXAMPLE_TAAL, miscela negativa/neutra/positiva rispetto all'atmosfera, casi dry e diagnostiche del D19. Nessuna Action avviata dal cambio.

### D.10 Correzione ulteriore A02: le proprietà della miscela appartengono a state (5 novembre 2026)

**Vincolo architetturale esplicitato dall'utente:** non basta spostare il calcolo di densita e gravita ridotta da Lie a `physics`; **le proprieta dello stato devono essere calcolate dai moduli `src/state/`**. L'intervento D.9 era quindi incompleto quanto a ownership.

**Nuova distribuzione:**
- `src/state/state_physical.f90` espone `evaluate_density_from_qp(qp,rho_m,red_grav,accepted)`: chiusura della composizione, ramo dry, densita della miscela tramite il kernel `state_mixture:evaluate_mixture_density`, gravita ridotta. Non costruisce `physical_state_real`; non dipende da `physics` ne da `numerics`. Il contratto `accepted` rifiuta composizioni non ammissibili.
- `src/physics/source_terms_2d.f90:evaluate_buoyancy_removal` chiede quelle **proprieta** a `state_physical` e applica solo il **criterio fisico di rimozione** `rho_mix<rho_a_amb`. Non valuta autonomamente composizione, densita, termodinamica o gravita ridotta.
- `src/numerics/lie_source_executor_2d.f90` invoca il processo fisico, interpreta `accepted`/errore 15, registra le diagnostiche e pubblica l'eventuale rimozione conservativa. Nessuna proprieta di stato calcolata in Lie.

**Distinzione vincolante per il resto del piano:** proprieta della miscela e derivate di `q`/`qp` in **state**; *leggi/processi e criteri fisici* in **physics**; integratori, accettazione, commit e rollback in **numerics**. Per A03 le leggi di iniezione possono coordinare il calcolo una volta per sweep, ma devono ottenere densita, cp e gravita ridotta dai kernel di `state` senza duplicarne le formule.

**Verifica pendente sul nuovo SHA:** compilazione, TAAL, neutrali e wet/dry, confronti di `red_grav` con il costruttore fisico completo e diagnostiche D19. Nessuna Action avviata automaticamente.

### D.11 Ownership A02 definitiva: criterio di densita come proprieta dello stato (5 novembre 2026)

**Decisione architetturale dell'utente:** anche la classificazione `rho_mix<rho_a_amb` appartiene a `state`. Non deve esistere un wrapper in `physics` che contenga solo questo confronto.

**Modifica effettuata:** `state_physical:evaluate_density_from_qp(qp,rho_m,red_grav,accepted,below_ambient)` restituisce, tramite output opzionale `below_ambient`, la classificazione di stato, inizializzata a `.FALSE.` e valorizzata solo dopo chiusura composizionale riuscita. La stessa routine continua a produrre la densita e la gravita ridotta con i kernel canonici, senza costruire un physical state completo. Eliminata interamente `source_terms_2d:evaluate_buoyancy_removal`. `lie_source_executor_2d` chiama direttamente `state_physical`, controlla `accepted`, registra `red_grav` e applica l'azzeramento del conservativo quando `below_ambient` e vero. Il criterio resta **stretto**: `rho_mix<rho_a_amb`; le miscele neutrali sono mantenute. Non esiste una seconda formula in `physics` o `numerics`.

**Confine definitivo:** `state` possiede le proprieta derivate e i predicati di classificazione dello stato; `physics` possiede modelli/processi/rate fisici (senza contenitori per un semplice confronto di stato); `numerics` orchestra i gruppi, decide quando applicare gli aggiornamenti e gestisce transazioni. Nessun accesso a `j,k` in `state_physical` o `physics`. D.9-D.10 rimangono cronologia e sono superati da questa decisione. Build/test/TAAL sul nuovo SHA pendenti, Actions manuali.

### D.12 Interventi A04 e A06 e analisi preliminare A05 (5 novembre 2026)

**Baseline build confermata dall'utente:** `3dd1bfc`. Gli interventi seguenti sono modifiche successive e **non ancora compilate/verificate dall'utente**. Nessuna Action o benchmark avviato.

**A04 — IMPLEMENTATO, post-spatial senza conversione ridondante:** `numerics/time_integration_2d.f90:advance_lie` continua a costruire e verificare `q_trial=q^n-dt*spatial_rate`, senza limiter sulle rate. Superato il controllo di ammissibilita e prima dell'eventuale conversione wet, se **l'intero vettore** `q_trial(:,j,k)==state%q(:,j,k)` (uguaglianza esatta, nessuna tolleranza), ora copia `qp_pre(:,j,k)` in `qp_trial(:,j,k)`. Il `qp_pre` corrisponde allo stesso `q^n` sia sul percorso CFL/HP-PCCU, sia nel preflight esplicito che aggiorna tutte le celle della reconstruction worklist. Quando `q_trial` cambia, rimane intatta `conservative_to_primitive(...,numerical_reconstruction=.TRUE.)`; i rami dry, sea sink, accettazione/rollback e diagnostiche non vengono alterati. **Risparmio:** fino a una conversione completa per cella wet invariata sul passo spaziale; il costo del confronto e della copia va misurato. Controllare celle wet, near-dry, sorgenti con rate spaziale nullo e presenza/assenza della preflight CFL.

**A06 — IMPLEMENTATO, Newton solo termico a temperatura scalare REAL/COMPLEX:** `state/state_physical.f90` espone `evaluate_temperature_from_q` con due procedure tipate REAL/COMPLEX per ricavare la temperatura dal conservativo usando unicamente il calore specifico e le frazioni massiche. Il ramo REAL riproduce la formula e i test su mass del costruttore REAL; il COMPLEX mantiene la chiusura composizionale e la temperatura perturbata del costruttore COMPLEX, essenziali per il complex-step. `physics/thermal_laws/thermal_cooling_2d.f90` espone due evaluator scalari per lo stesso rate; per evitare duplicazioni, la formula in `thermal_cooling_rate.inc` usa una variabile locale `temperature` sia nei percorsi esistenti con stato fisico completo sia nei due nuovi percorsi scalari. `thermal_laws_2d.f90` coordina i due evaluator con una dispatch per ID del modello. `physics/source_terms_2d.f90` usa il ramo minimo solo se `group%thermal_active` e **non** `group%rheology_implicit`, sia per residuo REAL sia per Jacobiano COMPLEX. I gruppi combinati mantengono il costruttore completo e la formula condivisa; non esiste cache cross-perturbazione. **Risparmio:** una costruzione dello stato fisico completo evitata per chiamata al rate termico-only, lasciando necessario il calcolo di T e `cp_m`. Confrontare residui, Jacobiani complex-step, convergenza e risultati dei casi thermal-only, lava e gruppi misti.

**A05 — ANALISI, NON IMPLEMENTATO:** `IMEX_SfloW2D.f90` imposta `pre_spatial_primitives_ready` soltanto per CFL abilitata, nessun forcing stocastico e nessuna delle sorgenti di iniezione abilitate. Il flag comporta il riuso non solo di `qp`, ma delle tracce ricostruite e delle velocita di faccia in `spatial_operator_2d`/`hyperbolic_2d`. Non basta constatare che l'iniezione e applicata dopo lo spatial: bisogna provare che condizioni al bordo, forcing temporale, worklist, topografia e selezione di faccia coincidono fra CFL e spatial. **Nessun allargamento del flag** finche la verifica non e completata.

**Residui C5:** build e regression tests sul nuovo HEAD; confronto `EXAMPLE_TAAL`, caso thermal-only con Jacobiano, scenario con sorgenti ed effetti wet/dry; contatori conversione per cella e valutazioni termodinamiche REAL/COMPLEX; profiling OpenMP e diagnostiche scientifiche. Non dichiarare completato C5 in assenza di misure.

### D.13 A05 — confronto CFL/spatial ed estensione selettiva alle iniezioni conservative (5 novembre 2026)

**Baseline verificata dall'utente:** build, `EXAMPLE_TAAL` ed `EXAMPLE_ETNA_LAVA` sul commit `565fdf8`. **Nuova modifica in attesa di build e regressioni.**

**Audit di identita del dato lungo il call graph.**

1. `IMEX_SfloW2D.f90` richiama `domain%check_solve` e `spatial_operator%prepare_hp_workset` **prima** di CFL. `compute_spatial_timestep` converte tutte le celle dell'halo di ricostruzione, ricostruisce le primitive e le eta HP, valuta le velocita e i gamma di faccia. Il passo `dt` puo essere ridotto per l'output o il passo precedente, ma questo **non modifica** la ricostruzione, dipendente dal tempo di inizio passo e dalle tracce, non da `dt`.
2. `simulation%advance_timestep` chiama `set_source_time_coefficients` tra CFL e Lie: modifica solo i coefficienti di attivazione/velocita delle **sorgenti conservative**; le leggi circolare, radiale e fissurale vengono valutate nella successiva transazione Lie, **dopo** `spatial_operator%evaluate`. I coefficienti di finestra temporale non sono consumati dalla `reconstruction_2d`. Per la sorgente radiale, la geometria `source_cell` e i vettori di bordo sono gia configurati prima della CFL, e `reconstruct` usa le stesse tracce conservative (inclusa la condizione della sorgente radiale) in entrambe le valutazioni.
3. `stochastic%prepare_timestep` puo mutare la componente `q(idx_stoch)` e `qp(idx_stoch)` **se** `stoch_transport_flag` e attivo; anche quando modifica solo i campi OU esterni, la compatibilita non e stata testata per ogni configurazione. Si mantiene quindi prudentemente l'esclusione `stochastic_flag`.
4. `lateral_source_flag`: `compute_spatial_timestep` riceve un override esplicito `boundary_source_enabled` ottenuto dal registry dei gruppi espliciti, mentre la ricostruzione spatial fresca, senza override, usa il fallback `lateral_source_flag`. I due valori potrebbero differire in configurazioni particolari; il riuso viene **conservativamente disabilitato** in presenza di questa sorgente. Questo non implica che le tracce siano necessariamente diverse, ma evita un cambiamento non dimostrato delle condizioni al bordo.
5. `cfl==-1`: `compute_spatial_timestep` ritorna prima di convertire, ricostruire e generare le velocita: nulla e riusabile. Resta il percorso completo.

**Modifica A05 implementata:** in `src/IMEX_SfloW2D.f90`, il gate `pre_spatial_primitives_ready` richiede ora soltanto **CFL attivo, `stochastic_flag=.FALSE.`, `lateral_source_flag=.FALSE.`**. Non esclude piu `radial_source_flag`, `bottom_radial_source_flag` e `bottom_fissural_source_flag`: questi modelli non invalidano `q/qp`, tracce di faccia o gamma tra le due visite. Il percorso opt-in continua a riusare `state%qp`, `reconstruction_workspace`, velocita e coefficienti idrostatici di faccia, **senza nuovi array, flag di versione o cache globali**. Nessuna modifica ai kernel delle sorgenti, alla fisica, alla geometria, all'integratore o all'ordine Lie. La validazione numerica sugli input con sorgenti e ancora **pendente**.

**Regressioni necessarie prima di chiudere A05:** build; esempio con iniezione radiale e fissurale/circolare, almeno una sorgente attiva in celle wet e dry; confronto conservativi, tracce, CFL, termine spaziale e diagnostics con e senza cache; casi con `lateral_source_flag`, `stochastic_flag` e `cfl=-1` per verificare che i fallback restino invariati. Misurare numero di ricostruzioni, conversioni e tempi sulle configurazioni con iniezione. Questa e una **estensione selettiva implementata**, non una prova di equivalenza completa per tutte le sorgenti.

### D.14 Nuovo splitting stocastico — trasporto prima di Ornstein-Uhlenbeck (5 novembre 2026)

**Decisione utente:** il trasporto conservativo di `q_Z=MZ` viene eseguito nello spatial, **prima** del fractional-step OU. L'evoluzione OU viene quindi applicata sul candidato post-spatial e precede i gruppi fisici Lie, cosi la reologia consuma i valori stocastici aggiornati. Questa e una **modifica intenzionale allo splitting**: non e un semplice speedup bitwise-equivalent rispetto all'ordine precedente.

**Nuovo ordine operativo verificabile dai sorgenti:**
1. `IMEX_SfloW2D.f90` calcola CFL, primitive, tracce e velocita dal conservativo pubblicato `q^n`. Nessun OU muta il conservativo tra CFL e spatial. A05 puo attivare `pre_spatial_primitives_ready` anche con `stochastic_flag` e `stoch_transport_flag` (rimangono esclusi soltanto `cfl==-1` e `lateral_source_flag`).
2. `time_integration_2d:advance_lie` costruisce `q_trial,qp_trial` post-spatial, verifica la composizione, gestisce il sink marino e conserva la proprieta del rollback.
3. Nello stesso integratore, **prima dei gruppi Lie**, `stochastic%prepare_timestep(state,dt,q_trial,qp_trial,domain)` aggiorna OU leggendo gli stati post-spatial **solo nelle celle `domain%solve_mask`**, mentre per quelle inattive legge lo stato pubblicato. L'evoluzione OU continua anche nelle celle asciutte e conserva la generazione seriale del rumore, l'eventuale convoluzione spaziale e il kernel di intensita Froude/Lucas. `Z` e prima sincronizzato da `q_trial(idx_stoch)/q_trial(idx_totMassEqn)` per le celle wet attive, poi evolve con Euler-Maruyama e l'intensita ottenuta dalla **velocita post-spatial**. Dopo OU si scrivono `q_trial(idx_stoch)=q_trial(mass)*Z` e `qp_trial(idx_stoch)=Z` solo sulle celle attive: il conservativo pubblicato non e toccato.
4. `execute_lie_source_groups` riceve `stochastic%effective_Z` aggiornato, mantenendo i processi reologici nel medesimo passo; gruppo 0 resta fedele al suo snapshot pre-spatial `q_original`, mentre la componente stocastica appartiene al candidato aggiornato dopo spatial.
5. In caso di rifiuto delle sorgenti, `time_integration_2d` ripristina `stochastic%Z` ed `effective_Z` da due buffer persistenti allocati solo quando necessari, oltre a conservativo autorevole, substrato e geometria. La sequenza random draw non viene riavvolta: il fallimento produttivo esistente termina comunque in `fatal_error`, e non e implementato un retry. In caso di successo il commit pubblica `q,qp` dopo i gruppi.
6. `simulation_2d:advance_simulation_timestep` non modifica piu lo stato prima della chiamata all'integratore; l'evento `TIMESTEP_EVENT_STOCHASTIC_PREPARED` viene emesso tramite callback **dentro l'integratore subito dopo lo stage OU** (anche con stochastic disabilitato, per mantenere l'evento di fase), seguito da `TIMESTEP_EVENT_LIE_COMPLETE` e dal clock.

**Ownership:** `stochastic_module` conserva l'evoluzione OU e il calcolo della sua intensita; `state` rimane owner dei campi conservativi/primitivi e l'integratore possiede i candidati; `physics` consuma soltanto `effective_Z` nella reologia. Non vengono aggiunti calcoli fisici nel Lie executor, né indici di griglia in `physics`. Cambiamento del segnale `sigma` per Froude/Lucas: e ora valutato sulla velocita dopo spatial, non prima.

**Verifiche obbligatorie (non eseguite in questo intervento):** build, `EXAMPLE_TAAL`, `EXAMPLE_ETNA_LAVA`, simulazione realmente stocastica con e senza `stoch_transport_flag`, CFL attivo e `cfl=-1`, seed fisso con OpenMP 1/2/4/8, norm e bilancio di `q_Z` prima/dopo spatial/OU, dipendenza di `sigma` dalla velocita, wet/dry, corrispondenza `effective_Z`–reologia, rollback D20 e restart binario. Le realizzazioni non vanno attese identiche alla versione con OU prima di spatial; confrontare proprieta statistiche e conservazione anziche richiedere identita bit a bit.

**Stato:** implementazione del nuovo ordine **IN CORSO / CODE COMPLETE, NOT TESTED**, A05 riuso CFL/spatial stocastico **DA VERIFICARE**. L'ultimo SHA confermato dall'utente era `565fdf8`; `bd73c30` ha esteso il riuso alle iniezioni conservative, non ancora validato. Non dichiarare scientificamente equivalente il nuovo splitting senza test.

### D.15 A08 — riuso selettivo del coefficiente idrostatico delle facce nei path interni (5 novembre 2026)

**Scopo ristretto su decisione utente:** implementare soltanto A08. Gli interventi A07, A09, A10, A11 e A12 non vengono eseguiti.

**Audit endpoint HP-PCCU.** Nel calcolo di CFL e nella valutazione delle velocita di faccia `hyperbolic_2d:eval_speeds`, `eval_interface_speed_bounds` memorizza il coefficiente idrostatico `gamma_interface_xL/xR/yB/yT`, valutato dagli endpoint `qp_adv_interface*` dopo averne sostituito la profondita con `h_selected_interface*`. Nel path interno `eval_cell_hydrostatic_paths` si riparte dalle **stesse tracce primitive** `qp_adv_interfaceR(:,j,k)`, `qp_adv_interfaceL(:,j+1,k)` in x e `qp_adv_interfaceT(:,j,k)`, `qp_adv_interfaceB(:,j,k+1)` in y, ma si sostituisce la profondita con `MAX(eta_cell-B_face,0)`. Composizione e temperatura **non vengono mutate** in questo passaggio; solo altezza (e, se dry, velocita). Il kernel termodinamico `eval_spatial_buoyancy` usa la profondita unicamente per il ramo `h<=EPSILON(1.0_wp)`; quando entrambi gli endpoint sono wet, il valore di Gamma e quindi lo stesso.

**Ottimizzazione implementata:** `src/spatial/hyperbolic_2d.f90` passa al path interno i coefficienti di faccia gia calcolati, insieme alle rispettive profondita selezionate. `eval_cell_hydrostatic_path` riusa un coefficiente **solo se sia `qp_path(1)` sia `h_selected_interface` sono `>EPSILON(1.0_wp)`**. Negli altri casi chiama ancora `eval_hydrostatic_coefficient`. Resta inalterata la successiva `regularize_dry_gamma_pair`, fondamentale per evitare discontinuita spurie alle rive; non sono modificati la quadratura, gli endpoint dei path, i flussi, le grandezze geometriche o le interfacce con la fisica. Nessun nuovo array persistente e nessuna cache indipendente dal lavoro gia svolto per CFL. Il risparmio previsto e fino a **quattro calcoli di densita/composizione per cella attiva 2D** nei path interni quando le tracce sono wet; il numero reale dipende dalla porzione asciutta e dall'attivazione delle facce.

**Verifica pendente:** build del nuovo HEAD, EXAMPLE_TAAL, EXAMPLE_ETNA_LAVA e test lake-at-rest/well-balanced; confronto dei `P_cell_x/P_cell_y`, `Gamma`, flux, conservativi e costi computazionali rispetto al baseline. Attenzione ai casi wet/dry dove `h_selected` e `h_path` sono su lati opposti della soglia EPSILON: il fallback deve mantenere il comportamento precedente. Nessuna Action avviata.

### D.16 A08 ripristinato per regressione in EXAMPLE_TAAL (5 novembre 2026)

**Segnalazione utente sul commit `3e2b02d`:** TAAL esegue fino a circa `t=27.3632` e poi termina con `FATAL ERROR: Lie timestep rolled back: D19: inadmissible source candidate`. La build precedente e i due esempi avevano avuto esito positivo sul baseline `c8605e4`. L'unico cambiamento ai kernel numerici del commit sospetto era il riuso di `Gamma` tra path interni HP-PCCU e face speeds in `src/spatial/hyperbolic_2d.f90`.

**Azione correttiva immediata:** ripristinato **integralmente** `src/spatial/hyperbolic_2d.f90` dalla revisione `c8605e429fee93440be1097779a34cfe0fe8d920`, antecedente ad A08. Sono dunque ripristinate entrambe le valutazioni originarie di `eval_hydrostatic_coefficient` nei path di cella, senza valori `Gamma` cached. Gli altri file produttivi di A05 e dello splitting stocastico non sono modificati. **A08 = ROLLED BACK / NOT ACCEPTED**, da non riapplicare sulla sola ipotesi di equivalenza dei coefficienti.

Il messaggio D19 deriva dal controllo dell'ammissibilita del candidato sorgente: **non dimostra da solo** quale variazione locale del percorso idrostatico abbia causato il rifiuto; il rollback A08 e una misura diagnostica e cautelativa per isolare la regressione. `IEEE_UNDERFLOW_FLAG` non e prova di causalita. Per riprendere A08 servono un confronto diretto di Gamma cached/fresh, path, flux e invarianti e una prova lake-at-rest/wet-dry.

**Gate:** far rieseguire build ed EXAMPLE_TAAL sul commit di rollback fino oltre `t=27.3632`; idealmente EXAMPLE_ETNA_LAVA. Se TAAL fallisce ancora, analizzare D19 con log di gruppo, cella e candidati, senza attribuire automaticamente l'errore alla sola A08.

### D.17 Rettifica diagnostica TAAL t circa 27.36 — errore D19 preesistente al rollback A08 (5 novembre 2026)

**Nuova evidenza dall'utente:** `EXAMPLE_TAAL` si arresta nello stesso punto anche su `b373f6a`, dove `hyperbolic_2d.f90` e stato ripristinato alla versione precedente ad A08. L'utente precisa che i run precedenti, ritenuti riusciti, potrebbero essere stati interrotti prima del tempo critico. **Non esiste ancora una baseline verificata oltre t circa 27.36.** E pertanto scorretto classificare il difetto come regressione dimostrata di A08. A08 rimane revertito per prudenza, non giudicato come causa.

**Punto di fallimento nel codice:** `lie_source_executor_2d.f90`, `cell_error=6`, `conservative_mass_admissible(q_candidate(1:n_eqns))` falso durante un gruppo esplicito; `q_candidate=q_local+dt*total_rate`. I controlli e i limiter dei rate non vanno alterati prima di identificare la causa. I dati essenziali sono gia stampati con `IMEX_LIE_DIAGNOSTICS=1`: `D19_BUDGET cell/group/dt/recipient/scaled_rate/candidate/reconstruction_match`, oltre alla spiegazione del fallimento. Non serve una nuova patch strumentale.

**Piano diagnostico:** rieseguire `EXAMPLE_TAAL` oltre 27 s con `IMEX_LIE_DIAGNOSTICS=1`, acquisire tutti i record `D19_BUDGET` e la riga `FATAL ERROR`, individuare la sorgente `group`, la componente negativa/non fisica e il bilancio del rate. Solo dopo risalire al termine sorgente di quella componente. Distinguere bug del modello da instabilita condizionale da timestep; nessuna correzione numerica o rollback generalizzato senza evidenza.

### D.18 TAAL D19 — perdita di ammissibilita per arrotondamento subnormale (5 novembre 2026)

**Diagnostica dell'utente sul baseline senza A08:** `t=27.3632`, `dt=1.9974649388436672E-01`, cella `(237,394)`, gruppo esplicito `4` (il numero gruppo da solo non dimostra quale sorgente del registry vi appartiene). `D19_BUDGET recipient(5)=1.0280442421462131E-313`, `scaled_rate(5)=-5.1467448673525914E-313`; con `q_candidate=q_local+dt*total_rate` risulta `candidate(5)=-4.9406564584124654E-324`, un solo quanto subnormale negativo. Il `reconstruction_match=T` conferma che i dati stampati ricostruiscono esattamente il candidato. Il problema e nel limite numerico del rate sorgente, non e una prova contro A08 o lo splitting stocastico.

**Causa numerica identificata:** `source_terms_2d:limit_explicit_source_rates` gia limita le sorgenti usando gli inventari fisici e un buffer relativo `4096*EPSILON`. Tuttavia `4096*EPSILON*1E-313` e minore della minima quantita rappresentabile, dunque il buffer si annulla per arrotondamento; una differenza di un ulp subnormale tra il candidato calcolato nel limitatore e quello costruito dal Lie executor puo produrre una componente negativa rigorosamente inammissibile. Il secondo controllo con al piu 256 passi `NEAREST` del limiter non garantisce la riserva necessaria a questa scala.

**Fix minimo:** SOLO nella routine gia esistente `limit_explicit_source_rates`, prima del controllo finale del candidato, per ciascuna componente conservativa vincolata (massa, solidi, gas aggiuntivi) con inventario positivo subnormale e rate negativo, riservare `MIN(q_reference,32*SPACING(q_reference))` prima di costruire il bound `limiter <= (q_reference-reserve)/(-dt*rate)`. Il margine e espresso in **quanti rappresentabili di quella componente**, senza threshold fisico arbitrario, e si applica soltanto a specie davvero subnormali soggette a consumo. La modifica riduce soltanto il fattore comune al gruppo esplicito; non taglia né modifica `q`, non allenta `conservative_mass_admissible` e non introduce limiter sullo spatial. Tutta l'operazione di sorgente (comprese le componenti di exchange) mantiene lo stesso fattore comune.

**Test pendenti:** ricompilazione; TAAL oltre `t=27.3632` e fino a fine esempio, con e senza `IMEX_LIE_DIAGNOSTICS=1`; ETNA_LAVA; confronto dei budget su celle normali. Se il problema persiste o migra a un'altra cella, raccogliere nuovamente `D19_BUDGET` e valutare la coerenza della costruzione del candidato nel limitatore e nel Lie executor. Non e dimostrata ancora la risoluzione fino a un run completo.

### D.19 Anticipo del sink di galleggiamento e condizione di neutralita (5 novembre 2026)

**Decisione scientifica dell'utente:** la miscela va eliminata quando `rho_mix <= rho_a_amb`, comprendendo uguaglianza e stati numericamente near-dry per cui `evaluate_density_from_qp` restituisce il valore atmosferico. Il criterio di rimozione appartiene a `state_physical`, mentre la posizione temporale del sink e responsabilita della transazione in `numerics`.

**Implementazione:** in `src/state/state_physical.f90` il predicato restituito da `evaluate_density_from_qp(...,below_ambient=)` e stato modificato da `rho_m<rho_a_amb` a `rho_m<=rho_a_amb`. In `src/numerics/lie_source_executor_2d.f90`, prima di ogni gruppo con ID diverso da zero, viene esaminato lo stato locale `q_local,qp_local` **gia accettato dal gruppo precedente**: se massa positiva e la densita soddisfa la condizione di rimozione, il candidato viene azzerato, `qp_local` viene inizializzato come dry e si interrompe la catena dei gruppi per quella cella, prima di generare un nuovo candidato soggetto a controllo di ammissibilita. Il gruppo 0 non e anticipatamente saltato, cosi rimane possibile l'iniezione in celle inizialmente dry. Il sink finale e mantenuto per verificare la rimozione anche dopo l'ultimo gruppo. Nessun candidato conservativo negativo viene convertito per calcolarne la densita. Le modifiche dei depositi o della topografia gia accettate restano nella transazione; non vengono anticipati/ripetuti i termini di exchange.

**Motivazione del caso TAAL:** il log `D19_BUDGET` mostra massa totale `5.7006955782118619E-16` e classe solida 1 `1.0280442421462131E-313` prima del gruppo 4; il suo aggiornamento produce `-4.9406564584124654E-324`. Con il nuovo ordine il controllo di buoyancy dello stato precedente gia valido puo eliminare la massa prima del consumo, qualora `qp_local` classifichi la miscela come neutra o galleggiante. Il criterio `<=` e esplicitamente voluto e differisce dal criterio storico `<`; modifica la fisica anche per miscele esattamente neutrali.

**Verifica ancora necessaria:** compilazione sul nuovo SHA, TAAL oltre 27.36 s e fino a fine simulazione, ETNA_LAVA, casi d'iniezione in dry, depositi e topografia in transazione, confronto mass budget e D19 diagnostics. Il precedente intervento subnormal-rate nel limitatore resta presente ma distinto; non viene rivendicata l'efficacia sul caso TAAL fino al test.

### D.20 Conferma utente: regressioni deterministiche complete superate (5 novembre 2026)

L'utente ha eseguito fino al **tempo finale** tutti e tre i casi `EXAMPLE_TALL`, `EXAMPLE_ETNA` ed `EXAMPLE_ETNA_LAVA`, senza errori, dopo il commit `de873faf8b2cc167c43f7aed36f4fdbd7681c241` che anticipa la rimozione delle miscele neutrali/galleggianti tra gruppi Lie. In precedenza aveva confermato risolto anche il blocco `D19` a `t≈27.36` del caso denominato `EXAMPLE_TAAL`. I nomi `TALL` e `TAAL` compaiono nelle comunicazioni dell'utente: verificarne la corrispondenza nei file di test se necessaria, senza presumere che identifichino esperimenti diversi.

**Stato:** regressioni deterministiche riferite dall'utente **PASS FINO A FINE RUN**; la validazione specifica con `stochastic_flag` / `stoch_transport_flag` attivi, test lake-at-rest, OpenMP riproducibilita e prestazioni restano **PENDING**. A08 resta revertito e non accettato; la conferma del run non ne dimostra l'equivalenza. Nessuna ulteriore modifica alla fisica introdotta in questo aggiornamento documentale.

### D.21 A08 reintrodotto dopo risoluzione D19 — validazione richiesta (5 novembre 2026)

**Decisione utente:** reintrodurre A08 ora che `EXAMPLE_TAAL/TALL`, `EXAMPLE_ETNA` e `EXAMPLE_ETNA_LAVA` sono stati eseguiti completamente con il fix di buoyancy `de873fa`.

**Modifica effettuata:** `src/spatial/hyperbolic_2d.f90` ripristinato **esattamente** dal commit originario A08 `3e2b02debd5c6cdc3c0b8e54d1c62ce07b4f9082`. Prima della sostituzione era identico al file baseline `c8605e429fee93440be1097779a34cfe0fe8d920`; non sono state perse altre modifiche spaziali. Il path idrostatico interno riusa `gamma_interface_xR/xL/yT/yB` gia calcolato sulle tracce di faccia, soltanto quando profondita del path **e** profondita selezionata sulla faccia sono entrambe maggiori di `EPSILON(1.0_wp)`. Negli altri casi mantiene la valutazione originale `eval_hydrostatic_coefficient` e la `regularize_dry_gamma_pair`. Nessuna nuova cache, nessun cambiamento agli endpoint, nessuna modifica a `state` o ai gruppi Lie.

**Rettifica storica:** il precedente blocco `D19` a `t≈27.36` non era attribuibile ad A08 in modo causale, dato che si ripresentava anche dopo il rollback del file spaziale. Il problema di TAAL e stato poi risolto con la rimozione anticipata delle miscele neutrali/galleggianti (`rho_mix<=rho_a_amb`), con conferma di run deterministici completi. A08 ora e **REINTRODOTTO / DA VALIDARE**, non ancora PASS.

**Test richiesti:** ricompilazione, completamento di `EXAMPLE_TAAL/TALL`, `EXAMPLE_ETNA`, `EXAMPLE_ETNA_LAVA`, confronto quantitativo con la baseline immediatamente precedente, e test lake-at-rest, wet/dry e well-balancing. Solo a seguire stimare il beneficio di performance (numero di valutazioni termodinamiche saltate, wall-clock e scaling). Lo stocastico rimane fuori dal gate di validazione immediata, come concordato.

### D.22 A08 confermato dall'utente: build e regressioni esempi OK (5 novembre 2026)

L'utente ha confermato **build e EXAMPLES OK** sul commit `365cc236e9344f79420b22bb9d388a17a68c7240`, che reintroduce il riuso di `Gamma` nei path idrostatici interni HP-PCCU. Nel contesto del test richiesto gli esempi sono `EXAMPLE_TAAL/TALL`, `EXAMPLE_ETNA` e `EXAMPLE_ETNA_LAVA`; questa conferma segue i precedenti run deterministici completi sulla baseline senza A08.

**Stato A08: PASS build + regressioni esempi (conferma utente)**. Restano aperte verifiche mirate lake-at-rest/well-balancing, confronto quantitativo dei campi e misure del beneficio prestazionale. Non confondere tali verifiche con il risultato positivo dei tre esempi. Validazione stocastica ancora rimandata.

### E. Evidenze richieste e stati dei blocchi

Per ciascun blocco riportare esplicitamente **DA FARE / IN CORSO / CODICE COMPLETATO / VERIFICATO**; separare l'evidenza statica dalla misura dinamica. Il commit architetturale `d8f4f57`, la rinomina `c9e576d` e la pulizia del workflow `c5d6eb2` non certificano un miglioramento prestazionale. Una Action Build verde certifica compilabilita sullo SHA effettivamente eseguito, non numero minimo di conversioni o parita scientifica.

Checklist finale:
- [ ] Inventario del cammino produttivo con frequenza effettiva delle conversioni e dei kernel per versione di cella/faccia e per iterazione Newton.
- [ ] Nessuna duplicazione dei cp/densita o altre proprieta nello stesso costruttore e nessuna copia allocabile non necessaria nel percorso caldo.
- [ ] Conversioni Lie motivate da una transizione effettiva di stato, con riuso di `qp` CFL e invalidazione corretta.
- [ ] Sorgenti e Newton riusano gli intermedi compatibili, mai uno scratch obsoleto; REAL/COMPLEX e complex-step restano coerenti.
- [ ] Proprietà minime di faccia calcolate sulla traccia finale e riusate, senza alterare HP-PCCU.
- [ ] Nessun nuovo accesso di `physics/` agli indici/campi globali, nessuna dipendenza ciclica o cache su griglia.
- [ ] Build, suite numerica/scientifica, OpenMP e confronto dei costi eseguiti **manualmente** e associati allo stesso SHA.
- [ ] Roadmap generale aggiornata con conclusione computazionale solo dopo tutte le evidenze.

**Prossimo intervento concreto: C0 + C1.** Prima raccogliere la mappa dei prodotti e degli input su cui si basa il riuso; poi rimuovere il ricalcolo della capacita termica nei costruttori da `q`, mantenendo identici REAL/COMPLEX e il comportamento numerico. Non riscrivere di nuovo la struttura di `src/state/`.

---

## 0. Sintesi delle decisioni vincolanti

- [ ] Unico stato fisico **generale, completo e coerente**, con versioni REAL e COMPLEX identiche per contenuto fisico; nessuno stato specializzato per reologia o modello ML.
- [ ] `q` resta la rappresentazione conservativa autorevole; `qp` resta la rappresentazione primitiva numerica; le proprietà derivate non sono una seconda fonte autorevole della griglia.
- [ ] Distinguere i *prodotti richiesti* (conversione numerica, proprietà di faccia, stato fisico completo), mai differenziare lo stato in base alla legge fisica.
- [ ] Operazioni in `physics/` strettamente **cell-local**; nessun indice `j,k`, accesso alla griglia o cache globale mutabile.
- [ ] Tutti i kernel numerici REAL/COMPLEX condividono le formule fisiche e preservano il complex-step; nessuna conversione COMPLEX che scarti perturbazioni.
- [ ] Le formule condivise REAL/COMPLEX hanno una sola implementazione sorgente: gli include `.inc` o un meccanismo equivalente sono mantenuti per evitare duplicazione manuale dei kernel tipati.
- [ ] Nessun wrapper o alias legacy solo per far compilare i chiamanti storici; migrazione diretta degli utilizzatori e rimozione del codice sostituito.
- [ ] Conservarne **quattro** interfacce per ogni termine sorgente: EXPLICIT, SEMI_IMPLICIT, IMPLICIT, ANALYTIC; i metodi non disponibili devono generare errore esplicito.
- [ ] Preservare Lie splitting, ammissibilità dello stato, conservazione, well-balanced HP-PCCU, gestione wet/dry e semantica delle sorgenti senza cambiare la fisica durante la migrazione.

**Criterio guida:** ottimizzare la quantità di lavoro svolto per valutazione, non semplicemente ridurre il numero di righe o di moduli.

## 1. Ambito, non-obiettivi e criteri di successo

Questa roadmap definisce il progetto di destinazione per la gestione degli stati, non la suddivisione dei file attuali né una sequenza di adattamento incrementale dell’architettura esistente. La struttura finale degli stati e delle relative API deriva da ruoli, dipendenze dei dati, costo per valutazione e contratti numerici. Comprende conversioni conservative/primitive, proprietà della miscela, costruzione e aggiornamento dello stato fisico, workspace e l’integrazione necessaria con flussi PCCU, termini sorgente e Newton/complex-step.

La riorganizzazione deve ottimizzare il lavoro complessivo dell’algoritmo di soluzione: ogni intermedio già calcolato deve poter essere passato al consumatore successivo, senza ricostruzioni implicite o rivalutazioni dello stesso contenuto. La semplicità richiesta è semplicità del grafo computazionale e dei contratti di stato, non conservazione dei moduli o dei nomi attuali.

Fuori ambito nella prima migrazione: modifica delle leggi di drag, delle reologie, della discretizzazione PCCU, del Lie splitting, delle soglie wet/dry e dei modelli di pressione di poro. Ogni eventuale revisione fisica va trattata in un progetto separato.

| Criterio | Definizione verificabile |
| --- | --- |
| Correttezza | Conservazione, limiti wet/dry, proprietà, Jacobiani e risultati di regressione equivalenti entro tolleranze giustificate. |
| Efficienza | Nessuna rivalutazione di un intermedio già disponibile nella stessa valutazione; allocazioni assenti nei loop caldi. |
| Architettura | Dipendenze unidirezionali, interfacce minimali, responsabilità singole, nessun wrapper pass-through. |
| Parallelismo | Kernel senza stato mutabile condiviso; scratch per worker o per cella; design device-friendly. |
| Estensibilità | Nuove reologie/ML leggono lo stato generale senza cambiare costruttori o introdurre campi specifici. |


## 2. Domanda reale dei consumatori

| Consumatore | Input richiesto | Output/proprietà richiesti | Vincolo dominante |
| --- | --- | --- | --- |
| Time integration / CFL | q oppure qp | qp, spessore, velocità, controllo ammissibilità | Costo per cella e semantica wet/dry |
| Ricostruzione e PCCU | Tracce primitive di faccia | Densità, g ridotta, gamma, energia, onde | Un solo calcolo per endpoint; well-balanced |
| Sorgenti espliciti/analitici | q o qp + geometria cella | Stato REAL completo | Condivisione fra processi e gruppi |
| Newton / complex-step | q perturbato + geometria | Stato COMPLEX completo | Coerenza delle derivate e assenza allocazioni |
| Output / inizializzazione | q, qp, composizione | Conversioni e diagnostica | API esplicite e controlli |

Non tutte le richieste devono transitare per un oggetto fisico completo: la numerica di faccia non deve pagare il costo di Richardson, curvature o quantità termodinamiche non consumate. Viceversa ogni kernel che restituisce un physical_state deve produrlo completo e coerente, senza profili “parziali” dipendenti dalla reologia.


## 3. Modello di stato e ownership


### 3.1 Stato conservativo autorevole

q è la sola fonte autorevole delle quantità conservate nel time integration: M = rho_m h, quantità di moto, energia termica M cp_m T, masse componenti e variabili trasportate. La rappresentazione qp è una vista numerica derivata o ricostruita, non una seconda fonte indipendente di verità. Le regole di sincronizzazione fra q e qp sono esplicite a ogni stage e durante le transazioni Lie.


### 3.2 Stato primitivo numerico

qp contiene h, u, v, T, frazioni massiche e variabili trasportate secondo il layout canonico. Non memorizza densità o proprietà reologiche ricostruibili; le tracce di faccia richiedono nuove proprietà valutate dalla composizione finale ricostruita, evitando di ricostruire indipendentemente quantità derivate.


### 3.3 Stato fisico generale

| Classe | Campi/semantica |
| --- | --- |
| Fondamentali | h, u, v, T; variabili trasportate Zs e pressione di poro in eccesso. |
| Miscela e carrier | rho_m, inv_rhom, rho_c, inv_rho_c, cp_m, cp_c. |
| Composizione | xs, xg, alphas, alphag, incluse frazioni zero quando la fase è assente. |
| Galleggiamento | red_grav e Ri (con soglia e saturazione originali). |
| Cinematica / topografia | w, mod_vel_hor, mod_vel2, mod_vel, curvature_acceleration. |

I campi sono proprietà generali. Viscosità specifiche, tensione di snervamento, coefficienti di drag, velocità di settling per specie e rate dei processi non appartengono allo stato comune; appartengono ai moduli costitutivi. REAL e COMPLEX hanno parità di campi e semantica. Nessun modello specifico decide quali campi vengano inizializzati.

Gli array `xs`, `xg`, `alphas` e `alphag` appartengono allo stato fisico e sono la rappresentazione persistente delle composizioni e delle frazioni volumetriche della valutazione corrente. Lo stato è proprietario della loro memoria e il costruttore dello stato li scrive direttamente; non devono esistere copie equivalenti con la stessa responsabilità nel workspace. La loro dimensione è determinata dalla configurazione del modello e viene allocata o verificata al momento dell’inizializzazione del workspace, mai durante la valutazione di una cella.

Il workspace possiede uno stato REAL o COMPLEX riutilizzabile e ne garantisce il ciclo di vita. Non è una seconda rappresentazione dello stato: contiene soltanto il risultato corrente e, se necessario, scratch distinti con una responsabilità computazionale esplicita. Un eventuale scratch aggiuntivo per composizione può essere introdotto solo se evita una copia o una rivalutazione misurabile, con ownership, durata, dimensionamento e percorso di aggiornamento documentati. In assenza di tale necessità, le formule devono scrivere direttamente negli array dello stato.


### 3.4 Coerenza, immutabilità e invalidazione

Il tipo fisico viene costruito completo a partire da un input valido. Le routine costitutive lo ricevono preferibilmente con INTENT(IN). Se una legge analitica produce uno stato successivo con u, v, T o composizione modificati, la modifica deve avvenire mediante un costruttore/operazione di aggiornamento che ricomputi tutte le dipendenze, oppure ricostruendo lo stato. Una semplice scrittura diretta a state%u lasciando Ri e mod_vel invariati è vietata.

In particolare: variazione di u/v → aggiornare w, mod_vel_hor, mod_vel2, mod_vel, curvature_acceleration, Ri; variazione di T o composizione → aggiornare cp, rho, frazioni volumetriche, red_grav e Ri; variazione di h → aggiornare Ri e ogni derivata che dipende dalla profondità; variazione della topografia → aggiornare solo quantità cinematiche geometriche. La topografia non va assunta costante tra gli stage: i dati geometrici vengono passati esplicitamente alla valutazione.

La proprietà degli array di composizione non cambia durante una valutazione: il workspace mantiene la stessa memoria tra chiamate e ogni nuovo stato sovrascrive tutti gli elementi pertinenti. Non sono ammesse copie implicite fra array di workspace e array dello stato, né assegnamenti di derived type che producano deep copy nel percorso caldo.


## 4. Kernel fondamentali e contratti

La struttura logica è composta da quattro kernel principali, più helper matematici necessari e distinti. Le firme riportate sono illustrative: non devono essere interpretate come obbligo di mantenere gli argomenti esatti dell’implementazione corrente.

| Kernel | Contratto | Lavoro escluso |
| --- | --- | --- |
| conservative_to_primitive | q → qp + esito ammissibilità + pressione dinamica opzionale | Ri, curvature, calori del carrier non necessari, stato fisico completo |
| primitive_to_conservative | qp → q con chiusura e proprietà minime | Stato fisico, Ri, curvature |
| evaluate_mixture | T + composizione → densità, capacità termiche, frazioni volumetriche e inversi | Topografia, velocità, reologia |
| build_physical_state_from_q / build_physical_state_from_qp | q o qp + geometria → physical_state completo REAL/COMPLEX | Accesso a celle adiacenti, coordinamento dei gruppi Lie |


### 4.1 Conversione q → qp: pipeline a passaggio singolo

1. Esaminare ammissibilità prima delle divisioni; trattare M ≤ soglia con la convenzione dry esistente.
1. Calcolare inv_M una volta, e da esso xs, xg e scalari trasportati; applicare esattamente la chiusura con la semantica di accettazione/fallimento.
1. Calcolare cp_m una volta; ricavare T = E/(M cp_m) rispettando fallback di temperatura e soglia eps_sing.
1. Calcolare rho_m/inv_rhom tramite la formula canonica, ma non costruire alphas/alphag se l’implementazione specializzata può evitarlo senza duplicare la fisica.
1. Calcolare h = M/rho_m e velocità usando la regolarizzazione corrente; produrre qp e l’eventuale p_dyn. Non costruire un physical_state.
Le duplicazioni oggi presenti nella validazione di qc_to_qp e nella decodifica r_phys_var vanno eliminate, ma senza anticipare normalizzazioni che cambierebbero energia o massa. Le formule e l’ordine delle operazioni devono essere confrontati con l’oracolo attuale.


### 4.2 Conversione qp → q

Applicare la chiusura delle frazioni una sola volta, calcolare la composizione carrier, cp_m e rho_m, quindi costruire M, momenti ed energia con intermedi riutilizzati. Preservare l’interpretazione di energia puramente termica (senza energia cinetica) e i casi con fase assente. Le variabili trasportate sono moltiplicate per M una sola volta. Le variabili geometriche non entrano nel kernel.


### 4.3 Valutazione termodinamica

Un unico punto di definizione per le formule della miscela, con specializzazioni REAL e COMPLEX coerenti. Le formule comuni devono essere scritte una sola volta, preferibilmente negli `.inc` già adottati dal progetto, e incluse da wrapper tipati minimi REAL/COMPLEX. Non duplicare manualmente i kernel in due procedure parallele. Evitare kernel “universali” che calcolino tutte le proprietà su ogni faccia: prevedere una valutazione compatta per density/energy/pressure e una valutazione completa usata dal costruttore dello stato. Le due viste condividono le formule elementari e gli intermedi, non due implementazioni indipendenti. L’assenza del carrier o dei gas aggiuntivi deve avere una semantica definita, compresi i limiti di frazione totale e le configurazioni liquide/gassose.


### 4.4 Costruttori dello stato fisico

Da qp: leggere direttamente le primitive; chiudere la composizione una volta; valutare le proprietà termodinamiche; calcolare red_grav, Richardson e cinematica dalla geometria corrente. Da q: decodificare una volta composizione/energia/velocità e usare gli stessi intermedi; evitare il percorso q → qp → physical_state che rifarebbe la chiusura. Il COMPLEX usa lo stesso grafo di calcolo e conserva le perturbazioni immaginarie secondo i criteri di branching sulle parti reali, senza operazioni non analitiche aggiuntive.


## 5. Valutazione selettiva: dove ammessa e dove vietata

La distinzione ammessa riguarda il prodotto richiesto, non il modello sorgente. I percorsi numerici producono primitive o proprietà di faccia senza creare uno stato fisico generale; ogni percorso che produce physical_state deve invece inizializzare tutti i campi, in REAL e COMPLEX. Non introdurre flag “requires_cp”, “is_MuI”, “state_partial” o registri per proprietà lette dalle singole reologie. Una reologia futura deve potersi collegare senza modificare la conversione degli stati.

Per Newton il guadagno viene da una decodifica COMPLEX compatta, scratch preallocato e riuso degli intermedi nella singola valutazione, non da stati incompleti per modello. La versione COMPLEX deve eseguire le stesse formule della versione REAL tramite gli stessi include `.inc` o sorgenti condivisi, preservando le perturbazioni immaginarie. L’eventuale scelta di valutare un intero gruppo sorgente dipende dal registro fisico, ma non cambia il contratto dello stato che il gruppo riceve.


## 6. Layout dati, memoria e parallelismo


### 6.1 Scelta di layout

Gli array di griglia q(component,x,y) e qp(component,x,y) restano orientati alla coalescenza delle componenti per cella, coerentemente con l’accesso corrente e con i chiamanti. Una modifica globale del layout richiederebbe benchmark e una fase distinta. Lo stato fisico è un oggetto cell-local transitorio, non un array allocato su tutta la griglia.


### 6.2 Workspace e allocazioni

Allocazione degli array di composizione appartenenti allo stato e degli eventuali scratch distinti al setup del worker o della configurazione; nessun ciclo ripetuto ALLOCATED/SIZE/DEALLOCATE durante i loop caldi. Il workspace è proprietario dello stato riutilizzabile e ne garantisce dimensione, inizializzazione e durata. La costruzione accetta workspace già dimensionati e scrive direttamente negli array dello stato; il controllo del layout viene eseguito al confine, una sola volta per configurazione. Evitare assegnamenti di derived types con componenti allocabili quando producono deep copy non necessarie. Documentare per ogni buffer ownership, durata, accesso e motivo della sua esistenza.


### 6.3 OpenMP / MPI / GPU

- OpenMP: nessuna scrittura condivisa nei kernel; uno scratch per worker o oggetti locali privati, senza lock né allocazioni ripetute.
- MPI: halo exchange e distribuzione della griglia a livello spaziale; il kernel riceve soltanto dati locali, mai j/k o riferimenti al dominio.
- GPU: separare host setup e device evaluation; niente ALLOCATE, I/O o fatal_error dentro il kernel offload. In futuro esplorare array a dimensione fissa/layout specializzato, senza imporre oggi un numero di specie hard-coded.
- Predisporre interfacce che possano essere compilate senza dipendere da OpenMP o runtime esterni all’interno della fisica, e verificare vectorization/inlining con report del compilatore.

## 7. Integrazione con la numerica spaziale e PCCU

Le proprietà di faccia devono essere calcolate dalla traccia finale dopo ricostruzione, clipping, HP correction e ogni selezione wet/dry. Non ricostruire rho_m, red_grav o gamma come variabili indipendenti se questo viola il bilanciamento idrostatico. Un endpoint workspace piccolo contiene almeno rho_m, red_grav, gamma = rho_m red_grav e cp_m se serve la conversione energetica.

1. Identificare l’esatta sequenza di costruzione delle tracce PCCU e dei loro stati di riferimento.
1. Preparare una sola volta le proprietà termodinamiche per ciascun endpoint effettivamente utilizzato.
1. Passare gli intermedi direttamente a velocità caratteristiche, coefficiente idrostatico e conversione in conservative; evitare ulteriori chiamate a mixt_var e qp_to_qc che ripetano le proprietà.
1. Preservare ordinamento e formule del path non conservativo, dry/wet handling, coefficienti geometrici e well-balanced; non anticipare il calcolo a monte di un limiter.
1. Confrontare a livello di singola faccia flussi, gamma, lambda min/max e contributi del path prima dei benchmark globali.

## 8. Integrazione con sorgenti, Lie e Newton


### 8.1 Termini sorgente

La numerica costruisce lo stato a partire da q o qp locali e dalla geometria disponibile al gruppo; i moduli in physics ricevono lo stato pronto con INTENT(IN). Nessuna legge sorgente ricostruisce densità, calori specifici o Richardson da input già presenti. Rate e aggiornamento dello stato restano responsabilità dei processi e dell’integratore rispettivamente. Le quattro interfacce di trattamento restano obbligatorie, con errore esplicito per quelle non implementate.


### 8.2 Lie splitting

Non conservare uno stato fisico tra gruppi che modificano q: il gruppo successivo deve ricevere lo stato ricostruito a partire dal q autorevole aggiornato. Il riuso è valido solo per più processi sullo stesso input congelato, non fra stadi temporali con contenuto differente. Nessun limiter dei rate dopo l’avanzamento puramente spaziale: verificare esclusivamente l’ammissibilità dello stato come nel contratto Lie stabilito.


### 8.3 Newton / complex-step

Il residuo e ciascuna perturbazione COMPLEX devono costruire uno stato coerente da q perturbato, in maniera deterministica e senza allocazioni. Conservare le condizioni dry e i branch sulla parte reale già adottati quando necessari alla differenziabilità complex-step. Le uscite REAL e COMPLEX devono coincidere sullo stato base entro tolleranze strette. Il Jacobiano complex-step è il test più sensibile della correttezza del nuovo percorso.


### 8.4 Reologie analitiche sequenziali

Le routine che modificano velocità dopo ogni legge devono usare un aggiornamento coerente di tutti i campi dipendenti. Il nuovo contratto proibisce modifiche dirette a u/v se poi lo stato viene riutilizzato senza aggiornare Ri, w, mod_vel e curvature. Preferire un aggiornamento kinematics_from_velocity completo o ricostruzione dello stato, preservando il fatto che T/composizione e geometria possano rimanere invariati.


## 9. Dipendenze logiche e organizzazione del codice

I nomi dei file si decidono solo dopo la pulizia. La seguente è una gerarchia logica di responsabilità; può tradursi in pochi moduli con dipendenze senza cicli, senza duplicare le formule e senza wrapper di compatibilità.

| Livello | Possesso logico | Dipendenze consentite |
| --- | --- | --- |
| A. Layout e controlli | Indici canonici, soglie, validazione del modello, chiusura composizione | Configurazione e diagnostica |
| B. Kernel materiali | Densità, cp, proprietà delle fasi, helper matematici REAL/COMPLEX | A + parametri fisici |
| C. Conversioni | q ↔ qp, proprietà di faccia, costruzione fisica completa | A + B + geometria per il solo costruttore |
| D. Consumatori | Spaziale, Lie, Newton, processi fisici, I/O | API di C, mai implementazioni private |

Non mantenere pass-through o alias di compatibilità per conservare l’organizzazione storica. Due procedure possono coesistere soltanto quando rappresentano contratti realmente distinti, per esempio conversione con verifica di ammissibilità versus conversione garantita per input già validato. La differenza deve essere esplicita nel contratto, misurabile nei costi e testata; non deve dipendere dal nome precedente della routine.


## 10. Piano operativo e deliverable

**Legenda stato:** `[ ]` non iniziato, `[~]` iniziato/in corso, `[x]` completato con evidenza verificabile. Ogni fase registra commit SHA, test, benchmark e rischi aperti. Un commit intermedio può non compilare; nessuna fase è verificata solo perché la modifica è stata committata.

- [~] **F0 — Baseline e grafo dei consumatori.** Inventariare i prodotti realmente richiesti dall’algoritmo, i campi utilizzati, la frequenza, gli intermedi riutilizzabili e gli invarianti; salvare riferimenti CI e risultati numerici sul commit iniziale. L’inventario non deve imporre la struttura finale dei moduli.
- [~] **F1 — Contratti e invarianti di destinazione.** Progettare ownership, semantica dry/wet, composizione, energia, velocità, Richardson, geometria e REAL/COMPLEX a partire dal grafo computazionale; definire API nuove, senza alias storici, e il confine degli include condivisi.
- [~] **F2 — Termodinamica condivisa.** Realizzare un solo punto di definizione delle formule e kernel mirati (minimo per la numerica; completo per i sorgenti); verificare casi solidi/gas/liquidi e chiusure. Implementati `evaluate_mixture` REAL/COMPLEX e gli include condivisi; la suite isolata verifica il caso gas/solido, mentre restano casi completi e audit delle chiusure.
- [~] **F3 — Conversioni conservative/primitive.** Riscrivere `q → qp` e `qp → q` in passaggio singolo con controllo di ammissibilità e intermedi riutilizzati; testare near-dry e conservazione. Implementate le quattro procedure tipate; la suite isolata verifica il round-trip wet, mentre restano near-dry, stati inadmissibili e confronto esteso col backend attuale.
- [~] **F4 — Stato fisico generale REAL/COMPLEX.** Implementare costruzione completa da `q` e `qp`, conversioni e aggiornamento coerente dopo mutazioni; preparare scratch configurati una volta. Implementati i tipi, l’inizializzazione dei workspace, i costruttori `build_physical_state_from_q` e `build_physical_state_from_qp` REAL/COMPLEX e `update_physical_state_kinematics`; la suite verifica parità di densità, Richardson, cinematica e dry state, mentre restano mutazioni e audit complesso completo.
- [~] **F5 — Migrazione dei sorgenti.** Reologie, mass exchange, Darcy, termica, iniezioni, dispatcher fisico e scratch Lie/Newton utilizzano i tipi fisici canonici `state_new` REAL/COMPLEX; i contributi impliciti condivisi usano `.inc` anche nell'assemblaggio esplicito REAL. I costruttori analitici e lo scratch dei gruppi espliciti combinati sono riutilizzati senza nuove validazioni locali. Restano: audit delle conversioni legacy ancora nei chiamanti numerici, test di equivalenza dei rate (incluse normalizzazione/overfill), regressioni Newton/complex-step, CI sul nuovo commit e verifica dei costi hot-path. Il riassetto del rollback Lie e differito a dopo F5.

**F5 test-client migration (5 novembre 2026):** test reologie cell-local, Froude, Darcy, stress mass exchange, D17 e source registry aggiornati ai tipi canonici REAL/COMPLEX. I test di costruzione esplicita inizializzano gli array una sola volta; il test mass exchange confronta la promozione numerica dei campi scalari senza dipendere dalla vecchia routine di estrazione. Esecuzione CI ancora richiesta.

**F5 validation workflow (5 novembre 2026):** `state-architecture-validation.yml` estende la matrice dei test esistenti ai runner `TEST_RHEOLOGY_CELL_LOCAL` (incluso il nuovo test di equivalenza REAL/COMPLEX e Jacobiani), `TEST_RHEOLOGY_FROUDE` e `TEST_SOURCE_REGISTRY`, oltre alle regressioni precedenti e ai gate state_new/old-new. La workflow e' manuale (`workflow_dispatch`), non e' stata eseguita su questo SHA. Accettazione: tutti i job verdi sullo stesso commit, log conservati, quindi regressioni Lie/OpenMP e confronto numerico multi-caso prima di chiudere F5.

**F5 test shared implicit (5 novembre 2026):** aggiunto `TESTS/TEST_RHEOLOGY_CELL_LOCAL/test_shared_implicit_kernels.f90` al runner strict/optimized: sette modelli, verifica EXPLICIT = SEMI_REAL + REAL(IMPLICIT_COMPLEX), Jacobiani complex-step contro differenze finite del contributo `EXPLICIT-SEMI` e stato a velocita zero. Test scritto e collegato al runner; esecuzione CI ancora PENDENTE.

**F5 rheology analytic dispatch (5 novembre 2026):** l'esecutore analitico del singolo modello valuta ora direttamente lo stato REAL congelato tramite un dispatcher di legge unico, evitando la copia profonda di composizioni `current_state=state`. La copia resta solo nel percorso multi-modello sequenziale, che aggiorna cinematiche/Richardson dopo ogni sottopasso. Contratti e formule dei modelli non modificati; test numerici ancora necessari.

**F5 Lie analytic scratch (5 novembre 2026):** l'API analitica accetta lo scratch REAL di worker gia inizializzato. L'esecutore Lie lo riusa per tutte le valutazioni analitiche locali evitando allocate/deallocate di quattro array nella fase analitica; chi usa l'API senza scratch ha un fallback esplicito. La trasformazione dell'operatore analitico in rate conservativo e invariata. Test e benchmark da eseguire.

**F5 build dependency migration (5 novembre 2026):** aggiornate le dipendenze esplicite dei moduli Fortran che importano i nuovi tipi e i costruttori (reologie, termica, Darcy, exchange, injection, source assembly, Lie/Newton) per evitare race di compilazione parallela sulla generazione dei file `.mod`. I test shell che compilano in ordine di `src/Makefile.am` sfruttano l'ordine esistente; i nuovi `.inc` sono cercati nella cartella dei sorgenti dai compilatori supportati. Rimane obbligatoria la CI di integrazione.

**F5 consumer migration (5 novembre 2026):** termica, Darcy, mass exchange e injection usano ora direttamente i tipi `state_new`. `source_terms_2d` costruisce lo stato dai costruttori canonici q/qp; Lie e Newton hanno scratch REAL/COMPLEX di questi tipi. La memoria dei worker viene preparata prima del loop OpenMP, con fallback locali esplicitamente inizializzati. Il modulo `state_new_types` espone inizializzatori diretti dei tipi fisici e ne mantiene gli initializer workspace deleganti. **Non verificato:** compilazione, Action, regressioni numeriche, Jacobiani complex-step, esiti `accepted` dei costruttori, e costi di allocazione nei fallback analitici. F5 resta IN CORSO.

**Nota sul registro storico:** i quattro paragrafi seguenti (rifiuto Newton, propagazione `accepted`, chiusura composizionale e relativi test) descrivono change set successivamente **annullati** dal commit `9952cd2`; non rappresentano le API correnti. Vengono conservati qui unicamente come traccia di audit.

**F5 rifiuto dei trial Newton (5 novembre 2026):** il dispatcher `eval_source_implicit` ora espone `accepted` opzionale e restituisce un rifiuto senza chiamare `fatal_error` su stati fisici conservativi non ammissibili. Il risolutore lo propaga con `eval_f` e `eval_jacobian`: per un punto iniziale o Jacobiano non ammissibile termina il gruppo con `linear_info=-9`; il primo trial full-step fallito porta al rollback, mentre la line search dei passi successivi dimezza il passo prima di riprovare. Non si accetta un rate zero come soluzione di un trial rifiutato. Commit `a93196b`, `0789860`; test CI/Newton dinamici ancora DA ESEGUIRE. F5 IN CORSO.

**F5 propagazione ammissibilita nel dispatcher sorgenti (5 novembre 2026):** `source_terms_2d` verifica l'esito `accepted` dei costruttori REAL/COMPLEX nei percorsi esplicito, mass exchange, semi-implicito e Newton; uno stato non ammissibile non viene passato alle leggi fisiche. Il percorso analitico restituisce `status=-9` e costruisce lo stato una sola volta prima del dispatch, preservando il no-op Darcy su massa nulla; rimossa l'annotazione `RECURSIVE` non necessaria. Negli altri percorsi, che non espongono uno stato di rifiuto, l'errore e esplicito. Resta da evolvere il contratto dei rifiuti Newton per gestire eventuali trial non ammissibili senza terminazione e da validare con CI/test numerici. F5 IN CORSO.

**F5 chiusura al limite e complex-step (5 novembre 2026):** `close_composition_complex` applica la stessa soglia di annullamento del valore REAL sulle frazioni prossime a zero, preservandone la parte immaginaria e quindi la sensibilita complex-step. `TEST_STATE_NEW` ora verifica anche una specie conservativa negativa entro il roundoff e una perturbazione immaginaria non nulla. I test della nuova serie di commit sono ancora da lanciare in CI.

**F5 chiusura composizionale conservativa (5 novembre 2026):** il costruttore `build_physical_state_from_q` REAL applica ora `close_composition` prima dei kernel di proprieta, come gia avviene nel corrispondente COMPLEX. Una composizione non ammissibile restituisce `accepted=.FALSE.`; un sovrariempimento entro la tolleranza viene chiuso senza divergenza REAL/COMPLEX. `TEST_STATE_NEW` aggiunge tre casi conservativi (sovrariempimento da arrotondamento, sovrariempimento reale, specie negativa), con confronto di densita e calore. Test CI del nuovo commit ancora da eseguire; F5 IN CORSO.

**F5 thermal shared REAL/COMPLEX rate (5 novembre 2026):** `thermal_cooling_2d` usa un solo `thermal_cooling_rate.inc` per il contributo implicito sia REAL (assemblaggio esplicito) sia COMPLEX (Newton). L'esplicito somma il semi-implicito REAL e il contributo implicito REAL senza promozione di stato; aggiunti test di parità REAL/COMPLEX e derivata complex-step in `test_lava_cooling.f90`. Registrata la dipendenza `EXTRA_DIST`; compilazione e Action del nuovo commit ancora da verificare. F5 resta IN CORSO.

**F5 thermal analytic (5 novembre 2026):** il kernel esatto di raffreddamento riceve ora direttamente la temperatura iniziale (unico dato dello stato consumato). Il dispatcher avanza lo scalare fra modelli senza `current_state=state` e senza copia profonda degli array di composizione. Aggiornato `TEST_RHEOLOGY_CELL_LOCAL/test_lava_cooling.f90` per il nuovo contratto. Restano da eseguire i test CI e la migrazione dei tipi termici REAL/COMPLEX e dei chiamanti Lie/Newton alla nuova API; F5 non e conclusa.
- [~] **F6 — Numerica spaziale e facce HP-PCCU.** Migrati termodinamica e conversione conservativa delle tracce, chiusura numerica condivisa, e conversione conservativo/primitivo nel CFL. Restano riuso delle proprieta tra endpoint, test scientifici HP-PCCU e profiling; non modificare Lie/rollback durante F6.
- [~] **F7 — Eliminazione delle API precedenti.** Rimossi 26 prerequisiti legacy nel grafo Make; migrate conversioni del driver diagnostico e dell'inizializzazione. Restano integrazione Lie/time-step e I/O, eliminazione API legacy realmente inutilizzate, audit dei cicli.
- [ ] **F8 — Gate numerici e scientifici.** Strict-debug, suite di conversione, complex-step/Newton, Lie/rollback, well-balanced, test per modelli e confronti di output.
- [ ] **F9 — Profiling, scaling e chiusura.** Benchmark ripetuti 1/2/4/8 thread dove disponibili, profili granulari per valutazione e per intermedio, confronto multi-esempio e aggiornamento di documenti e roadmap generale.

### Registro operativo sintetico

| Fase | Stato | Evidenza (SHA / Action / test) | Blocco aperto |
| --- | --- | --- | --- |
| F0 | IN CORSO | Analisi consumatori e roadmap architetturale; baseline numerica da consolidare | Grafo completo dei consumatori, riferimenti CI e risultati numerici iniziali |
| F1 | IN CORSO | Contratti, ownership degli array e API concettuali documentati; tipi/workspace implementati | Invarianti formali e contratto completo del workspace da chiudere |
| F2 | IN CORSO | `src/state_new/state_new_mixture.f90`; `mixture_heat_capacity.inc`; `mixture_from_mass.inc`; `tests/TEST_STATE_NEW`; compilazione e caratterizzazione isolate superate | Casi liquido/gas/solido, chiusure limite e confronto con il backend esistente |
| F3 | IN CORSO | `src/state_new/state_new_conversions.f90`; round-trip wet in `tests/TEST_STATE_NEW` superato | Near-dry, ammissibilità e audit delle operazioni ripetute |
| F4 | IN CORSO | `src/state_new/state_new_types.f90`, `src/state_new/state_new_physical.f90`, `kinematics.inc`, `richardson.inc`; test REAL/COMPLEX e dry state superati | Mutazioni coerenti, copertura COMPLEX perturbata e integrazione |
| F5 | IN CORSO | Reologie, exchange, iniezioni, termica, Darcy e dispatcher sui tipi `state_new`; scratch per worker, `.inc` REAL/COMPLEX, test di equivalenza; baseline CI verde `9952cd2` | Validare i commit successivi, audit conversioni legacy nei chiamanti, regressioni Newton/complex-step e costi hot-path. Rollback globale: DOPO F5 |
| F6 | IN CORSO | Proprieta idrostatiche, facce PCCU e CFL migrati a `state_new` con mode numerico legacy | Validare nuovo commit, ottimizzare riuso di endpoint, regressioni HP-PCCU e profiling |
| F7 | IN CORSO | Eliminati 26 archi Make legacy; driver diagnostico, CFL e inizializzazione migrati | Migrare integratore e I/O nel change set dedicato, audit delle API effettivamente inutilizzate, validare CI |
| F8 | DA FARE | — | Prove scientifiche e Actions sul commit finale |
| F9 | DA FARE | — | Benchmark e chiusura documentata |


**Regola di esecuzione:** l'ordine F0–F9 è una sequenza logica di responsabilità; F2/F3/F4 possono essere sviluppate in un unico change set consistente e validato alla sua conclusione. Non introdurre wrapper per mantenere compilabili i commit intermedi.

### Registro di migrazione F5 — primo change set (5 novembre 2026)

**Contratti consolidati in questo change set:**
- [x] Il setup valida le dimensioni di **tutti e quattro** gli array `xs`, `xg`, `alphas`, `alphag`, sia REAL sia COMPLEX.
- [x] `update_physical_state_kinematics` aggiorna anche Richardson quando cambiano `u,v`; i costruttori e `refresh_physical_state` non ne duplicano il calcolo.
- [x] `TEST_STATE_NEW` include un test dell'aggiornamento di velocita e della parita REAL/COMPLEX del Richardson.

**F5 — estensione a empirical/lahar (5 novembre 2026):** anche i contributi impliciti empirical e lahar sono compilati direttamente da un singolo `.inc` ciascuno in REAL (assemblaggio esplicito) e COMPLEX (Newton); nessuna conversione dello stato e nessuna formula duplicata. La viscosita lahar utilizza la funzione generica REAL/COMPLEX gia disponibile. Complessivamente sette reologie migrate con questa strategia. Test e CI ancora DA ESEGUIRE.

**F5 — formule implicite condivise REAL/COMPLEX (5 novembre 2026):** per Voellmy, mu(I)-Voellmy, lava, plastic e quadratic le interfacce esplicite sommano ora il contributo semi-implicito REAL e il contributo della legge implicita valutata in REAL; le interfacce Newton COMPLEX includono la medesima formula tramite quattro file `.inc`. Voellmy e mu(I)-Voellmy condividono lo stesso kernel turbolento. Eliminate le promozioni REAL→COMPLEX per queste cinque leggi, mantenendo quattro interfacce e semantica delle soglie. Verifica statica del change set: effettuata; compilazione/CI/equivalenza numerica/complex-step: DA ESEGUIRE. Gli altri modelli e la migrazione dei chiamanti restano aperti.

**F5 — cambio contratti delle reologie (5 novembre 2026):** tutti i 12 moduli reologici (dispatcher e leggi) ora dichiarano direttamente `physical_state_real`/`physical_state_complex` da `state_new_types`. La promozione REAL → COMPLEX usata dalle valutazioni esplicite composte è un'operazione numerica esplicita in `state_new_types`, con buffer COMPLEX preinizializzato: non è un forwarding wrapper. Il vecchio costruttore del `state_conversion_2d` non viene più usato dalle reologie. La nuova promozione richiede aggiornamento dei chiamanti e della gestione degli scratch; compilazione e test restano NON VERIFICATI fino alla migrazione coordinata di `source_terms_2d` e dell'esecutore Lie/Newton. Non marcare F5 completata.

**Inventario storico F5 (snapshot `1ca1fe6`, SUPERATO):** i checkbox seguenti registrano lo stato di un audit precedente e NON descrivono il master attuale. La migrazione delle leggi fisiche a `state_new` e dello scratch Lie/Newton e gia implementata; l'audit residuo riguarda conversioni numeriche e validazione finale, non la reintroduzione di `accepted` nei dispatcher.
- [ ] Reologie: dispatcher `rheology_2d`, Voellmy, plastic, lava, lahar, empirical, quadratic, Froude, Lucas, mu(I), mu(I)-Voellmy, newmodel. Migrare i tipi, trasferimenti REAL/COMPLEX e la copia profonda `current_state=state` nel percorso analitico.
- [ ] Mass exchange: `mass_exchange_2d`, deposition, erosion, entrainment, vaporization, liquid loss e lava deposition. Preservare budget, rate e processi locali.
- [~] Darcy, termica e iniezioni: tre iniezioni volumetriche (circolare, fissurale, radiale) migrano ai kernel mirati `evaluate_mixture_density`/`evaluate_mixture_heat_capacity` senza costruire `qp_source`; resta migrare i tipi di stato, Darcy, termica e altri consumatori.
- [ ] `source_terms_2d`: unificare la costruzione dello stato per valutazione e propagare `accepted` senza introdurre wrapper storici.
- [ ] `lie_source_executor_2d` e `nonlinear_solver_2d`: inizializzare workspace `state_new` per worker, preservare Lie, rollback, integrita dei rate, Newton e complex-step.
- [ ] F6/F7 separati: `hyperbolic_2d`, `reconstruction_2d`, `spatial_operator_2d`, `equation_terms_2d`, `time_integration_2d`, `init_2d`, `inpout_2d` ancora importano le vecchie API.

**Vincolo F6:** gli operatori in `spatial/` devono restare indipendenti dalle formule costitutive. Proprietà di faccia solo dalla traccia finale; `pccu_2d` continua a ricevere coefficienti gia valutati, senza `USE` di moduli fisici.

**Secondo change set F5 (iniezioni e build):** sostituita `mixt_var` nelle tre sorgenti volumetriche con due kernel della nuova miscela (densità/calore) e gravità ridotta locale; rimossi i `qp_source` artificiali. La composizione del materiale sorgente viene chiusa su piccoli array automatici locali con `state_new_composition`, mantenendo il controllo di ammissibilità prima dei kernel; per composizioni sovrariempite oltre la tolleranza la nuova API rifiuta anziché rinormalizzare come la legacy (caso limite da caratterizzare nei test). `src/Makefile.am` elenca i sei moduli `state_new`, gli include e le dipendenze; `TEST_STATE_EQUIVALENCE` confronta ora proprietà del materiale iniettato con il backend storico per gas/solido, liquido/solido e casi limite. Le quattro interfacce e i rate restano invariati. **Non ancora migrati** i tipi fisici in `physics/`, `source_terms_2d`, Lie e Newton; validare build e test prima del prossimo blocco.

**Validazione:** l'ultima Action verde precedente e `State architecture validation` su `3df7017`; il test di equivalenza introdotto successivamente va rieseguito e ampliato. Il nuovo change set modifica solo `state_new`, test e roadmap; la migrazione dei consumatori di produzione resta aperta.

## 11. Matrice dei test e accettazione

| Area | Casi obbligatori | Esito richiesto |
| --- | --- | --- |
| Stato e composizione | Secco, near-dry, massa nulla/negativa, composizione sovrariempita, componenti zero/negative per roundoff. | Semantica identicale; errori espliciti sugli stati inadmissibili. |
| REAL vs COMPLEX | Base reale, perturbazioni sulle masse, momenti, temperatura ed energia, specie e scalari. | Parte reale equivalente; derivate complex-step coerenti. |
| Conservazione | Round-trip q↔qp, massa, specie, energia termica e quantità di moto. | Tolleranze motivate, senza regressioni sistematiche. |
| PCCU / HP | Lake-at-rest, wet/dry, slope, roundoff idrostatico, tracce corrette. | Well-balanced e flussi invariati entro tolleranze strette. |
| Sorgenti | Reologia esplicita/analitica/Newton, sedimentazione, erosione, entrainment, Darcy e gas loss. | Rate e stati finali coerenti con baseline. |
| Parallelismo | 1/2/4/8 thread, seed deterministico, lifecycle scratch, sanitizers e memcheck. | Nessuna race o allocazione nel percorso caldo. |
| Prestazioni | Microbenchmark q→qp, qp→q, build REAL/COMPLEX, facce e TAAL. | Costo quantificato per sottosistema e nel tempo totale. |

È essenziale non dichiarare “bitwise equivalent” come obiettivo universale se cambia l’ordine delle operazioni floating-point: occorrono test bitwise dove realisticamente preservabile e tolleranze fisicamente/numericamente motivate altrove. Registrare mass error, energia, runout, campi per cella, fallimenti di Newton e iterazioni, non solo il tempo totale.


## 12. Strategia di profiling e decisioni prestazionali

Misurare non solo il wall time, ma anche numero di valutazioni e costo medio per valutazione. Profilare: `conservative_to_primitive`, `primitive_to_conservative`, `build_physical_state_from_q`, `build_physical_state_from_qp`, `evaluate_mixture`, stato COMPLEX per Jacobiano, proprietà di faccia e source group. Distinguere costo delle formule, controlli, copie e allocazioni. Benchmarkare configurazioni gas/solid, liquid/solid, gas aggiuntivi, differenti numeri di specie e regimi wet/dry.

Prima ipotesi qualitativa: il riuso sulle facce PCCU può offrire il beneficio maggiore; la conversione q→qp è una seconda opportunità; la costruzione fisica COMPLEX pesa soprattutto nei casi Newton. Non attribuire percentuali prima delle misure. Il criterio di accettazione non è solo velocità: ogni accelerazione deve preservare la fisica e non peggiorare in modo significativo scaling, footprint di memoria o manutenibilità.


## 13. Rischi e contromisure

| Rischio | Contromisura |
| --- | --- |
| Alterare soglie o ordini di arrotondamento wet/dry | Golden tests focalizzati e confronto ramo per ramo sui casi-limite. |
| Diversità delle derivate COMPLEX | Formula condivisa e test complex-step rispetto a differenze/derivate attese. |
| Ri o cinematica non aggiornati dopo cambi velocità | Operazioni di aggiornamento atomiche per dipendenza e test invarianti dello stato. |
| Rottura well-balanced in PCCU | Calcolare proprietà dopo la traccia HP definitiva; test di equilibrio per faccia. |
| Overengineering di kernel e tipi | Massimo quattro kernel primari; helper solo con responsabilità distinta e chiamanti reali. |
| Allocabili e copie nascoste nei loop | Preallocazione, INTENT espliciti, profiling allocazioni e report compilatore. |
| Integrazione ML che inquina lo stato | Proprietà generali nello stato; feature/preprocessing e parametri del modello nei processi ML. |


## 14. Definition of done

- Matrice dei consumatori verificata e nessun uso delle vecchie API residue nei sorgenti di produzione.
- Quattro kernel primari con contratti chiari, formule REAL/COMPLEX condivise e test che coprono ogni ramo essenziale.
- Stato generale REAL/COMPLEX sempre completo e coerente, senza dipendenza dalla reologia o da flag di proprietà disponibili.
- Nessuna allocazione/deallocazione, controllo di dimensione ridondante o copia profonda evitabile nei loop caldi.
- Nessuna duplicazione di calcoli termodinamici nelle chiamate sequenziali alla stessa cella o allo stesso endpoint di faccia.
- Separazione numerica/fisica confermata: nessun j/k in physics, nessuna cache globale e nessuna dipendenza circolare dei moduli.
- Accettazione CI e test scientifici, inclusi path idrostatici, casi di iniezione, TAAL, Newton/complex-step e scaling.
- Profiling prima/dopo pubblicato con misure per sottosistema e valutazione numerica documentata.
- Documentazione aggiornata: ownership, invarianti, mappa delle dipendenze, interfacce e prescrizioni per aggiungere una reologia o un modello ML.

## 15. Registro di correzione e verifica P0–P3

Questo registro è parte del piano corrente di lavoro. Ogni voce deve risultare `NON INIZIATA`, `IN CORSO` o `COMPLETATA` con file e test associati; le voci non possono essere considerate concluse sulla sola base della compilazione.

| Voce | Stato | Evidenza corrente | Lavoro residuo |
| --- | --- | --- | --- |
| P0.1 Pressione dinamica | COMPLETATA | Formula corretta in `state_new_conversions.f90`; test isolato superato | Aggiungere confronto esplicito REAL/COMPLEX nel test |
| P0.2 Celle asciutte | CODICE COMPLETATO | Costruttori REAL/COMPLEX da `q` ora seguono lo stesso ramo per ogni massa positiva; la regolarizzazione near-dry è applicata solo alle velocità, mentre massa nulla resta dry | Validare formalmente densità/proprietà dry contro il backend canonico |
| P0.3 Regolarizzazione near-dry | COMPLETATA | Costruttori da `q` usano `eps_sing`/`eps_sing4`; conversioni già regolarizzate | Aggiungere casi near-dry dedicati alla suite |
| P0.4 Richardson REAL/COMPLEX | CODICE COMPLETATO | Formula REAL/COMPLEX condivisa in `richardson.inc`; cap verificato prima della divisione per evitare overflow; velocità nulla e valori negativi preservati | Estendere la copertura numerica ai test dedicati del backend di riferimento |
| P0.5 Stato COMPLEX completo | CODICE COMPLETATO | Reset espliciti REAL/COMPLEX; `Zs`, `exc_pore_pres` e `T` inizializzati e trasferiti da `q`/`qp` anche nei percorsi COMPLEX e dry | Aggiungere test di riuso con processi abilitati/disabilitati |
| P0.6 Chiusura composizione | CODICE COMPLETATO | Unico modulo `state_new_composition.f90` usato da conversioni, costruzione e refresh; residui entro tolleranza normalizzati, stati inadmissibili rifiutati | Confrontare la normalizzazione con il backend canonico nella fase di test |
| P0.7 Fallimenti riconoscibili | CODICE COMPLETATO | `accepted` opzionale nelle conversioni, nei costruttori e nel refresh; dry e inadmissibile hanno esiti distinti | Mantenere il contratto durante l’integrazione con i consumatori |
| P1.1–P1.5 Duplicazioni e costo | CODICE COMPLETATO | Kernel minimi per calori e densità; costruttori REAL/COMPLEX da `q` diretti; conversioni riusano `inv_rho_m`; chiusura composizionale condivisa | Misurare il costo e verificare i percorsi finali |
| P1.4 Kernel minimi | CODICE COMPLETATO | `density_from_mass.inc` separa densità/reciproco dalle proprietà complete; wrapper REAL/COMPLEX condivisi | Verificare tutte le conversioni e impedire proprietà non richieste nei percorsi minimi |
| P1.6 Stato generale | CODICE COMPLETATO | Tipi unificati REAL/COMPLEX senza campi per reologia; `refresh_physical_state` aggiorna le dipendenze dopo mutazioni | Verificare tutti i consumatori e i campi obbligatori |
| P2.1 Verifiche ALLOCATED/SIZE | COMPLETATA | `require_*_storage` separa i controlli di allocazione dalle chiamate a `SIZE` | Mantenere il controllo al confine di setup |
| P2.2 Layout dei quattro array | COMPLETATA | Inizializzazione workspace verifica e prepara `xs`, `xg`, `alphas`, `alphag` | Aggiungere test di layout parzialmente allocato |
| P2.3–P2.6 Memoria e workspace | CODICE COMPLETATO | `layout_valid` impostato al setup; i costruttori verificano un solo flag; reset espliciti impediscono valori obsoleti; nessun array duplicato nel workspace | Verificare indipendenza worker e assenza di copie profonde nei consumatori |
| P3.1–P3.7 Test e benchmark | IN CORSO | Aggiunto `TESTS/TEST_STATE_EQUIVALENCE` per il confronto diretto legacy/`state_new`; la workflow manuale esegue anche il confronto; caso locale superato | Eseguire la Action sul repository, estendere le perturbazioni COMPLEX e aggiungere benchmark |

## Appendice A — API concettuale di destinazione

```fortran
! Vettori locali; nessun indice di griglia nel nucleo fisico.

call conservative_to_primitive(q, qp, accepted, p_dyn)
call primitive_to_conservative(qp, q, accepted)

call evaluate_mixture(T, xs, xg, mixture)
call build_physical_state_from_q(q, geometry, state_real)

call build_physical_state_from_qp(qp, geometry, state_real)
call build_physical_state_from_q(q_complex, geometry, state_complex)

call update_physical_state_kinematics(state_real, geometry, u_new, v_new)
! Proprietà di faccia dalla traccia finale ricostruita:

call prepare_face_properties(qp_face, face_properties)
```

Le due forme REAL da q e da qp dovranno avere nomi diversi o un discriminante di tipo, perché in Fortran non possono essere distinte automaticamente da un’interfaccia generica quando gli argomenti hanno lo stesso tipo e rango. Le firme sono quindi esempi del contratto, non sorgente pronto alla compilazione.


## Appendice B — Vincoli scientifici da non alterare

- Energia conservativa: M cp_m T, non energia cinetica totale.
- Richardson: h g’/(u²+v²), convenzione a velocità nulla e saturazione a 10^15, con protezioni prima della divisione.
- Composizione: chiusura e trattamento dei residui di arrotondamento coerenti con gli attuali controlli.
- Regolarizzazione delle velocità near-dry e politiche separate di accettazione/ammissione.
- Galleggiamento e gamma ricalcolati sulla traccia finale; preservazione del well-balanced.
- Geometria fornita per stage/cella, senza precomputazione assunta immutabile nel tempo.
- Complex-step: branching sulle parti reali dove richiesto, perturbazioni immaginarie preservate.
- Interfacce fisiche complete, fallimento esplicito dei trattamenti non implementati.
Riferimento tecnico: repository demichie/IMEX_SfloW2D_dev, master; analisi dei chiamanti e delle implementazioni di state_conversion_2d, source_terms_2d, equation_terms_2d, nonlinear_solver_2d, mass_exchange, rheology e moduli spaziali al 5 novembre 2026.

## 15. Regole di manutenzione del documento

- Modificare questo file insieme ai change set che avanzano le fasi F0–F9; aggiornare **stato, commit, test e rischi** senza marcare completato il lavoro sulla sola base di un audit statico.
- Tenere le decisioni definitive in cima; nel registro operativo annotare le prove più recenti, evitando cronologie prolisse che rendono il piano illeggibile.
- Riportare nella [roadmap generale](IMEX_SfloW2D_Development_Roadmap.md) il progresso aggregato, senza duplicare i dettagli tecnici.
- Se un'ottimizzazione modifica l'ordine di operazioni floating-point, registrarlo e usare tolleranze scientificamente giustificate; non dichiarare equivalenza bitwise non verificata.
- Le regole di formattazione Fortran sono in [FORTRAN_STYLE.md](FORTRAN_STYLE.md); i test e le evidenze di accettazione rimangono in `TESTS/` e `docs/acceptance/`.

## Audit architetturale finale - Blocco 1 (5 novembre 2026)

Stato: **IN CORSO**. Ricognizione statica del ramo `master` al commit `a84f6e093539` prima delle modifiche di questo blocco.

- Le famiglie di termini sorgente esaminate espongono le interfacce `explicit`, `semi_implicit`, `implicit` e `analytic`; i metodi non implementati devono continuare a terminare con errore esplicito. Presenza delle firme non equivale a validazione dei comportamenti.
- Nelle implementazioni `src/physics/**/*.f90` controllate non sono emersi accessi diretti a `q(j,k,...)` o `qp(j,k,...)`; le dipendenze da `geometry_2d` osservate nei kernel di iniezione, pore pressure e coefficienti sono limitate a utility geometriche (`pi_g`, `bottom_source_area_ratio`), da mantenere distinte dai campi globali.
- Il `src/Makefile.am` distribuiva gia `rheology_newmodel_linear.inc` e `thermal_cooling_rate.inc` tramite `EXTRA_DIST`, ma mancavano i due archi di dipendenza per la ricompilazione dei rispettivi oggetti: aggiunti in questo blocco.
- Restano da verificare con controlli riproducibili: dipendenze `USE` vs regole Makefile, assenza di cicli sull'intero grafo, duplicazioni e wrapper realmente inutilizzati, correttezza dei quattro trattamenti per ciascun modello, e build/test CI del ramo di audit.
- Nessuna rimozione di interfacce obbligatorie o modifica delle formule fisiche autorizzata dal solo audit statico.

### Blocco 1 - verifica sistematica delle dipendenze (5 novembre 2026)

- Sono stati verificati i `USE` dei 68 file `.f90` di produzione elencati in `src/` (esclusi i diagnostici Fortran sotto `python_diagnostics`). Il confronto tra i moduli interni e le regole esplicite del Makefile ha individuato **14 archi diretti mancanti**, ora aggiunti a `src/Makefile.am`.
- Il grafo delle dipendenze di compilazione, integrato con questi 14 archi, non presenta cicli (verifica statica sul grafo delle regole); il risultato non sostituisce un build parallelo.
- Le implementazioni esaminate sotto le cinque famiglie di sorgenti (`rheology`, `mass_exchange`, `mass_injection`, `pore_pressure`, `thermal`) hanno tutte le quattro firme previste. I metodi non supportati devono continuare a segnalare errore; la validazione dinamica dei singoli metodi non e inclusa nella verifica delle firme.
- Aggiunto `.github/scripts/check_fortran_dependencies.py` per controllare su ogni PR: copertura dei sorgenti nel Makefile, `USE` diretti, `.inc`, cicli e presenza delle quattro interfacce; nuova workflow `.github/workflows/fortran-dependency-audit.yml`.
- Wrapper e duplicazioni non sono stati rimossi: servono prima evidenze sui chiamanti e regressioni; l'interfaccia pubblica di `state.f90` resta per ora invariata.
- **Da validare:** esecuzione della nuova Action sul commit, build parallelo e regressioni scientifiche. Stato Blocco 1: `IN CORSO` fino a tali prove.

### Blocco 1 - esito del primo controllo CI

- La Action `Fortran architecture dependency audit`, run [42115367142](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/42115367142), ha restituito **success** sul commit `a51f4fc6c04a` (verifica statica automatica, non compilazione).
- Aggiunto un secondo job per una build `gfortran` parallela (`make -j4`), per confermare il comportamento delle dipendenze Fortran su ambiente CI. La build resta da osservare dopo il commit.

### Blocco 1 - verifica build e contratti delle sorgenti

- La build completa `gfortran` con `make -j4` e il controllo statico sono **PASS** nella [Action 42115531813](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/42115531813). La dipendenza incrociata fra moduli e Makefile e quindi verificata anche da un compilatore in parallelo.
- Audit dei modelli reologici: gli 11 modelli hanno le quattro interfacce; nei modelli con contributi reali/complessi comuni, i file `.inc` eliminano la duplicazione dei kernel senza introdurre conversioni improprie. Gli `explicit` che assemblano il semi-implicito e l'implicito reale costituiscono composizione prevista, non wrapper di compatibilita eliminabili.
- `rheology_froude`, `rheology_lucas`, `rheology_mui`: la routine implicita restituisce zero, ma i tre termini sono catalogati come `has_newton=.FALSE.` e `default_method=EXPLICIT`; il test `TEST_RHEOLOGY_FROUDE` controlla anche il valore zero della sua interfaccia complessa. Non modificare il comportamento senza riprogettare esplicitamente il contratto di partizione semi-implicito/implicito.
- `rheology_newmodel`: `explicit` e `analytic` falliscono esplicitamente e sono marcati non disponibili nel catalogo; il `semi_implicit` restituisce zero mentre il contributo fisico e nell'`implicit`. Questo non e un metodo pubblicizzato come implementato per i gruppi non consentiti.
- Le varianti di dispatch per nome della reologia sono ancora chiamate dai test di accettazione; la loro rimozione non rientra in una pulizia dimostrabilmente sicura. `eval_source_analytic_state` e un helper chiamato da `eval_source_analytic` per la mappa di stato esatta, non un wrapper di compatibilita inattivo.
- Aggiunta una matrice CI con sei suite (`SOURCE_REGISTRY`, `SOURCE_MASS`, `MASS_EXCHANGE`, `RHEOLOGY_CELL_LOCAL`, `RHEOLOGY_FROUDE`, `STATE_CONVERSION`) per verificare il comportamento dopo l'audit. **Da osservare gli esiti sul nuovo commit.**
- Resta una verifica esaustiva degli entry point pubblici contro tutti i consumatori prima di poter attestare l'assenza assoluta di wrapper inutilizzati; la sola ricerca testuale non costituisce prova sufficiente.

### Blocco 1 - allineamento dei runner di regressione

- La prima matrice di regression CI ha segnalato fallimenti in piu suite, pur con build parallelo e audit statico superati. Il confronto con il Makefile ha rilevato che cinque script di test compilavano i sorgenti senza gli include path delle reologie e del cooling; allineate le opzioni `-I` a quelle del Makefile.
- `TEST_SOURCE_MASS` richiede l'eseguibile, non lo compila autonomamente: la workflow ora costruisce il programma prima di eseguire quella suite.
- I risultati precedenti delle regressioni non sono considerati verifiche di fisica; attendere l'esito dei nuovi job prima di attribuire eventuali errori ai kernel.

### Blocco 1 - gate finale delle regressioni (Action 42116271701)

- Sul commit `03502fe3a811`, i job **Fortran static audit** e **parallel build** sono `PASS`.
- I job `TEST_STATE_CONVERSION`, `TEST_SOURCE_MASS`, `TEST_SOURCE_REGISTRY`, `TEST_MASS_EXCHANGE` e `TEST_RHEOLOGY_CELL_LOCAL` risultano `FAIL` al passo `Run regression`, anche dopo la correzione dei percorsi di include e la compilazione dell'eseguibile per il test integrato. Il job Froude non aveva ancora un esito definitivo al momento della rilevazione.
- La cause di questi fallimenti non sono state accertate: lo stato dei job non contiene i messaggi del compilatore o le asserzioni che permetterebbero di distinguere problemi del runner da regressioni del codice. Non modificare la fisica sulla base dei soli exit code.
- La PR #5 e il **Blocco 1 restano aperti**. Prima del merge occorre ispezionare i log/artifact della Action, distinguere cause infrastrutturali da difetti reali, correggere i problemi riproducibili e riottenere un gate di regressione verde.
