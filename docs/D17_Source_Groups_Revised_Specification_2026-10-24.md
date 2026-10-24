# D17 — Organizzazione unificata dei termini sorgente e dei gruppi di integrazione

**Specifica di revisione autonoma | 24 ottobre 2026 | stato: IMPLEMENTAZIONE IN CORSO (gate numerici aperti)**

> Questo documento riguarda **solo D17**. Sostituisce, per gli aspetti in conflitto, la precedente impostazione D17 che consentiva più solve Newton indipendenti all'interno dello stesso gruppo. Non costituisce una dichiarazione che il codice sia già allineato né un'autorizzazione a cambiare la fisica o gli schemi RK-IMEX. Le indicazioni «attuale» riportano il comportamento osservato nel repository `demichie/IMEX_SfloW2D_dev` sul branch `master` consultato il 24 ottobre 2026 e nella documentazione D17 allegata alla roadmap. Le indicazioni «richiesto» definiscono invece il contratto da implementare.

> **Audit consolidato al 24 ottobre 2026 — HEAD validata `35be9f0a`:** entrambe le Actions verdi confermate. Catalogo di 24 sorgenti e 96 procedure concrete, registro senza sotto-solver, assemblatori in `source_terms_2d`, `equation_terms_2d` limitato a cinque procedure, correzione della quantità di moto e finalizzazione post-RK/halo in `state_2d`. Restano aperti dispatch scientifico, Newton accoppiato, equivalenza legacy/nativo, `rheology_newmodel`, regressioni D17.20 e gate D17.21.

## D17.0 — Obiettivo e principio di semplificazione

D17 ha un solo obiettivo funzionale: per **ogni termine sorgente fisico attivo** decidere **(a) come integrarlo**, **(b) quando integrarlo** e **(c) a quale gruppo di integrazione appartiene**. Il raggruppamento serve a stabilire quali equazioni trattare *insieme*, in particolare quando è richiesto un solve Newton accoppiato. Nessun livello addizionale di raggruppamento o scheduling deve essere introdotto senza una necessità concreta.

Le tre proprietà sono:

1. **Metodo**: `EXPLICIT`, `NEWTON` (trattamento semi-implicito più implicito esistente), oppure `ANALYTIC` (aggiornamento con soluzione analitica disponibile).
2. **Cadenza**: `STAGE` (all'interno dell'integrazione RK-IMEX, a ciascuno stadio previsto dal tableau) oppure `STEP` (dopo il completamento dell'intero passo RK).
3. **Gruppo**: insieme di termini e relative equazioni trattati come una sola unità numerica, con metodo e cadenza uniformi.

**Regola vincolante: un gruppo = un metodo + una cadenza + una unità di soluzione.** Nei gruppi `NEWTON` questo implica **un unico sistema non lineare con tutte le equazioni target del gruppo**, mai una sequenza di solve Newton separati dentro lo stesso gruppo. Se due problemi devono essere risolti separatamente, appartengono a gruppi differenti. Un gruppo `ANALYTIC` contiene **un solo aggiornamento/solve analitico**. Un gruppo `EXPLICIT` assembla i contributi espliciti, senza un solve implicito.

**Esempio vincolante.** Velocità `NEWTON` e temperatura `ANALYTIC`, entrambe `STAGE`, **non** possono condividere il gruppo: servono due gruppi. Velocità e temperatura entrambe `NEWTON` e `STAGE` possono condividere un gruppo: allora Newton risolve **congiuntamente** il sistema per `(q_2,q_3,q_4)` (salvo layout specifico del modello). Se devono essere risolte una dopo l'altra, servono due gruppi Newton ordinati. Il raggruppamento non autorizza a trascurare derivate incrociate realmente presenti nel residuo.


**Stato di avanzamento D17.0 — PARZIALE — architettura implementata; gate numerici aperti.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.0a — Contratto uniforme delle leggi fisiche

**Ogni famiglia di termini sorgente è implementata nel proprio modulo fisico, con lo stesso contratto di quattro interfacce: `explicit`, `semi_implicit`, `implicit`, `analytic`.** Le interfacce esistono anche quando una particolare famiglia non supporta quel trattamento. Il catalogo dichiara esclusivamente il supporto dei tre metodi `EXPLICIT`, `NEWTON`, `ANALYTIC`; le interfacce dei metodi non disponibili devono, se chiamate erroneamente, produrre un errore diagnostico (mai un rate nullo fittizio). Il parser/validatore impedisce già prima della simulazione di selezionarle. `NEWTON` richiede **entrambe** le interfacce `semi_implicit` e `implicit`: la prima può restituire un contributo nullo se non occorre un trattamento semi-implicito. Il catalogo non deve descrivere come la legge ripartisce internamente i contributi. **Le quattro interfacce fisiche non introducono un quarto metodo numerico**.

Il requisito vale indistintamente per reologia, raffreddamento/termica, pressione dei pori, trasferimenti di massa e ogni altra legge sorgente presente o futura. `source_terms_2d` invoca solo le interfacce pertinenti alla configurazione validata e **assembla tutti i contributi**, senza eccezioni o canali laterali per una famiglia fisica. Le routine esistenti possono essere mantenute come implementazione interna delle interfacce: evitare duplicazione delle formule e wrapper inutilmente stratificati.


**Stato di avanzamento D17.0a — IMPLEMENTATO / CI VERDE — quattro interfacce concrete per 24 sorgenti (96 procedure); dispatch numerico da verificare.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.1 — Confini di responsabilità

| Componente | Responsabilità ammessa | Non deve fare |
|---|---|---|
| `source_catalog_2d` | Identità stabile e capacità/compatibilità dei termini fisici | Calcolare sorgenti, scegliere `dt`, eseguire solve |
| Configurazione e validazione (`source_input_2d`, `source_configuration_2d`, registro) | Tradurre input legacy/nativo in gruppi semplici, verificare integralmente il contratto **una sola volta prima della simulazione** | Valutare rate per cella, modificare la fisica, eseguire iterazioni Newton |
| `source_terms_2d` | **Unico punto di assemblaggio fisico** dei contributi espliciti, semi-impliciti, impliciti, analitici e dei contributi di scambio di massa richiesti dal gruppo | Determinare la cadenza, scegliere il metodo, gestire RK o i cicli temporali |
| Moduli fisici (`rheology_2d`, `thermal_laws_2d`, `pore_pressure_laws_2d`, `mass_exchange/mass_exchange_laws_2d` e moduli dei processi) | Leggi fisiche dietro quattro interfacce uniformi (`explicit`, `semi_implicit`, `implicit`, `analytic`), con disponibilità dichiarata per ciascuna famiglia; kernel REAL/COMPLEX solo dove necessari | Conoscere gruppi, tableaux RK o scheduling |
| `time_integration_2d` e `simulation_2d` | Eseguire gruppi secondo cadenza e ordine, applicando correttamente i coefficienti RK e il fractional step post-RK | Implementare formule costitutive o reinterpretare configurazioni non valide |
| `nonlinear_solver_2d` | Risolvere **un sistema Newton per gruppo**, con residuo/Jacobiana ricavati dalle API di `source_terms_2d`; richiamare nei punti prescritti le correzioni dello stato tramite `state_2d` | Contenere dispatch della fisica, implementare direttamente le correzioni dello stato o più solve indipendenti interni al gruppo |
| `state_2d` | Implementare verifiche, correzioni di roundoff, proiezioni delle celle asciutte e correzioni numeriche sui blocchi di variabili, indipendentemente dall’integratore chiamante | Conoscere tableaux RK, cadenzamento dei gruppi o scegliere leggi fisiche |
| `state_conversion_2d` | Conservare le trasformazioni conservativo/primitivo e i relativi kernel esistenti, riutilizzati da `state_2d` dove opportuno | Introdurre scheduling delle sorgenti o integrare temporalmente le equazioni |
| Esecutore generico dei gruppi `STEP`, geometria | Applicare gli incrementi prodotti dai gruppi e coordinare aggiornamenti di inventari/terreno quando previsti | Introdurre un percorso temporale speciale per mass-exchange o riassemblare la fisica |

**Rinomina obbligatoria (D17):** l’attuale modulo `equation_terms_2d` deve diventare `source_terms_2d` (file `src/physics/source_terms_2d.f90`), perché il suo compito è assemblare i termini sorgente. La rinomina non deve tuttavia trasferire automaticamente tutte le API attuali nel nuovo modulo: le routine di **flusso inerziale**, **coefficienti idrostatici**, **velocità caratteristiche** e gli altri operatori non-sorgente devono essere assegnati ai moduli già responsabili del trasporto/spazio o a un modulo mirato soltanto se necessario. Non modificare le formule o la discretizzazione HP-PCCU durante lo spostamento. Le sorgenti al bordo producono uno stato imposto, non un rate di cella: mantenere la valutazione fisica nell’interfaccia appropriata di `source_terms_2d` e l’applicazione nel relativo hook di bordo. Aggiornare tutti i `USE`, il grafo delle dipendenze Fortran, `Makefile.am`, test, documentazione e gli audit statici; evitare wrapper permanenti con il vecchio nome. Prima del cambio verificare il call graph completo per non lasciare API non-sorgente senza proprietario.


**Stato di avanzamento D17.1 — PARZIALE — separazione catalogo, sorgenti e stato avanzata; equivalenza e confini residui aperti.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.2 — Catalogo minimo: nomi, target e tre capabilities (vincolante)

**Ogni termine sorgente configurabile si registra una sola volta mediante il nome canonico univoco, il target e tre capacità:** `has_explicit`, `has_newton`, `has_analytic`. Questa è l'intera registrazione funzionale minima; evitare metadati duplicativi quali ID numerici fisici, `family`, `model_name`, `bundle_id`, `role`, `complex_rate`, `has_semi_implicit` e `has_implicit` nel descrittore della sorgente. I tre booleani rappresentano *metodi di integrazione disponibili*, non quattro modalità del solver. `NEWTON` **implica** che le due routine `semi_implicit` e `implicit` siano entrambe disponibili; `semi_implicit` può legittimamente non contribuire e restituire zero. Le interfacce non supportate devono fallire esplicitamente quando invocate per errore, mentre l'input incompatibile viene rifiutato prima di iniziare.

Schema concettuale Fortran, nomi/target da normalizzare all'input:

```fortran
TYPE :: source_term_type
  CHARACTER(LEN=64) :: name = ''
  CHARACTER(LEN=32) :: target = ''
  LOGICAL :: has_explicit = .FALSE.
  LOGICAL :: has_newton = .FALSE.
  LOGICAL :: has_analytic = .FALSE.
END TYPE source_term_type
```

**Non esistono ID numerici permanenti** né per i termini sorgente né per le identità dei modelli fisici, comprese le reologie. Sostituire `SOURCE_TERM_*` numerici, `rheology_id_*`, `rheology_model_id`, `model_id`, `term_id` e conversioni nome–ID nel percorso fisico/configurativo con nomi canonici; sono ammessi gli interi come indici temporanei di array o come indici delle equazioni, mai come identità alternative. Il catalogo rimane statico e una ricerca lineare per nome è sufficiente; non introdurre factory, tabelle dinamiche o procedure pointer. Normalizzare una volta i nomi (per esempio minuscolo, trim), verificarne l'unicità e rifiutare nomi non riconosciuti.

**Granularità:** un termine del catalogo è una *sorgente fisica completa configurabile*, non ciascuna addenda della sua formula. Per esempio, `rheology_voellmy` è **una sola sorgente** con target `momentum`: Coulomb e turbulent drag sono componenti interne valutate dalle routine `eval_voellmy_semi_implicit_momentum_rate` e `eval_voellmy_implicit_momentum_rate`, non due sorgenti registrate, schedulate o validate separatamente. Analogamente le componenti viscosa/turbolenta/yield degli altri modelli reologici restano dettagli interni della legge. La rappresentazione interna non può essere ripartita in gruppi di scheduling distinti; il metodo `NEWTON` combina SI/IM rispettando il contratto numerico esistente. L'accorpamento non autorizza a contare due volte un rate.

La selezione del metodo e della cadenza rimane esclusivamente **nel gruppo**. Il catalogo descrive le capacità del metodo; la validazione generale applica le cadenze supportate dall'infrastruttura (`STAGE` per Newton/analytic, `STAGE` o `STEP` solo quando il percorso explicit e gli aggiornamenti dello stato sono correttamente implementati). Le restrizioni specifiche e imprescindibili di aggiornamenti atomici (per esempio inventari e terreno per mass-exchange) vengono garantite dal relativo assemblatore/esecutore, senza aggiungere una seconda registrazione della stessa informazione. Non generalizzare automaticamente `EXPLICIT/STAGE` per mass-exchange: preservare per ora il trattamento `STEP/EXPLICIT` supportato.

La selezione del target è nominale (per esempio `momentum`, `temperature`, `pore_pressure`) e viene risolta una sola volta in indici di equazioni tramite metadata delle equazioni; i solver continuano a lavorare solo con tali indici. I contributi dipendenti alle altre equazioni, per esempio quelli prodotti dal mass-exchange, vengono assemblati dalla legge/coordinatore fisico senza registrare una voce sorgente distinta per ciascun termine algebrico. Per la migrazione dei vecchi input nativi, eventuali nomi storici vanno convertiti o rifiutati esplicitamente durante la lettura, senza creare sorgenti attive duplicate.


**Stato di avanzamento D17.2 — IMPLEMENTATO / CI VERDE — catalogo canonico e ID reologici rimossi; rheology_newmodel aperto.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.3 — Semantica dei metodi

L’insieme dei metodi di integrazione **non deve essere considerato definitivamente chiuso** a `EXPLICIT`, `NEWTON` e `ANALYTIC`. Questi sono i tre metodi implementati in D17, ma il meccanismo di selezione e validazione deve consentire l’aggiunta successiva di altri metodi numerici, inclusi eventuali integratori basati su modelli appresi. Questa estensibilità è un requisito architetturale, non richiede di implementare ulteriori metodi in D17.

### D17.3.1 `EXPLICIT`

`source_terms_2d` assembla il vettore di rate dei termini selezionati dal gruppo, `S_G^EX(q)`. A cadenza `STAGE`, il vettore entra nella parte esplicita del tableau RK. A cadenza `STEP`, se il trattamento è supportato, il gruppo viene eseguito dopo RK una volta sull'intervallo `dt`, nel proprio ordine di esecuzione. La valutazione e l'aggiornamento degli stati al bordo, quando presenti, restano separati dal rate di cella.

### D17.3.2 `NEWTON`

`NEWTON` indica **un unico trattamento numerico di gruppo**, che preserva la distinzione già esistente tra componenti **semi-implicite** (eventuale predictor/correzione) e **implicite** (residuo e Jacobiana). `source_terms_2d` assembla `S_G^SI` e `S_G^IM` e ne limita la selezione ai termini del gruppo. `nonlinear_solver_2d` raccoglie l'unione delle equazioni target del gruppo e costruisce **un unico vettore di incognite e un unico residuo congiunto**. La Jacobiana è dimensionata all'intero gruppo e comprende gli accoppiamenti realmente presenti, valutati lungo il medesimo stato delle incognite.

Per uno stadio diagonale RK, l'intervallo implicito pertinente è `a_ii^IM * dt`, come nel percorso esistente; la rappresentazione esatta del predictor SI e del residuo deve mantenere la semantica numerica validata, senza rifattorizzazioni matematiche arbitrarie. `NEWTON` a cadenza `STEP` è un obiettivo di configurabilità **solo se** viene specificato e validato un operatore di step coerente: non è attualmente consentito dal registro; non dichiararne il supporto prima dell'implementazione e della validazione. Nel frattempo l'input `STEP+NEWTON` deve fallire con diagnostica esplicita.

### D17.3.3 `ANALYTIC`

`source_terms_2d` offre l'operatore analitico `q_after = Phi_G(q_before, delta_t)` per **un solo problema analitico**. L'operatore restituisce lo stato aggiornato per il target, non un rate o un incremento arbitrario. A `STAGE` usa l'intervallo indicato dal tableau (nella configurazione attuale `a_ii^IM*dt`); a `STEP`, solo qualora quel particolare solve abbia una semantica e un'implementazione supportate, riceve `dt` e agisce dopo RK. Fino ad allora `STEP+ANALYTIC` deve essere rifiutato all'avvio. **Non** accumulare più solve analitici indipendenti nello stesso gruppo: usare gruppi distinti.


**Stato di avanzamento D17.3 — PARZIALE — capacità e interfacce SI/IM presenti; comportamento quantitativo per gruppo da verificare.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.4 — Cadenza e ordine

- `STAGE`: esecuzione secondo il tableau RK-IMEX, non un semplice aggiornamento Euler ripetuto. I contributi espliciti e impliciti vengono combinati secondo i coefficienti e i pesi degli stadi. La valutazione ripetuta del rate in Newton (residuo/Jacobiana/line search) non equivale a più applicazioni temporali della sorgente.
- `STEP`: esecuzione **una sola volta** dopo l'assemblaggio RK e prima dell'avanzamento dell'orologio. Un gruppo di trasferimento di massa `STEP/EXPLICIT` conserva l’atomicità e gli aggiornamenti coerenti di stato, inventari e terreno **attraverso il medesimo esecutore dei gruppi STEP**.
- L'ordine è **tra gruppi** compatibili con la fase di esecuzione. Nessuna numerazione di solve indipendenti dentro un gruppo. L'ordine `STAGE` deve rispettare la separazione dei contributi espliciti e impliciti del tableau: un unico ordinamento numerico non deve essere usato per pretendere di cambiare l'ordine matematico del metodo RK.
- Identità di gruppo: combinazione `(cadence, group_id)` oppure ID globalmente univoco, secondo il contratto di input scelto; usare una regola unica e documentata. L'ordine non costituisce l'identità del gruppo. Conservare il namespace stage/step indipendente senza ambiguità.


**Stato di avanzamento D17.4 — PARZIALE — cadenza/ordine disponibili; dispatch STEP generale e confronto temporale ancora aperti.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.5 — Struttura minima del gruppo

Contratto concettuale (non è prescrizione di una firma Fortran identica):

```fortran
TYPE :: source_group_type
  INTEGER :: id
  INTEGER :: cadence         ! STAGE / STEP
  INTEGER :: method          ! EXPLICIT / NEWTON / ANALYTIC
  INTEGER :: execution_order
  CHARACTER(LEN=64), ALLOCATABLE :: term_names(:)
  INTEGER, ALLOCATABLE :: equation_ids(:)
  ! Per ANALYTIC: l'identità dell'unico operatore analitico, se necessaria.
END TYPE
```

**Eliminare concettualmente** dal modello del gruppo le liste `implicit_solver_ids(:)`, `analytic_solver_ids(:)` e i `source_solver_type` intesi come livello di scheduling con più solve per gruppo. Se per ragioni di migrazione si conserva temporaneamente un descrittore helper di solve, deve esserci **esattamente un descrittore per gruppo implicito** e non deve poter definire un sottopartizionamento diverso dalle equazioni del gruppo. Preferire la rappresentazione minima che riduce duplicazioni e stati inconsistenti. `legacy_id` può sopravvivere temporaneamente solo nell'adattatore o per tracciabilità, senza influenzare la nuova semantica.

Gli `equation_ids` del gruppo derivano dai termini e dal modello fisico. Per Newton identificano esattamente il sistema accoppiato (ad esempio, quantità di moto x/y e temperatura se condividono il gruppo). Per `ANALYTIC` devono coincidere con il target dell'unico solve disponibile. Per `EXPLICIT` definiscono il footprint delle sorgenti, non un sistema lineare da risolvere.


**Aggiornamento 24 ottobre — commit `1d916bcd`:** eliminate le argomentazioni `solvers(:)` dalle API di orchestrazione in `nonlinear_solver_2d` e la copia `source_solvers` in `time_integration_2d`. I gruppi e i relativi `equation_ids` determinano la chiamata al Newton e all'operatore analitico; `source_solver_type`, `solver_entries` e le liste degli ID di solve permangono nel registro e nel parser e devono essere rimossi con una successiva modifica coordinata. Build e test sul nuovo commit non ancora eseguiti.

**Aggiornamento — rimozione strutturale (24 ottobre 2026):** preparata la sostituzione del registro basato sui descrittori dei sotto-solver con il registro fondato unicamente sui gruppi. Eliminati `source_solver_type`, `solver_entries`, `implicit_solver_ids`, `analytic_solver_ids` e le relative API; adattati parser nativo, adapter legacy e test. Questa modifica non è ancora stata verificata mediante le Actions e non dimostra equivalenza numerica; lo stato resta IN CORSO.

**Stato di avanzamento D17.5 — COMPLETATO (STRUTTURA) / CI VERDE — nessun sotto-solver o term_ids, singolo metodo e target di gruppo.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.6 — Gruppi `ANALYTIC`

- Un solo metodo: `ANALYTIC`; una sola cadenza; **un solo solve**.
- La soluzione analitica deve essere effettivamente disponibile per tutti i termini del gruppo considerati nel loro insieme, non supposta dalla sola presenza di un singolo termine analitico.
- Non creare automaticamente un solve analitico composito concatenando più operatori individuali.
- Se velocità, temperatura e pressione dei pori sono integrate analiticamente in modo indipendente, usare gruppi differenti (anche se il vecchio input le riferiva allo stesso `IMEX_group`).
- L'operatore aggiorna soltanto le equazioni target, lasciando invariato il restante stato conservativo, a meno di una trasformazione dichiarata e validata.


**Stato di avanzamento D17.6 — PARZIALE — operatore analytic unico implementato; confronto con baseline aperto.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.7 — Gruppi `NEWTON` (requisito corretto e vincolante)

**Il gruppo Newton corrisponde a un solo sistema non lineare.** Tutte le equazioni del gruppo entrano nello stesso vettore di incognite `U_G`, con residuo `R_G(U_G)=0` e una Jacobiana della dimensione di `U_G`. Questo vale anche se alcune derivate incrociate sono zero in un modello particolare.

Esempi:

| Configurazione | Gruppi necessari | Risultato |
|---|---|---|
| Velocità `NEWTON/STAGE`, temperatura `NEWTON/STAGE`, **accoppiate** | Uno | Un unico Newton sulle variabili `(q2,q3,q4)` |
| Velocità e temperatura `NEWTON/STAGE`, **separate** | Due | Newton `(q2,q3)` e Newton `(q4)` in gruppi distinti |
| Velocità `NEWTON/STAGE`, temperatura `ANALYTIC/STAGE` | Due | Newton per velocità; un distinto aggiornamento analitico della temperatura |
| Velocità `NEWTON/STAGE`, temperatura `NEWTON/STEP` | Due | Cadenze distinte; la seconda configurazione è attualmente non supportata e va rifiutata |

**Obsoleto:** il requisito precedente «un gruppo Newton può contenere più solve Newton indipendenti eseguiti in ordine». Le sezioni D17.7, D17.10, D17.19 e i vecchi report che lo descrivono documentano un'implementazione da correggere, non il nuovo obiettivo. Il vecchio blocco Newton accoppiato deve continuare a essere risolto congiuntamente; il nuovo contratto deve estendere tale regola ai gruppi nativi. Quando due gruppi condividono dipendenze fisiche, l'esecuzione sequenziale può differire dal solve accoppiato: non prometterne equivalenza senza verifica.


**Stato di avanzamento D17.7 — PARZIALE — singolo Newton per gruppo; residui/Jacobiane accoppiati da verificare.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.8 — Modulo `source_terms_2d` e API di assemblaggio

**Conservare le quattro interfacce di assemblaggio**, alimentate in modo uniforme dalle corrispondenti quattro interfacce dei moduli fisici, adattandole a un *descrittore di gruppo* già validato e rimuovendo il passaggio di un selettore di sotto-solve distinto dal gruppo quando non più necessario:

- `eval_expl_terms(..., group)` → rate esplicito delle equazioni del gruppo;
- `eval_semi_implicit_terms(..., group)` → contributo SI di un gruppo Newton;
- `eval_implicit_terms(..., group)` → contributo IM dello stesso gruppo, in REAL e COMPLEX per residuo/Jacobiana;
- `eval_analytic_terms(group, ..., interval, ...)` → stato aggiornato dall'unico operatore analitico del gruppo.

L'elenco `term_names` seleziona i contributi e il modulo assembla i vettori sulle equazioni appropriate; non deve scegliere autonomamente la cadenza o il metodo. Non introdurre un dispatcher universale con tipo di output ambiguo: rate SI/IM/EX e aggiornamento analitico hanno contratti matematici diversi. **Regola generale senza eccezioni:** ogni contributo sorgente delle equazioni, **inclusi tutti i trasferimenti di massa**, deve essere assemblato in `source_terms_2d`. In particolare erosione, deposizione, perdite della fase continua, entrainment, vaporizzazione e contributi accoppiati a quantità di moto, energia, componenti e pressione dei pori non devono essere assemblati una seconda volta nell'esecutore temporale. L'assemblatore deve produrre un risultato coerente che comprenda sia il vettore delle sorgenti delle equazioni conservative sia i rate accessori necessari per gli inventari dei materiali e per la variazione del letto. Gli eventuali termini derivati vengono inclusi automaticamente nelle equazioni pertinenti, senza essere schedulati come sorgenti indipendenti. Tenere distinta la *valutazione/assemblaggio dei rate* dall'*applicazione atomica* dello step.


**Stato di avanzamento D17.8 — IMPLEMENTATO (ARCHITETTURA) / CI VERDE — quattro assemblatori comuni in source_terms_2d; regressioni aperte.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.9 — Rate espliciti e sorgenti al bordo

**Base già implementata:** routing dei rate di cella `STAGE/EXPLICIT` filtrato dal gruppo e hook al bordo separato. Nel riallineamento non duplicare le sorgenti al bordo come rate di cella e non alterare l'operatore spaziale, la curvatura o le finestre di attivazione. Una sorgente di bordo produce uno stato imposto: la sua partecipazione alla configurazione deve rispettare tale rappresentazione. Conservare i test già validi e aggiungere verifiche mirate alla struttura semplificata.


**Stato di avanzamento D17.9 — PARZIALE — routing esplicito e hook al bordo migrati; rate e BC da confrontare.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.10 — SI/IM, Newton e Jacobiana congiunta

**Correzione necessaria:** oggi `solve_implicit_groups` scorre `implicit_solver_ids` ed esegue `solve_newton_group` più volte all'interno di un gruppo. Sostituire questo comportamento con **una sola chiamata Newton per gruppo** e una selezione SI/IM che includa **tutti i termini** del gruppo.

Azioni:

1. Costruire un unico insieme ordinato di equazioni target, senza duplicati, associato al gruppo.
2. Preparare il contributo SI pertinente all'insieme e conservarne la correzione/limitazione legacy senza doppio conteggio.
3. Costruire residuo e Jacobiana sull'intero vettore incognito; ogni valutazione usa lo stesso stato di prova di tutte le incognite e i rate REAL/COMPLEX coerenti.
4. Eseguire un solo solve Newton, con line search e diagnostica riferite al gruppo.
5. Se il solve fallisce, propagare il fallimento dello stadio secondo il contratto attuale; nessun fallback silenzioso a solve separati.
6. Non modificare le regole fisiche di attrito, raffreddamento o pressione dei pori, e non alterare tableaux o tolleranze senza motivo documentato.

Per evitare regressioni, verificare almeno un sistema accoppiato velocità+temperatura (con dimensione e Jacobiana realmente controllate), i casi separati in due gruppi, i modelli legacy con solve congiunto e **le correzioni sia dopo SI sia dopo Newton**, con riferimento e ordine corretti; vedere D17.10.1.

### D17.10.1 — Correzioni numeriche associate ai blocchi di equazioni (vincolante)

La correzione attualmente chiamata `correct_friction_momentum` **non si seleziona in base alla famiglia o al modello reologico né alla presenza di uno specifico ID di termine sorgente**. È una correzione numerica associata alle **variabili/equazioni risolte**. Il suo nome storico non deve indurre a registrarla come proprietà della legge reologica. Un blocco può essere vettoriale, come le due componenti della quantità di moto, o scalare, come temperatura o pressione dei pori. La descrizione del blocco usa soltanto **indici di equazione** e, dove necessario, la dimensione del blocco: non deve generare sotto-solve Newton.

**Sequenza effettiva nella `dev` da preservare nel refactoring:**

1. Per il gruppo Newton, identificare tramite i metadati delle equazioni se il gruppo contiene il **blocco completo** al quale si applica una correzione; per l'attuale quantità di moto sono le due equazioni `q2` e `q3`. Non basta trovare una sola componente. La validazione preventiva deve rifiutare la separazione di un blocco che richiede un trattamento congiunto.
2. Se il percorso semi-implicito è attivo, assemblare `S_G^SI`, costruire il predictor `q_SI` a partire da `q_known` e applicare **una prima correzione** al blocco pertinente: `correct_momentum(q_known(momentum_ids), q_SI(momentum_ids))` (notazione illustrativa). Ricavare poi il contributo semi-implicito *effettivo* dallo stato corretto, come nel codice attuale, evitando di contabilizzare due volte l'incremento.
3. Costruire il residuo Newton sull'**intero insieme di equazioni del gruppo** usando lo stato post-SI corretto come riferimento appropriato, quindi eseguire **un solo Newton congiunto** sul gruppo. Immediatamente prima di Newton, salvare i valori delle componenti cui si applicherà la seconda correzione.
4. **Dopo il successo di Newton**, applicare **una seconda correzione** sul blocco di equazioni pertinente, usando come riferimento lo stato salvato *prima di Newton* e come stato da correggere il risultato Newton. Questa correzione post-Newton può essere pertinente **anche quando il percorso SI non è stato applicato**: non deve essere condizionata semplicemente dall'esistenza di un rate SI. In caso di fallimento Newton, conservare la gestione di errore e rollback esistente, senza correzioni su uno stato non convergente.
5. Le due correzioni hanno **riferimenti distinti** e non vanno fuse: la prima confronta `q_known -> q_SI`; la seconda `q_before_Newton -> q_after_Newton`. Conservare inoltre l'ordine corrente delle eventuali conversioni di stato, limitazioni e aggiornamenti dei contributi effettivi.

**Separazione delle responsabilità (decisione aggiornata):** `source_terms_2d` continua a occuparsi **soltanto dell’assemblaggio fisico** di `explicit`, `semi_implicit`, `implicit`, `analytic`. L’**implementazione** delle correzioni numeriche delle variabili viene collocata in `state_2d`. `nonlinear_solver_2d` conserva la responsabilità di **quando** chiamarle: dopo il predictor semi-implicito e dopo Newton, con i due riferimenti distinti sopra descritti. La selezione del blocco avviene tramite gli *indici delle equazioni* definiti nei metadati (non in base alla reologia); il solver non deve contenere condizioni fisiche come `group_has_momentum` con indici letterali disseminati nel codice. La routine storica `correct_friction_momentum` va migrata in `state_2d`, rinominandola in modo descrittivo e indipendente dalla reologia, ad esempio `correct_momentum_direction(old_momentum,new_momentum)`; opera sulle **due componenti conservative** e riceve lo stato prima e dopo. Non cambiare la sua formula né costruire un framework/registro di correzioni. Questa routine **non** deve essere eseguita automaticamente nella finalizzazione generica dello stato: richiede due stati e va chiamata solo nei punti numerici previsti.

**Portata verificata:** la descrizione dei due punti di correzione sopra deriva da `src/numerics/nonlinear_solver_2d.f90`, `solve_implicit_groups`, nel repository `dev` (consultazione 24/10/2026). Il criterio più generale per nuovi blocchi scalari/vettoriali è un requisito progettuale, non una capacità dichiarata già implementata. Prima di qualsiasi spostamento delle chiamate, caratterizzare con test la funzione esatta di `correct_friction_momentum` e assicurare la stessa sequenza e i medesimi risultati numerici nei casi legacy.

### D17.10.2 — Conservare il Newton indicizzato di v2/dev (vincolante)

**Non riscrivere il metodo Newton per renderlo genericamente indicizzato: lo è già nella sua parte numerica.** In `IMEX_SfloW2D_v2`, `implicit_map(:)` selezionava righe e colonne della Jacobiana, e in `IMEX_SfloW2D_dev` la medesima idea è già implementata mediante `group_equations(:)` / `active_equations(:)`. Il refactoring D17 deve **riutilizzare il solver della dev**, non ripristinare integralmente `solver_2d.f90` della v2 né creare un nuovo solver. La modifica consiste nel passare, per ciascun gruppo Newton, l'**intero insieme di indici del gruppo** una sola volta.

Conservare il calcolo del residuo sui soli indici attivi, la Jacobiana complex-step con tutte le derivate incrociate tra equazioni del gruppo, l'estrazione della sottomatrice e i risolutori lineari già disponibili (percorsi specializzati `2x2` e `3x3`, percorso denso generale per altre dimensioni), insieme a scaling, line search, criteri di convergenza e gestione dei fallimenti. La semplificazione elimina **solo** il ciclo di solve indipendenti *interni* a un gruppo e i relativi metadati ridondanti. Non reintrodurre indici fisici hard-coded nell'algoritmo Newton.

**Blocchi di correzione e gruppi non coincidono necessariamente.** Un gruppo Newton può comprendere un blocco vettoriale, ad esempio i due indici della quantità di moto, e una o più equazioni scalari: la correzione agisce sul blocco pertinente, ma **Newton risolve congiuntamente tutte le equazioni del gruppo**. Le correzioni non definiscono sottosolve, non cambiano il metodo e non impongono gruppi aggiuntivi. Le due applicazioni delle correzioni descritte in D17.10.1 restano distinte e nelle stesse posizioni temporali.

**Tracciabilità:** `IMEX_SfloW2D_v2/src/solver_2d.f90` (`implicit_map`, `solve_rk_step`, `eval_jacobian`) e `IMEX_SfloW2D_dev/src/numerics/nonlinear_solver_2d.f90` (`solve_newton_group`, `eval_f`, `eval_jacobian`).


**Aggiornamento 24 ottobre — commit `3b54bc01`:** nel runtime Newton il footprint numerico viene passato dall'insieme delle equazioni del gruppo; il descrittore legacy è controllato per corrispondenza esatta (inclusi i casi di liste permutate), ma non determina più la lista passata al solve. Il descrittore non è ancora stato eliminato, e servono test compile/registry e accoppiamento multi-target.

**Aggiornamento 24 ottobre — commit `075d3c95`:** gli assemblatori `eval_semi_implicit_terms` e `eval_implicit_terms` accettano ora il solo `group`, senza parametro `solver`; la scelta delle sorgenti dipende direttamente dai termini del gruppo. Il runtime passa il gruppo nelle chiamate SI, IM REAL e IM COMPLEX. Restano la migrazione delle firme interne del Newton, la rimozione del registro dei sotto-solver e i test scientifici/complex-step.

**Aggiornamento — commit `ab515430`:** tolto `source_solver_type` dalle firme interne di `solve_newton_group`, residuo e Jacobiana REAL/COMPLEX; il gruppo è l'unica selezione trasmessa nelle valutazioni del rate. Il livello di orchestration conserva ancora `solvers(:)` per convalida e deve essere migrato insieme al registro.

**Stato di avanzamento D17.10 — PARZIALE — Newton indicizzato e due correzioni preservati; prove 1x1/2x2/3x3 aperte.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.11 — Kernel REAL/COMPLEX

**Base già realizzata nel perimetro verificato:** riutilizzo del kernel COMPLEX della correlazione lahar dal wrapper REAL. Conservare la compatibilità complex-step e prevenire la duplicazione delle formule. Dove un termine selezionabile entra in Newton, verificare la disponibilità delle sue valutazioni necessarie al residuo e alla Jacobiana. Non imporre file `.inc` dove il riuso del kernel esistente è sufficiente.


**Stato di avanzamento D17.11 — PARZIALE — kernel REAL/COMPLEX mantenuti; Jacobiane da validare per legge.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.12 — Solve analitici

**Base già implementata:** interfaccia analitica in `source_terms_2d`, con target velocità, temperatura o pressione dei pori. Semplificare il contratto rimuovendo i riferimenti a un descrittore di sotto-solve indipendente quando la scelta è già univocamente determinata dal gruppo; mantenere l'intervallo fisico e l'aggiornamento del solo target. Rifiutare gruppi che richiederebbero due operatori analitici distinti.


**Stato di avanzamento D17.12 — PARZIALE — analitico e thermal flag indipendente implementati; bilanci quantitativi aperti.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.


**Aggiornamento 24 ottobre — commit `e4dda5ed`:** l'assemblatore analitico `eval_analytic_terms` riceve ora soltanto il gruppo e riconosce l'unico operatore da `group%equation_ids` (momentum 2/3, temperatura 4, pressione dei pori); non richiede il descrittore `source_solver_type`. Il runtime mantiene ancora il descrittore per i controlli del registro: la rimozione completa deve migrare anche SI/IM e Newton. Verifica del nuovo commit ancora da eseguire.
## D17.12a — Attivazione indipendente del raffreddamento termico (nuova decisione vincolante)

**Obiettivo.** L'attivazione della legge di raffreddamento della lava non deve dipendere dalla selezione della reologia lava. La reologia descrive la risposta meccanica; l'attivazione della sorgente termica è una scelta fisica distinta, configurata nei parametri della temperatura. Non introdurre nuove vie di integrazione o moduli dedicati: la legge termica conserva le quattro interfacce uniformi e i suoi contributi vengono assemblati da `source_terms_2d` nei gruppi ordinari.

**Situazione verificata sul `master` dev (da correggere):** `inpout_2d.f90` inizializza e valida i coefficienti di raffreddamento nel ramo `rheology_model == "lava"`; `equation_terms_2d.f90::init_problem_param` rende implicita l'equazione di energia attraverso la stessa condizione; `source_configuration_2d.f90` usa in almeno un percorso nativo l'identità della reologia lava per registrare i termini termici e il target energia; `thermal_laws_2d.f90::eval_rheology_cooling_energy_rate` dipende dal `model_id` reologico. Anche i fallback senza descrittore di gruppo in `equation_terms_2d` richiedono audit. La vecchia dipendenza è **obsoleta come regola di attivazione**.

**Contratto nuovo.**

1. Aggiungere un flag logico autonomo, proposto `thermal_cooling_flag`, in `&TEMPERATURE_PARAMETERS`, dichiarato nel modulo che possiede i parametri termici e letto/validato con gli altri parametri della temperatura. Il nome esatto può essere adattato alle convenzioni del codice, ma non cambiare il significato.
2. `thermal_cooling_flag = .TRUE.` attiva la legge termica e quindi i suoi termini sorgente nel catalogo/registro, **indipendentemente** da `rheology_model`; `FALSE` li disattiva anche quando la reologia lava è selezionata. La scelta di metodo (`EXPLICIT`, `NEWTON`, `ANALYTIC`) e cadenza (`STAGE`, `STEP`) resta esclusivamente nei gruppi, con i controlli di capacità già previsti in D17.
3. Sostituire con il flag termico i test di attivazione basati sulla reologia in: preparazione/validazione dei coefficienti, layout e attivazione dell'equazione termica nella partizione legacy, selezione termini e target nel costruttore nativo, assemblaggio e bridge analitico. **Non** eliminare selezioni reologiche dove servono davvero a definire la legge meccanica o le sue proprietà dipendenti da temperatura.
4. Rimuovere la dipendenza dal modello reologico dal dispatch del raffreddamento (`eval_rheology_cooling_energy_rate`): esporre una chiamata termica indipendente che riutilizzi la formula già presente per raffreddamento radiativo/convettivo; evitare copie della fisica. Il nome storico può essere mantenuto transitoriamente solo come wrapper di migrazione non più usato dai percorsi nuovi.
5. Validare subito dopo l'input che, quando la legge termica è attiva, esistano le equazioni/variabili energetiche richieste, che i coefficienti necessari siano fisicamente ammissibili, che i termini termici attivi siano registrati una sola volta e che il metodo/cadenza richiesti siano disponibili. Quando è disattivata, la reologia lava **non** deve attivare indirettamente sorgenti di raffreddamento.
6. Conservare il trattamento termico separato dagli altri termini energetici: disabilitare questa legge non deve disabilitare trasporto dell'energia, termini termici di iniezione o altri processi energetici indipendenti.

**Default e compatibilità:** definire e documentare esplicitamente la semantica dell'assenza del nuovo parametro prima dell'implementazione. Per non cambiare inavvertitamente i vecchi file di input, nel percorso di compatibilità legacy l'assenza del flag dovrebbe ricostruire l'attivazione storica (raffreddamento attivo con reologia lava); quando il flag è **esplicitamente fornito**, il suo valore prevale sempre sulla reologia, anche nel percorso legacy. Per nuovi input nativi, preferire una dichiarazione esplicita del flag o un default documentato e univoco. Non usare un default che silenziosamente riaccoppi le due scelte quando il flag è stato specificato.

**Casi minimi di regressione:** (a) reologia lava + flag attivo: confronto con baseline; (b) reologia lava + flag disattivo: nessuna sorgente radiativa/convettiva; (c) altra reologia + flag attivo: sorgente termica selezionata e integrata col metodo disponibile; (d) altra reologia + flag disattivo: nessuna sorgente termica; (e) assenza del flag in vecchi input; (f) gruppi `NEWTON/STAGE` e `ANALYTIC/STAGE` termici indipendenti dalla reologia, con verifica su residuo/Jacobiana o bridge analitico; (g) errori precoci per configurazioni termiche non supportate. Preservare la correttezza dell'input termico e le unità dei coefficienti.


**Stato di avanzamento D17.12a — IMPLEMENTATO (CONFIGURAZIONE) / CI VERDE — thermal_cooling_flag e matrice di fallback testati; bilanci aperti.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.13 — Default e compatibilità legacy

**Requisito obbligatorio:** in assenza di `SOURCE_GROUPS` il comportamento supportato deve coincidere con l'attuale baseline legacy (stessi termini attivi, trattamento SI/IM, ordine effettivo, finestre, rate, aggiornamenti e bilanci). Il costruttore legacy deve produrre gruppi che rispettino la nuova regola **un gruppo / un solve**:

- gruppi Newton legacy con target accoppiati → **un gruppo Newton con un solo sistema congiunto**;
- gruppi Newton legacy separati → gruppi Newton distinti;
- più target analitici legacy sotto una medesima label → **gruppi analitici distinti**, nell'ordine compatibile documentato;
- rate espliciti `STAGE` → gruppi espliciti secondo la baseline;
- batch di mass-exchange `STEP` → gruppo esplicito unico, conservando atomicità e posizione temporale.

La tabella dei default deve includere la nuova attivazione termica indipendente e il comportamento in assenza del relativo flag secondo D17.12a. Deve indicare condizioni di attivazione, metodo, cadenza, gruppo, equazioni target e precedenze per ogni combinazione fisica supportata. Confrontare sia il registro normalizzato sia l'effettiva evoluzione numerica. Non riprodurre una struttura di sotto-solve solo per imitare la nomenclatura legacy.


**Stato di avanzamento D17.13 — PARZIALE — mapping combinatorio legacy coperto dalla suite; confronto output aperto.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.14 — Trasferimenti di massa: applicazione del contratto fisico uniforme

**Stato verificato nel repository (non il risultato desiderato):** i trasferimenti di massa sono oggi espliciti e applicati dopo RK. `equation_terms_2d::eval_mass_exchange_terms` assembla rate conservativi, erosione/deposizione e modifica del letto; `mass_exchange_step_2d` gestisce un percorso post-RK specifico che valuta le sorgenti e applica l’incremento. Questo **percorso speciale è da eliminare come modalità di integrazione separata**. Devono essere conservate solo le operazioni di aggiornamento realmente necessarie, integrate nell'esecuzione generica del gruppo.

### D17.14.1 — Quattro interfacce fisiche uniformi

Come **ogni altro modulo di leggi sorgente**, `mass_exchange_laws_2d` deve offrire le **quattro interfacce uniformi** (nomi definitivi coerenti con il progetto). È un caso particolare del contratto generale, non una struttura speciale:

| Interfaccia di legge | Implementazione iniziale | Instradamento verso `source_terms_2d` |
|---|---|---|
| **Explicit** | Attiva: utilizza le formule esistenti | Assemblatore `eval_expl_terms` |
| **Semi-implicit** | Non implementata: errore esplicito se chiamata | Assemblatore `eval_semi_implicit_terms` |
| **Implicit** | Non implementata: errore esplicito se chiamata | Assemblatore `eval_implicit_terms` |
| **Analytic** | Non implementata: errore esplicito se chiamata | Assemblatore `eval_analytic_terms` |

La medesima regola vale per ogni modulo fisico. Quattro interfacce di valutazione **non significano quattro metodi configurabili**: restano i tre metodi `EXPLICIT`, `NEWTON` (che può utilizzare SI+IM) e `ANALYTIC`. Le interfacce non implementate **non** devono restituire zero né assumere che la sorgente sia inattiva; devono produrre un errore diagnostico se raggiunte. Il controllo primario, obbligatorio, è comunque **in fase di configurazione**: una scelta non disponibile viene rifiutata prima della simulazione, cosicché le interfacce non implementate non siano invocate nel normale percorso numerico.

### D17.14.1a — Organizzazione per processi fisici in `src/physics/mass_exchange/` (decisione vincolante)

Riorganizzare il mass-exchange in una **sottocartella specifica**, con **un modulo Fortran per processo fisico** e un unico modulo di coordinamento interno della famiglia. Non suddividere in file distinti le singole formule o le singole componenti conservative: la granularità è **il processo fisico**. La struttura proposta è:

```text
src/physics/mass_exchange/
  deposition_2d.f90          # deposizione dei solidi + perdita della fase continua indotta dalla deposizione
  erosion_2d.f90             # erosione e acquisizione dal substrato
  liquid_loss_2d.f90         # perdita della fase liquida
  darcy_loss_2d.f90          # perdita della fase continua legata alla pressione dei pori
  entrainment_2d.f90         # entrainment dell'aria
  vaporization_2d.f90       # vaporizzazione dipendente dalla deposizione
  mass_exchange_laws_2d.f90 # coordinamento e assemblaggio interno dei trasferimenti di massa
```

I nomi dei singoli file sono una proposta concreta, adattabile alle convenzioni Fortran esistenti senza modificarne le responsabilità. Migrare le routine oggi presenti in `src/physics/mass_exchange_laws_2d.f90`, quelle in `equation_terms_2d::eval_mass_exchange_terms` e le formule di vaporizzazione oggi in `thermal_laws_2d`, senza duplicare kernel e senza cambiare le formule. Aggiornare in modo coordinato dipendenze di build (`Makefile.am` e dipendenze moduli), import, test e manifest. Evitare dipendenze circolari: i moduli di processo calcolano rate e non richiamano l'assemblatore della famiglia né `source_terms_2d`.

**Deposizione e perdita della fase continua ad essa associata sono un solo processo fisico**, quindi un solo termine sorgente configurabile: `mass_exchange_deposition`. Il contributo della perdita indotta è calcolato nella stessa valutazione, ad esempio `L_dep = coeff_porosity * SUM(deposition_rate(:))`. `mass_exchange_carrier_deposition_loss` non deve più essere un termine fisico schedulabile separato: migrare il catalogo, le configurazioni legacy/native, i controlli di completezza e i test, gestendo eventuali nomi storici soltanto durante la compatibilità in lettura e senza creare una seconda sorgente attiva. Non confondere l'accorpamento fisico con la somma finale di perdite che condividono vincoli.

Il modulo **`mass_exchange_laws_2d`**, residente nella sottocartella, fornisce le quattro interfacce fisiche uniformi della famiglia richieste da D17. Inizialmente `explicit` è l'unica implementata: le altre rifiutano esplicitamente le chiamate e le configurazioni non supportate vengono respinte all'avvio. L'interfaccia `explicit`:

1. seleziona **solo i processi attivi appartenenti al gruppo richiesto**, utilizzando il catalogo/registro già esistente, senza creare un secondo registro;
2. richiede ai moduli dei processi i rate fisici necessari e ne conserva le dipendenze (deposizione → perdita indotta e vaporizzazione; erosione → inventario disponibile);
3. applica i **vincoli condivisi** nel punto corretto, in particolare **una sola limitazione complessiva** della perdita di fase continua richiesta da deposizione, perdita liquida e Darcy, senza limitatori autonomi incoerenti; documenta come l'eventuale saturazione viene attribuita alle componenti senza introdurre risultati ambigui;
4. costruisce **un solo vettore complessivo dei rate conservativi** e i corrispondenti rate accessori di inventari/topografia a partire dalle **stesse valutazioni**, applicando coerentemente maschere geometriche, coefficienti e limitatori;
5. non applica direttamente incrementi temporali, non modifica lo stato globale e non introduce uno scheduler speciale.

**Due livelli di assemblaggio, responsabilità diverse:** `mass_exchange_laws_2d` coordina i processi *interni al mass-exchange* e i loro vincoli condivisi; `source_terms_2d` è l'**unico assemblatore generale** di tutti i termini sorgente (inclusi quelli prodotti dall'interfaccia `explicit` del mass-exchange). Il secondo non deve ripetere le formule o i limitatori del primo. L'esecutore ordinario dei gruppi `STEP` applica atomicamente l'incremento conservativo e gli effetti accessori; l'eventuale codice di supporto per inventari e geometria non costituisce un percorso temporale separato.

**Configurazione e divisibilità:** l'esistenza di file distinti non autorizza automaticamente gruppi temporali indipendenti. Nella prima implementazione, i processi che condividono limitatori, rate dipendenti, inventari o aggiornamenti accoppiati restano nel **batch unico `STEP/EXPLICIT` legacy**, e la validazione rifiuta suddivisioni che ne compromettano la coerenza. Una futura separazione in gruppi richiede una semantica verificata per i rate condivisi, per la disponibilità di materiale e per l'ordine di esecuzione; non aggiungerla implicitamente come effetto del riordino dei file. Rimangono invariati i limiti dipendenti da `dt` (depositi disponibili, `erodible_inventory/dt`, fase continua, packing), finché un'estensione `STAGE` non sia progettata e validata separatamente.

### D17.14.2 — Un unico assemblaggio in `source_terms_2d`

`source_terms_2d::eval_expl_terms` deve includere i rate espliciti del mass-exchange **quando il termine appartiene al gruppo esplicito richiesto**, così come già fa per ogni altra sorgente. Non mantenere `eval_mass_exchange_terms` come canale di assemblaggio **parallelo** a `eval_expl_terms`: riutilizzare il suo contenuto fisico e i kernel esistenti per realizzare la normale interfaccia esplicita. Analogamente, i contributi SI, IM e analitici futuri confluiranno nei tre rispettivi assemblatori già previsti, senza un secondo dispatcher mass-exchange.

L'interfaccia esplicita deve restituire il vettore dei rate conservativi **e i dati accessori coordinati** indispensabili per erosione/deposizione, inventari e variazione del letto. Usare un risultato/contesto appropriato, disponibile ai gruppi che lo richiedono, senza forzare array di inventario inutili nelle firme di tutte le sorgenti. Il modulo fisico calcola i rate; il modulo sorgenti li assembla; l'esecutore del gruppo applica l'incremento e coordina gli aggiornamenti consistenti. Le maschere geometriche e i limitatori condivisi devono essere applicati in modo coerente a tutti gli output dello stesso gruppo; la loro collocazione esatta va verificata prima di rifattorizzare, mantenendo i bilanci fisici.

### D17.14.3 — Riproduzione del comportamento legacy tramite configurazione ordinaria

La configurazione di default che riproduce il codice attuale è:

- **Metodo:** `EXPLICIT`.
- **Cadenza:** `STEP`.
- **Gruppo:** un unico gruppo contenente il batch di termini mass-exchange fisicamente inseparabili.
- **Ordine:** **ultimo fra i gruppi `STEP`** eventualmente presenti (prima dell'avanzamento del clock).

Il driver deve eseguire i gruppi `STEP` in ordine, attraverso **lo stesso percorso generico usato per gli altri gruppi `STEP`**, senza una chiamata mass-exchange aggiuntiva o indipendente. Quando incontra il gruppo mass-exchange, ottiene i rate mediante `source_terms_2d::eval_expl_terms` e applica una volta `dt * rate` allo stato conservativo, insieme alle variazioni corrispondenti di inventari e topografia. I dettagli di conversione di stato, disponibilità di materiale, limitazione e aggiornamento geometrico rimangono quelli convalidati dal percorso storico, ma non giustificano una modalità temporale distinta.

**Equivalenza richiesta:** in assenza di altri gruppi `STEP`, il nuovo `STEP/EXPLICIT` del mass-exchange deve riprodurre la sequenza e i bilanci odierni; in presenza di altri gruppi `STEP`, quelli precedenti devono essere applicati prima, e il mass-exchange deve utilizzare lo stato risultante. Nessun doppio conteggio fra assemblatore esplicito degli stadi RK e assemblatore dei gruppi `STEP`: la stessa routine di assemblaggio è chiamata nei due contesti, ma valuta **soltanto i termini appartenenti al gruppo corrente**.

### D17.14.4 — Limiti deliberati

Al momento il mass-exchange supporta **solo `EXPLICIT/STEP`** come combinazione verificata. I metodi SI/IM/ANALYTIC non sono disponibili; anche `EXPLICIT/STAGE` non va accettato senza aver definito e validato come trattare coerentemente gli aggiornamenti accessori a ogni stadio RK. Si evita così di confondere la presenza di un'interfaccia `explicit` con la disponibilità automatica di ogni cadenza.

Il batch indivisibile e i limitatori condivisi devono continuare a essere rispettati. Non introdurre un metodo implicito o analitico fittizio per rendere l'interfaccia uniforme. La semplificazione richiesta riguarda il **percorso di assemblaggio e scheduling**, non una modifica delle leggi fisiche.


**Stato di avanzamento D17.14 — PARZIALE — processi mass exchange e limiter comune; equivalenza STEP aperta.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.15 — Sequenza del timestep

Sequenza production da preservare:

1. Il driver determina il `dt` effettivo applicando i vincoli esistenti.
2. Calcolo delle medie/fattori delle finestre temporali delle sorgenti.
3. Preparazione OU condizionale.
4. Avanzamento IMEX-RK con i gruppi `STAGE`.
5. Esecuzione una sola volta dei gruppi `STEP` in ordine, con il gruppo `EXPLICIT` mass-exchange **ultimo** nella configurazione legacy, senza percorso temporale speciale.
6. Avanzamento `runtime%t`.

È necessario testare le sorgenti circolari, fissurali, radiali e al bordo, inclusi periodi inattivi, rampe, finestre parziali e wrap, nonché `dt`, output e restart nei casi pertinenti. La verifica già svolta delle cinque fasi source-window/OU/RK/step/clock rimane evidenza storica nel suo perimetro.


**Stato di avanzamento D17.15 — PARZIALE — ordine RK/STEP presente e stato estratto; regressioni temporali aperte.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.15a — `time_integration_2d` indipendente dalla fisica e operazioni sullo stato in `state_2d` (vincolante)

**Obiettivo:** `time_integration_2d` deve coordinare tableaux, stadi, contributi spaziali, gruppi `STAGE` e `STEP` e composizione della soluzione, **senza incorporare controlli sulle specifiche variabili fisiche**. Il suo codice attuale (~750 righe) contiene ancora controlli sulla massa totale, sulle componenti solide, sulla temperatura e sulle celle asciutte. D17 deve centralizzare queste operazioni sullo stato senza cambiare formule o ordine dei calcoli.

**Destinazione e confini:**

- Centralizzare in **`src/core/state_2d.f90`** le operazioni numeriche di verifica e modifica dello stato. Introdurre o riutilizzare poche routine con nomi riferiti all'**operazione compiuta**, mai all'algoritmo chiamante. Nomi candidati: `validate_cell_state` (solo verifica e diagnostica), `correct_cell_roundoff` (correzioni di roundoff), `finalize_dry_cell_state` (stato asciutto, azzeramenti, proiezioni), `enforce_component_mass_bounds` (limiti e compatibilità delle masse). La scelta definitiva dei nomi dipende dall'inventario delle routine già esistenti e dall'esatta loro semantica. **Non introdurre `finalize_rk_cell_state` o altri nomi dipendenti da RK/IMEX.**
- **`state_conversion_2d` mantiene i kernel di conversione** conservativo/primitivo e di proiezione già implementati. Le nuove routine in `state_2d` li richiamano ove necessario; non duplicare `qc_to_qp` o le trasformazioni termodinamiche.
- Centralizzare in **`state_2d`** anche l'**inizializzazione dello stato primitivo asciutto**, inclusi gli stati delle celle dell'**halo di ricostruzione**. Fornire una routine dal nome indipendente dall'integratore, ad esempio `initialize_dry_primitive_state(qp_cell)`, che definisca per intero lo stato primitivo asciutto, inclusa la temperatura di riferimento attualmente pari a `T_ambient`, senza indici fisici o valori di riferimento codificati in `time_integration_2d`. Utilizzare la stessa routine sia per inizializzare gli stati primitivi asciutti delle celle sia per l'halo; preservare valori, condizioni di chiamata e numero/ordine delle inizializzazioni nel percorso numerico attuale.
- `time_integration_2d` richiama le routine di `state_2d` **negli stessi punti** in cui oggi esegue i controlli e le correzioni; non sostituire indiscriminatamente controlli di stadio e di fine step con un unico passaggio globale. Se è utile una routine di coordinamento dello stato, usare un nome neutro come `enforce_cell_state_consistency` e preservarne rigorosamente il contratto.
- Le maschere di dominio `B_nodata` sono informazioni geometriche, non leggi costitutive: mantenere il coordinamento al livello giusto, delegando l'azzeramento dello stato se opportuno senza introdurre nuove dipendenze inverse.

**Sequenza e invarianti da preservare esattamente:** registrazione dello stato grezzo tramite observer prima delle correzioni; controllo NaN; clipping di roundoff delle masse; controllo e azzeramento della massa totale quasi nulla (e errore per valori negativi non recuperabili); positività delle masse solide; conversione conservativo/primitivo; proiezione della quantità di moto nelle celle asciutte e sincronizzazione delle primitive; controllo della soglia della temperatura (attualmente 273 K); verifica e correzione ammessa del caso somma-masse-solide maggiore della massa totale con accesso allo stato precedente; trattamento delle celle `B_nodata`. Preservare anche le **soglie**, i messaggi/percorsi d'errore e l'assenza di cambiamenti nel trasporto delle masse o dell'energia durante le proiezioni. Valutare se alcune primitive debbano essere risincronizzate dopo correzioni, ma non modificarne implicitamente la sequenza durante il puro refactoring.

**Correzione direzione/verso della quantità di moto:** portare l'attuale `correct_friction_momentum` in **`state_2d`**, sotto un nome descrittivo dell'operazione (candidato `correct_momentum_direction`). Passare **blocco conservativo prima e dopo**, ad esempio due vettori di due componenti. La routine non deve conoscere Newton, SI, RK o il modello reologico; confronta e corregge i due stati secondo la formula esistente. Gli indici dei blocchi provengono da metadati e non da numeri codificati nel solver. `nonlinear_solver_2d` la chiama (1) **dopo SI**, prima del solve Newton, con `q_known` e predictor SI; (2) **dopo Newton riuscito**, con stato salvato subito prima di Newton e risultato Newton. Conservare il ricalcolo del contributo SI effettivo dopo la prima correzione. La seconda chiamata resta possibile anche quando SI non è attivo, se previsto dal percorso attuale. **Non** invocarla automaticamente in `validate_cell_state`, `finalize_dry_cell_state` o in una finalizzazione universale, perché il suo significato dipende esplicitamente da due stati consecutivi.

**Non confondere pulizia del modulo con eliminazione dei controlli:** tutte le verifiche numeriche a runtime restano obbligatorie. L'integratore può mantenere `n_eqns`, indici forniti come metadati, coefficienti RK, stato e operatori numerici; ciò che va rimosso sono i dettagli specifici di temperatura, solidi, quantità di moto e nomi delle leggi fisiche. L'attività è un **refactoring a comportamento invariato**, non una revisione della formulazione numerica.


**Stato di avanzamento D17.15a — IMPLEMENTATO (MIGRAZIONE PRINCIPALE) / CI VERDE — finalizzazione post-RK e dry halo/stage in state_2d; test stato/observer mancanti.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.16 — Input e validazione preventiva e unica

### D17.16.1 Regola temporale

**I parametri di configurazione dei termini sorgente devono essere controllati subito dopo la lettura e l'interpretazione completa dell'input, prima dell'inizializzazione del solver e di qualsiasi timestep.** Possono essere necessari i parametri fisici per determinare termini attivi e layout delle equazioni; non serve invece alcuno stato di cella o avanzamento temporale.

Sequenza richiesta: `read_param` → risoluzione parametri/modello ed elenco termini attivi → costruzione dei gruppi (nativo oppure adapter legacy) → **validazione completa** → congelamento del registro → inizializzazione del workspace numerico → simulazione.

### D17.16.2 Controlli obbligatori

1. Parser: namelist assente/valida, sintassi, dimensioni e limiti, nomi sconosciuti, duplicati e conflitti fra input nativo e `IMEX_*` legacy esplicitamente specificati.
2. Ogni termine **attivo e schedulabile** appartiene **esattamente a un gruppo**; nessuna omissione, doppia assegnazione o termine derivato indipendente.
3. Ogni gruppo ha ID, ordine, metodo e cadenza validi; nessuna collisione/ambiguità tra gruppi e fasi.
4. Un solo metodo e una sola cadenza per gruppo; nessun raggruppamento di termini con proprietà incompatibili.
5. Ogni termine supporta davvero la coppia `method/cadence` assegnata; in particolare, le coppie non ancora implementate (`STEP+NEWTON`, `STEP+ANALYTIC`) **non** devono essere accettate nel codice corrente.
6. Target/equazioni esistenti nel layout attivo, senza ripetizioni; raggruppamenti compatibili con le dipendenze e le leggi fisiche selezionate.
7. `NEWTON`: **un solo sistema** per gruppo, con tutti i target del gruppo; vietata qualunque lista di solve Newton indipendenti. Disponibilità di residuo, Jacobiana e contributi SI/IM richiesti.
8. `ANALYTIC`: **un solo operatore analitico** per gruppo, target e combinazione di termini supportati; vietata la concatenazione automatica di solve analitici distinti.
9. `EXPLICIT`: nessun solve Newton/analitico; contratto rate/stato al bordo e cadenza coerente.
10. Operazioni fisiche indivisibili, limitatori condivisi, accoppiamenti e batch di mass-exchange non spezzati; preservare integrità di inventari, componenti e topografia.
11. Conservare l'ordine legacy per i default nei percorsi che non devono cambiare; rifiutare input nativo incompleto anziché disattivare implicitamente una sorgente.
12. Messaggi di errore specifici: nome del termine, gruppo, metodo, cadenza e ragione del rifiuto.

Dopo la validazione, l'integratore **non** deve rivalidare continuamente gli stessi attributi nei loop di cella/stadio/Newton. Possono restare assert/guardie di sicurezza in punti di confine, ma non copie, parsing o ricostruzione del registro in hot loop. Restano invece obbligatori a runtime i controlli **fisici e numerici**: positività, finitezza, condizioni wet/dry, convergenza Newton, linear solver, ammissibilità e gestione dei fallimenti.

### D17.16.3 Input nativo

Riutilizzare la namelist opzionale `SOURCE_GROUPS` dove possibile; **semplificarla** eliminando i campi `source_group_solver_count`, `source_solver_order`, `source_solver_target_count`, `source_solver_targets` quando i target si ricavano in maniera univoca dai termini e dai dati del modello. Se alcuni target non sono ricavabili in modo non ambiguo (per esempio solve di velocità senza rate fisici attivi), mantenere **un singolo campo target a livello di gruppo**, non una lista di solve interni. Il formato esatto è una scelta implementativa da formalizzare prima di modificare il parser; i vecchi input nativi già distribuiti devono avere un percorso di migrazione esplicito o essere rifiutati con diagnostica chiara, non reinterpretati silenziosamente.

Esempio *concettuale* (non sintassi Fortran già accettata dal parser):

```text
GRUPPO 1: STAGE, NEWTON, termini reologici + termini termici compatibili, target VELOCITY+TEMPERATURE
GRUPPO 2: STAGE, ANALYTIC, termine di pressione dei pori, target PORE_PRESSURE
GRUPPO 3: STEP, EXPLICIT, batch mass-exchange completo
```

Non dichiarare questo esempio valido per qualunque reologia: la compatibilità di ogni termine con il metodo scelto deve risultare dal catalogo e dai test.


**Stato di avanzamento D17.16 — PARZIALE — parser e validazioni nella CI; equivalenza input numerica aperta.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

**Aggiornamento successivo — commit `04b2bc07`:** i quattro campi del sotto-solver sono stati rimossi dalla namelist `SOURCE_GROUPS` e dal parser nativo; gli input nativi storici che li contengono ricevono errore di namelist, mentre il percorso legacy senza `SOURCE_GROUPS` rimane separato. I target ora si specificano con `source_group_target_count` e `source_group_targets`. Il livello interno `source_solver_type` e il catalogo a ID numerici restano da migrare. La modifica richiede ancora verifica di build e test.

## D17.17 — Separazione fisica/numerica

Conservare l'isolamento dei moduli `*_laws_2d` rispetto a RK, gruppi e infrastruttura spaziale; l'assemblaggio dei contributi nelle equazioni rimane in `source_terms_2d`. La selezione fisica deve essere guidata dall'elenco termini del **gruppo**, non da identificatori di sotto-solve non necessari. Il solver legge residuo/rate dal bridge equazionale, non importa direttamente leggi fisiche. Completare audit statico e call graph, distinguendo le API di sorgente dalle API fisiche non-sorgente dell'operatore spaziale.


**Stato di avanzamento D17.17 — PARZIALE — equation_terms ripulito; assemblatori in source_terms e finalizzazione in state_2d; audit generale API aperto.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.18 — Adattamento delle routine esistenti

Minimizzare le modifiche alle formule: **uniformare tutti i moduli di leggi sorgente alle quattro interfacce** e adattare la selezione nei quattro assemblatori (`eval_expl_terms`, `eval_semi_implicit_terms`, `eval_implicit_terms`, `eval_analytic_terms`) alla nuova struttura semplice. **Integrare il mass-exchange in questi stessi assemblatori:** la sola interfaccia di legge attiva è ora `explicit`, che alimenta `eval_expl_terms`; le altre tre interfacce di `mass_exchange_laws_2d` rifiutano esplicitamente le chiamate. **Non conservare `eval_mass_exchange_terms` come percorso di integrazione distinto**, né `mass_exchange_step_2d` come orchestratore speciale: riutilizzare soltanto le operazioni di aggiornamento necessarie nell’esecutore ordinario dei gruppi `STEP`, senza spezzare il batch. Rimuovere eventuali chiamate duplicate, rispettare finestre e maschere, nessuna integrazione doppia fra `STAGE` e `STEP`, nessun cambiamento involontario di segno, conservazione, stato o Jacobiana.

**Riorganizzazione del modulo:** prima classificare tutte le routine pubbliche dell’attuale `equation_terms_2d` come *sorgente* o *non-sorgente*. Spostare soltanto le seconde verso il proprietario appropriato, mantenendo signature e comportamento ove possibile; rinominare quindi modulo e file in `source_terms_2d` e aggiornare import, compilazione, test e documentazione. Separare i commit/logical steps di spostamento/rinomina da quelli che modificano il comportamento dei gruppi, così da rendere le regressioni attribuibili. Non creare un ulteriore dispatcher generico.


**Stato di avanzamento D17.18 — PARZIALE — quattro interfacce per 24 sorgenti; dispatch ed equivalenza da validare.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## D17.18a — Procedura minima per aggiungere una sorgente fisica (vincolante)

### D17.18a.1 Principio e percorso di chiamata

Una sorgente che utilizza equazioni, target e metodi già esistenti deve richiedere modifiche soltanto **alla sua legge, a una voce del catalogo, al dispatcher della famiglia fisica, alla lettura degli eventuali parametri, al build e ai test**. `source_terms_2d` non deve avere `SELECT CASE` per i singoli modelli reologici: richiama il dispatcher reologico con il nome della sorgente; `rheology_2d` seleziona il modello e le sue interfacce. `source_configuration_2d` risolve automaticamente, dove appropriato, il nome `rheology_` + nome reologia; il registro deriva target, metodi ammessi e appartenenza al gruppo dal catalogo senza rami dedicati al modello. Il collegamento nel dispatcher è **statico, esplicito e diretto**: i nomi canonici delle reologie selezionano tramite `SELECT CASE` le corrispondenti routine Fortran. Non introdurre procedure pointer, factory, registri dinamici, generatori automatici del dispatcher o passaggi di generazione del codice nel build. La modifica manuale dei `CASE` pertinenti in `rheology_2d.f90` è parte intenzionale della procedura di aggiunta del modello.

```text
INPUT: rheology_model='newmodel'
  -> nome canonico 'rheology_newmodel'
  -> SOURCE_CATALOG: name, target='momentum', capabilities
  -> SOURCE_GROUPS: metodo NEWTON, cadenza STAGE, gruppo ordinato
  -> source_terms_2d: assemblatore generico, nessun CASE newmodel
  -> rheology_2d: dispatcher sul nome 'newmodel'
  -> rheology_newmodel_2d: due routine SI + IM per NEWTON
  -> nonlinear_solver_2d: un solo Newton sulle equazioni target del gruppo
```

### D17.18a.2 Esempio completo: `rheology_newmodel`

**A. Nuovo file** `src/physics/rheology_laws/rheology_newmodel_2d.f90`: contiene la fisica e le quattro entrypoint uniformi, con firme compatibili con il contratto della famiglia. In un modello Newton-only:

```fortran
MODULE rheology_newmodel_2d
  ! USE dei parametri fisici necessari
CONTAINS
  SUBROUTINE eval_newmodel_explicit(...)
    ! Non supportata: errore diagnostico, mai zero silenzioso.
  END SUBROUTINE

  SUBROUTINE eval_newmodel_semi_implicit(...)
    ! Richiesta sempre se has_newton=.TRUE.
    ! Restituisce il contributo SI o un contributo nullo intenzionale.
  END SUBROUTINE

  SUBROUTINE eval_newmodel_implicit(...)
    ! Richiesta sempre se has_newton=.TRUE.
    ! Rate implicito e valutazione compatibile con la Jacobiana complex-step.
  END SUBROUTINE

  SUBROUTINE solve_newmodel_analytic(...)
    ! Non supportata: errore diagnostico.
  END SUBROUTINE
END MODULE
```

Le firme e i nomi di esempio sono *illustrativi*, non sorgente compilabile; vanno derivati dalle interfacce uniformi consolidate in D17. Se le due rappresentazioni SI/IM corrispondono a termini della stessa legge, devono ricostruire la sorgente completa senza doppio conteggio. La correzione della direzione della quantità di moto resta in `state_2d`, non nella legge.

**B. Un'unica voce** in `src/numerics/source_catalog_2d.f90`:

```fortran
source_term_type( &
  name='rheology_newmodel', target='momentum', &
  has_explicit=.FALSE., has_newton=.TRUE., has_analytic=.FALSE.)
```

`has_newton` implica la disponibilità di **entrambe** le interfacce SI e IM. Il catalogo non contiene nomi di sub-termini, flag SI/IM separati, dettagli della Jacobiana o nomi delle formule.

**C. Dispatcher statico** `src/physics/rheology_laws/rheology_2d.f90`: aggiungere manualmente i `USE` del nuovo modulo e i `CASE` richiesti nelle routine di dispatch effettivamente supportate. Il dispatcher usa esclusivamente **nomi canonici** e inoltra le chiamate alla legge, senza valutare formule fisiche, decidere la disponibilità del metodo o eseguire validazioni specifiche del modello. La disponibilità di `EXPLICIT`, `NEWTON`, `ANALYTIC` è determinata esclusivamente dalle tre capabilities nel catalogo e verificata prima della simulazione. Le due entrypoint SI e IM fanno entrambe parte del percorso `NEWTON` (SI può restituire zero). Le routine di dispatch possono mantenere `SELECT CASE` distinti per rappresentazione, perché le firme REAL, COMPLEX e analitiche possono essere differenti: **non** introdurre un'interfaccia universale con argomenti opzionali o un flag di trattamento soltanto per avere un singolo `SELECT CASE`.

Esempio *illustrativo* del ramo aggiunto al dispatcher implicito esistente:

```fortran
SELECT CASE (TRIM(source_name))
CASE ('rheology_voellmy')
  CALL eval_voellmy_implicit_momentum_rate(...)
CASE ('rheology_newmodel')
  CALL eval_newmodel_implicit_momentum_rate(...)
CASE DEFAULT
  CALL fatal_error('Unknown rheology source in implicit dispatch')
END SELECT
```

Aggiungere analogamente il caso nel dispatcher semi-implicito di Newton; per interfacce non supportate, la validazione ne proibisce la selezione e un'eventuale invocazione errata fallisce esplicitamente. I nomi/firme nel frammento sono schematizzati: non sostituire le interfacce esistenti senza necessità. **Nessun dispatch reologico** in `source_terms_2d`, nessuna conversione nome–ID, nessuna seconda registrazione delle capabilities in `rheology_2d`.

**D. Input** `src/io/inpout_2d.f90`: leggere i nuovi coefficienti, fissarne i default, normalizzare il nome e richiamare dove appropriato un controllo dei vincoli fisici nel modulo della legge. La namelist Fortran rimane statica; evitare una piattaforma di input dinamico. L'attivazione fisica della reologia resta quella dell'input, non viene inferita dal semplice fatto che una voce esista nel catalogo.

**E. Build e test:** aggiungere il file e le dipendenze a `src/Makefile.am`; testare valori e limiti dei rate, SI nullo/non nullo, derivate complex-step quando richieste, risposta NEWTON e rigetto di EXPLICIT/ANALYTIC se indisponibili; eseguire test di compatibilità e regressione. Esempio nella cartella `EXAMPLES/` quando utile.

### D17.18a.3 Esempio di migrazione Voellmy

**Prima:** due voci distinte per gli identificativi numerici storici `SOURCE_TERM_RHEOLOGY_VOELLMY_COULOMB` e `SOURCE_TERM_RHEOLOGY_VOELLMY_TURBULENT_DRAG`, con ruoli semi-implicit/implicit nel catalogo. **Dopo:** una sola voce nominale `rheology_voellmy`, target `momentum`, con i tre booleani del metodo secondo le capacità *verificate*. Coulomb e turbulent drag rimangono nelle attuali routine `eval_voellmy_semi_implicit_momentum_rate` e `eval_voellmy_implicit_momentum_rate`, senza diventare sorgenti schedulabili; `NEWTON` usa entrambe in ordine e con il residuo corretto. Non attribuire capacità EXPLICIT o ANALYTIC senza test di disponibilità, né alterare le formule/limiti/correzioni esistenti. Applicare lo stesso principio ai componenti viscous/turbulent/yield di altre reologie: **una legge reologica completa = una sorgente configurabile**.

### D17.18a.4 File che devono restare invariati per una nuova reologia ordinaria

A refactoring D17 completato, l'aggiunta di una reologia sulle equazioni esistenti **non richiede modifiche** a `src/physics/source_terms_2d.f90`, `src/numerics/source_registry_2d.f90`, `src/physics/source_configuration_2d.f90`, `src/numerics/nonlinear_solver_2d.f90`, `src/numerics/time_integration_2d.f90`, `src/core/state_2d.f90` e `state_conversion_2d.f90`. Eccezioni solo se si aggiungono **nuove equazioni/variabili, nuovi target, nuovi metodi** o un vero nuovo contratto di aggiornamento fisico (per esempio inventari/topografia): non sono una semplice aggiunta di legge.

### D17.18a.5 Vincoli della validazione generale

Il validatore legge nomi, target e tre capabilities, controlla unicità, copertura dei termini attivi, numero/ordine dei gruppi, metodo-cadenza implementati, equazioni e un unico solve per gruppo Newton. **Non contiene `CASE` su singole sorgenti reologiche o su componenti Coulomb/turbulent.** Un gruppo `NEWTON` è ammesso soltanto se tutti i termini hanno `has_newton=.TRUE.`; ogni legge deve fornire SI e IM anche se SI produce zero. Un gruppo `ANALYTIC` richiede capacità e un unico solve analitico dell'insieme: **il supporto individuale non garantisce la composabilità analitica di più leggi**. Un gruppo `EXPLICIT` richiede `has_explicit=.TRUE.` e una cadenza attuabile. Qualsiasi vecchia configurazione che ripartiva le componenti interne della stessa reologia fra gruppi deve essere migrata verso la singola sorgente oppure rifiutata con diagnostica; non reinterpretarla silenziosamente.

**Test del dispatcher:** verificare che i nomi canonici inoltrino alle routine corrette in ciascuna rappresentazione supportata, che le capabilities non siano ricontrollate o replicate nel dispatcher e che eventuali chiamate a interfacce non supportate falliscano esplicitamente. La registrazione nel catalogo non determina automaticamente un indirizzo di procedura: il collegamento statico rimane esplicito e ispezionabile.

**Test di estensibilità obbligatorio:** aggiungere in un branch di test una reologia Newton-only fittizia seguendo A–E, verificare che compilazione e validazione funzionino modificando il **solo dispatcher statico** fra i moduli eseguibili preesistenti (oltre a catalogo/input/build), *senza alcuna modifica* ai moduli elencati in D17.18a.4, e verificare il rifiuto delle combinazioni non supportate. Tale esercizio rende misurabile la promessa di semplificazione.


**Stato di avanzamento D17.18a — DA INIZIARE (PROVA NUOVO MODELLO) — rheology_newmodel e rimozione ID reologici interni ancora aperti.** **CI:** due Actions verdi fino a `35be9f0a` confermate; test scientifici D17.20–21 ancora aperti.

## Guida di lettura e priorità di migrazione

**Fonte normativa:** D17.0–D17.18a definiscono il comportamento finale e non costituiscono una sequenza di commit. **Unica sequenza operativa:** D17.19, integrata dalla dipendenza vincolante di seguito. D17.20 definisce i test, D17.21 i gate finali, D17.22 la cronologia. Lo stato riportato in coda a ciascuna sezione tecnica deve essere aggiornato dopo i commit pertinenti: `DA INIZIARE`, `IN CORSO` oppure `COMPLETATA` (quest'ultima solo con evidenze di verifica).

**La migrazione è una sostituzione, non una manutenzione perpetua della struttura transitoria.** Il catalogo finale contiene esattamente `name`, `target`, `has_explicit`, `has_newton`, `has_analytic`. Un'unica sorgente rappresenta l'intera legge fisica (per esempio `rheology_voellmy`), non i contributi Coulomb e turbolento. Non restano `SOURCE_TERM_*` numerici come identità, `rheology_model_id`, `role`, `bundle_id`, `model_id`, `family`, conversioni nome-ID, `term_ids` in parallelo a `term_names` o un registro permanente dei sotto-solver. Sono invece consentiti gli indici numerici dei gruppi e delle equazioni.

**Vincolo sul gruppo:** `source_group_type` rappresenta direttamente il metodo, la cadenza, l'ordine, l'insieme dei nomi fisici e gli indici delle equazioni. Un gruppo Newton ha un unico sistema accoppiato; un gruppo analytic ha un unico operatore. `implicit_solver_ids(:)`, `analytic_solver_ids(:)`, `source_solver_type` e gli array `source_solver_*` dell'input sono **solo impalcature di transizione da rimuovere**, non API finali da cristallizzare con cardinalità uno. Il residuo e la Jacobiana Newton usano `group_equations(:)` e includono tutti i rate SI/IM del gruppo; conservare le due correzioni direzionali distinte ai punti prescritti in D17.10.1.

**Come eseguire i passi 3, 4, 5, 6, 9 e 10c di D17.19 senza doppie architetture permanenti:**

1. **Contratto e baseline:** preservare l'inventario legacy, i test dei default e la baseline numerica. Il gate GNU compile+registry raggiunto sul commit `decee60b` non vale come regressione scientifica end-to-end.
2. **Metadati di gruppo come proprietario unico:** rendere l'appartenenza delle sorgenti e gli indici delle equazioni autoritativi a livello del gruppo. Nei passaggi intermedi verificare l'uguaglianza *esatta* fra gruppo e vecchio descrittore helper; non aggiungere nuove funzioni o stati persistenti che dipendano dalle liste di sotto-solver.
3. **Catalogo nominale e accorpamento delle leggi:** introdurre le voci fisiche complete; migrare insieme configurazione, validazione, dispatch statico per nome, assemblatori e test. Conservare le formule nei moduli fisici. Eliminare gli ID e i metadata ridondanti una volta migrati tutti i chiamanti, senza sostituirli con un nuovo strato di lookup per ID.
4. **Input `SOURCE_GROUPS` minimale:** rimuovere solver count/order/targets per-solver, derivare il target dal catalogo; quando serve il caso senza sorgenti di momentum attive, prevedere soltanto un target opzionale *di gruppo*. Migrare o rifiutare esplicitamente la vecchia sintassi nativa; nessuna reinterpretazione tacita. Preservare il legacy input senza `SOURCE_GROUPS`.
5. **Eliminazione del livello sotto-solver:** dopo la migrazione dei chiamanti, eliminare `source_solver_type`, `solver_entries`, `implicit_solver_ids`, `analytic_solver_ids`, il passaggio `solver=` negli assemblatori e gli array duplicati di metadati; il solver indicizzato deve ricevere il solo gruppo e la sua lista di equazioni.
6. **Gate di accettazione della semplificazione:** test di gruppo Newton 1x1, 2x2 e 3x3, Jacobiana complex-step con derivate incrociate, gruppi distinti, analitico, input nativo/legacy, test negativo per sotto-solve, assenza di termini mancanti o doppi e prova `rheology_newmodel` di D17.18a. Nessuna cancellazione della vecchia rappresentazione prima di aver migrato ogni chiamante; nessun rinvio indefinito dell'eliminazione dopo il gate.

**Politica dei commit:** usare unità coerenti e verificabili, registrare gli SHA nelle sezioni interessate e in D17.22. Avviare le due Actions manualmente ai gate significativi, mai su ogni push. Un esito verde di compilazione o registro non certifica automaticamente l'equivalenza numerica.

## D17.19 — Piano implementativo aggiornato (ordine consigliato)

**Audit aggiornato a `d5e3c932` (due Actions verdi):** assemblatori in `source_terms_2d`; `equation_terms_2d` contiene cinque procedure e non riesporta la curvatura; correzione momentum, finalizzazione post-RK e inizializzazione asciutta di halo/stage sono in `state_2d`; coefficienti di iniezione in `mass_injection_coefficients_2d`. Le regressioni D17.20–21 rimangono aperte.

| Passo | Azione | Tipo | Accettazione mirata | Stato |
|---|---|---|---|---|
| 1 | Congelare inventario termini e baseline dei casi legacy; registrare test storici e SHA | Conservare | Catalogo e comparatori disponibili | PARZIALE — inventario CI verificato; baseline numeriche incomplete |
| 2 | Formalizzare contratto **un gruppo/un metodo/una cadenza/un solve**; eliminare la semantica dei multi-solve Newton | **Correggere** | Specifica, test negativi per multi-solve | COMPLETATO (STRUTTURA) / CI VERDE — un metodo e un solve per gruppo; test coupled Newton separati |
| 3 | Semplificare `source_group_type` e il registro, con un unico insieme di target e nessun scheduling interno | **Correggere** | Build, test di registry e immutabilità | COMPLETATO (STRUTTURA) / CI VERDE — no term_ids o sotto-solver |
| 4 | Semplificare il parser `SOURCE_GROUPS` e adapter legacy; validare subito dopo input | **Correggere** | Input validi, conflitti e messaggi negativi | IMPLEMENTATO (STRUTTURA) / CI VERDE — parser/legacy nominali; confronto end-to-end aperto |
| 5 | Adattare l'assemblaggio SI/IM a **tutti i termini del gruppo** | **Correggere** | Test dei rate per gruppo, REAL/COMPLEX | PARZIALE / CI VERDE — SI/IM in source_terms; test rate accoppiati aperti |
| 6 | **Riutilizzare il Newton indicizzato già presente nella dev** (`group_equations(:)`), senza riscriverlo: una sola chiamata Newton per gruppo con residuo/Jacobiana congiunti; preservare le **due correzioni numeriche** dopo SI e dopo Newton, riconoscendo i blocchi tramite indici di equazione (D17.10.1) | **Correggere** | Test sistemi 2x2, 1x1, 3x3, accoppiamento e correzioni nei due punti | PARZIALE / CI VERDE — Newton indicizzato e due correzioni; test quantitativi mancanti |
| 7 | Adattare bridge analitico all'operatore unico di gruppo | Adeguare | Un solo analytic; confronto legacy | PARZIALE / CI VERDE — bridge analytic presente, baseline numeriche aperte |
| 8 | Preservare routing esplicito `STAGE` e hook al bordo | Conservare/adeguare | Rate e BC senza duplicazione | IMPLEMENTATO (ROUTING) / CI VERDE — STAGE e boundary hook; test rate aperti |
| 8a | Separare il raffreddamento dalla selezione reologica: aggiungere `thermal_cooling_flag` in `TEMPERATURE_PARAMETERS`, migrare i controlli su layout, sorgenti, coefficienti e kernel, senza cambiare la legge fisica; preservare compatibilità dei vecchi input quando il flag è assente | **Correggere** | Matrice 2x2 reologia lava/non lava e raffreddamento on/off, parser e bilanci energetici | IMPLEMENTATO (CONFIGURAZIONE) / CI VERDE — flag e fallback, bilanci aperti |
| 9 | Uniformare **tutte** le famiglie di leggi sorgente al contratto fisico `explicit`/`semi_implicit`/`implicit`/`analytic`, conservando i kernel esistenti; registrare per ciascuna sorgente **un solo nome, un target, tre capabilities** (`EXPLICIT`/`NEWTON`/`ANALYTIC`), con SI+IM obbligatorie per NEWTON | **Correggere** | Test per ogni famiglia: dispatch coerente, disponibilità e rifiuto di interfacce non implementate | IMPLEMENTATO (96 INTERFACCE) / CI VERDE — dispatch quantitativo da verificare |
| 10 | Riorganizzare `src/physics/mass_exchange/` per processi, accorpare deposizione + perdita indotta, migrare i kernel esistenti e introdurre il coordinamento interno in `mass_exchange_laws_2d` senza duplicare i limitatori condivisi; integrare la sola interfaccia `explicit` in `source_terms_2d::eval_expl_terms`, eliminare il routing temporale separato e preservare il batch legacy come **ultimo gruppo `STEP/EXPLICIT`** | **Correggere** | Bilanci, test dei singoli processi, limitatore comune e dipendenze, conservazione dell’inventario, nessun doppio aggiornamento, corretto ordine STEP e build aggiornato | PARZIALE / CI VERDE — processi mass exchange e limitatore condiviso; confronto STEP mancante |
| 10a | Estrarre da `time_integration_2d` i controlli/correzioni fisico-specifiche in `state_2d` con nomi neutrali (`finalize_dry_cell_state`, controlli su masse/temperatura, roundoff), riusando `state_conversion_2d` e conservando sequenza/soglie; centralizzare anche l'inizializzazione delle primitive asciutte nell'halo | **Correggere** | Audit delle chiamate e test di equivalenza sugli stati intermedi/finali, celle asciutte, errori, observer e halo di ricostruzione | IMPLEMENTATO (MIGRAZIONE PRINCIPALE) / CI VERDE — finalizzazione post-RK e dry halo/stage ora in state_2d; test equivalenza stati, observer e casi asciutti ancora da eseguire |
| 10b | Migrare `correct_friction_momentum` in `state_2d` come correzione della direzione della quantità di moto a due stati; conservarne le due chiamate dal solver dopo SI e dopo Newton, tramite metadati di indici | **Correggere** | Test before/after di entrambi i punti, SI attivo/inattivo, sistema accoppiato e assenza di chiamate automatiche in finalizzazione generica | IMPLEMENTATO (MIGRAZIONE) / CI VERDE — correct_friction_momentum in state_2d; verifiche before/after aperte |
| 10c | Migrare catalogo, gruppi, input e dispatcher da ID numerici fisici a nomi canonici; accorpare le componenti delle reologie sotto un solo termine configurabile e rendere il validatore indipendente dai modelli | **Correggere** | Esempio `rheology_newmodel` D17.18a, test compatibilità e assenza di ID fisici | COMPLETATO (RIMOZIONE ID) / CI VERDE — nomi canonici, test rheology_newmodel aperto |
| 11 | Eliminare verifiche di configurazione ridondanti dagli hot loop, non i controlli numerici | Semplificare | Audit e prestazioni senza regressioni | IMPLEMENTATO (PULIZIA) / CI NUOVA PENDENTE — audit dei quattro moduli ed eliminazione doppie validazioni |
| 12 | Eseguire regressioni proporzionate al rischio su legacy/nativo, finestre, RK, OU, output, restart, ifx/CI pertinenti | Verificare | Evidenze e baseline approvate | PARZIALE / CI VERDE — build e registry, regressioni numeriche/restart/ifx da completare |
| 13 | Chiudere D17 solo quando D17.20–D17.21 sono soddisfatti | Gate | Decisione documentata | NON SODDISFATTO — gate scientifici D17.20–21 aperti |

La voce D17.2 richiede ora una **migrazione strutturale del catalogo e dei nomi**; le voci D17.9 routing esplicito, D17.11 riuso kernel, D17.12 bridge analitico e parti D17.15 sono **basi già realizzate**; il batch D17.14 esiste ma il suo routing speciale deve essere sostituito; non ripianificare da zero ciò che è corretto. Il cambiamento prioritario riguarda il significato dei gruppi Newton, il parser/registro e il residuo congiunto. L'implementazione di metodi ancora non supportati alla cadenza `STEP` non è una precondizione per chiudere il refactoring di compatibilità: serve una decisione esplicita di scope, altrimenti restano combinazioni vietate e validate come tali.

## D17.20 — Verifiche obbligatorie e matrice minima

### Contratto e input

- Migrazione nominale: assenza di ID numerici persistenti di sorgenti e reologie, unicità/canonicalizzazione dei nomi, lookup del catalogo statico e messaggi per nomi sconosciuti.
- Voellmy e altre reologie: una sola sorgente per legge completa; Coulomb/turbulent/viscous/yield rimangono componenti interne, con stessa fisica e stesso risultato dei casi legacy.
- Per ogni voce `has_newton`: presenza di entrambe le interfacce SI/IM, test esplicito del caso SI nullo, complessità complex-step e assenza di doppio conteggio.
- Test end-to-end `rheology_newmodel` senza modifiche ai cinque moduli numerici/assemblativi genericamente riutilizzabili (D17.18a.4), comprese configurazione e validazione preventiva.

- Input legacy senza `SOURCE_GROUPS`: riproduzione di gruppi e comportamento supportati.
- Flag `thermal_cooling_flag` in `TEMPERATURE_PARAMETERS`: on/off indipendente da reologia lava/non lava, con comportamento storicamente compatibile in assenza del flag e rifiuto anticipato di configurazioni termiche non supportate.
- Input nativo valido con più gruppi `NEWTON/STAGE` distinti.
- Gruppo Newton con **velocità e temperatura nello stesso solve**, dimensione 3 e Jacobiana completa, contro riferimento indipendente appropriato.
- Due gruppi Newton separati per velocità e temperatura; ordine verificato senza confonderli con il caso monolitico.
- Gruppi distinti `NEWTON/STAGE` velocità e `ANALYTIC/STAGE` temperatura.
- Un unico `ANALYTIC` per gruppo; rifiuto di più solve analitici nello stesso gruppo.
- Rifiuto di sotto-solve Newton multipli nello stesso gruppo e di duplicati/termini omessi.
- Rifiuto di metodo/cadenza non supportati, target incoerenti, processi fisici indivisibili spezzati, termini derivati schedulati, conflitti nativo/legacy e ordine ambiguo.
- Gruppo di velocità senza rate fisici attivi, se supportato dall'input legacy e dal layout, senza inventare un termine fittizio.

### Matematica/numerica

- Verificare separatamente **correzione dopo SI** (riferimento `q_known`) e **correzione dopo Newton** (riferimento post-SI/pre-Newton), anche con SI disabilitato dove pertinente; confrontare stati prima/dopo, condizioni di attivazione, ordine e assenza di doppio conteggio. Le correzioni dipendono dai blocchi di equazioni, **non dal modello reologico**.
- Raffreddamento attivo con reologia non lava e inattivo con reologia lava; confronto dei termini radiativo/convettivo, della partizione energetica, dei metodi disponibili e dei bilanci; nessuna attivazione implicita nascosta dipendente dal modello reologico.
- Confronto REAL/COMPLEX e complex-step sui sistemi accoppiati; controllare esplicitamente le derivate incrociate, la selezione esclusivamente tramite indici (`implicit_map` nella v2; `group_equations` nella dev), i percorsi lineari 2x2/3x3/denso e l’assenza di sottosolve Newton interni al gruppo.
- Stadi RK `N_RK=2,3,4` nei casi pertinenti, coefficienti/ordine invariati e non regressione di gruppi misti.
- Fonte al bordo unica, sorgente di curvatura, finestre circolari/fissurali/radiali (impulso parziale, rampa, wrap, inattività), OU quando attivo.
- Trasferimento di massa post-RK atomico: inventari, disponibilità, positività, conservazione, aggiornamento Q1 del letto e restart.
- Nuova struttura `mass_exchange/`: build dei moduli per processo; equivalenza dei rate di erosione, deposizione, perdita liquida, Darcy, entrainment e vaporizzazione rispetto alla baseline.
- Deposizione e perdita indotta come **un solo termine fisico**: nessuna seconda sorgente schedulabile, corretti contributi di massa/momento/energia/specie e mantenimento delle dipendenze verso vaporizzazione.
- Perdita complessiva di fase continua: confrontare il **limitatore unico post-somma** anche con deposizione, perdita liquida e Darcy simultaneamente attivi; verificare scenari in cui il limite è saturo e non saturo.
- Validazione del batch inseparabile: respingere separazioni fra gruppi che renderebbero incoerenti i rate, i limiti `dt` o gli aggiornamenti accessori; mantenere `STEP/EXPLICIT` come sola combinazione inizialmente accettata.
- Confronto default/nativo sull'intera matrice di combinazioni **effettivamente supportate**; `dt`, output, restart e bilanci nei test end-to-end selezionati.

- `time_integration_2d`: nessuna correzione hard-coded di temperature, solidi, quantità di moto o celle asciutte; inizializzazione dello stato primitivo asciutto e dell'halo tramite `state_2d`, con i valori e le convenzioni legacy invariati; le routine `state_2d` conservano observer, ordine delle operazioni, soglie, comportamento d’errore e stati Q/QP. Test specifici per cella asciutta, quasi asciutta, massa negativa per roundoff, massa solida superiore alla totale e temperatura sotto soglia.
- `correct_momentum_direction` in `state_2d`: verificare con due stati in ingresso la stessa correzione legacy di direzione/verso; preservare separatamente le due chiamate post-SI e post-Newton e i loro stati di riferimento, evitando l’attivazione automatica in una correzione generica.

### Gates e tracciabilità

Applicare il principio generale di regressione **risk-based**: T1 per contratto/validazione e gruppi, T2 per casi end-to-end rappresentativi, T3 su modifiche architetturali e gate, T4 solo ove richiesto dalla portata scientifica. Registrare SHA, toolchain, profilo strict/optimized, numero thread, casi, tolleranze, risultati ed esclusioni motivate. I test locali storici non costituiscono automaticamente accettazione del nuovo HEAD. Prima della chiusura verificare le piattaforme/CI pertinenti e approvare le baseline; non richiedere campagne scientifiche indiscriminate per ogni micro-modifica.

## D17.21 — Criteri di chiusura

**Gate aggiuntivo — estensibilità:** nessuna conversione nome–ID per sorgenti/modelli fisici, catalogo limitato a nome/target/capabilities, una voce per reologia completa, e prova `rheology_newmodel` descritta in D17.18a con il solo collegamento manuale nel dispatcher statico della famiglia, senza modifiche ad assemblatore, configurazione, registro o solver e senza procedure pointer, factory o generazione automatica del dispatch.


D17 è **COMPLETATO** soltanto quando tutti i punti seguenti risultano verificati nel perimetro approvato:

1. Ogni termine fisico attivo ha un nome canonico univoco e appartiene esattamente a un gruppo valido; i quattro contributi derivati non diventano sorgenti indipendenti.
2. Ogni gruppo ha metodo unico, cadenza unica, ordine definito e insieme di equazioni consistente.
3. **Ogni gruppo Newton è risolto come un unico sistema accoppiato contenente tutte le sue equazioni**, selezionate mediante indici come già nella v2 e nella dev; sono preservati la Jacobiana complex-step e i risolutori lineari esistenti; non esistono più multi-solve Newton indipendenti nello stesso gruppo.
4. **Ogni gruppo analitico contiene un unico solve analitico**; metodi diversi richiedono gruppi diversi anche a cadenza uguale.
5. I metodi/cadenze non supportati vengono rifiutati **subito all'avvio**; i default legacy restano equivalenti nei casi supportati.
6. **Ogni modulo di leggi sorgente** espone `explicit`, `semi_implicit`, `implicit`, `analytic`, con le capacità realmente disponibili dichiarate e validate all’avvio; le interfacce mancanti non simulano un contributo nullo. `source_terms_2d` assembla **tutti** i contributi attraverso i quattro assemblatori comuni. Mass-exchange usa oggi solo `explicit`, che alimenta `eval_expl_terms` e produce rate conservativi e accessori; nessun assemblaggio fisico alternativo post-RK.
7. I trasferimenti di massa sono organizzati in `src/physics/mass_exchange/` per processo fisico, con deposizione e perdita indotta unificate in un solo termine schedulabile; `mass_exchange_laws_2d` coordina rate e limitatori condivisi senza duplicare l’assemblaggio generale di `source_terms_2d`. Nessuna duplicazione delle sorgenti attive o dei rate e nessuna separazione di gruppi fisicamente incompatibile.
8. Le sorgenti `STAGE` seguono i tableaux; tutti i gruppi `STEP` seguono un **unico esecutore ordinato**, con mass-exchange `STEP/EXPLICIT` ultimo nella configurazione legacy e aggiornamento atomico coerente, senza doppie applicazioni.
9. I controlli sulla configurazione si svolgono una volta dopo l'input; i controlli fisico-numerici rimangono attivi a runtime. Le correzioni numeriche dei blocchi di equazioni (attualmente quantità di moto a due componenti) sono selezionate tramite indici, non tramite leggi reologiche; **sono conservati entrambi i punti applicativi, dopo SI e dopo Newton**, con riferimenti temporali distinti come in D17.10.1.
10. `time_integration_2d` non implementa direttamente correzioni su temperatura, masse solide, quantità di moto o celle asciutte; l'inizializzazione dello stato primitivo asciutto, compreso l'halo, è delegata a `state_2d`. Le operazioni sono in `state_2d` e richiamate nei punti previsti, con nomi che descrivono l’operazione e senza riferimenti a RK; `state_conversion_2d` conserva le conversioni. La correzione direzione/verso a due stati è in `state_2d`, chiamata da `nonlinear_solver_2d` nei due momenti previsti, senza variazioni numeriche.
11. I test mirati, le regressioni end-to-end, la compatibilità REAL/COMPLEX, i bilanci e le evidenze CI/platform necessarie al gate sono superati e documentati.
12. La roadmap e le evidenze distinguono capacità supportate, limitazioni deliberate e funzionalità future; il precedente requisito di multi-solve per gruppo è dichiarato **OBSOLETO**.

### Riferimenti di codice e documentali per l'attuazione

- `src/numerics/source_catalog_2d.f90`
- `src/numerics/source_registry_2d.f90`
- `src/numerics/source_input_2d.f90`
- `src/physics/source_configuration_2d.f90`
- `src/physics/source_terms_2d.f90` (nuovo nome; attualmente `src/physics/equation_terms_2d.f90`)
- `src/numerics/nonlinear_solver_2d.f90`
- `src/numerics/time_integration_2d.f90`
- `src/core/state_2d.f90` (routine di validazione, correzione e proiezione dello stato; verificare il percorso effettivo prima di implementare)
- `src/physics/state_conversion_2d.f90` (posizione indicativa da verificare nel repository; riutilizzare le API di conversione già esistenti)
- `src/physics/mass_exchange/` (nuovi moduli per processo e `mass_exchange_laws_2d.f90`; migrazione dall’attuale `src/physics/mass_exchange_laws_2d.f90`)
- `src/numerics/mass_exchange_step_2d.f90` (supporto storico da assorbire nell’esecutore ordinario STEP)
- `src/core/simulation_2d.f90`
- `docs/acceptance/D17_SOURCE_CONTRACT_2026_10_23.md` e `D17_SOURCE_TERMS_MANIFEST_2026_10_23.json`
- `docs/acceptance/D17_10_GROUP_FILTERED_IMPLICIT_RATES.md`, `D17_13_DEFAULT_CONFIGURATION.md`, `D17_14_STEP_ROUTING.md`, `D17_15_TIMESTEP_ORDER.md`, `D17_16_VALIDATION.md`, `D17_17_STEP_SOURCE_SEPARATION.md` (e altri report D17 della stessa cartella)

**Nota di supersessione:** i report di accettazione storici non vanno cancellati né riscritti retroattivamente; il nuovo documento diventa la fonte normativa per la semantica dei gruppi. I report precedenti sono evidenza storica del codice al momento del test, non approvazione del comportamento multi-solve ora dichiarato obsoleto.

---

## D17.22 — Registro di attuazione e ripresa del lavoro (aggiunta, 24 ottobre 2026)

> **NOTA DOCUMENTALE:** le sezioni D17.0–D17.21 che precedono questa aggiunta sono state conservate **integralmente e senza riscritture** dal file di riferimento `D17_Source_Groups_Revised_Specification_2026-10-24(1).md` fornito dall'utente. Questa sezione aggiuntiva registra fatti di implementazione e non modifica requisiti o priorità della specifica. In caso di divergenza, prevalgono i requisiti vincolanti D17.0–D17.21. Non usare la vecchia sezione D17 della roadmap generale come sostituto di questa specifica aggiornata.

### D17.22.1 — Repository, ramo, criteri di tracciabilità

- **Repository:** [demichie/IMEX_SfloW2D_dev](https://github.com/demichie/IMEX_SfloW2D_dev).
- **Branch esclusivo di implementazione:** `feature/d17-source-refactor`; non modificare `master` durante il lavoro.
- **Stato complessivo:** **IN CORSO / NON COMPLETATO**. Entrambe le Actions su `bebbf3a` confermate verdi dall'utente; la nuova estrazione dell'erosione richiede CI. «Commit pubblicato», «compilazione superata», «test del registro superati» e «equivalenza numerica dimostrata» sono stati distinti.
- **Regola per le nuove chat:** leggere prima questo documento per intero; poi controllare HEAD del branch, gli ultimi commit e gli ultimi run di GitHub Actions. Aggiornare *questo stesso file* dopo ogni passo significativo con attività effettuate, SHA, test ed eventuali problemi. Non cancellare né riassumere il piano originale.
- **Non introdurre architetture alternative:** il doppio storage `term_ids`/`term_names` è temporaneo; l'obiettivo rimane una sola rappresentazione nominale. Non aggiungere nuovi adapter, metadati o registri permanenti non presenti nei requisiti.

### D17.22.2 — Attività effettivamente svolte sul branch

| Commit | Intervento | Stato verifiche |
|---|---|---|
| [`fcb3d373`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/fcb3d373056d95b84bd0adcb977b09fe4599c9f2) | Lookup di catalogo tramite nome, controllo di unicità; `group_contains_name`. | Integrato; regressione del registro superata successivamente. |
| [`bd2d1cbc`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/bd2d1cbc58fd2c8040e5785ee60f300038ccdd85) | Introdotti `term_names` accanto ai vecchi `term_ids`. | **Solo transitorio**; da eliminare nel refactoring finale. |
| [`922cb5e2`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/922cb5e2dc363fb5e086d152f6a5d293bf19c875) | Compatibilità della validazione con descrittori legacy. | Transitorio. |
| [`d4ffc17a`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/d4ffc17a7dc2e3977cd3ab225b92edd0ed825c6c) | Test su nomi canonici, lookup, gruppi e solver. | Test registro CI successivamente superati. |
| [`27758c8f`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/27758c8fd8f159c884de2b91c3285f2f8fcc8947) | Capabilities `has_explicit`, `has_newton`, `has_analytic` e verifica metodi basata su tali campi. | Mappa transitoria costruita dalle vecchie informazioni; da riconciliare con il catalogo minimo definitivo. |
| [`bd3cd731`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/bd3cd731925d7d5b2ba34872f596d0a81a4e43cb) | Test matrice delle capabilities. | Test registro CI successivamente superati. |
| [`533cbe9d`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/533cbe9d8e5f769affbdb3a5f6bdbaafb52760f3) | Recupero per nome del descrittore nella validazione dei gruppi. | Compila e supera la regressione mirata. |
| [`447ce9a5`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/447ce9a57cca335c187c3cdf2cbd9bf7c40c441a) | Test round-trip nome/capabilities. | Test registro CI successivamente superati. |
| [`66ea646e`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/66ea646e6e94dccf89e7330170348ac5991677aa) | Selezione nominale di quattro sorgenti esplicite in `equation_terms_2d.f90`. | Compila; equivalenza end-to-end non ancora verificata. |
| [`9c908505`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/9c908505805391e09a2fa4baca8050d1f025b1d5) | Verifica appartenenza al gruppo Newton e lookup per nome nel percorso compatibile. | Primo build CI fallito per un import mancante, corretto di seguito. |
| [`fa1820a5`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/fa1820a53467c5db6fa846ca1b29bcdb1d517789) | Workflow `.github/workflows/d17-compile.yml`: sola compilazione e artifact dei log. | Il primo run ha rilevato il problema di import. |
| [`e3a9f743`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/e3a9f743faa4683c1804c4b81fdad08131d71078) | Import della funzione `source_term_descriptor_by_name` in `equation_terms_2d`. | **Compilazione GNU CI superata**. |
| [`61b9b707`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/61b9b7073c491b6b0150d4a1aca141508bece7b2) e aggiornamenti documentali successivi | Creato registro transitorio `D17_SOURCE_REFACTOR_HANDOFF.md`, aggiornato con evidenze. | Documento ausiliario, **non** sostitutivo di questa specifica. |
| [`788bbc70`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/788bbc704eec6c021ba16573710ee7e8f8245ba0) | Workflow `d17-source-registry.yml` per `TESTS/TEST_SOURCE_REGISTRY/run_test.sh` in `strict` e `optimized`. | Primo avvio fallito per permessi sullo script. |
| [`8b8bb486`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/8b8bb486a191ed6628162ea63a8f1f6c66125f11) | Workflow corretto con `sh TESTS/TEST_SOURCE_REGISTRY/run_test.sh`. | **Test registro CI superati**. |
| [`fa2a3187`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/fa2a3187fa3939b5ee04eeda708585b14ae9c3af) e commit documentali successivi | Aggiornato il registro con build e test superati. | Documentazione; nessuna accettazione D17. |
| [`e27bd9f2`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/e27bd9f22322ac3cb446f368f73d6daac8f5906c) e [`66e24375`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/66e2437530ec72ac3468cfd55b4cdd8298120992) | Copia storica D17 della vecchia roadmap e nota sulla differenza con l'handoff. | **Non usare** quella copia storica come specifica normativa; la base normativa è il testo caricato dall'utente, qui preservato. |

### D17.22.3 — Verifiche svolte ed evidenze

| Verifica | Stato | Evidenza e limiti |
|---|---|---|
| Compilazione GNU su Ubuntu 24.04 | **SUPERATA** | [Action run 40084588414](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/40084588414), commit `e3a9f743`, dopo correzione dell'import. Anche [run 40084676958](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/40084676958) passato sul commit documentale `e7a27f5`. |
| Test del registro delle sorgenti | **SUPERATI** | [Run 40085400556](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/40085400556), [run 40085421672](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/40085421672). Script esistente eseguito in `strict` e `optimized`. |
| Primo errore CI registrato | **RISOLTO** | [Run 40084188010](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/40084188010): simbolo `source_term_descriptor_by_name` non importato in `equation_terms_2d.f90`. |
| Primo errore CI del registro | **RISOLTO** | [Run 40085043424](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/40085043424): `Permission denied` dello script, prima di eseguire qualunque test. |
| Test numerici end-to-end, equivalenza legacy, Jacobiane REAL/COMPLEX, bilanci, combinazioni di raffreddamento | **NON VERIFICATI sul refactoring finale** | Non attribuire a CI di compilazione o test metadata proprietà che non ha verificato. |

### D17.22.4 — Stato per sezione della specifica

| Sezione normativa | Stato attuale | Da completare |
|---|---|---|
| D17.0–D17.1 e D17.0a | **PARZIALE** | Contratto uniforme e separazione completa responsabilità. |
| D17.2 | **PARZIALE** | Catalogo davvero minimo: nomi, target e 3 capabilities; eliminazione ID, ruoli, bundle e voci per sottocomponente. |
| D17.3–D17.7 | **PARZIALE** | Identità nominale, cadenze e gruppi definitivi; **un solo Newton accoppiato per gruppo** e un solve analitico per gruppo. |
| D17.8–D17.10.2 | **PARZIALE** | Rename `source_terms_2d`, quattro assemblatori, quattro metodi fisici; correzione direzionale nei due punti e mantenimento Newton indicizzato. |
| D17.11–D17.12 | **PARZIALE** | Revisione copertura kernel REAL/COMPLEX e analitici nel contratto definitivo. |
| D17.12a | **NON IMPLEMENTATO** | Raffreddamento termico attivabile indipendentemente dalla reologia, con default storico conservato. |
| D17.13 | **PARZIALE** | Configurazioni di default e confronto end-to-end con baseline. |
| D17.14–D17.14.4 | **PARZIALE** | Moduli `mass_exchange/` per processo; transazione atomica e limiter condiviso; schedulazione solo consentita. |
| D17.15–D17.15a | **PARZIALE** | Esecutore STEP comune, gestione dello stato e celle secche centralizzata in `state_2d`. |
| D17.16 | **PARZIALE** | Validazione preventiva una sola volta, messaggi d'errore per configurazioni non supportate. |
| D17.17–D17.18a | **NON COMPLETATO** | Moduli con confini definitivi; dimostrazione di nuova reologia senza cambiare generici. |
| D17.19 | **IN CORSO** | Seguirà l'ordine implementativo originale, aggiornando questo registro con commit. |
| D17.20 | **PARZIALE** | Compilazione e test registro superati; matrice completa scientifica/numerica non ancora eseguita. |
| D17.21 | **NON SODDISFATTO** | Non dichiarare chiuso D17 finché tutti i gate non sono documentati. |

### D17.22.5 — Prossime azioni e gestione delle versioni

1. Considerare **questa specifica integrale** il documento di riferimento: non tornare a strutture della roadmap generale del 15 ottobre quando incompatibili.
2. Eseguire la regressione numerica già presente sul branch e registrare commit/test/baseline prima del cambio strutturale.
3. Preparare un **cambiamento coordinato e compilabile** che unifichi le sorgenti fisiche reologiche, normalizzi il catalogo, sostituisca gli ID con nomi in input/registry/assembly/solver e aggiorni i test; evitare stratificazioni temporanee indefinite.
4. Proseguire nell'ordine dettagliato di **D17.19** con la semantica delle sezioni D17.0–D17.18a, verifiche D17.20 e gate D17.21.
5. Dopo ciascun commit importante, **aggiornare questa stessa D17.22** con SHA, modifiche, CI, regressioni e attività residue, **senza alterare il testo preesistente** salvo revisione esplicitamente richiesta.



### D17.22.6 — Politica di esecuzione manuale delle Actions (24 ottobre 2026)

Le due Actions D17 dedicate, `d17-compile.yml` e `d17-source-registry.yml`, devono utilizzare **esclusivamente `workflow_dispatch`**: nessuna esecuzione automatica ai commit o ai push. Compilazione e test vengono avviati intenzionalmente soltanto quando necessari a verificare un gruppo coerente di modifiche o un gate D17. Registrare qui l'esito dei run eseguiti; non interpretare l'assenza di nuove esecuzioni come validazione dei commit successivi.


### D17.22.7 — Avanzamento D17.19, passo 2 (24 ottobre 2026)

Commit [`079225a0`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/079225a0a8be366c5816f4d896b9b2ac34441bd6): il validatore del registro richiede ora esattamente **un solo solver Newton per gruppo**. Il test che prima accettava due solver indipendenti nello stesso gruppo è stato trasformato in test negativo. Non sono state cambiate formule fisiche, implementazione di Newton o sequenze RK.

**Stato:** prima correzione relativa al contratto D17.19 passo 2 pubblicata, **non ancora verificata con CI**. Restano da adeguare parser, struttura dei gruppi e runtime; più gruppi Newton distinti restano consentiti. Le Actions rimangono manuali.


### D17.22.8 — Validazione precoce dei solve nativi (24 ottobre 2026)

In `source_configuration_2d::initialize_native_source_configuration`, il parsing ora rifiuta immediatamente i gruppi `NEWTON` o `ANALYTIC` con numero di solver diverso da uno, senza attendere la validazione del registro. Il test `TEST_SOURCE_REGISTRY` aggiunge il caso di due solver Newton nel gruppo nativo e verifica il rifiuto transazionale. Non sono state modificate fisica, integratore o leggi. **Stato: commit pubblicato, Actions da avviare manualmente solo quando richiesto; test della modifica non ancora eseguiti.**

 
### D17.22.9 — Esecutore Newton a solve unico per gruppo (24 ottobre 2026)

Commit [`db750a5b`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/db750a5b5430b6e0e37c3b6d2013d79171ba2627): in `nonlinear_solver_2d.f90`, la gestione di un gruppo `NEWTON` verifica **esattamente un** `implicit_solver_id` e invoca una sola volta `solve_newton_group`, eliminando lo sweep `DO solver_index` ora vietato. Il solver mantiene la lista indicizzata delle equazioni, la valutazione SI, la correzione direzionale prima del Newton e quella successiva, nonché il solver lineare e la Jacobiana esistenti. Le modifiche sono limitate al controllo di cardinalità, alla scelta del primo/unico ID e alla rimozione del loop.

**Stato:** commit pubblicato; non sono ancora state avviate Actions, per rispettare la politica manuale D17.22.6. Resta da semplificare la struttura del gruppo e da completare la migrazione nominale. Prima del prossimo refactoring strutturale eseguire compilation e test pertinenti su richiesta, oltre alle regressioni numericali proporzionate al cambiamento.


### D17.22.10 — Unificazione della mappa delle equazioni Newton (24 ottobre 2026)

Commit [`521f752c`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/521f752c4b99e850f4d77e5ec45b5e5d740be29c): il percorso `NEWTON` di `nonlinear_solver_2d.f90` ora richiede che il descrittore dell'unico solver abbia esattamente tutte le equazioni del gruppo, controlla la copertura senza richiedere lo stesso ordine, usa `group_equations(:)` negli aggiornamenti SI e passa questa stessa lista al Newton indicizzato. Non modifica le routine della Jacobiana, le leggi fisiche né i tableaux. **Stato: implementato, non ancora verificato tramite compilazione/test dopo questo commit.** Restano da completare l'assemblaggio congiunto e la riorganizzazione dei descrittori secondo D17.19.

**Politica Actions:** compilazione e registro possono essere eseguiti manualmente al prossimo gate di verifica del gruppo di modifiche; nessun trigger `push` deve essere reintrodotto.


### D17.22.11 — Gate manuale compilazione e registro superato (24 ottobre 2026)

Le Actions sono state avviate manualmente tramite `workflow_dispatch` sul commit [`decee60b`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/decee60b7b48cc563e77e9569c9de0760a6348bf). **Entrambe terminate con successo**:

- [D17 compile diagnostics — run 40088178193](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/40088178193): compilazione GNU riuscita.
- [D17 source registry tests — run 40088314870](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/40088314870): regressione del registro `strict` e `optimized` riuscita.

Questo gate copre anche gli ultimi cambiamenti ai gruppi Newton fino al commit indicato (validatore, input, singola esecuzione e selezione delle equazioni). **Non dimostra ancora l'equivalenza numerica end-to-end né la correttezza della Jacobiana per gruppi Newton con più target**. I gate scientifici della sezione D17.20 restano aperti. Le Actions restano manuali, senza trigger automatico sui push.


### D17.22.12 — Copertura dei rate SI/IM di gruppo (24 ottobre 2026)

Commit [`67aab829`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/67aab8293ae288fc70877826fbfb4fd0f00e2712): il validatore dell'assemblatore `equation_terms_2d::validate_newton_selection` richiede ora che l'unico solver Newton includa **tutte e sole** le sorgenti fisiche presenti nel gruppo. In assenza di corrispondenza esatta, la selezione viene rifiutata prima di assemblare i rate SI/IM. I controlli restano compatibili con i descrittori legacy basati su `term_ids` durante la migrazione nominale. Nessuna formula fisica o Jacobiana modificata.

**Stato:** commit pubblicato; questa modifica **non è stata ancora compilata o testata**. Completare i test specifici di selezione delle sorgenti, poi avviare manualmente compile e registry Actions al prossimo gate coerente. Il successivo intervento dovrà continuare con il gruppo nominale unico previsto da D17.19, senza conservare permanentemente gli ID fisici.


### D17.22.13 — Chiarimento architetturale e ordine di migrazione (24 ottobre 2026)

Inserita la **Guida di lettura e priorità di migrazione** immediatamente prima di D17.19. D17.19 resta l'unico piano sequenziale; i passi 3/4/5/6/9/10c si devono eseguire come **una migrazione coordinata** verso il catalogo e il gruppo minimali, non come accumulo di strutture legacy con controlli aggiuntivi. La specifica tecnica D17.0–D17.18a è conservata integralmente, inclusi gli stati di avanzamento. Il gate compile+registry precedente rimane valido per lo SHA verificato; le nuove modifiche richiederanno un nuovo controllo.


### D17.22.14 — Prima semplificazione strutturale del validatore (24 ottobre 2026)

Commit [`8fc3e17a`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/8fc3e17abd6f45f2ef17672299cb9f61ff8cb3d6): nel validatore `source_registry_2d` rimossi il vettore `collected` e la maschera `equation_seen` che erano necessari ad aggregare più sotto-solve indipendenti. Ora l'unico descrittore implicito deve avere **esattamente lo stesso insieme di termini ed equazioni del gruppo**, verificato tramite `same_set`. Restano temporaneamente `source_solver_type` e gli ID nel formato del registro: saranno eliminati nell'unità coordinata di migrazione nominale definita nella guida prima di D17.19; non sono il design finale.

**Verifica:** implementato, non ancora verificato tramite le Actions sul nuovo HEAD. La baseline compile+registry verde del precedente SHA rimane documentata in D17.22.11; le regressioni numeriche restano aperte. Nessuna modifica a formule fisiche o schemi RK.


### D17.22.15 — Validazione del singolo solve senza ciclo interno (24 ottobre 2026)

Commit [`0d9d444f`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/0d9d444f2dd9d505591b13fcd38701071a0b143d): eliminata dal validatore del registro l'iterazione sui sotto-solver del gruppo, sostituita da una singola validazione per gruppo non esplicito. Rimossi anche un controllo di inclusione delle equazioni ormai ridondante rispetto alla successiva verifica di uguaglianza degli insiemi. Non sono cambiate le formule fisiche né il calcolo Newton.

**Stato:** D17.5 ancora IN CORSO. Il registro conserva temporaneamente `source_solver_type` e gli ID di solve, che la specifica obbliga a rimuovere dopo la migrazione coordinata dei chiamanti. Il nuovo commit non ha ancora un esito compile/registry: i test verdi precedenti non vengono attribuiti a questo HEAD. I workflow restano manuali.


### D17.22.16 — Validazione diretta del descrittore unico (24 ottobre 2026)

Commit [`5c01d13d`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/5c01d13d49d562a395e5b805e2d4ffc4321bd09b): in `source_registry_2d::validate_registry` eliminati l'array locale aggregato `solve_ids` e il controllo di inclusione dei termini già garantito dall'uguaglianza esatta dei due insiemi. Il descrittore unico è selezionato direttamente dall'ID Newton oppure analytic del gruppo; l'ordine del descrittore deve essere 1. Si tratta di rimozione di codice di transizione nel validatore, non dell'eliminazione definitiva degli ID e di `source_solver_type`, ancora utilizzati dai chiamanti. Nessuna formula o integrazione numerica modificata.

**Stato:** D17.5 IN CORSO; questo commit non è ancora stato verificato da CI. I gate manuali precedenti si riferiscono a uno SHA antecedente. La rimozione definitiva delle strutture duplicate richiede ancora la migrazione coordinata di parser, runtime e catalogo prevista prima di D17.19.


### D17.22.17 — Costruzione diretta dei gruppi impliciti nativi (24 ottobre 2026)

Commit [`7ea4969c`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/7ea4969c2053d79a202242c93706b8a5357951e2): in `source_configuration_2d` rimossi i due cicli che scandivano i sotto-solver nei gruppi impliciti nativi, sia durante la derivazione del partizionamento delle equazioni sia durante la costruzione del descrittore. Poiché la validazione preliminare richiede già esattamente un solver, l'implementazione usa direttamente il primo e unico elemento, con ordine 1, e assegna al gruppo l'intera lista delle sue equazioni. Il formato `SOURCE_GROUPS` e il percorso legacy restano per ora immutati; la rimozione dei campi solver-specifici dall'input è ancora prevista nella migrazione coordinata.

**Verifiche:** nuovo commit non ancora testato con le due Actions. D17.16 rimane IN CORSO; l'eliminazione finale degli ID e di `source_solver_type` non è ancora realizzata. Nessun cambiamento alle formule o ai tableaux RK.


### D17.22.18 — Stato operativo nella tabella D17.19 (24 ottobre 2026)

Aggiunta la colonna **Stato** nella tabella D17.19, con categorie DA INIZIARE / IN CORSO / COMPLETATA. Nessun passo è dichiarato COMPLETATO senza tutti i controlli di accettazione. Le attività di rimozione dei campi sotto-solver dall'input, eliminazione di `source_solver_type` e conversione del catalogo alle sole identità nominali restano un **unico programma coordinato**: aggiornare parser, adapter legacy, registry, solver, assemblatori, test e build prima di dichiarare conclusa la migrazione. Ogni avanzamento successivo deve aggiornare la colonna e le sezioni tecniche pertinenti.


### D17.22.19 — Primo passaggio coordinato verso l'input a target di gruppo (24 ottobre 2026)

Commit [`a4f18652`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/a4f186524ea09413ea79680a2b91de75d266f70c): `source_input_2d` ammette i nuovi campi `source_group_target_count` e `source_group_targets`; `source_configuration_2d` traduce questi valori nell'attuale descrittore di gruppo/solver e rifiuta esplicitamente la mescolanza di target di gruppo con `source_group_solver_count`. Il formato nativo precedente rimane un percorso legacy esplicito in questa fase, non la soluzione finale. `TEST_SOURCE_REGISTRY` ora verifica una configurazione con target di gruppo e il conflitto dei due formati.

**Stato reale:** D17.16 e D17.19 passo 4 restano IN CORSO; la rimozione dei quattro vecchi campi dalla namelist, l'eliminazione di `source_solver_type` e la migrazione del catalogo agli identificativi esclusivamente nominali **non sono ancora state implementate**. Non sono stati modificati kernel o schemi RK. **Verifiche:** nuovo codice non ancora compilato/testato; le Actions rimangono manuali. Questa sezione non costituisce un gate di accettazione della migrazione completa.


### D17.22.20 — Controlli input nativo (24 ottobre 2026)

Commit [`b3549e79`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/b3549e79a131dfd7d497b09a25ecbfe6109bbab0): il parser ora rifiuta metadati legacy di sotto-solver residui quando si usano target di gruppo, anche con conteggio solver zero; rifiuta target incompleti e orfani. Aggiunti test negativi. D17.16 resta IN CORSO: i vecchi campi, source_solver_type e gli ID fisici devono ancora essere rimossi. Nessuna modifica alla fisica. Verifiche di compilazione e test sul nuovo HEAD ancora non eseguite.


### D17.22.21 — Controlli sui metadati residui della sintassi legacy nativa (24 ottobre 2026)

Commit [`5349d89c`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/5349d89c9faf2ace73f69d524c445dbc78c099e9): nel parser `source_configuration_2d` vengono ora rifiutati gli attributi `source_solver_order`, `source_solver_target_count` e `source_solver_targets` quando manca il relativo solver legacy, nonché gli attributi del secondo sotto-solver ormai vietato. Aggiunti test negativi in `TEST_SOURCE_REGISTRY`. La transizione a target di gruppo non modifica le formule fisiche.

**Stato:** D17.16 IN CORSO. Il formato legacy nativo rimane temporaneamente supportato per migrazione; rimozione definitiva dei campi, di `source_solver_type` e degli ID fisici ancora da completare. **Verifica:** commit non ancora coperto da nuove Actions, che rimangono manuali.


### D17.22.22 — Rimozione definitiva dei campi sotto-solver da SOURCE_GROUPS (24 ottobre 2026)

Commit [`04b2bc07`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/04b2bc074007e6947826a1ea1b3607de9d8dd26b): rimossi dalla namelist e da `source_input_2d` i quattro campi `source_group_solver_count`, `source_solver_order`, `source_solver_target_count`, `source_solver_targets`, nonché il limite `SOURCE_INPUT_MAX_SOLVERS`. Il parser `source_configuration_2d` deriva direttamente l'unico descrittore interno dai target del gruppo. Aggiornati i test di inizializzazione, di configurazione e della namelist, con casi negativi per target mancanti, sconosciuti o duplicati. La vecchia sintassi **nativa** con campi sotto-solver è ora intenzionalmente non valida; il precedente input IMEX legacy, che non usa `SOURCE_GROUPS`, mantiene il suo percorso distinto. I target possono ancora essere elencati esplicitamente a livello del gruppo; l'inferenza automatica dagli attributi del catalogo verrà completata con la migrazione nominale.

**Stato:** la rimozione dei campi dall'input è **implementata ma non ancora verificata con le due Actions**. D17.16 e D17.19 passo 4 restano IN CORSO finché parser, migrazione nominale e verifiche non soddisfano l'intero contratto. `source_solver_type` rimane temporaneamente una struttura interna e non è stato dichiarato eliminato. Nessuna modifica alle formule delle leggi fisiche o agli schemi RK.


### D17.22.23 — Il gruppo determina le incognite del Newton (24 ottobre 2026)

Commit [`3b54bc01`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/3b54bc01ba3de2ae2d71144aea05852e62b6936d): il runtime `nonlinear_solver_2d` passa a `solve_newton_group` gli indici delle equazioni del gruppo, non la copia nel descrittore del solver. Prima della chiamata richiede l'uguaglianza esatta degli insiemi e considera il momentum dal gruppo. Il controllo è transitorio e impedisce che un sotto-solver definisca un sistema numerico differente. **Non** è la rimozione definitiva di `source_solver_type`: il descrittore resta necessario alle attuali interfacce SI/IM, Newton e analytic, da migrare insieme. D17.5 e D17.10 restano IN CORSO. Nuovo HEAD non ancora compilato o testato con le Actions manuali; nessuna formula fisica o RK modificata.


### D17.22.24 — Analytic guidato dal gruppo (24 ottobre 2026)

Commit [`e4dda5ed`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/e4dda5ed5e4f0ffde6e38038b1c737bd0faf457e): `eval_analytic_terms` non riceve più `source_solver_type`; l'operatore viene selezionato dalla dimensione e dal contenuto di `group%equation_ids`. Aggiornata la chiamata in `nonlinear_solver_2d`. Le formule analitiche rimangono invariate. **Non** è ancora stata eliminata la struttura globale `source_solver_type`, necessaria alle altre interfacce. D17.12 e D17.5 sono IN CORSO per la migrazione coordinata. Compile/registry e regressioni analitiche sul nuovo SHA sono ancora da eseguire.


### D17.22.25 — Assemblatori SI/IM senza selettore del sotto-solver (24 ottobre 2026)

Commit [`075d3c95`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/075d3c95d265abc179ee35eadf7bdd387a1bf069): in `equation_terms_2d` sono state eliminate le dipendenze da `source_solver_type` per `eval_semi_implicit_terms` ed `eval_implicit_terms`, insieme ai relativi helper di selezione; i termini di gruppo sostituiscono il doppio filtro gruppo/solver. Adeguati i tre punti di invocazione in `nonlinear_solver_2d`. I kernel delle leggi fisiche non sono modificati. **Non è la rimozione definitiva di `source_solver_type`**: rimangono firma interna Newton, registry, adapter legacy e configurazione nativa. **Verifiche:** commit non ancora compilato né testato sulle Actions manuali; D17.5, D17.8, D17.10 ancora IN CORSO.


### D17.22.26 — Firme interne Newton senza descrittore solver (24 ottobre 2026)

Commit [`ab515430`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/ab51543037946bf92fe0b69b5617fa87539764e3): eliminate le argomentazioni `source_solver` dalle routine interne di Newton, residuo e Jacobiana in `nonlinear_solver_2d`. Le equazioni rimangono definite dal gruppo e il percorso REAL/COMPLEX passa ai rate IM con il solo gruppo. Il vettore di descrittori `solvers(:)` rimane ancora nella API esterna e nei controlli di transizione: **non è stata ancora eliminata la definizione `source_solver_type` dal registro**. Non sono state cambiate formule fisiche o metodi RK. **Verifica:** nuova compilazione e test Actions ancora da eseguire sul commit; D17.5 e D17.10 IN CORSO.


### D17.22.27 — Rimozione di solvers(:) dal percorso numerico esterno (24 ottobre 2026)

Commit [`1d916bcd`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/1d916bcdad552f7c08e23e13a2679626d3235c46): eliminate le argomentazioni `solvers(:)` da `solve_imex_stage`, `solve_implicit_stage`, `solve_rk_step` e `solve_implicit_groups`, e la copia `source_solvers=source_registry%solvers()` nell'integratore temporale. I due percorsi Newton e analytic consumano gli indici delle equazioni dei gruppi validati. I kernel fisici, le formule SI/IM e le due correzioni di quantità di moto non sono stati modificati. **Stato:** D17.5 IN CORSO: rimangono `source_solver_type`, `solver_entries`, `implicit_solver_ids` e `analytic_solver_ids` in registro/configurazione e i relativi test da adeguare. **Verifica:** Actions e regressioni ancora da eseguire sul nuovo HEAD; non attribuire a questo commit i test verdi del precedente SHA `19380916`.


### D17.22.28 — Eliminazione dei descrittori di sotto-solver dal registro (24 ottobre 2026)

**Intervento strutturale:** `source_registry_2d` non definisce più `source_solver_type`, né memorizza `solver_entries` o gli array `implicit_solver_ids` / `analytic_solver_ids` nei gruppi. `initialize`, `initialize_legacy` e `make_group` costruiscono e validano direttamente gruppi con i propri `equation_ids` e `term_ids` temporanei. La validazione della forma dell'operatore analitico e dell'accorpamento della legge reologica ora usa solo il gruppo. In `source_configuration_2d` non si costruiscono più descrittori `solvers`; `TEST_SOURCE_REGISTRY` controlla l'assegnazione delle equazioni e i casi negativi direttamente sui gruppi, sostituendo i precedenti test sugli ID dei sotto-solver.

**Limiti:** l'eliminazione degli ID numerici delle *sorgenti fisiche* e la migrazione nominale del catalogo sono attività separate ancora aperte in D17.2 / D17.19 passo 10c. Il contratto di scheduling è già a livello di gruppo ma l'equivalenza scientifica multi-target Newton e REAL/COMPLEX richiede i gate D17.20–D17.21. **Verifica corrente:** nessuna nuova Action eseguita; la baseline verde precedente resta legata allo SHA `c5dc9e9a`. D17.5 e D17.19 passi 3–4 sono IN CORSO in attesa di compilazione e test.


### D17.22.29 — Selezione nominale delle sorgenti negli assemblatori SI/IM (24 ottobre 2026)

Nell'assemblatore `equation_terms_2d` la selezione dei contributi impliciti per ruolo/modello usa prioritariamente `group%term_names` e `source_term_descriptor_by_name`, mentre il controllo della presenza del singolo termine risolve il suo nome canonico e interroga `group_contains_name`. Il fallback a `term_ids` è confinato ai descrittori legacy non ancora normalizzati. Questa modifica non elimina ancora gli ID numerici dal catalogo, né cambia i kernel fisici e i metodi numerici; è un passo preparatorio della migrazione D17.2/10c. **Stato:** IN CORSO; commit da verificare con compilazione e test sul nuovo HEAD.


### D17.22.30 — Membership nominale anche per le API di compatibilità (24 ottobre 2026)

`group_contains(group,term_id)` risolve prima l'identità nominale dal catalogo, quindi delega a `group_contains_name(group,name)`: un ID numerico non è più usato come chiave per verificare l'appartenenza a un gruppo normalizzato. Rimane il fallback ai `term_ids` soltanto per gruppi grezzi pre-normalizzazione. Il catalogo conserva temporaneamente gli ID e le suddivisioni delle leggi reologiche; **D17.2 e D17.19/10c non sono completati**. Nessun cambiamento ai kernel fisici. Verifica con le Actions necessaria sul nuovo commit.


### D17.22.31 — Validazione dei nomi canonici nel registro (24 ottobre 2026)

Il validatore `source_registry_2d` rifiuta ora nomi canonici vuoti o duplicati all'interno di un gruppo, oltre a verificare la coerenza con la rappresentazione numerica transitoria. Aggiunte due regressioni negative al test del registro. Questo controllo prepara la rimozione degli ID numerici: **non** costituisce ancora la migrazione completa del catalogo o l'accorpamento dei contributi delle leggi reologiche. D17.2 e D17.19 passo 10c rimangono IN CORSO. Non ancora verificate compile/registry Actions sul nuovo HEAD.


### D17.22.32 — Identità canonica univoca delle leggi reologiche (24 ottobre 2026)

Aggiunta al catalogo `rheology_source_name(model_id)` che associa ciascuna delle dieci famiglie reologiche supportate a un'unica identità nominale (`rheology_voellmy`, `rheology_lahar`, ecc.), indipendente dal numero di addendi della relativa legge. I test del registro verificano che i nomi esistano, siano distinti e corrispondano a leggi legacy con contributi fisici. Questa è la definizione dell'identità per la successiva migrazione: **non** unifica ancora le voci schedulabili nel catalogo né altera le liste legacy `rheology_source_terms` o la fisica. `source_term_type` rimane transitorio, D17.2 / D17.19 passo 10c IN CORSO. Build e Actions sul nuovo commit ancora da verificare.


### D17.22.33 — Nomi reologici canonici nell'input SOURCE_GROUPS (24 ottobre 2026)

`source_configuration_2d` accetta ora nei gruppi nativi una singola dichiarazione del nome canonico della legge reologica configurata (ad esempio `rheology_lava`) e la espande nella lista completa dei contributi legacy della stessa legge. Il conteggio delle dichiarazioni di input è distinto dal conteggio delle componenti espanse, così la validazione globale dei contributi attivi rimane invariata. Una dichiarazione mista del nome completo con un suo addendo viene respinta come duplicazione. Aggunti test di equivalenza e rifiuto in `TEST_SOURCE_REGISTRY` per la legge lava. Non cambiano le formule SI/IM, la scelta del gruppo, gli schemi RK o la conservazione. **Stato:** D17.16/input nominale in avanzamento; catalogo ancora a 34 voci e `source_term_type` transitorio. La fusione interna delle voci reologiche e il catalogo minimale D17.2/10c non sono conclusi. Compilazione e regressioni sul nuovo HEAD da verificare.


### D17.22.34 — Atomicità della legge reologica in Newton e analytic (24 ottobre 2026)

Il registro rifiuta ora ogni gruppo implicito contenente soltanto un sottoinsieme degli addendi legacy appartenenti alla stessa legge reologica. La verifica di completezza usa `rheology_source_terms(family)` e impedisce di considerare schedulabile una singola formula fisica di Voellmy, lava, lahar ecc. Un test negativo aggiuntivo copre esplicitamente la legge lava incompleta in Newton; il controllo analitico di completezza preesistente rimane attivo. La misura prepara l'unificazione delle voci senza cambiare il calcolo SI/IM né il comportamento legacy con leggi complete.

**Stato:** passo intermedio, non ancora catalogo nominale definitivo. `source_term_type`, `term_ids` e il catalogo dei 34 addendi sono ancora presenti e devono essere sostituiti con una voce per legge in una migrazione coordinata di registry, input e dispatcher. Build e test del nuovo HEAD non ancora eseguiti.


### D17.22.35 — Catalogo separato delle identità reologiche canoniche (24 ottobre 2026)

Introdotto in `source_catalog_2d` l'array `RHEOLOGY_LAWS` con una sola voce nominale per ciascuna delle dieci leggi reologiche complete, distinta dall'array legacy `SOURCE_CATALOG` dei contributi matematici. La funzione `rheology_source_name` ora risolve il nome attraverso il catalogo canonico; `canonical_rheology_catalog_unique` ne valida nomi e model ID univoci, con test dedicato. Il parser `SOURCE_GROUPS` continua a espandere il nome della legge nelle componenti fisiche esistenti. Questa separazione costituisce la base per l'eliminazione successiva degli ID delle componenti dai gruppi e dal dispatcher; **non** rappresenta ancora la rimozione del catalogo legacy né il completamento del `source_term_type` minimale. Le Actions devono ancora verificare il nuovo commit.


### D17.22.36 — Vista canonica delle sorgenti schedulabili per gruppo (24 ottobre 2026)

`source_registry_2d` esporta `group_source_names(group)`, che espone ogni legge reologica come una sola sorgente nominale e deduplica gli addendi matematici appartenenti alla stessa famiglia, mantenendo distinti gli altri termini fisici. La rappresentazione interna `term_ids` non viene ancora eliminata: rimane necessaria agli assemblatori REAL/COMPLEX e alla compatibilità legacy. Inserita una regressione nella matrice delle configurazioni legacy che verifica l'identità canonica completa per i gruppi reologici. **Questo passo non equivale ancora a un registro esclusivamente nominale né alla rimozione dei 34 contributi del catalogo transitorio.** La compilazione e le Actions del nuovo HEAD sono ancora da verificare.


### D17.22.37 — Identità canoniche persistenti nel gruppo (24 ottobre 2026)

`source_group_type` contiene ora `source_names(:)`, popolato una sola volta durante `initialize` insieme alla rappresentazione transitoria degli addendi (`term_ids` e `term_names`). `source_names` espone una sola identità per legge reologica completa e una voce per ogni altra sorgente; la validazione verifica che il campo coincida con il raggruppamento canonico calcolato dal catalogo. La matrice di regressione legacy controlla la disponibilità e coerenza di questo campo per ogni gruppo di stage. **Non sono stati eliminati `term_ids` o `term_names`**, tuttora richiesti dai kernel; la migrazione del dispatcher e la rimozione degli addendi dal catalogo schedulabile restano in corso. Nuovo HEAD ancora da verificare con le Actions.


### D17.22.38 — Selezione nominale effettiva nel dispatcher SI/IM (24 ottobre 2026)

In `equation_terms_2d`, `group_has_role` verifica l'appartenenza della legge reologica attraverso `group%source_names` prima del controllo interno della presenza del ruolo e degli addendi. In `source_registry_2d`, `group_contains_name` risolve prioritariamente i nomi canonici `source_names`, mantenendo l'accesso transitorio agli addendi `term_names` per la compatibilità con il dispatcher corrente. Estesa la matrice dei test legacy alla risoluzione nominale delle leggi complete. Non sono mutate le formulazioni fisiche SI/IM, il Newton o gli schemi RK. **La migrazione non è completata**: `term_ids` e il catalogo legacy restano necessari alla valutazione fisica. Il nuovo HEAD richiede verifica delle Actions e delle regressioni numeriche.


### D17.22.39 — Eliminazione della dipendenza da term_ids nel dispatcher fisico (24 ottobre 2026)

`src/physics/equation_terms_2d.f90` non legge più `group%term_ids`: i controlli d'ingresso per i gruppi espliciti e Newton richiedono `source_names` allocato e la ricerca di ruoli SI/IM risolve gli addendi dai nomi `term_names`, mantenendo il filtro sulla legge canonica `source_names`. Un controllo nel test del registro verifica la disponibilità delle due rappresentazioni nominali e la selezione nominale della legge lava. **Limite:** `term_ids` rimane nell'oggetto gruppo e in registry/configurazione/test, perciò non è stata completata la sua eliminazione globale. Questa modifica non cambia le formule fisiche; la verifica del nuovo HEAD attraverso le Actions è ancora richiesta.


### D17.22.40 — Rimozione effettiva di `term_ids` da `source_group_type` (24 ottobre 2026)

Il campo allocabile `term_ids(:)` viene **eliminato** da `source_group_type`. I gruppi di stage/step conservano `term_names(:)` per la compatibilità temporanea degli addendi fisici e `source_names(:)` per le identità delle sorgenti schedulabili; `equation_ids(:)` continua a definire il blocco numerico. `initialize_legacy` converte gli ID d'ingresso in nomi durante la costruzione di ciascun gruppo, mentre `source_configuration_2d` converte allo stesso modo le dichiarazioni native. `validate_registry` risolve ogni nome mediante `source_term_descriptor_by_name`, verifica unicità, schedulabilità, target, duplicazioni globali e leggi complete; le verifiche dei bundle e dei gruppi analitici ricavano gli ID soltanto localmente dal catalogo, senza memorizzarli nel gruppo. `equation_terms_2d` resta nominale come dal passo D17.22.39. La suite `TEST_SOURCE_REGISTRY` è stata migrata ai confronti nominali e conserva i casi negativi di sorgenti sconosciute, derivate, duplicate e reologie incomplete. **Stato:** rimozione del campo implementata, in attesa dei gate di compilazione e test sul nuovo SHA. La successiva eliminazione di `term_names` e del catalogo numerico transitorio resta distinta; non dichiarare ancora conclusa D17.2/10c né l'equivalenza numerica.


### D17.22.41 — Fix importazione del dispatcher dopo la rimozione degli ID (24 ottobre 2026)

Il gate `D17 compile diagnostics` sul commit `91b16562` ha evidenziato un errore di compilazione in `equation_terms_2d.f90`: `rheology_source_name` veniva chiamata da `group_has_role` senza essere importata nella clausola `USE source_catalog_2d, ONLY`. Aggiunta l'importazione esplicita. La procedura del catalogo era già `PURE`: il secondo errore del compilatore era una conseguenza della risoluzione fallita del simbolo. Nessun cambiamento di fisica o scheduling. La compilazione e i test devono essere rieseguiti sul nuovo commit.

### D17.22.42 — Chiusura CI della rimozione di `term_ids` dai gruppi (24 ottobre 2026)

Dopo il passo D17.22.41, la compilazione ha superato l'errore di importazione del simbolo `rheology_source_name`. Il successivo errore di build nel modulo `mass_exchange_step_2d.f90` è stato risolto nel commit [`609f293`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/609f2931d937cf78ae08a9788bc4e8373e6f17a6): il dispatcher post-RK usa `group%term_names` e `source_term_descriptor_by_name`, senza accedere a `%term_ids`.

I test del registro hanno poi richiesto interventi mirati: [`07bd36d`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/07bd36d0650facb97091c68189064d27d08592f0) evita il costruttore dinamico di stringhe di lunghezza non uniforme in `group_source_names`; [`acac97e`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/acac97ee22368e887d844c1417cc65ea9c24bc44) aggiunge diagnostica stderr per le combinazioni legacy rifiutate; [`c02ba1f`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/c02ba1fbcd12082188dae01008e06a21fbcd95c2) azzera il contesto `family` per ogni gruppo durante la validazione; [`f8420e5`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/f8420e5e5f22f03b8bf04b3e9d325ef3ae2491c1) adegua il test negativo dei bundle reologici al rifiuto più specifico `complete law bundle`.

**Esito comunicato dall'utente:** entrambe le GitHub Actions (compilazione e test del registro) sono passate dopo questi interventi, sullo stato del branch culminato nel commit `f8420e5`. Non sono stati forniti in questa conversazione gli identificatori dei due run: la conferma è riportata come evidenza comunicata dall'utente, non come audit indipendente dei log CI. I precedenti errori della sequenza sono quindi risolti rispetto a questi gate; ciò non certifica l'equivalenza numerica complessiva D17.20–D17.21.

### D17.22.43 — Audit del campo eliminato e inizializzazione esplicita di `make_group` (24 ottobre 2026)

Sul commit `f8420e5` sono stati esaminati tutti i 50 file Fortran di `src/`: non risultano più riferimenti al campo `term_ids` né accessi `%term_ids` nei sorgenti di produzione. **La rimozione di `term_ids` da `source_group_type` e dai relativi consumatori è completata a livello di sorgente e sostenuta dai gate di build/test comunicati dall'utente.**

Nel commit [`d3cfc58`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/d3cfc58fc1d01164d3669a19fd3db5f0889c6542), `make_group` in `source_registry_2d.f90` alloca esplicitamente `term_names` ed `equation_ids` prima delle assegnazioni, per affrontare i warning GNU `-Wmaybe-uninitialized` sugli allocabili. **La scomparsa dei warning e i nuovi gate CI dopo questo commit devono ancora essere verificati.**

**Limite architetturale vincolante:** non confondere il risultato precedente con il completamento di D17.2 e D17.19/10c. Il catalogo corrente contiene ancora i 34 descrittori numerici dei contributi, `source_term_descriptor(id)`, la conversione `source_names_for_ids`, `group_term_ids` e la rappresentazione compatibile `term_names`; `source_names` rappresenta già le identità canoniche delle leggi complete, ma non è ancora l'unica struttura delle sorgenti. La specifica finale rimane il **catalogo minimale nominale, una sorgente per legge fisica**, senza ID numerici delle sorgenti né addendi schedulabili individualmente. Servono modifiche coordinate a catalogo, configurazione, registro, dispatcher e test, con verifica delle regressioni numeriche D17.20–D17.21.

**Stato:** rimozione del campo `term_ids` dai gruppi **COMPLETATA**; ulteriore migrazione verso il catalogo minimale **IN CORSO**; accettazione D17 complessiva **NON COMPLETATA**.

### D17.22.44 — Documento normativo unico per i successivi aggiornamenti (24 ottobre 2026)

Per richiesta esplicita dell'utente, **questo file** `docs/D17_Source_Groups_Revised_Specification_2026-10-24.md` è il solo documento D17 da mantenere aggiornato. Le precedenti menzioni di altri handoff, documenti storici o roadmap in questa specifica descrivono il contesto della migrazione e non costituiscono un obbligo di aggiornarli. **Ogni successiva modifica di stato o registrazione di un'implementazione D17 va riportata qui**, preservando la specifica integrale e distinguendo tra codice modificato, compilazione riuscita, test del registro e regressioni numeriche scientifiche.

### D17.22.45 — Confronto nominale dei bundle reologici nel validatore (24 ottobre 2026)

Nel modulo `source_registry_2d` il controllo di completezza della legge reologica nei gruppi Newton e analytic confronta ora direttamente i nomi degli addendi presenti in `group%term_names` con quelli attesi dal catalogo. Eliminata la funzione inversa `group_term_ids`, che ricostruiva temporaneamente una lista di ID numerici dalla rappresentazione nominale del gruppo; il confronto usa `name_membership` e la conversione di sola andata `source_names_for_ids`. I controlli di completezza e cardinalità del gruppo analytic restano invariati nella loro intenzione.

**Stato:** modifica al codice pubblicata, **compilazione e test da verificare sul nuovo commit**. Questo intervento elimina una conversione nome→ID dal registro ma **non** conclude la migrazione D17.2/D17.19-10c: il catalogo numerico dei 34 contributi, `source_names_for_ids`, e `term_names` rimangono transitori e richiedono un intervento successivo coordinato. Nessuna modifica prevista ai kernel fisici, ai tableau o al solutore Newton.

### D17.22.46 — Regola di esecuzione dei gate CI e successiva semplificazione del costruttore (24 ottobre 2026)

**Regola operativa vincolante concordata con l'utente:** **non eseguire o richiedere GitHub Actions dopo ogni modifica o commit.** Gli interventi incrementali appartenenti a uno stesso blocco architetturale D17 possono essere accumulati in più commit coerenti, accompagnati da revisione statica del codice e delle interfacce. Eseguire le Actions **solo ai gate sostanziali**, in particolare dopo il completamento di una migrazione strutturale, la modifica di un contratto numerico o quando occorre verificare esplicitamente una proprietà di comportamento, compilazione o regressione. Un risultato verde vale per lo SHA verificato, non per commit successivi. Le regressioni di equivalenza fisica/numerica restano gate distinti dai test di build e registry. Questa regola disciplina anche tutti i passi successivi di D17.19–D17.21.

**Gate superato:** l'utente conferma che **entrambe le Actions di compilazione e test sono passate** dopo il commit [`2f52b28`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/2f52b2866cdd29ccc7dfc26713b6643c9ff745c2). La conferma non equivale al completamento delle regressioni scientifiche D17.20 né è estensibile automaticamente ai commit successivi.

**Passo incrementale successivo:** in `source_registry_2d`, `make_group` riceve ora i nomi degli addendi fisici `names(:)` anziché ID numerici `terms(:)`. `initialize_legacy_registry` converte gli ID legacy con `source_names_for_ids` **al confine di costruzione dei gruppi**, mentre il gruppo e il suo costruttore operano direttamente con nomi. La rappresentazione temporanea `term_names`, il catalogo degli addendi e la conversione legacy in ingresso rimangono; non è stato ancora raggiunto il catalogo minimale finale D17.2 / D17.19-10c. Nessuna modifica intenzionale ai kernel o all'ordine RK.

**Verifica:** questa modifica incrementale è **PENDING AL PROSSIMO GATE SOSTANZIALE**; non richiede un'Action dedicata immediata in base alla regola sopra.

### D17.22.47 — Aggiornamento sistematico dello stato D17.x e selezione nominale dei target (24 ottobre 2026)

**Regola documentale vincolante:** l'intero testo delle sezioni D17.0–D17.21 resta normativo e non deve essere accorciato o sostituito con una sintesi. Dopo ogni commit, aggiornare gli **stati locali dopo ciascun punto D17.x coinvolto** e mantenere sincronizzata la **tabella operativa D17.19**, oltre al registro cronologico D17.22. Distinguere rigorosamente implementato, verificato da build/registry e verificato scientificamente. I resoconti storici non equivalgono allo stato attuale; non rimuovere dettagli o requisiti non ancora implementati.

**Stato CI comunicato:** compilazione e test del registro superati su `2f52b28`; i commit successivi `3d2f326` e questo passo incrementale restano in attesa del prossimo gate sostanziale, senza richiedere Actions per ciascun commit.

**Implementazione:** in `source_registry_2d`, `names_for_targets` sostituisce `terms_for_targets` come selettore dei termini legacy per target, restituendo direttamente i nomi degli addendi nell'ordine storico. `initialize_legacy_registry` passa questi nomi al costruttore nominale `make_group`; scompare il passaggio intermedio che costruiva la lista selezionata di ID numerici per poi riconvertirla in nomi. I metadati del catalogo legacy restano ancora necessari per interpretare i vecchi input; non viene modificata la fisica e non sono state rimosse le singole voci di contributo reologico.

**Gate futuro:** dopo ulteriori commit coordinati, compilazione e strict/optimized registry; prima della chiusura D17, matrice numerica D17.20, controlli D17.21 e prova di estensibilità D17.18a.

### D17.22.48 — Normalizzazione unica dei nomi attivi prima del raggruppamento (24 ottobre 2026)

Commit [`308979d`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/308979d85b4220305574bd15d258c58c1caf174e): `initialize_legacy_registry` conserva gli ID necessari all'interfaccia di input legacy e alla validazione iniziale, ma ne ricava una sola volta `active_names=source_names_for_ids(active_terms)`. La selezione successiva `names_for_targets` riceve i nomi e consulta i metadati con `source_term_descriptor_by_name`. I gruppi continuano a contenere soltanto i nomi degli addendi e delle sorgenti schedulabili; non viene ripristinato `term_ids`. Ordine dei termini, cadenze, gruppi, schemi RK e kernel fisici non sono stati intenzionalmente modificati.

**Stato:** D17.2, D17.5, D17.16 e D17.19 passi 3, 4, 10c **IN CORSO**. Il catalogo finale minimale con una sorgente per legge fisica resta da realizzare. **Verifica:** commit pubblicato, ma non ancora coperto dalle Actions; queste verranno eseguite solo al prossimo gate sostanziale secondo D17.22.46. Questo aggiornamento non elimina, riassume né sostituisce i requisiti tecnici originali.

### D17.22.49 — Validazione nominale anche al confine legacy (24 ottobre 2026)

`initialize_legacy_registry` esegue `active_names=source_names_for_ids(active_terms)` **prima** dei controlli di schedulabilità e compatibilità dei target; questi controlli utilizzano ora `source_term_descriptor_by_name(active_names(i))`. Rimangono immutati i messaggi di rifiuto, le verifiche della partizione implicita, la selezione dei gruppi e il successivo instradamento per nomi. È stato rimosso il precedente secondo punto di conversione, ora ridondante. L'interfaccia legacy accetta ancora ID numerici: l'intervento non introduce il catalogo finale minimale e non modifica kernel o schemi RK.

**Stati delle sezioni pertinenti:** D17.2, D17.5 e D17.16 **IN CORSO**; D17.19 passi 3, 4 e 10c **IN CORSO**, con righe della tabella aggiornate. La documentazione normativa esistente rimane integralmente conservata. **Verifica:** modifica incrementale, in attesa di un gate CI sostanziale secondo D17.22.46; nessuna equivalenza numerica complessiva dichiarata.

### D17.22.50 — Atomicità dei bundle verificata per nomi, senza mappa ID→gruppo (24 ottobre 2026)

In `source_registry_2d::validate_registry` viene eliminato l'array transitorio `term_group(SOURCE_TERM_COUNT)`, precedentemente usato per registrare la posizione numerica di ogni addendo nel gruppo. La validazione di un batch atomico utilizza ora `same_batch(names,groups)`, che verifica direttamente `group%term_names` e accetta soltanto due configurazioni: batch totalmente assente, oppure tutti i suoi membri presenti nello stesso gruppo. La protezione preesistente `seen(id)` contro le duplicazioni globali resta invariata; sono mantenuti sia il messaggio di errore `D17: coupled evaluation batch split or incomplete` sia la lettura temporanea delle membership dal catalogo legacy.

**Stato locale:** D17.2 **IN CORSO**, D17.5 **IN CORSO**, D17.16 **IN CORSO**; tabella D17.19 aggiornata ai passi 3, 5 e 10c, tutti **IN CORSO**. Questo è un passo intermedio: il catalogo numerico e le voci separate degli addendi reologici non sono ancora stati eliminati. Il testo normativo della specifica precedente resta integralmente invariato. **Verifica:** modifica pubblicata, in attesa del prossimo gate CI sostanziale; nessuna conferma di equivalenza fisica/numerica scientifica.

### D17.22.51 — Parser nativo nominale: eliminazione della matrice di ID dei gruppi (24 ottobre 2026)

Commit [`2ee7d5e`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/2ee7d5e71da809f0061f0f52fe9e6a2e87e9b471). `initialize_native_source_configuration` non conserva più i membri dei gruppi nella matrice numerica `parsed_terms`: utilizza invece `parsed_names` con nomi canonici. La costruzione di `group%term_names` copia direttamente tali nomi, senza ricostruire i descrittori dagli ID alla fine del parsing. Eliminato anche il lookup ausiliario `native_term_id`; il riconoscimento dei termini di input utilizza `source_term_descriptor_by_name`. La dichiarazione di una legge reologica completa continua temporaneamente a espandersi negli addendi del vecchio catalogo, e la deduplicazione globale conserva un indicizzamento legacy: questo passo **non** introduce ancora il catalogo minimale definitivo.

**Stato dei punti interessati:** D17.2 **IN CORSO**, D17.5 **IN CORSO**, D17.16 **IN CORSO**; tabella D17.19 passi 3, 4 e 10c aggiornata, tutti **IN CORSO**. **Verifica:** commit pubblicato, nuovo gate di build/registry ancora da eseguire; regressioni numeriche D17.20–D17.21 non certificate. È stata mantenuta integralmente la specifica normativa precedente, senza eliminazioni o sintesi.

### D17.22.52 — Unicità e copertura del parsing nativo per nomi (24 ottobre 2026)

Commit [`41e20a4`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/41e20a495c8fd8320c951ba4aecbee0145f55fed). Nel parser nativo `source_configuration_2d`, la mappa temporanea `seen_terms(SOURCE_TERM_COUNT)` degli ID attivati è stata sostituita da `seen_names` con confronto diretto dei nomi canonici. La validazione rifiuta gli stessi casi di ripetizione globale e mantiene il controllo di copertura: ogni contributo attivo della configurazione legacy deve essere rappresentato una volta nel calendario nativo. Non è cambiato il significato dei gruppi, né la classificazione delle sorgenti; il catalogo legacy numerico è ancora impiegato per costruire i nomi attivi e per espandere le leggi reologiche complete.

**Stato locale:** D17.2 **IN CORSO**, D17.5 **IN CORSO**, D17.16 **IN CORSO**; D17.19 passi 3, 4 e 10c aggiornati e **IN CORSO**. **Verifica:** commit incrementale pubblicato, gate CI da eseguire al termine del blocco strutturale; D17.20–D17.21 ancora aperti. Nessun requisito o dettaglio tecnico preesistente della specifica è stato eliminato.

### D17.22.53 — Elenco delle sorgenti native attive per nome (24 ottobre 2026)

In `source_configuration_2d`, l'insieme `active_terms` che descriveva per ID le sorgenti fisiche attive della configurazione nativa è sostituito da `active_names`. Le attivazioni esplicite e i controlli finali di copertura utilizzano direttamente i nomi catalogati, mentre i nomi delle componenti delle leggi reologiche sono ancora ricavati dagli ID del catalogo transitorio durante la costruzione dell'elenco. `parsed_names` e `seen_names` restano gli unici elenchi di membership del percorso nativo. Le condizioni di attivazione e i messaggi di errore rimangono invariati.

**Stati:** D17.2 e D17.16 **IN CORSO**; tabella D17.19 passi 4 e 10c aggiornata e **IN CORSO**. Il passo non rimuove ancora le costanti `SOURCE_TERM_*` dal resto del repository, i sottotermini del catalogo né i metadati legacy. **Verifica:** attesa al prossimo gate sostanziale, senza affermare equivalenza numerica; tutti i requisiti originali sono preservati.

### D17.22.54 — Espansione nominale delle leggi nel parser e disaccoppiamento dalla cardinalità legacy (24 ottobre 2026)

Commit [`249f07a`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/249f07a7ec6a47a984ccba3a70c5c7654a0f3ec4). In `initialize_native_source_configuration` le dichiarazioni delle reologie complete vengono espanse temporaneamente dai componenti del catalogo in `law_names`, prima della verifica di unicità e della memorizzazione in `parsed_names`. Per gli altri nomi nativi, il descrittore catalogato fornisce direttamente il nome; i controlli sul gruppo non usano più un `term_id` numerico locale. La capacità di `seen_names` è legata al limite dei termini dichiarabili nei gruppi (`SOURCE_INPUT_MAX_GROUPS*SOURCE_INPUT_MAX_TERMS`), non più a `SOURCE_TERM_COUNT`. I vecchi addendi reologici rimangono una compatibilità temporanea, da eliminare con la migrazione simultanea del catalogo e del dispatcher.

**Stato dopo questo intervento:** D17.2, D17.5, D17.16 e D17.19 passi 3, 4, 10c **IN CORSO**. Questo non modifica deliberatamente i kernel fisici, non introduce un nuovo catalogo permanente e non sostituisce ancora il tipo minimale a cinque campi. **Verifica:** pubblicato ma non sottoposto a nuove Actions; il prossimo gate sostanziale resta necessario prima di attribuire successo di compilazione e test al nuovo HEAD. I requisiti tecnici precedenti sono conservati integralmente.

### D17.22.55 — API nominale di adattamento reologico per il parser nativo (24 ottobre 2026)

Nel modulo `source_catalog_2d` è introdotta la funzione temporanea `rheology_component_names(model_id)`, che restituisce **direttamente i nomi** degli addendi storici di una legge reologica in ordine di catalogo, senza costruire un array di ID fisici. Nel modulo `source_configuration_2d` il percorso `SOURCE_GROUPS` elimina la variabile `law_terms`: sia l'espansione della legge nativa sia l'elenco delle sorgenti attive usano `law_names` e nomi canonici. La selezione del modello legacy resta numerica all'ingresso e non si modificano le valutazioni fisiche.

**Limite dichiarato:** questa API è esclusivamente un adattatore transitorio, non un secondo catalogo permanente. Gli addendi storici e le 34 voci catalogate non sono ancora sostituiti da un'unica sorgente per legge, e il `source_term_type` a cinque campi resta da implementare in coordinamento con il dispatcher. **Stati:** D17.2, D17.16 e tabella D17.19 passi 4, 9, 10c **IN CORSO**. **Verifica:** nuovo commit da confermare al prossimo gate CI sostanziale; preservata integralmente la specifica e le evidenze precedenti.

### D17.22.56 — Completezza reologica verificata direttamente per nomi nel registro (24 ottobre 2026)

Il validatore `source_registry_2d` utilizza ora `rheology_component_names(family)` per ottenere direttamente la membership storica completa delle leggi reologiche nei gruppi NEWTON e ANALYTIC. Nei due controlli sono state eliminate le liste intermedie `terms=rheology_source_terms(family)` e la successiva conversione `source_names_for_ids(terms)`. Sono preservati confronto di membership, cardinalità degli analytic velocity e diagnostiche. Restano temporanei `family`, i metadati del catalogo dei 34 contributi e l'espansione degli addendi: la vera identità finale rimane `rheology_<model>` e dovrà essere l'unica voce catalogata.

**Stati aggiornati:** D17.2, D17.5, D17.16 **IN CORSO**; tabella D17.19 passi 3, 9 e 10c **IN CORSO**. **Verifica:** modifica ancora da validare con Actions nel prossimo gate sostanziale; regressioni D17.20–D17.21 aperte. Nessun requisito della specifica preesistente è stato rimosso o sintetizzato.

### D17.22.57 — Unicità globale delle sorgenti del registro per nome (24 ottobre 2026)

In `source_registry_2d::validate_registry` la mappa booleana `seen(SOURCE_TERM_COUNT)`, indicizzata dalla posizione numerica del vecchio catalogo, è stata sostituita da `seen_names`: una lista nominale dei contributi già incontrati nei gruppi. L'unicità globale viene verificata con `ANY(seen_names==term%name)`, preservando il rifiuto dei duplicati e i relativi messaggi. La lista viene costruita soltanto durante la validazione iniziale, fuori dai cicli temporali e di cella. Il controllo temporaneo `term%id` per le voci legacy rimane intenzionalmente fino alla sostituzione coordinata del catalogo.

**Stati dopo il commit:** D17.2, D17.5, D17.16 **IN CORSO**; tabella D17.19 passi 3, 4, 10c aggiornati e **IN CORSO**. `source_term_type` definitivo a cinque campi e la rimozione delle costanti `SOURCE_TERM_*` sono ancora da realizzare. **Verifica:** modifica pubblicata, non ancora verificata con compilazione o Actions; si applica la regola dei gate sostanziali D17.22.46. Tutti i requisiti e dettagli normativi antecedenti sono mantenuti.

### D17.22.58 — Validazione degli sconosciuti senza identità numeriche nel registro (24 ottobre 2026)

In `source_registry_2d::validate_registry` il controllo di validità dei membri dei gruppi non usa più `term%id` e `SOURCE_TERM_COUNT`. Il descrittore richiesto tramite `source_term_descriptor_by_name` è valido solo se `LEN_TRIM(term%name)>0`, seguito dal controllo preesistente `term%schedulable`; gli altri controlli di capability, target, unicità nominale e legge reologica rimangono invariati. L'eliminazione dei campi legacy dal tipo richiederà comunque migrazione coordinata del catalogo e dei consumatori rimasti. Gli indici di equazione e l'indice locale `id` usato per equazioni e scansioni del catalogo sono concetti distinti dalle identità fisiche.

**Stato:** D17.2, D17.5 e D17.16 **IN CORSO**, tabella D17.19 passi 3, 9 e 10c aggiornata. Le costanti `SOURCE_TERM_*` e il catalogo con 34 voci restano da rimuovere. **Gate:** modifica non ancora testata dalle Actions, da verificare al prossimo gate sostanziale insieme agli altri commit successivi a `2f52b28`. I requisiti precedenti della specifica rimangono integralmente preservati.

### D17.22.59 — Verifica nominale dei batch atomici senza liste ID intermedie (24 ottobre 2026)

Nel modulo `source_catalog_2d` viene introdotto l'adattatore transitorio `source_bundle_names(bundle_id)`, che restituisce direttamente, nell'ordine storico, i nomi degli addendi di un batch atomico legacy. In `source_registry_2d::validate_registry` sono state eliminate le liste numeriche temporanee `bundle_terms`: sia l'identificazione del primo componente del batch sia la verifica `same_batch` usano nomi. La semantica rimane «batch completamente assente oppure integralmente presente in un unico gruppo» e non modifica alcun kernel fisico o solutore. Il catalogo dei 34 addendi, il metadato `bundle_id` e la scansione di compatibilità degli ID del catalogo non sono ancora rimossi; l'adattatore è una fase intermedia da eliminare nella migrazione finale D17.2.

**Stati:** D17.2, D17.5, D17.16 **IN CORSO**, tabella D17.19 passi 3, 9, 10c aggiornata e **IN CORSO**. **Verifica:** commit pubblicato, build/registry da eseguire al prossimo gate sostanziale; D17.20–D17.21 restano aperti. Nessun requisito precedente eliminato o sintetizzato.

### D17.22.60 — Verifica dell'atomicità guidata dai nomi effettivi dei gruppi (24 ottobre 2026)

**Blocco architetturale (un unico commit):** `source_registry_2d::validate_registry` non esegue più `DO id=1,SOURCE_TERM_COUNT` per trovare i batch fisici atomici. Dopo aver validato le identità dei membri, itera sulle **sorgenti effettivamente presenti** in ciascun `group%term_names`, ottiene il descrittore tramite lookup nominale e, soltanto per un membro con `bundle_id` non nullo, confronta il batch completo con `source_bundle_names` e `same_batch`. Si preserva il rifiuto di batch incompleti o distribuiti tra gruppi, anche quando il primo componente del batch è assente. Nessuna formula reologica o altro kernel numerico è modificata.

**Limite e gate:** la scansione numerica del catalogo è eliminata dal validatore dei batch, ma `bundle_id`, `SOURCE_TERM_*`, `rheology_model_id` e la registrazione dei 34 addendi rimangono transitori. La validazione preventiva conserva i controlli sui nomi e sulle capabilities. Se un batch dichiarato nel catalogo è malformato ma del tutto assente dai gruppi, non viene più ispezionato da questo controllo del registro: l'invariante delle definizioni del catalogo deve essere verificato autonomamente al momento della migrazione finale, senza caricare il registro di responsabilità estranee.

**Stati:** D17.2, D17.5, D17.16 **IN CORSO**; D17.19 passi 3, 4, 9, 10c **IN CORSO** e aggiornati direttamente nella tabella. **Verifica:** CI e regressioni da eseguire a un gate sostanziale; non dichiarare completata la migrazione a cinque campi. Specifica precedente interamente conservata.

### D17.22.61 — Migrazione coordinata dell'interfaccia legacy configurazione–registro (24 ottobre 2026)

In un unico blocco sono aggiornati sia `source_configuration_2d::initialize_configured_sources` sia `source_registry_2d::initialize_legacy_registry`. La configurazione continua a rilevare le sorgenti attive con gli stessi flag e nello stesso ordine, ma compone `active_names` con i nomi degli addendi storici; per la reologia utilizza `rheology_component_names`. Il registro accetta ora `active_names(:)` direttamente e controlla schedulabilità, target, completezza e raggruppamento tramite lookup nominale. Eliminata da `source_registry_2d` la conversione `source_names_for_ids`, non più necessaria a questo confine. Non si introducono ID sostitutivi né nuovi campi nei gruppi.

**Resta da realizzare:** l'unificazione delle voci reologiche, il catalogo fisico definitivo con cinque campi, il target nominale e la rimozione generale delle costanti numeriche `SOURCE_TERM_*`. Il parser nativo e l'adattamento legacy conservano temporaneamente i nomi degli addendi per non alterare il dispatcher fisico ancora attuale. **Stati:** D17.2, D17.5, D17.13, D17.16 e D17.19 passi 3, 4, 10c **IN CORSO**. **Verifica:** gate di compilazione/registro ancora necessario prima dell'accorpamento fisico; nessuna regressione numerica attestata. Tutto il testo normativo precedente resta preservato.

### D17.22.62 — Correzione gate di compilazione GNU dopo migrazione nominale legacy (24 ottobre 2026)

**Evidenza del gate sul commit `50a81f57ce7b27ae8045c22fa959cf7590e760e8`:** i log forniti dall'utente mostrano successo di `autoreconf` e `./configure` e fallimento della compilazione con GNU Fortran 13.3.0 in `src/physics/source_configuration_2d.f90`. Alle righe originali 69, 72 e 74 le espressioni che aggiungono rispettivamente `SOURCE_TERM_MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS`, `SOURCE_TERM_MASS_EXCHANGE_CARRIER_DARCY_LOSS` e `SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION` superavano 132 colonne, generando `-Werror=line-truncation` e `Unterminated character constant`. La compilazione non è arrivata al linking: nessuna conclusione sui test runtime né sulla validità numerica è consentita.

**Correzione:** suddivise le tre istruzioni usando il carattere di continuazione Fortran `&` preservando esattamente i nomi, le condizioni di attivazione e l'ordine di inserimento in `active_names`. Effettuato controllo statico: nessuna riga di `source_configuration_2d.f90` supera ora il limite di 132 caratteri. Nessun'altra modifica al comportamento o alle API.

**Stati aggiornati:** D17.2, D17.5, D17.13, D17.16 **IN CORSO**, tabella D17.19 passi 3, 4, 10c e 12 **IN CORSO**. La correzione richiede una nuova build sul commit che la contiene, senza anticipare l'esito delle Actions. Il catalogo finale a cinque campi resta da realizzare. La specifica tecnica preesistente è conservata integralmente.

### D17.22.63 — Riparazione test del registro dopo migrazione dell'API legacy nominale (24 ottobre 2026)

**Evidenze:** la build principale sul commit `528ce0f56e91c80748f49e748433329bbde96ee0` è stata confermata riuscita dall'utente. Il log `test-run(6).log` mostra sei errori di compilazione in `TESTS/TEST_SOURCE_REGISTRY/test_source_registry.f90`, alle righe 119, 306, 319, 341, 344, 346: `initialize_legacy` richiede ora `active_names(:)` di tipo `CHARACTER(*)`, mentre i test passavano `INTEGER(4)` (in tre casi vettore di ID, nei restanti tre costruttore vuoto `[INTEGER::]`). I test quindi non sono stati eseguiti: non si tratta di un fallimento di asserzioni.

**Correzione coordinata:** tutti i punti di chiamata del test utilizzano l'interfaccia nominale. Gli scenari che costruiscono `terms(:)` come ID di catalogo per copertura combinatoria li convertono **soltanto nel test**, tramite `catalog_names_for_test(terms)` prima di chiamare il registro; tre casi negativi con input vuoto passano `[CHARACTER(LEN=64)::]`. La conversione non viene reintrodotta nel codice del solver. Preservati il numero e i contenuti delle 864 combinazioni legacy, i test dei batch incompleti e i casi di errore.

**Stati:** D17.5, D17.13, D17.19 passi 3, 4, 10c, 12 **IN CORSO**. La build principale precedente è verificata, ma il test aggiornato e le sue asserzioni attendono un nuovo gate CI sul commit correttivo. La specifica D17 originale, incluso l'obiettivo del catalogo a cinque campi, resta invariata e integrale.

### D17.22.64 — Riparazione lunghezze CHARACTER nella verifica nominale del registro (24 ottobre 2026)

**Evidenza:** il log `test-run(7).log` della CI mostra che il test del registro supera la fase di compilazione ma si interrompe a runtime in `source_registry_2d.f90:353`, chiamato da `test_source_registry.f90:119`: `Different CHARACTER lengths (30/64) in array constructor`. La lista allocabile `seen_names` ha elementi `CHARACTER(LEN=64)`, ma il compilatore tratta nell'array constructor `[seen_names,term%name]` il lato sinistro con una diversa lunghezza effettiva durante la verifica runtime. La causa non è un errore di soluzione numerica: l'esecuzione non arriva alle asserzioni di regressione successive.

**Correzione:** sostituita l'unica costruzione `seen_names=[seen_names,term%name]` nel registro con `seen_names=[CHARACTER(LEN=64) :: seen_names,term%name]`, mantenendo invariati membership e diagnostiche. Il tipo esplicito forza una lunghezza coerente dei componenti della lista. Non sono introdotti identificatori numerici né modificati kernel fisici. Il gate successivo deve accertare anche che non ci siano altri problemi runtime nascosti.

**Stati:** D17.2, D17.5, D17.13 e D17.19 passi 3, 4, 10c, 12 **IN CORSO**. Build e compilazione del test confermate dal gate precedente; risultato runtime ancora pendente. Requisiti originali D17 preservati.

### D17.22.65 — Correzione uniforme dei costruttori di nomi reologici nella configurazione (24 ottobre 2026)

**Evidenza:** `test-run(8).log` del gate sul commit `faa69c2` segnala un errore runtime `Different CHARACTER lengths (0/64) in array constructor` in `src/physics/source_configuration_2d.f90:49`, chiamato durante `test_configured_default_schedule`. La concatenazione `active_names=[active_names,law_names]` non dichiara la lunghezza dei caratteri nel costruttore, pur essendo `active_names` allocato inizialmente vuoto e `law_names` prodotto dall'adattatore nominale. Il controllo statico del sorgente rivela una seconda istruzione identica nel percorso nativo, che potrebbe determinare lo stesso problema.

**Correzione coordinata:** entrambe le occorrenze, nei percorsi legacy e nativo, diventano `active_names=[CHARACTER(LEN=64) :: active_names,law_names]`. Tutte le altre concatenazioni nominali del modulo erano già dotate di type-spec esplicita. Nessuna condizione di attivazione, identità fisica o formula è modificata. Il fix precedente per `seen_names` nel registro rimane valido.

**Stati:** D17.2, D17.13, D17.16 e D17.19 passi 4, 10c, 12 **IN CORSO**. La build supera la compilazione, il gate runtime `faa69c2` è fallito, il risultato sul nuovo commit resta da verificare; non dichiarare conclusa la suite. Requisiti tecnici precedenti preservati.

### D17.22.66 — Identità fisica canonica risolta centralmente dal catalogo (24 ottobre 2026)

**Gate precedente:** come confermato dall'utente, build del blocco precedente e Action Source Registry sul commit `985f6d7` sono entrambe verdi; questo non attesta le regressioni numeriche D17.20–D17.21.

**Blocco architetturale:** viene introdotta l'interfaccia nominale `canonical_source_name(component_name)` nel catalogo, che risolve il nome del contributo matematico storico nel nome della sorgente fisica completa. Ogni componente reologico, anche quando internamente ancora rappresentato da più addendi, produce una sola identità `rheology_<law>`; tutte le altre sorgenti schedulabili mantengono la loro identità catalogata, e i nomi sconosciuti producono la stringa vuota. La funzione concentra la traduzione legacy `model_id -> rheology_source_name` all'interno del catalogo; `group_source_names` nel registro non interroga più direttamente `term%model_id`.

**Regressioni:** aggiunta a `TESTS/TEST_SOURCE_REGISTRY/test_source_registry.f90` la verifica per ciascuna voce del catalogo: i componenti reologici devono ricondursi alla stessa legge completa, le altre voci conservano il proprio nome, e un nome sconosciuto non è accettato. L'estensione non altera le 864 configurazioni legacy e non modifica kernel o schedule.

**Limiti espliciti:** l'API non rende definitivo il vecchio catalogo a 34 addendi, né rimuove ancora `model_id` da tutti i consumatori, `SOURCE_TERM_*` o i campi superflui di `source_term_type`. L'accorpamento fisico vero richiede ancora una migrazione coordinata del dispatcher. **Stati:** D17.2, D17.5, D17.16 e D17.19 passi 3, 9, 10c **IN CORSO**. La nuova HEAD deve ancora superare la CI.

### D17.22.67 — Selezione della legge completa nel dispatcher IM/SI (24 ottobre 2026)

**Blocco coordinato catalogo, dispatcher, test:** `equation_terms_2d::group_has_role` non percorre più `group%term_names` per selezionare gli addendi di una legge. La membership temporale si decide esclusivamente tramite `group_contains_name(group,rheology_source_name(model_id))`: la sorgente schedulata è la **legge reologica completa**. Il catalogo espone `rheology_has_role(model_id,role)` per risolvere se l'evaluatore interno di quella legge include contributi `SOURCE_ROLE_IMPLICIT` o `SOURCE_ROLE_SEMI_IMPLICIT`. Non vengono creati due schedule indipendenti e resta intatta la valutazione numerica e l'ordine dei kernel IM e SI. La verifica distingue correttamente la presenza nel gruppo dalla disponibilità interna di un ruolo.

**Regressioni:** il test del registro controlla per ogni contributo legacy che il ruolo della sua legge sia riconosciuto da `rheology_has_role` e che il modello nullo non restituisca ruoli. Restano preservate tutte le 864 configurazioni legacy e le verifiche negative precedenti.

**Limiti:** la tabella legacy a 34 addendi rimane fonte temporanea dei metadati interni al catalogo; la sua sostituzione con una voce per legge e la definizione definitiva del tipo a cinque campi richiedono ancora la migrazione del codice fisico e dei controlli `role`, `model_id` e `SOURCE_TERM_*` rimasti. **Stati:** D17.2, D17.5, D17.16 e tabella D17.19 passi 3, 9, 10c **IN CORSO**. **Verifica:** questo commit non è stato ancora compilato né testato dalla CI; il gate precedente sul commit `985f6d7` risulta verde da conferma dell'utente. Nessun requisito normativo originale eliminato.

### D17.22.68 — Rimozione dei lookup numerici per le sorgenti termiche e di porosità nel dispatcher (24 ottobre 2026)

Il percorso `equation_terms_2d` non chiama più `group_has_term(group,term_id)` per determinare la presenza del raffreddamento radiativo, del raffreddamento convettivo e della dissipazione della pressione dei pori. Le tre selezioni utilizzano direttamente `group_contains_name` con i nomi attuali delle sorgenti. Vengono eliminate dal dispatcher la funzione `group_has_term`, la conversione locale `source_term_descriptor(term_id)` e le importazioni inutilizzate per queste tre costanti numeriche. Resta invariato il controllo di accoppiamento radiativo/convettivo, la selezione delle due vie real/complex e il comportamento fisico. La selezione reologica IM/SI continua a usare l'identità della legge completa introdotta in D17.22.67.

**Obiettivo aperto:** non è ancora stato sostituito il catalogo legacy delle 34 componenti con il `source_term_type` a cinque campi, né sono stati rimossi tutti i nomi di addendi dalla configurazione. La modifica non aggiunge tipi, adattatori o altre strutture. **Gate:** il blocco richiede verifica di build e test dopo il commit. D17.2, D17.5, D17.12a, D17.16 e D17.19 passi 9,10c restano **IN CORSO**; il testo normativo è stato conservato.

### D17.22.69 — Eliminazione delle API pubbliche di enumerazione numerica delle sorgenti (24 ottobre 2026)

**Blocco coordinato senza nuovi adattatori di produzione.** `source_catalog_2d` elimina le due funzioni `rheology_source_terms(model_id)` e `source_bundle_terms(bundle_id)` e le relative esportazioni `PUBLIC`. Entrambe restituivano elenchi di identificatori numerici di addendi fisici storici. Il codice di produzione, inclusi configurazione legacy/nativa e validazione del registro, era già migrato a `rheology_component_names` e `source_bundle_names`: pertanto le due API numeriche non sono più necessarie. I due corpi eliminati non vengono sostituiti nel solver.

**Test:** le enumerazioni dei vecchi indici rimangono soltanto nella suite di compatibilità, mediante `test_rheology_ids` e `test_bundle_ids`, funzioni locali al programma di test che interrogano i descrittori per mantenere i 864 scenari, i casi negativi di legge incompleta e batch atomico e i controlli di equivalenza. Sono esplicitamente strumenti di regressione, **non** nuove strutture o interfacce di produzione. Questo commit non modifica il dispatcher né i kernel fisici.

**Limite D17 ancora aperto:** le voci storiche dei 34 addendi, `source_term_type` con campi legacy e alcune costanti `SOURCE_TERM_*` permangono e devono essere eliminate nella conversione completa a una sorgente per legge. Non dichiarare conclusi D17.2 o 10c. **Stati:** D17.2, D17.5, D17.16 e tabella D17.19 passi 3,9,10c **IN CORSO**. **Gate:** compilazione e test della nuova HEAD non ancora eseguiti. Conservati i requisiti normativi D17 precedenti.

### D17.22.70 — Ripristino delle interfacce nominali accidentalmente eliminate (24 ottobre 2026)

Audit del commit `9627b2e`: la cancellazione per intervalli del vecchio catalogo ha rimosso involontariamente tre routine che sono ancora invocate dal codice di produzione e dai test: `rheology_source_name(model_id)`, `canonical_rheology_catalog_unique()` e `rheology_component_names(model_id)`. Questa condizione causa un guasto di compilazione indipendente dalle regressioni numeriche. Il ripristino delle tre definizioni originali conserva le due API numeriche `rheology_source_terms` e `source_bundle_terms` effettivamente eliminate; non vengono aggiunti nuovi adattatori e non viene alterata la fisica.

**Stato:** recupero della coerenza di compilazione da verificare alla prossima Actions. Il catalogo a 34 contributi, i campi legacy e le costanti numeriche restano **da eliminare**, non sono stati dichiarati completati. Gli stati D17.2, D17.5, D17.16 e la tabella D17.19 rimangono **IN CORSO**. Nessun requisito originale è stato rimosso.

### D17.22.71 — Schedulazione diretta delle leggi reologiche complete (24 ottobre 2026)

Modificati insieme `source_catalog_2d`, `source_configuration_2d` e `source_registry_2d`. I due percorsi di configurazione legacy e nativo ora registrano direttamente `rheology_<law>`, senza espansione in contributi matematici. Il catalogo riconosce il nome canonico della legge come singola identità schedulabile. Nel registro sono eliminati i controlli di completezza basati sulla presenza di tutti gli addendi reologici: viene richiesta direttamente la sorgente reologica completa. Le dichiarazioni native dei singoli addendi vengono rifiutate. Le formule fisiche implicite e semi-implicite non cambiano.

**Limiti e verifiche ancora necessarie:** le costanti numeriche e le 34 voci del catalogo storico devono ancora essere rimosse, insieme ai test legacy che costruiscono gruppi con addendi singoli. Questo passaggio è una migrazione funzionale del percorso di configurazione e registro, **non** il completamento dell'obiettivo a cinque campi. L'Action di compilazione e la suite del registro devono essere eseguite e i loro eventuali errori risolti prima di considerare stabile il nuovo blocco. Stati D17.2, D17.5, D17.16 e D17.19 **IN CORSO**.

### D17.22.72 — Adeguamento delle regressioni alle sorgenti reologiche complete (24 ottobre 2026)

Dopo D17.22.71 i gruppi del solver non contengono più addendi reologici schedulabili. È quindi aggiornata la suite `TEST_SOURCE_REGISTRY` affinché `catalog_names_for_test` trasformi le vecchie sequenze numeriche dei casi di regressione in nomi fisici canonici **deduplicati**, preservando i 864 scenari legacy. Gli assert relativi alla configurazione standard e al parsing nativo verificano esplicitamente che l'intera reologia lava sia una sola voce `rheology_lava`, non due componenti. Gli input legacy contenenti una sotto-selezione di addendi vengono normalizzati al nome dell'intera legge nel test, mentre gli input native che dichiarano il singolo addendo restano vietati. I controlli di batch termici continuano separatamente.

**Limitazioni:** il catalogo dei 34 addendi e il `source_term_type` esteso non sono ancora stati rimossi: questa è la sincronizzazione necessaria della suite alla nuova identità schedulabile, non la chiusura dell'intero D17. Le eventuali regressioni fisiche e le Actions sul nuovo commit sono ancora da verificare. D17.2, D17.5, D17.16, D17.19 **IN CORSO**.

### D17.22.73 — Distinzione dei test di input legacy e dei gruppi diretti (24 ottobre 2026)

Precisati i tre assert negativi nel test del registro: quando i vecchi nomi dei sottocomponenti vengono inseriti **direttamente** in `group%term_names` o nei SOURCE_GROUPS nativi, la nuova validazione deve respingerli. La sola conversione nominale per i 864 scenari legacy avviene in `catalog_names_for_test`. Nessuna modifica alla fisica. L'Action completa resta necessaria per verificare il comportamento. D17 resta IN CORSO.

### D17.22.74 — Correzione degli input nativi nei test di integrazione (24 ottobre 2026)

Aggiornata la routine **esistente e solo di test** `set_group`: quando genera il namelist nativo, riconduce tutti gli ID dei componenti di una stessa legge al nome `rheology_<law>`, deduplicando le occorrenze e mantenendo coerente `source_group_term_count`. L'insieme delle sorgenti non reologiche e degli input di test espliciti resta immutato. In questo modo i test del percorso nativo verificano la vera API a legge completa e non dichiarano più due addendi che il parser correttamente respinge. Il codice di produzione non acquisisce nuove strutture.

**Verifica pendente:** l'integrazione D17.22.71–74 richiede la compilazione e l'esecuzione della suite; non è ancora implementata la rimozione definitiva del catalogo a 34 voci.

### D17.22.75 — Isolamento dei test negativi dei target nativi e rimozione codice inutilizzato (24 ottobre 2026)

Nel log `test-run(9).log` il test del registro viene compilato ed eseguito, ma fallisce l'assert «missing implicit group targets rejected». La verifica precedente aveva inserito due componenti reologici legacy, deliberatamente non ammessi dal parser. La verifica dei target proseguiva senza ripristinare il nome canonico valido e intercettava quindi prima l'errore di componente reologico. **Correzione:** prima dei casi relativi a target assenti, sconosciuti e duplicati, ripristinati `source_group_term_count(3)=1`, `source_group_terms(:,3)=''` e `source_group_terms(1,3)=rheology_source_name(rheology_id_lava)`. Questo mantiene ogni caso negativo indipendente e consente di verificare il diagnostico corretto.

Dallo stesso log eliminate la variabile locale `law_names` inutilizzata da `source_registry_2d`, la funzione locale `name_membership` non più chiamata e la variabile `law_names` inutilizzata nel percorso legacy di `source_configuration_2d`. Nessun nuovo tipo, adattatore o livello di catalogazione introdotto; la fisica resta invariata.

**Gate:** nuova HEAD ancora da validare con Source Registry Tests; il fatto che il test sia arrivato al caso 646 dimostra compilazione ed esecuzione parziale ma non la conclusione della suite. La sostituzione definitiva del catalogo a cinque campi resta **IN CORSO**.

### D17.22.76 — Esclusione degli addendi reologici dai test di atomicità dei batch (24 ottobre 2026)

Il log `test-run(10).log` riporta `FAIL: legacy rheology batch normalizes to a complete source` in `test_source_registry.f90:332`. Il test combinava due contratti incompatibili: costruiva volutamente un sottoinsieme dei vecchi addendi reologici (in particolare il batch della legge lahar) e pretendeva che fosse un batch schedulabile al pari di più sorgenti fisiche distinte. In D17 la reologia ha una sola identità per legge: non va verificata con il test di integrità dei vecchi sotto-addendi.

Nella suite, il ciclo delle verifiche di *atomic batch* esclude ora gli elementi con `model_id/=0`, conservando la verifica negativa dei batch non reologici, incluso il raffreddamento termico. Rimangono separatamente i test delle identità reologiche complete e del rifiuto dei componenti nativi. Nessun codice di produzione o algoritmo fisico modificato. **Gate runtime sul nuovo commit pendente**, D17.19 e migrazione del catalogo a cinque campi ancora aperti.

### D17.22.77 — Rimozione diretta dell'espansione reologica nel parser nativo (24 ottobre 2026)

Il parsing di `SOURCE_GROUPS` non alloca più una lista temporanea di componenti per ogni sorgente e non esegue più un secondo ciclo sui sottotermini. Ogni voce immessa è validata e registrata una sola volta, direttamente con l'identità fisica canonica; i vecchi nomi di addendi reologici non sono schedulabili. Eliminata dal catalogo l'interfaccia `rheology_component_names`, non più utilizzata dalla configurazione o dal registro. La selezione dei metodi Newton/analitici, il rilevamento dei duplicati e la validazione dei target restano nel percorso esistente.

**Questo non completa D17:** il catalogo a 34 voci, gli identificatori storici e i metadati extra di `source_term_type` devono ancora essere eliminati con la revisione delle verifiche e dei test. Non si dichiara valida la nuova HEAD prima della CI. Nessuna struttura o funzione nuova è stata aggiunta al solver.

### D17.22.78 — Eliminazione dell'API numerica di appartenenza ai gruppi (24 ottobre 2026)

Rimossa da `source_registry_2d` la funzione pubblica `group_contains(group,term_id)`, che riceveva identificatori numerici del vecchio catalogo per convertirli in nomi. Il driver principale `IMEX_SfloW2D`, l'integrazione temporale `time_integration_2d` e i test della registry usano ora direttamente `group_contains_name(group,name)` per controllare le sorgenti fisiche di bordo e di scambio di massa; eliminati in questi due chiamanti gli import degli ID numerici corrispondenti. Non viene aggiunto alcun adattatore o livello di compatibilità. I kernel e le formule restano invariati.

**Stato:** è un'eliminazione effettiva di un'API numerica obsoleta, ma il catalogo a 34 contributi e il tipo a cinque campi devono ancora essere sostituiti. CI sulla nuova HEAD pendente; D17.19 resta IN CORSO.

### D17.22.79 — Eliminazione della doppia rappresentazione dei nomi nei source groups (24 ottobre 2026)

Refactoring coordinato di `source_registry_2d.f90`, `source_configuration_2d.f90`, `equation_terms_2d.f90` e dell'intera suite `TEST_SOURCE_REGISTRY`: `source_group_type` contiene **un solo array** `source_names(:)` per i nomi canonici delle sorgenti fisiche, senza `term_names(:)`. Il costruttore dei gruppi legacy e il parser nativo popolano direttamente `source_names`. Eliminati il popolamento secondario `populate_registry_term_names` e la conversione ripetuta dei nomi nella funzione `group_source_names`, che ora restituisce semplicemente l'array del gruppo. `group_contains_name` legge soltanto l'array canonico; nessuna tolleranza per il vecchio array di componenti matematici. Test e validazione operano sull'unica rappresentazione. Eliminato anche dall'import di `equation_terms_2d` il simbolo `group_contains` già rimosso nella milestone precedente. La fisica e la schedulazione non cambiano.

**Stato:** eliminazione definitiva di `term_names` da `source_group_type` implementata. Rimane **DA COMPLETARE** la sostituzione del catalogo legacy a 34 contributi e di `source_term_type` a cinque campi. Nessun nuovo adattatore o array di compatibilità. Gate: Actions compile diagnostics e source registry tests sulla nuova HEAD, ancora da verificare.

### D17.22.80 — Correzione compilazione post-RK dopo la rimozione di term_names (24 ottobre 2026)

La Action `D17 compile diagnostics` su `245ef1671cfdc119bf924ad8513b8947684eb8fc` fallisce durante la compilazione di `src/numerics/mass_exchange_step_2d.f90`: alle righe 167, 169, 170 rimanevano tre accessi a `step_groups(group_index)%term_names`, eliminato da `source_group_type` in D17.22.79. Il messaggio successivo `Expecting END IF statement` è secondario a questi errori di tipo. Correzione diretta: i tre accessi sono ora `%source_names` e la validazione continua a rifiutare ogni sorgente post-RK che non abbia `SOURCE_ROLE_MASS`. Nessun cambiamento alle formule o alla logica degli scambi di massa, nessun array sostitutivo. Compilation e Source Registry Actions sulla nuova HEAD da verificare. Il catalogo dei 34 addendi e la riduzione a cinque campi restano D17 aperti.

### D17.22.81 — Eliminazione di role/model_id/bundle_id dalla validazione del registro (24 ottobre 2026)

Rifattorizzata `source_registry_2d` eliminando l'accesso ai campi legacy `term%role`, `term%model_id` e `term%bundle_id`: identificazione della legge reologica completa, classificazione delle sorgenti post-RK di scambio di massa e validazione dell'atomicità radiativa–convettiva si basano sulle identità nominali delle sorgenti già registrate. Restano invariati i controlli di capacità, cadenza, destinazione delle equazioni e completezza del gruppo; i kernel fisici non vengono modificati. Nessun catalogo parallelo o adattatore introdotto.

**Stato:** il registro è ora indipendente da tre campi da eliminare nella struttura finale a cinque campi. Il catalogo legacy a 34 contributi e i test storici devono ancora essere trasformati. Le due GitHub Actions devono validare questa nuova HEAD.

### D17.22.82 — Correzione test negativi Newton/analitico e rimozione same_batch (24 ottobre 2026)

Il log `test-run(20261024-171638).log` relativo alla milestone D17.22.81 interrompe la suite alla verifica `legacy component is rejected from a directly specified Newton group`. L'asserzione del test confrontava erroneamente `error==''` pur richiedendo esplicitamente il rifiuto dei componenti reologici legacy. Ripristinata la condizione negativa corretta `error/=''` e aggiunto il controllo `.NOT.trial%is_initialized()` nei due test Newton e analitico. Rimossa la funzione privata `same_batch`, non più utilizzata dopo la verifica diretta dell'atomicità del raffreddamento nel registro. Nessuna alterazione del comportamento del solver o delle formule fisiche. CI sul nuovo commit ancora pendente; catalogo definitivo a cinque campi **NON ANCORA COMPLETATO**.

### D17.22.83 — Rifiuto degli addendi reologici legacy nella validazione nativa (24 ottobre 2026)

Il test eseguito dopo D17.22.82 fallisce ancora su `legacy component is rejected from a directly specified Newton group`. La causa effettiva è nel percorso `validate_registry`: la ricerca nel catalogo restituiva un descrittore valido per un addendo `SOURCE_TERM_RHEOLOGY_*` e lo considerava schedulabile, anche se il gruppo deve registrare esclusivamente la legge fisica completa `rheology_*`. Aggiunto un rifiuto esplicito degli addendi `SOURCE_TERM_RHEOLOGY_*` dopo il lookup, prima della registrazione e del controllo delle capacità, per entrambi i metodi Newton e analitico. I test negativi rimangono invariati e non sono stati indeboliti. Questa guardia è destinata a scomparire insieme al catalogo degli addendi; non è un adattatore per accettare la vecchia rappresentazione. CI da verificare, catalogo definitivo a cinque campi ancora aperto.

### D17.22.84 — Sostituzione strutturale del catalogo (24 ottobre 2026)

**Implementata la riscrittura coordinata** del catalogo: da 34 addendi matematici/dipendenti a **24 sorgenti fisiche schedulabili**. Esattamente 10 leggi reologiche complete, una sorgente termica completa `thermal_cooling`, una dissipazione di pressione, cinque sorgenti esplicite e sette scambi di massa. `source_term_type` ha esattamente `name`, `target`, `has_explicit`, `has_newton`, `has_analytic`. Eliminati i vecchi ID `SOURCE_TERM_*`, i relativi descrittori dipendenti, i bundle e il catalogo parallelo reologico. Registro e configurazione usano target nominali; mass exchange controlla il target fisico; dispatcher usa nomi fisici. I ruoli matematici SI/IM della legge reologica rimangono privati alla logica di valutazione e non sono sorgenti schedulabili. Riscritta la suite del catalogo con verifica del contratto fisico, test di rigetto addendi e ciclo completo di 864 configurazioni legacy.

**Attenzione:** questo è un commit architetturale ampio intenzionalmente pubblicato PRIMA del gate CI. La suite nativa storica è stata sostituita e le verifiche puntuali in precedenza contenute devono essere ripristinate in forma coerente col nuovo catalogo; non asserire equivalenza numerica senza eseguire i test. Eventuali errori compiler/CI vanno corretti sul nuovo assetto, senza ripristinare cataloghi o adattatori.

### D17.22.85 — Identificatori numerici dei soli blocchi-equazione

I vecchi ID numerici delle singole sorgenti fisiche restano eliminati. Sono mantenute le costanti `SOURCE_TARGET_VELOCITY`, `SOURCE_TARGET_TEMPERATURE`, `SOURCE_TARGET_PORE_PRESSURE` per gli indici dei target della partizione IMEX: non sono voci del catalogo delle sorgenti né campi aggiuntivi in `source_term_type`. Questa distinzione evita modifiche non necessarie alle interfacce della partizione delle equazioni.

### D17.22.86 — Capacità degli input nativi collegate al nuovo catalogo (24 ottobre 2026)

Il log `compile(7).log` sulla HEAD `e14ad96` segnala in `source_input_2d.f90` un simbolo rimosso: `SOURCE_PHYSICAL_TERM_COUNT`. L'input nativo ora deriva `SOURCE_INPUT_MAX_GROUPS` e `SOURCE_INPUT_MAX_TERMS` da `SOURCE_TERM_COUNT`, che rappresenta le 24 sorgenti fisiche complete. Non si ripristina l'alias legacy. Gli errori successivi sulle dimensioni degli array e sulla dichiarazione NAMELIST sono conseguenze del simbolo non risolto. Compilazione successiva ancora da verificare.

### D17.22.87 — Copertura regressiva del catalogo fisico e parsing case-insensitive (24 ottobre 2026)

Ripristinati test della configurazione di default lava (gruppi analitici energia e pressione, Newton velocità, esplicito di bordo, post-RK deposizione), della lettura NAMELIST, del conflitto tra SOURCE_GROUPS e campi legacy IMEX, del rifiuto di duplicati, della costruzione di cinque gruppi nativi e dei casi negativi su addendi matematici, target sconosciuti e sorgenti mancanti/duplicate. La suite continua a verificare le 864 assegnazioni legacy. Corretto il parser della configurazione nativa: confronta i nomi delle sorgenti in modo case-insensitive mediante normalizzazione in minuscolo, senza alias legacy e con identità canoniche di catalogo. Gli input SOURCE_GROUPS contenenti nomi maiuscoli sono quindi confrontati con i nomi fisici minuscoli. Nessun cambio ai kernel numerici. Verifiche CI pendenti.

### D17.22.88 — Correzione duplicazione thermal_cooling nella configurazione legacy (24 ottobre 2026)

Il test CI del commit `3878eb1` ha scoperto un residuo della fusione radiative + convective: `initialize_configured_sources` inseriva due volte `thermal_cooling` nella lista legacy, facendo fallire la validazione di unicità. Ora registra una sola sorgente fisica termica; entrambe le componenti del raffreddamento restano valutate nel kernel. Eliminata anche la variabile inutilizzata `j` dal test. La CI successiva è da verificare.

### D17.22.89 — Audit sostanziale stati D17.0–D17.18a e matrice D17.19

Effettuata riconciliazione degli stati delle 23 sezioni D17.x con l'implementazione a 24 sorgenti e cinque campi; eliminati i paragrafi di stato obsoleti che menzionavano il catalogo a 34 addendi come ancora presente. Le righe della tabella D17.19 distinguono ora sottoblocchi **COMPLETATI — STRUTTURA**, lavori **PARZIALMENTE COMPLETATI**, lavori **DA INIZIARE** e gate di test ancora pendenti. Gli esiti CI verificati valgono fino alla HEAD `4403f05`; i cambi `3878eb1` e `055db6c` attendono nuovamente Actions. Le attività multidisciplinari non sono state promosse a completate senza evidenza end-to-end: D17.20/D17.21, flag termico indipendente, state_2d, mass-exchange e rheology_newmodel restano aperti.

### D17.22.90 — Attivazione termica indipendente dalla reologia (24 ottobre 2026)

Implementato `thermal_cooling_flag` opzionale in `TEMPERATURE_PARAMETERS`. Quando assente, è inizializzato in fase di lettura secondo il comportamento storico: attivo solo per lava; se presente, permette di disattivare il raffreddamento con lava o attivarlo con le altre reologie. Aggiornati layout IMEX energia (`init_problem_param`), catalogazione delle sorgenti attive dei percorsi legacy/nativo, selezione analitica, kernel implicito e validazione preventiva dei coefficienti termici. La legge fisica Stefan–Boltzmann/convezione rimane invariata; il dispatcher chiama ora direttamente la funzione di raffreddamento indipendentemente dal modello. Aggiunta regressione metadata 2×2 Voellmy/lava e off/on. Stato implementativo D17.12a e D17.19–8a aggiornati inline. **CI non ancora eseguita; equivalenza dei bilanci energetici non ancora dimostrata.**

### D17.22.92 — Correzione di quantità di moto in state_2d (24 ottobre 2026)

Realizzata la migrazione della routine pura `correct_friction_momentum` da `rheology_2d` a `state_2d`, senza cambiare né l'algoritmo di proiezione (raggio del vettore precedente, fattore limitato in [0,1]) né i due punti di chiamata del solver. La prima usa il riferimento `q_known(2:3)` immediatamente dopo SI, la seconda `momentum_before_newton` dopo Newton. Eliminati l'esportazione dalla reologia e il re-export transitorio da `equation_terms_2d`; `nonlinear_solver_2d` importa direttamente la funzione di correzione da `state_2d`. Aggiunti test isolati di direzione, saturazione, inversione e riferimento nullo in `TEST_STATE_INITIALIZATION`. Aggiornato l'avanzamento della riga D17.19–10b; la restante D17.15a è parziale. Verifica CI pendente.

### D17.22.93 — Audit degli stati dopo gate termico e migrazione momentum (24 ottobre 2026)

Riconciliati gli stati di tutte le 23 sezioni D17.x e le righe pertinenti di D17.19 con la HEAD `ecef5c5`. Il gate confermato dalle due Actions arriva a `61fee22`: comprende il nuovo flag termico e i relativi test 2×2. Il refactoring `ecef5c5` sposta la correzione di quantità di moto in `state_2d` e aggiunge test, ma manca conferma CI: D17.19–10b resta **IMPLEMENTATO / CI PENDENTE**, non chiuso scientificamente. D17.5 è completato strutturalmente; D17.2 e 10c rimangono parziali per identità numeriche reologiche e il requisito di estensione `rheology_newmodel`. Restano esplicitamente aperti D17.20–D17.21, mass-exchange, i quattro assemblatori e la finalizzazione dello stato.

### D17.22.94 — Estrazione coordinata dell'assemblatore post-RK di mass exchange

L'intero corpo di `eval_mass_exchange_terms` è stato trasferito senza alterare le formule da `equation_terms_2d` al nuovo modulo `mass_exchange_assembly_2d`. Il modulo `mass_exchange_step_2d` dipende ora direttamente dall'assemblatore fisico; aggiunte sorgente e dipendenze esplicite a `src/Makefile.am`. Conservati i limitatori esistenti, la gestione dell'inventario, le correzioni per le celle di sorgente e l'aggiornamento topografico. Aggiornati D17.14 e D17.19–10 come parzialmente completati: non si dichiara ancora la sostituzione dell'esecutore STEP né la scomposizione nei moduli di processo. Gate compilazione/CI ancora da eseguire.

### D17.22.95 — Gate verde su assemblatore e correzione linker (24 ottobre 2026)

L'utente conferma entrambe le Actions verdi su `ce18d8c`. Il linking iniziale di `98ab45e` falliva per la mancante importazione di `eval_pore_pressure_gas_loss`; il successivo `ce18d8c` ha ripristinato l'import. Questo gate non equivale ai confronti quantitativi D17.20.

### D17.22.96 — Estrazione del processo carrier-loss con limitatore condiviso (24 ottobre 2026)

Creato `src/physics/mass_exchange/mass_exchange_carrier_loss_2d.f90`. Trasferiti senza cambiare le espressioni: perdita liquida prescritta, perdita indotta dalla deposizione, gas-loss Darcy condizionata dal packing e **un solo limitatore dopo la somma**, dipendente anche dai rate di deposizione e da `dt`. L'assemblatore richiama la routine nello stesso punto e mantiene il batch indivisibile. Aggiornati `src/Makefile.am` e le dipendenze Fortran. **D17.14 e D17.19–10: PARZIALMENTE COMPLETATI**, con primo sottoprocesso estratto; restano gli altri processi, la convergenza verso `STEP/EXPLICIT` comune, i test di bilancio e i gate D17.20–21. **Nuova CI da eseguire, nessuna equivalenza numerica dichiarata.**

### D17.22.97 — Ripristino dipendenze build carrier-loss (24 ottobre 2026)

Il log `compile(9).log` per `e564b7c` segnala che `mass_exchange_carrier_loss_2d.mod` manca. La causa è il `Makefile.am`: la regola di dipendenza del modulo esisteva, ma mancavano sia la registrazione in `IMEX_SfloW2D_SOURCES` sia la dipendenza diretta dell'assemblatore. Corrette entrambe senza modificare i kernel numerici. D17.14 e D17.19 passo 10 **PARZIALMENTE COMPLETATI**; verifica della CI successivamente confermata verde dall'utente per entrambe le Actions sul commit `bebbf3a`; test numerici D17.20–21 ancora aperti.

### D17.22.98 — CI verde sul build carrier-loss e estrazione dell'erosione (24 ottobre 2026)

L'utente conferma **entrambe le Actions verdi** sulla HEAD `bebbf3a`. Estratto il blocco erosivo in `src/physics/mass_exchange/mass_exchange_erosion_2d.f90`, con la chiamata originale a `eval_erosion_rates`, la soglia su `erosion_coeff`, lo zero nei casi disattivi e il limite sull'inventario `erodible` determinato dalla legge esistente. L'assemblatore ora richiama il processo prima del calcolo di deposizione, senza spostare i controlli preliminari di cella asciutta e cutoff dei solidi. Registrati sorgente e dipendenze in `src/Makefile.am`. **D17.14 e D17.19 passo 10: PARZIALMENTE COMPLETATI**; carrier-loss ed erosione separati, deposizione e coordinamento STEP generale aperti. Nuova CI ancora pendente; confronti scientifici e inventari D17.20–21 non chiusi.

### D17.22.99 — Estrazione della deposizione solida e gate verde erosione (24 ottobre 2026)

L'utente conferma entrambe le Actions verdi su `051ca03`. Creato `src/physics/mass_exchange/mass_exchange_deposition_2d.f90`: la valutazione locale della viscosità Sutherland (quando attiva), la selezione delle proprietà di sedimentazione, il ciclo per classe, la limitazione tramite `eval_solid_deposition_rate` e il controllo del rate negativo sono trasferiti in un unico processo fisico. Costanti storiche di sedimentazione (0.6, 4.65) e formule invariate. In `mass_exchange_assembly_2d` resta la medesima sequenza: erosione, proprietà miscela, deposizione, carrier-loss; il limitatore carrier-loss è applicato dopo che tutte le deposizioni sono note. Build aggiornato con voce sorgente e dipendenza `.mod`. **D17.14 e D17.19–10: PARZIALMENTE COMPLETATI**; tre processi ora estratti, routing STEP/EXPLICIT e verifiche scientifiche D17.20–21 aperti. Nuovo gate CI pendente.

### D17.22.100 — Estrazione della vaporizzazione dopo gate verde deposizione (24 ottobre 2026)

L'utente conferma entrambe le Actions verdi su `d90411a`. Creato `src/physics/mass_exchange/mass_exchange_vaporization_2d.f90`, che conserva le chiamate a `eval_vaporization_mass_rate` e `eval_vaporization_energy_rate` e gli aggiornamenti di massa, energia e gas nel vettore conservativo, condizionati da `gas_flag`, `liquid_vaporization_flag` e `B_zone`. La chiamata rimane prima delle sorgenti stocastiche e di pressione dei pori, preservandone la dipendenza dal contributo totale di massa. Aggiornate lista sorgenti e dipendenze Fortran. **D17.14 e D17.19 passo 10 restano PARZIALMENTE COMPLETATI:** quattro blocchi fisici isolati, ma ancora da completare la riorganizzazione per processi e il routing STEP/EXPLICIT unificato. Nuova CI e verifiche quantitative D17.20–21 ancora pendenti.

### D17.22.101 — CI verde per vaporizzazione ed estrazione air entrainment (24 ottobre 2026)

L'utente conferma entrambe le Actions verdi sulla HEAD `f20e97c`. Il kernel del contributo di air entrainment è stato estratto in `src/physics/mass_exchange/mass_exchange_entrainment_2d.f90`; preservati l'attivazione mediante `entrainment_flag`, i controlli su spessore e Richardson, la chiamata a `eval_air_entrainment_rate` e il rate nullo nei casi inattivi. Il calcolo resta dopo erosione/deposizione/carrier-loss e prima dell'assemblaggio conservativo. Aggiornato il Makefile con sorgente e dipendenze esplicite. **D17.14 e D17.19 passo 10: PARZIALMENTE COMPLETATI**; i blocchi erosione, deposizione, carrier-loss, vaporizzazione ed entrainment sono separati, ma non sono ancora soddisfatti l'assemblaggio generico `STEP/EXPLICIT`, il coordinamento definitivo e i gate D17.20–21. La nuova CI è pendente.

### D17.22.102 — Centralizzazione del bilancio conservativo e topografico (24 ottobre 2026)

L'utente conferma entrambe le Actions verdi sul commit `6ba4fa6`. Estratto da `mass_exchange_assembly_2d` il blocco con la costruzione del vettore conservativo massa/momento/energia/massi solidi, l'addendo di vaporizzazione, le sorgenti derivate di trasporto stocastico e pressione dei pori, e il rate di variazione del letto. Ora risiede in `src/physics/mass_exchange/mass_exchange_balance_2d.f90`, chiamato **una sola volta** dopo il calcolo completo e limitato dei rate dei cinque processi. Sono preservati formule, ordine degli addendi e dipendenza delle due sorgenti derivate dal termine di massa dopo vaporizzazione. Aggiornate lista delle sorgenti e dipendenze Fortran in `src/Makefile.am`. **D17.14 e D17.19 passo 10: PARZIALMENTE COMPLETATI**, non chiusi: `mass_exchange_step_2d` è ancora esecutore dedicato e non esiste ancora routing uniforme `STEP/EXPLICIT`. **CI sul nuovo commit pendente** e verifiche scientifiche D17.20–21 non concluse.

### D17.22.103 — Primo passaggio all'esecutore ordinario STEP (24 ottobre 2026)

L'utente conferma entrambe le Actions verdi su `4924f18`. Spostato il corpo dell'esecuzione post-RK in `src/numerics/step_source_executor_2d.f90`; `simulation_2d` chiama ora `execute_step_source_groups` direttamente, dopo RK e prima dell'avanzamento del clock, preservando observer e sequenza di inventari e aggiornamento topografico. `mass_exchange_step_2d` resta un modulo di compatibilità che riesporta `apply_mass_exchange_step` e `apply_mass_exchange_increment` dal nuovo esecutore, evitando di rompere chiamanti storici. Aggiornato `Makefile.am` con ordine e dipendenze dei moduli. **D17.14/D17.19–10 PARZIALI:** il passo sposta il proprietario temporale ma NON implementa ancora l'assemblatore comune `source_terms_2d::eval_expl_terms`, né elimina il vincolo del batch mass-exchange o la vecchia API. Nuova CI da eseguire, bilanci D17.20–21 aperti.

### D17.22.104 — Bridge fisico STEP/EXPLICIT nel modulo equazionale (24 ottobre 2026)

L'utente conferma due Actions verdi sul commit `bf2513c`. Il modulo equazionale `equation_terms_2d` espone ora `eval_step_mass_exchange_terms`: interfaccia cell-local che restituisce rate conservativi, rate erosione/deposizione, perdita della fase continua e rate topografico senza applicare direttamente inventari o geometria. `step_source_executor_2d` la usa al posto dell'import diretto di `mass_exchange_assembly_2d`, preservando un solo punto di applicazione post-RK. Introdotto il controllo esplicito `SOURCE_METHOD_EXPLICIT` sui gruppi STEP. Aggiornate dipendenze Fortran in `Makefile.am`. **D17.14 / D17.19–10 rimangono PARZIALI**: questo è un bridge di migrazione, non l'assemblatore definitivo `source_terms_2d::eval_expl_terms`; quest'ultimo gestisce ancora soltanto gli stadi. Richiesti nuova CI e test quantitativi D17.20–21, senza dichiarare equivalenza end-to-end.

### D17.22.105 — Rimozione delle dipendenze dall'API STEP transitoria (24 ottobre 2026)

Confermate dall'utente entrambe le Actions verdi sulla HEAD `8f2ce06`. Aggiornati i due test di `TESTS/TEST_MASS_EXCHANGE` per usare direttamente `step_source_executor_2d::execute_step_source_groups`, mantenendo per chiarezza un alias locale al vecchio nome della routine. Il modulo `mass_exchange_step_2d` non e' piu' incluso nelle sorgenti di `src/Makefile.am` e sono eliminate le dipendenze di build verso il modulo transitorio; il file storico rimane nel repository ma non e' compilato. **D17.14 / D17.19 passo 10: PARZIALMENTE COMPLETATI**: eliminato l'utilizzo dell'adapter dal build e dai test, ma manca ancora il dispatch comune in `eval_expl_terms`. Test numerici D17.20–21 e nuova CI pendenti.

### D17.22.106 — API esplicita generica per STAGE e STEP (24 ottobre 2026)

L'utente conferma entrambe le Actions verdi sulla HEAD `af9d954`. In `equation_terms_2d` e' definita l'interfaccia pubblica Fortran `eval_expl_terms` con due `MODULE PROCEDURE`: `eval_stage_expl_terms`, ottenuta rinominando la precedente routine esplicita di stadio, ed `eval_step_expl_terms`, ottenuta rinominando il bridge di mass exchange. `step_source_executor_2d` richiama l'API generica anziche' un bridge pubblico separato. Le firme incompatibili vengono risolte staticamente dal compilatore e conservano la cadenza propria; la transazione atomica sugli inventari rimane una sola per cella post-RK. **D17.14 e D17.19 passo 10 restano PARZIALMENTE COMPLETATI**: l'unificazione nominale dell'assemblatore esplicito e' attuata, ma mancano ancora selezione STEP direttamente dal contratto di gruppo nella routine fisica, completamento del coordinamento dei processi e verifiche di equivalenza numerica D17.20–21. Nuova CI pendente.

### D17.22.107 — Audit della collocazione fisica mass exchange (24 ottobre 2026)

L'utente segnala correttamente che l'organizzazione richiesta in D17.14 non era completata: nonostante l'estrazione di cinque processi e bilancio, due moduli dedicati esclusivamente al mass exchange erano ancora nella radice `src/physics/`. Dopo le due Actions verdi su `f1cf8a1`, trasferiti **senza modifiche ai sorgenti** `mass_exchange_laws_2d.f90` e `mass_exchange_assembly_2d.f90` sotto `src/physics/mass_exchange/`; aggiornate in `src/Makefile.am` le voci sorgente e tutte le dipendenze corrispondenti. I nomi dei MODULE Fortran sono invariati, quindi le istruzioni `USE` restano valide. Restano nella radice solo i moduli fisici *condivisi*, come `thermal_laws_2d` e `pore_pressure_laws_2d`, e il modulo equazionale generale `equation_terms_2d`. **D17.14 e D17.19 passo 10 restano PARZIALMENTE COMPLETATI** per le parti funzionali non ancora chiuse; la *collocazione dei moduli dedicati* e' invece sistemata. Nuova CI da eseguire; confronti D17.20–21 aperti.

### D17.22.108 — Eliminazione delle leggi fisiche centralizzate e assegnazione ai processi (24 ottobre 2026)

L'utente chiarisce il requisito D17.14: non basta collocare file nella stessa directory, ma **ogni processo deve possedere direttamente i propri kernel fisici**, senza dipendenze verso un aggregatore `mass_exchange_laws_2d`. Audit sul sorgente: l'aggregatore esponeva otto procedure. Migrati `eval_erosion_rates` in `mass_exchange_erosion_2d`; `settling_velocity`, `eval_settling_fluid_properties`, `eval_solid_deposition_rate`, `eval_lahar_deposition_rate` e `eval_hindered_settling_deposition_rate` in `mass_exchange_deposition_2d`; `eval_continuous_phase_loss_rate` in `mass_exchange_carrier_loss_2d`; `eval_air_entrainment_rate` in `mass_exchange_entrainment_2d`. Preservati i corpi originali delle routine, senza riscrivere le espressioni fisiche e senza duplicare routine. Eliminato `mass_exchange_laws_2d.f90` e aggiornata la lista sorgenti e le dipendenze in `src/Makefile.am`; rimossi import inutilizzati in `equation_terms_2d`. Adeguati `TESTS/TEST_RHEOLOGY_CELL_LOCAL/test_source_modules.f90` e `TESTS/TEST_MASS_EXCHANGE/test_mass_exchange.f90` ai moduli proprietari. La dipendenza dalla legge `eval_lahar_fluid_viscosity` rimane nel modulo deposizione in quanto proprieta' reologica condivisa, non in un aggregatore mass exchange. **D17.14 / D17.19 passo 10 non ancora chiusi** per contratto temporale e verifiche quantitative. Due Actions su `f1cf8a1` verdi, nuova CI pendente.

### D17.22.109 — Fix linker dopo distribuzione delle leggi fisiche (24 ottobre 2026)

Il log utente `compile(10).log` sulla HEAD `13c608d` conferma che ogni sorgente Fortran e' compilato, ma il linking fallisce con `undefined reference to eval_pore_pressure_gas_loss_` generato da `mass_exchange_carrier_loss_2d`. Nell'estrazione delle leggi era scomparsa l'istruzione `USE pore_pressure_laws_2d, ONLY : eval_pore_pressure_gas_loss`, trasformando la chiamata a procedura di modulo in chiamata a procedura esterna implicita. Ripristinato l'import esplicito, senza modificare le formule fisiche. La dipendenza `pore_pressure_laws_2d.$(OBJEXT)` in `Makefile.am` era gia' corretta. **D17.14/D17.19 passo 10 ancora PARZIALI.** Nuova CI richiesta; test quantitativi D17.20–21 aperti.

### D17.22.110 — Quattro interfacce uniformi dei processi mass exchange (24 ottobre 2026)

Le due Actions sul commit `edd0b38` sono verdi. Introdotte nelle cinque famiglie `mass_exchange_erosion_2d`, `mass_exchange_deposition_2d`, `mass_exchange_carrier_loss_2d`, `mass_exchange_entrainment_2d`, `mass_exchange_vaporization_2d` quattro interfacce pubbliche per i trattamenti `explicit`, `semi_implicit`, `implicit`, `analytic`. Ogni interfaccia `explicit` e' un alias generico della routine fisica esistente, quindi non altera il calcolo o la firma dell'implementazione. Gli ingressi non supportati hanno le stesse firme del processo e generano `fatal_error`, mai un contributo zero silenzioso. Modificati `mass_exchange_assembly_2d` e `mass_exchange_balance_2d` affinche' il percorso operativo chiami gli ingressi `explicit`. Il batch rimane esclusivamente `STEP/EXPLICIT` e la selezione di metodi non supportati resta vietata dalla validazione. **D17.14 e D17.19 passo 10 sono ancora PARZIALI:** i contratti di processo sono stati uniformati, ma non sono conclusi i bilanci quantitativi e l'unificazione del coordinamento secondo il requisito completo D17.14. CI del nuovo commit ancora da eseguire.

### D17.22.111 — Nomi effettivi uniformi delle quattro interfacce (24 ottobre 2026)

L'utente osserva che l'intervento precedente aveva creato **alias pubblici** `eval_<processo>_explicit` verso routine con nomi non uniformi. Per erosione, deposizione, carrier-loss, entrainment e vaporizzazione sono state rinominate **le procedure implementative** rispettivamente in `eval_erosion_explicit`, `eval_deposition_explicit`, `eval_carrier_loss_explicit`, `eval_entrainment_explicit`, `eval_vaporization_explicit`, eliminando le corrispondenti dichiarazioni `INTERFACE ... MODULE PROCEDURE`. Le tre altre procedure in ciascun modulo seguono esattamente gli stessi prefissi e i suffissi `semi_implicit`, `implicit`, `analytic`. Resta lecito che le leggi sottostanti usino nomi descrittivi diversi (per esempio `settling_velocity`, `eval_erosion_rates`) poiche' sono kernel interni e non le quattro interfacce del trattamento. Le chiamate dell'assemblatore e del bilancio usano gia' i nuovi nomi e non richiedono ulteriori modifiche. **Convenzione mass exchange applicata; D17.14/D17.19 passo 10 ancora PARZIALI** per il contratto di esecuzione e le prove numeriche. CI su nuovo commit pendente.

### D17.22.112 — Audit di tutte le famiglie sorgente: quattro trattamenti e semantica eval/solve (24 ottobre 2026)

L'utente conferma entrambe le Actions verdi su `32a5e4d` e richiede uniformita' su **tutti** i processi fisici, non solo mass exchange. Sono stati inventariati i moduli in `src/physics`: le dieci leggi reologiche, il dispatcher `rheology_2d`, `thermal_laws_2d`, `pore_pressure_laws_2d`, i cinque processi mass exchange e i moduli di sola infrastruttura/assemblaggio. Convenzione semantica scelta: **`eval_<sorgente>_explicit`**, **`eval_<sorgente>_semi_implicit`**, **`eval_<sorgente>_implicit`** calcolano termini/rate; **`solve_<sorgente>_analytic`** evolve lo stato con una soluzione analitica. Un kernel costitutivo di supporto (es. viscosita', velocita' di sedimentazione, coefficiente di attrito) non e' automaticamente una sorgente e non necessita del contratto quattro-trattamenti. I dieci modelli reologici espongono ora tutti i quattro nomi uniformi; dove un metodo e' supportato, il nome pubblico e' per questa fase un'interfaccia generica che riusa la routine numerica storica senza alterarne la firma; dove non e' supportato, la procedura segnala `fatal_error`. Aggiunte interfacce equivalenti al dispatcher reologico, al raffreddamento lava e alla pressione dei pori. **LIMITI / LAVORO RESIDUO IMPORTANTE:** queste ultime interfacce non costituiscono ancora una ridenominazione delle procedure implementative e le firme fisiche dei trattamenti sono differenti; inoltre `thermal_laws_2d` e `pore_pressure_laws_2d` contengono piu' sorgenti configurabili e richiedono un audit per sorgente e non solo per modulo. Evitare di chiudere D17.0a o la riga D17.19 corrispondente fino all'eliminazione di questi scostamenti e a test del dispatcher. Nuova CI pendente; nessuna modifica intenzionale a formule o scheduling.

### D17.22.113 — Nomi concreti dei modelli reologici e sorgenti termiche/pore (24 ottobre 2026)

Due Actions verdi sul commit `02b2214`. Audit del `SOURCE_CATALOG`: dieci sorgenti reologiche, due dissipazioni termica/pore-pressure, cinque sorgenti esplicite generali (bottom circular/fissural, radial, geometry curvature, boundary injection) e sette sorgenti nominate di mass exchange. I dieci moduli reologici e il dispatcher contenevano ancora alias pubblici verso vecchi nomi: rinominate le routine concrete, eliminati gli alias `MODULE PROCEDURE` dedicati al trattamento e aggiornati i chiamanti di produzione e i quattro test dedicati. `eval_rheology_<modello>_semi_implicit`, `eval_rheology_<modello>_implicit` e `solve_rheology_<modello>_analytic` sono ora nomi procedurali concreti dove le leggi esistono. Analogamente le implementazioni per `thermal_cooling` e `pore_pressure_dissipation` hanno nomi concreti `eval_..._implicit` e `solve_..._analytic`. I trattamenti non supportati conservano il diagnostico esistente. **Audit rigoroso: D17.0a non completo.** Le sorgenti STAGE/EXPLICIT di iniezione e curvatura del catalogo non sono ancora suddivise in quattro interfacce a livello del modulo fisico; gli stub non supportati vanno armonizzati anche nelle firme (e SI/IM dei metodi NEWTON). Inoltre, le cinque sorgenti di mass exchange sono rappresentate dal catalogo con granularita' diversa dalle interfacce aggregate `carrier_loss`; il coordinamento e il contratto nominale vanno verificati. Anche gli ingressi analitici non supportati dei cinque processi mass exchange sono stati rinominati secondo la convenzione `solve_..._analytic`. Il nuovo commit richiede CI, test quantitativi D17.20–21 ancora aperti. Non dichiarare completata l'uniformazione generale.

### D17.22.114 — Interfacce concrete per tutte le 24 sorgenti configurabili (24 ottobre 2026)

L'utente conferma due Actions verdi su `02b2214`. Da questa base sono stati realizzati due interventi coordinati. Primo: rinomina delle procedure concrete nei dieci moduli reologici e nelle leggi termiche/pore (D17.22.113), con adattamento dei chiamanti. Secondo: audit completo del `SOURCE_CATALOG` e creazione di tutte le entrate fisiche per le **24 sorgenti** e i **quattro trattamenti** (`eval_<nome_catalogo>_explicit`, `eval_<nome_catalogo>_semi_implicit`, `eval_<nome_catalogo>_implicit`, `solve_<nome_catalogo>_analytic`). I tre contributi di iniezione al fondo circolare, da fissura e radiale sono stati estratti dai rispettivi blocchi di `eval_stage_expl_terms` in procedure pubbliche operative; l'assemblatore le invoca solo quando il gruppo contiene il nome, senza ripetere i contributi. Il kernel di curvatura e la routine di iniezione al bordo hanno nomi concreti conformi, con chiamanti aggiornati. In mass exchange, i tre rate richiesti di carrier-loss (liquid, deposition, Darcy) hanno interfacce specifiche e si sommano **prima dell'unico limitatore di disponibilita'**, mantenendo l'accoppiamento fisico. I trattamenti non disponibili dal catalogo producono `fatal_error`, non un rate zero fittizio. Per sorgenti con `has_newton=.TRUE.`, una meta' del trattamento SI/IM che non contribuisce produce invece **zero valido** e non un errore di metodo: corretta questa distinzione nelle sei reologie interessate, nel raffreddamento e nella dissipazione pore-pressure. Le procedure dei kernel analitici sono denominate `solve_*`, poiche' evolvono lo stato, quelle di rate `eval_*`. Aggiunto `TESTS/ACCEPTANCE/test_d17_physical_interfaces.py` e relativo passo nella Action `d17-source-registry.yml`, che controllano presenza, unicita', pubblicita' e disponibilita' di tutte le 96 procedure fisiche. **Stato: implementazione nominale completata, nuova CI non ancora eseguita e prove end-to-end aperte.** D17.14 e gate D17.20–21 rimangono parziali; non si dichiara equivalenza numerica senza evidenza.


### D17.22 — Gate strutturale sull'ordine delle interfacce (24 ottobre 2026)

Il controllo `TESTS/ACCEPTANCE/test_d17_physical_interfaces.py`, gia eseguito dalla Action `d17-source-registry.yml`, verifica ora per ciascuna sorgente canonica del catalogo non solo l'esistenza delle quattro procedure e il rifiuto dei metodi non supportati, ma anche che le quattro routine risiedano nello stesso modulo e che dichiarazioni `PUBLIC` e definizioni seguano l'ordine `explicit`, `semi_implicit`, `implicit`, `analytic`. Questo e un gate statico: non dimostra equivalenza fisica/numerica. Il routing finale di mass exchange verso il solo assemblatore comune resta aperto ai passi 10 e 12 di D17.19, cosi come D17.20 e D17.21. L'esito della nuova Action e da verificare dopo il commit.


### D17.22 — Consolidamento del confine STEP/physics (24 ottobre 2026)

Aggiornato `TESTS/ACCEPTANCE/test_d17_physical_interfaces.py` per verificare staticamente che `step_source_executor_2d` invochi una sola volta per cella `source_terms_2d::eval_step_source_explicit`, non importi i vecchi moduli di assemblaggio, e che l'entry point comune deleghi al batch fisico atomico. Il controllo blocca regressioni strutturali del percorso introdotto nel commit `cf5823c6` e corretto in `33d606de`. Nessuna modifica alle formule fisiche o alla sequenza RK/STEP: rimangono aperti il trasferimento definitivo degli assemblatori e i test numerici D17.20. Compilazione e test sul nuovo commit da verificare.


### D17.22 — Riorganizzazione delle sorgenti di massa in mass_source_terms (24 ottobre 2026)

I sette moduli precedentemente in `src/physics/mass_exchange/` sono stati collocati in `src/physics/mass_source_terms/`. Il modulo principale `mass_exchange_assembly_2d` e stato rinominato `mass_source_terms_2d`, preservando la procedura `eval_mass_exchange_terms` e gli stessi processi fisici, limiti e bilanci. Gli altri sei file rimangono separati per funzione (erosione, deposizione, perdita di fase continua, entrainment, vaporizzazione e assemblaggio conservativo), senza mutare i loro kernel. Il chiamante generale `source_terms_2d` e le dipendenze esplicite in `src/Makefile.am` sono stati aggiornati. Questa riorganizzazione e strutturale e non costituisce la chiusura dei passi D17.19 n. 10, D17.20 e D17.21: restano l'accorpamento semantico delle sorgenti schedulabili di deposizione/perdita indotta, il dispatch finale e le regressioni quantitative. Il nuovo commit richiede compilazione e test di registro.


### D17.22 — Separazione delle sorgenti di massa: scambi e iniezioni (24 ottobre 2026)

Correzione della precedente organizzazione: `src/physics/mass_source_terms/` e sostituita da due cartelle indipendenti. `src/physics/mass_exchange_terms/` contiene il coordinatore `mass_exchange_terms_2d.f90` e i moduli separati per erosione, deposizione, perdita carrier, entrainment, vaporizzazione e bilancio conservativo. `src/physics/mass_injection_terms/` contiene i quattro file `mass_boundary_injection_2d.f90`, `mass_bottom_circular_injection_2d.f90`, `mass_bottom_fissural_injection_2d.f90` e `mass_radial_injection_2d.f90`; ogni file conserva esattamente le quattro procedure di trattamento della propria sorgente. Il modulo `mass_injection_coefficients_2d.f90` conserva i quattro coefficienti di finestra temporale condivisi, evitando duplicazione dello stato RK. `source_terms_2d` rimane il punto di accesso comune; aggiornato `src/Makefile.am` con i nuovi moduli e il relativo ordine di dipendenza. Le routine fisiche sono ricollocate senza cambi di formula; nuova build e test sul commit sono ancora necessari. D17.20 e D17.21 restano aperti.


### D17.22 — Migrazione dell'assemblatore STAGE/EXPLICIT (24 ottobre 2026)

L'implementazione di `eval_stage_expl_terms` e stata trasferita in `source_terms_2d` con il nome `eval_stage_source_explicit`: selezione per nome dei termini di curvatura, iniezione circolare, iniezione fissurale e iniezione radiale, controllo della cadence/method e guardia di cella asciutta rimangono invariati. La vecchia `eval_stage_expl_terms` di `equation_terms_2d` e ora un sottile adattatore per mantenere compatibile l'interfaccia generica pubblica `eval_expl_terms`, usata dai chiamanti. I kernel di iniezione restano nella cartella `mass_injection_terms`; i coefficienti mediati sul passo rimangono nel modulo condiviso. Avanzamento architetturale parziale: gli assemblatori SEMI_IMPLICIT, IMPLICIT e ANALYTIC, la rimozione definitiva dell'adattatore legacy e l'accettazione numerica D17.20/D17.21 restano aperti. Compilazione e test sul commit da verificare.


### D17.22 — Migrazione assemblatori SEMI_IMPLICIT e IMPLICIT (24 ottobre 2026)

Le implementazioni degli assemblatori `eval_semi_implicit_terms` e `eval_implicit_terms` sono state ricollocate in `source_terms_2d` con i nomi comuni `eval_source_semi_implicit` ed `eval_source_implicit`. In `equation_terms_2d` permangono gli adattatori compatibili con gli argomenti originali. Spostati nel modulo comune anche i controlli `validate_newton_selection` e `group_has_role`; mantenuti i percorsi di valutazione real e complex per il Newton e il complex-step, le selezioni nominali per rheology/thermal/pore-pressure e le chiamate ai medesimi kernel fisici. Nessuna modifica intenzionale ai metodi o alle leggi. Migrazione architetturale parziale: assemblatore ANALYTIC, eliminazione dei bridge legacy e regressioni numeriche D17.20/D17.21 ancora aperti. Il commit corrente richiede compilazione e test CI.


### D17.22 — Migrazione dell'assemblatore ANALYTIC (24 ottobre 2026)

`eval_analytic_terms` e stato migrato in `source_terms_2d` come `eval_source_analytic`, insieme al tipo `analytic_term_context_type`. `equation_terms_2d` mantiene solo un adattatore con la firma originale e riesporta il tipo per compatibilita dei chiamanti. I tre percorsi di soluzione analitica (reologia/momento, raffreddamento/energia e dissipazione pressione interstiziale) conservano gli stessi kernel, condizioni di selezione e status di errore. Tutti e quattro i trattamenti principali sono ora implementati in `source_terms_2d` ma i bridge legacy, l'unificazione completa del dispatch D17 e le verifiche numeriche quantitative D17.20/D17.21 sono ancora aperti. Compilazione e test del nuovo commit da verificare.


### D17.22 — Eliminazione dei cinque adattatori legacy degli assemblatori (24 ottobre 2026)

Eliminati da `equation_terms_2d` l'overload generico `eval_expl_terms` e i relativi adattatori `eval_stage_expl_terms` e `eval_step_expl_terms`, insieme a `eval_semi_implicit_terms`, `eval_implicit_terms` ed `eval_analytic_terms`. `time_integration_2d` ora importa e chiama direttamente `source_terms_2d::eval_stage_source_explicit`, mentre `nonlinear_solver_2d` importa e chiama direttamente `eval_source_semi_implicit`, `eval_source_implicit`, `eval_source_analytic` e `analytic_term_context_type`. Rimossi i relativi export legacy; nessun wrapper sostitutivo. Conservata in `equation_terms_2d` la gestione temporanea dei coefficienti di attivazione, ancora chiamata dalla pianificazione temporale: essa richiede una migrazione autonoma prima di poter eliminare ogni dipendenza dal modulo delle equazioni. Prossima pulizia: riduzione degli import/ri-export fisici e audit delle dipendenze residue. Compilazione e test del presente commit ancora da verificare, e regressioni numeriche D17.20/D17.21 non completate.


### D17.22 — Pulizia dei residui equation/source e cache delle iniezioni (24 ottobre 2026)

`set_source_time_coefficients`, `window_mean`, `window_integral` e `interval_mean_velocity` sono stati trasferiti da `equation_terms_2d` a `mass_injection_coefficients_2d`, insieme alla cache dei coefficienti utilizzata dalle iniezioni. `simulation_2d` ora importa direttamente il setter dalla cartella `mass_injection_terms`; il modulo delle equazioni non deve piu gestire i coefficienti di attivazione. Rimossi da `equation_terms_2d` gli import e gli export legacy delle 20 interfacce fisiche gia ricollocate nei moduli sorgenti. Aggiornate le dipendenze Fortran in `src/Makefile.am` e la descrizione del modulo delle equazioni. Stato: separazione degli assemblatori piu netta, ma verifiche CI sul presente commit ancora necessarie. Rimangono controlli completi sugli import storici e regressioni quantitative D17.20/D17.21; nessuna equivalenza numerica dichiarata.


### D17.22 — Audit dei riferimenti residui in equation_terms (24 ottobre 2026)

Dopo il passaggio dei coefficienti temporali a `mass_injection_coefficients_2d`, audit dei soli import e delle dichiarazioni di `equation_terms_2d`: eliminati gli import delle leggi reologiche, termiche e pore-pressure non piu chiamate dal modulo, gli import del catalogo/registro delle sorgenti, le costanti ed indici `parameters_2d` non piu referenziati e un gruppo di commenti Doxygen orfani lasciati dalla rimozione degli assemblatori. Il modulo mantiene le cinque procedure legate alle equazioni (`init_problem_param`, `eval_inertial_flux`, `eval_hydrostatic_coefficient`, `eval_local_speeds_direction`, `limit_component_mass_flux`) e la riesportazione storica di `curvature_acceleration` per non cambiare i chiamanti senza audit separato. Nessuna formula modificata. Questo commit richiede nuovo controllo di compilazione e registro; D17.20/D17.21 rimangono aperti. L'audit non certifica l'assenza di codice inutilizzato in tutto il repository.


### D17.22 — Audit interfacce e rimozione riferimenti legacy (24 ottobre 2026)

Il commit successivo a `a974b5b1` aggiunge alla suite statica D17 un gate contro la ricomparsa degli assemblatori `eval_expl_terms`, `eval_stage_expl_terms`, `eval_step_expl_terms`, `eval_semi_implicit_terms`, `eval_implicit_terms`, `eval_analytic_terms` e degli import legacy delle sorgenti fisiche attraverso `equation_terms_2d`. La verifica CI di questo nuovo gate e ancora pendente; fino a `a974b5b1` compilazione e test registro sono confermati verdi. La chiusura D17.20/D17.21 resta esplicitamente aperta.


### D17.22 — Audit mirato degli export e dei commenti Doxygen (24 ottobre 2026)

Dopo conferma di entrambe le Actions verdi sul commit `78070634`, ripulita la documentazione di `source_terms_2d`: la descrizione ora copre i quattro metodi e non indica piu l'inizializzatore temporale come bridge legacy. In `mass_exchange_terms_2d` rimossa una doppia intestazione `\\brief` relativa a una routine inesistente e un riferimento storico di migrazione non piu utile. In `mass_injection_coefficients_2d` introdotti `PRIVATE` di default e `PUBLIC` espliciti per setter e quattro coefficienti condivisi, lasciando privati i kernel ausiliari `window_mean`, `window_integral` e `interval_mean_velocity`. Nessuna modifica alle formule o al dispatch. Questo nuovo commit richiede Actions; le regressioni scientifiche D17.20–D17.21 restano aperte. L'audit degli export delle altre famiglie fisiche non e ancora concluso.


### D17.22 — Pulizia Doxygen residuo in equation_terms_2d (24 ottobre 2026)

Audit puntuale di `equation_terms_2d.f90`: eliminati i due blocchi Doxygen orfani in coda al modulo (curvatura del fondo e coefficienti temporali delle sorgenti), gia migrati rispettivamente in `source_terms_2d` e `mass_injection_coefficients_2d`. Eliminata anche la doppia importazione locale di `pore_pressure_flag` e `thermal_cooling_flag` in `init_problem_param`, gia visibili a livello di modulo; ridotta la spaziatura superflua. Restano le cinque procedure del modulo equazioni e la riesportazione di `curvature_acceleration`, da conservare fino all'audit dei chiamanti. Nessun algoritmo modificato; verifica CI del nuovo commit pendente. I gate numerici D17.20/D17.21 restano aperti.


### D17.22 — Rimozione della riesportazione di curvatura da equation_terms (24 ottobre 2026)

Dopo Actions verdi confermate sul commit `dd6ce6ff`, rimossi da `equation_terms_2d` `USE rheology_2d, ONLY : curvature_acceleration` e il corrispondente `PUBLIC`. La funzione resta definita e pubblica nel proprietario `rheology_2d`, mentre gli assemblatori fisici la importano direttamente in `source_terms_2d`. Aggiunto un controllo statico in `test_d17_physical_interfaces.py` contro la riesportazione legacy. Nessuna formula o firma di procedura modificata. Il nuovo commit richiede conferma delle Actions; D17.20/D17.21 restano aperti.


### D17.22 — Audit Doxygen dei moduli mass-exchange (24 ottobre 2026)

Dopo la validazione CI del commit `21f5165b`, corretti i commenti introduttivi dei moduli `mass_exchange_terms_2d`, `mass_exchange_balance_2d`, `mass_exchange_deposition_2d` e `mass_exchange_carrier_loss_2d` per esplicitare ownership del batch, limitatore condiviso e distinzione tra rate fisici e aggiornamenti di inventario. Eliminato linguaggio storico ambiguo su `legacy fluid properties` nel modulo deposition senza modificare la legge. Nessuna interfaccia fisica cancellata: per le sorgenti configurabili resta il contratto a quattro metodi D17. Pulizia esclusivamente documentale in attesa di verifica CI; audit globale delle API e D17.20/D17.21 ancora aperti.


### D17.22 — D17.19 passo 10a: estrazione della finalizzazione post-RK (24 ottobre 2026)

Dopo commit `0da5c7cb`, estratto da `time_integration_2d::imex_RK_solver` il blocco finale di validazione delle celle attive in `state_2d::finalize_imex_cell_state`. Il solver mantiene l'assemblaggio conservativo e la chiamata observer sullo stato grezzo, poi delega a `state_2d` senza alterare intenzionalmente l'ordine delle operazioni: NaN, clipping delle masse negative da roundoff, spessore negativo (soglia -1e-7), masse solide negative, conversione conservativa-primitive, proiezione momentum asciutto, soglia temperatura 273 K, riparazione conservativa delle masse solide con conversione storica da q_previous, azzeramento nodata. Aggiunte dipendenze Makefile esplicite da state_conversion e diagnostics per state_2d e da state_2d per time_integration_2d. Questo completa la migrazione della finalizzazione delle celle attive per il passo 10a, ma NON il trattamento halo né la verifica di equivalenza su stati intermedi/finali. Due Actions richieste; validazione numerica D17.20–D17.21 aperta.


### D17.22 — D17.19 passo 10a: inizializzazione dry dell'halo e degli stage (24 ottobre 2026)

Dopo Actions verdi sul commit `b2f61de2`, estratta in `state_2d::initialize_dry_physical_cell` l'inizializzazione fisica asciutta (tutte le componenti zero, temperatura `T_ambient`) e sostituite le due sequenze duplicate in `time_integration_2d`: inizializzazione della reconstruction halo ad inizio RK e fallback dry dello stato fisico durante gli stage espliciti. Preservati l'ordine dei loop OpenMP, le condizioni sui target e i valori assegnati, senza modifiche a formule o controlli. Assieme a `finalize_imex_cell_state`, risultano quindi migrati in `state_2d` sia il blocco finale della cella attiva sia l'inizializzazione dry dell'halo; restano da eseguire test quantitativi di equivalenza e audit di ulteriori preparazioni dry in moduli spaziali. CI di questo commit pendente, D17.20/D17.21 non chiusi.


### D17.22 — Audit coordinato degli stati D17.x e della tabella D17.19 (24 ottobre 2026)

**Base verificata:** commit `35be9f0a`, entrambe le Actions confermate verdi dall'utente. Aggiornate le schede di stato di D17.0–D17.18a e tutte le righe operative di D17.19. In particolare, D17.19 passo 10a: finalizzazione post-RK, proiezione dry, controlli di massa/temperatura e inizializzazione asciutta halo/stage sono in `state_2d` con CI verde; mancano ancora test quantitativi degli stati intermedi/finali, observer e condizioni limite, pertanto il punto NON è dichiarato chiuso scientificamente. Passo 10b: migrazione della correzione momentum verificata dalla CI, prove before/after ancora da chiudere. Sono chiusi come **struttura**, non come equivalenza numerica, i passi 2 e 3. D17.20 e D17.21 rimangono aperti; nessun esito scientifico è dedotto dalle sole Actions. Le voci storiche di D17.22 restano una cronologia e non prevalgono su questo audit consolidato.


### D17.22 — Unificazione inizializzazione dry post-RK/STEP mass exchange (24 ottobre 2026)

Dopo audit coordinato su `35be9f0a` e conferma delle due Actions verdi, rimossa la duplicazione dello stato fisico asciutto nel modulo `step_source_executor_2d`: sia il fallback durante la valutazione dei rate STEP, sia il ramo dry dopo `apply_mass_exchange_increment`, ora invocano `state_2d::initialize_dry_physical_cell`. Il contratto resta invariato (tutti i campi zero, temperatura ambiente `T_ambient`), e gli aggiornamenti conservativi e di inventario rimangono nello stesso ordine. Aggiunta dipendenza Makefile `step_source_executor_2d -> state_2d` e guardia statica D17 sulle due chiamate. Questo avanza D17.15a/D17.17 e la verifica dei percorsi dry del passo 10a, senza concludere le regressioni quantitative D17.20–D17.21. Nuovo commit non ancora validato tramite CI.


### D17.22 — Gate strutturale sulla sequenza finale RK (24 ottobre 2026)

Dopo CI verde confermata sul commit `e4f6dfca`, aggiunto `TESTS/ACCEPTANCE/test_d17_state_finalization.py` e integrato nella Action manuale `D17 source registry tests`. Il gate ispeziona la chiamata dell'observer sullo stato conservativo grezzo prima della finalizzazione, l'ordine dei controlli NaN/clipping/spessore/masse/conversione/proiezione dry/temperatura/roundoff/nodata, le soglie storiche, il recupero da `q_previous`, e la centralizzazione dry nei percorsi RK halo/stage e mass-exchange STEP. Questo e un **test strutturale**, non un confronto numerico degli stati: D17.19 passo 10a resta implementato con gate quantitativi incompleti. Il commit che introduce il gate necessita verifica CI; D17.20–21 ancora aperti.


### D17.22 — Primo test numerico eseguibile della finalizzazione dry (24 ottobre 2026)

Aggiunto `TESTS/ACCEPTANCE/test_d17_state_numeric.f90` ed esecuzione nella Action `D17 compile diagnostics`, dopo la build completa: verifica valori concreti dell'inizializzatore fisico dry con `T_ambient`, finalizzazione conservativa della cella a massa nulla, proiezione delle componenti di velocità, azzeramento di uno spessore negativo entro la tolleranza storica (`-1e-8`) e maschera `nodata`. Il programma è collegato agli oggetti Fortran costruiti dal progetto ed eseguito con runtime bounds checks; nessuna formula fisica modificata. Il test non copre ancora casi wet con costituzione multicomponente, roundoff delle classi solide, errori fatali attesi, observer dinamico né confronto legacy/nativo. La relativa CI è pendente e D17.20–D17.21 rimangono aperti.


### D17.22 — Prima migrazione nominale delle capabilities reologiche (24 ottobre 2026)

Dopo entrambe le Actions verdi su `c4e5cb03`, il catalogo introduce `rheology_name_has_role(source_name, role)`: la selezione delle contribuzioni semi-implicite e implicite usa i nomi canonici delle dieci sorgenti reologiche. `rheology_has_role(model_id,role)` rimane temporaneamente come adattatore che converte l'ID nel nome e delega alla nuova funzione; ciò consente migrazione incrementale senza modifiche alle formule o al dispatch fisico. Inserito controllo statico per impedire che la selezione nominale contenga ID reologici. **D17.10c rimane PARZIALE**: `rheology_model_id`, `resolve_rheology_model`, il mapping ID-nome e gli `SELECT CASE(model_id)` fisici sono ancora da migrare; non dichiarare la rimozione degli ID completata. CI di questo commit pendente.


### D17.22 — Migrazione sostanziale del dispatch reologico a nomi canonici (25 ottobre 2026)

Dopo la conferma di entrambe le Actions verdi su `c79250ac`, convertiti insieme i tre dispatcher di `rheology_2d` (`eval_rheology_semi_implicit`, `eval_rheology_implicit`, `solve_rheology_analytic`) da `INTEGER model_id` e `SELECT CASE` di costanti numeriche a `CHARACTER(*) model_name` e `SELECT CASE(TRIM(model_name))` su nomi canonici. `source_terms_2d` trasmette i nomi ai tre dispatcher e utilizza `rheology_name_has_role` per selezionare le contribuzioni del gruppo. La scelta delle implementazioni fisiche REAL, COMPLEX e analitiche non dipende quindi più dagli ID. Aggiunto gate statico su tutti e tre i dispatcher/chiamanti. Il mapping temporaneo `rheology_model_id -> rheology_source_name` resta al confine legacy: gli input, gli altri chiamanti e le costanti persistenti devono essere migrati per chiudere D17.10c. Le formule non sono modificate, CI del presente commit ancora pendente, D17.20–D17.21 aperti.


### D17.22 — Dispatch reologico dai nomi dei gruppi (25 ottobre 2026)

Dopo Actions verdi confermate su `220522a2`, il modulo `source_terms_2d` seleziona per i gruppi D17 la legge attiva direttamente da `group%source_names`, tramite `group_rheology_name`, e ne determina i contributi SI/IM con `rheology_name_has_role`. Il solve analitico usa lo stesso nome del gruppo; i dispatcher REAL, COMPLEX e analitici ricevono `selected_model_name` invece della conversione numerica come scelta operativa. Solo le chiamate prive di un gruppo esplicito conservano transitoriamente il fallback `rheology_source_name(rheology_model_id)` per compatibilità legacy. Rafforzato l'audit delle interfacce per verificare il collegamento nominale. **D17.10c resta parziale:** persistono ID nel parser/modello di input e nel fallback senza gruppo, e rimangono da validare estensibilità e regressioni scientifiche. Questo commit necessita nuova CI.


### D17.19-10c — Rimozione degli identificativi numerici della reologia (25 ottobre 2026)

Implementata una migrazione coordinata dei moduli Fortran e dei test: `rheology_model` rimane il nome canonico della legge; rimossi da `parameters_2d` il campo persistente `rheology_model_id` e le costanti `rheology_id_*`; rimossi gli adattatori nome–ID nel catalogo e in `rheology_2d`. Input, registry, fallback delle sorgenti, leggi di deposizione, raffreddamento e stochastic usano nomi. Il parser normalizza l'input legacy testuale e valida il modello direttamente rispetto al catalogo. Aggiornati i test che impostavano ID numerici e aggiunto un gate statico per impedire il ritorno di questi simboli nel codice e nei test. **Stato: implementazione proposta; verifica di compilazione e due GitHub Actions non ancora confermate verdi.** L'assenza dei nomi legacy non costituisce prova di equivalenza numerica: D17.20–D17.21 restano aperti.


### D17.19 passo 11 — Audit hot loop e consolidamento della validazione preventiva (25 ottobre 2026)

Rimossi da `source_terms_2d` i controlli ripetuti di appartenenza al metodo `NEWTON` nelle valutazioni semi-implicita e implicita e il controllo del metodo/cadenza `EXPLICIT/STAGE` dalla valutazione esplicita per cella. Queste proprietà sono validate all'assemblaggio dei gruppi in `source_registry_2d`, prima della simulazione: il kernel riceve solo gruppi già ammessi, senza revalidarne la forma ad ogni valutazione di rate/Jacobiana. Sono mantenuti i controlli numerici per cella e gli errori relativi alle interfacce fisiche non supportate. Le verifiche di schedule eseguite una volta per RK-step/STEP-batch restano intenzionalmente in `time_integration_2d` e `step_source_executor_2d` perché non appartengono al ciclo interno sulle celle. **Stato: audit mirato e rimozione dei controlli di metodo dall'hot path implementati; CI e regressioni quantitative da verificare.** D17.20–D17.21 restano aperti.


### D17.19 passo 11 — Audit finale dei quattro moduli esecutivi (25 ottobre 2026)

Audit di `time_integration_2d`, `nonlinear_solver_2d`, `source_terms_2d` e `step_source_executor_2d`: i controlli numerici, la verifica del registro inizializzato, le invarianti degli indici e il controllo del batch atomico rimangono. Sono rimossi dall'esecutore STEP i controlli duplicati di metodo, cadenza e target dei termini, già imposti da `source_registry_2d` durante la costruzione/validazione preventiva dei gruppi. Eliminato un doppio lookup identico di `thermal_cooling` durante la valutazione implicita; la presenza della sorgente rimane selezionata per nome canonico. Eliminata anche la doppia guardia sul descrittore vuoto nel validatore del registro. La prima pulizia hot-loop di `source_terms_2d` è stata verificata con due Actions verdi su `d5e3c932`. **Questa revisione coordinata richiede nuove Actions; D17.20 e D17.21 rimangono aperti.**


### D17.18a — Fixture di estensibilita `rheology_newmodel` (25 ottobre 2026)

Aggiunto il modello dimostrativo Newton-only `rheology_newmodel` con legge di drag lineare `S_x=-0.1 rho u`, `S_y=-0.1 rho v`, SI nullo intenzionale e rate impliciti COMPLEX derivabili con complex-step. Le interfacce EXPLICIT e ANALYTIC non supportate falliscono esplicitamente; la capability NEWTON e registrata nel catalogo. Collegati solo il catalogo, il modulo di legge, il dispatcher della famiglia reologica, la build e i test: **nessuna modifica a `source_terms_2d`, `source_configuration_2d`, `source_registry_2d` o `nonlinear_solver_2d`**. La ricerca dei nomi reologici nel catalogo non dipende piu dalla cardinalita fissa di 10 elementi. Il test `test_model_names` verifica il nome canonico e le capability. Questa e una prova strutturale con legge numerica illustrativa, non un nuovo modello costitutivo scientificamente calibrato. **CI e integrazione end-to-end ancora da verificare; D17.20–21 aperti.**


### D17.20 — Prima suite numerica indipendente per estensione Newton (25 ottobre 2026)

Aggiunto `TESTS/TEST_SOURCE_REGISTRY/test_newmodel_numerics.f90`, eseguito dai profili strict e optimized di `run_test.sh` e quindi dalla Action `d17-source-registry.yml`. Verifica quattro stati fisici, compreso riposo, velocita quasi nulla e cambi di segno; capability nominali della nuova legge; contributo SI intenzionalmente nullo; rate implicito contro riferimento analitico indipendente; Jacobiana 2x2 complex-step con derivata non perturbativa; Newton sul residuo conservativo confrontato con la mappa esatta backward Euler `q_after=q_before/(1+0.1*dt)` e residuo finale. Questa suite testa il kernel di legge e una soluzione Newton di riferimento **indipendente**, non ancora il percorso integrato `nonlinear_solver_2d`, ne la Jacobiana accoppiata 3x3, i casi mass exchange e legacy/native. **Stato: test inseriti nei due profili, esito CI da verificare**. Gli altri casi della matrice D17.20 e D17.21 restano aperti.


### D17.20 — Verifica quantitativa attraverso il dispatcher di produzione (25 ottobre 2026)

Aggiunta la regressione Fortran `TESTS/ACCEPTANCE/test_d17_production_rheology.f90`, collegata contro gli oggetti reali della build Autotools nella Action `d17-compile.yml`. La prova esercita il **dispatcher `rheology_2d` di produzione** per `rheology_newmodel` in quattro condizioni, confrontando i rate con una formula indipendente, l'interfaccia semi-implicita nulla, due colonne della Jacobiana COMPLEX tramite complex-step e la consistenza di un passo implicito con soluzione esatta. Si aggiunge alla prova matematica 2x2 e ai controlli del catalogo precedenti. **Non equivale ancora a un test del solver `nonlinear_solver_2d` sul sistema accoppiato 3x3**, ne a una regressione legacy/native, RK, mass exchange o termica. **Stato: test e workflow implementati, nuova Action da eseguire; gate D17.20–21 aperti.**


### D17.20 — Esecuzione Newton di produzione per gruppi 2x2 e 3x3 (25 ottobre 2026)

Aggiunto `TESTS/ACCEPTANCE/test_d17_production_newton.f90` alla Action `d17-compile.yml`, compilato contro gli oggetti della build e chiamante `nonlinear_solver_2d::solve_rk_step`. La suite costruisce gruppi nominali NEWTON di dimensione 2 (momentum) e 3 (momentum e temperatura), controlla convergenza/diagnostica e confronta il drag lineare `rheology_newmodel` con la soluzione esatta del passo diagonale `q_m^{n+1}=q_m^n/(1+0.1*a_ii*dt)`; verifica anche che massa e temperatura rimangano invariati quando i coefficienti di raffreddamento sono nulli. Il sistema 3x3 esercita il percorso del solver ma **non testa ancora le derivate incrociate fisiche**, poiche la temperatura non genera contributi attivi nel fixture. Il test non e ancora validato da CI: eventuali difetti di configurazione dovranno essere corretti dopo i log. D17.20–D17.21 restano aperti.


### D17.20 — Newton 3x3 con accoppiamento fisico lava/raffreddamento (25 ottobre 2026)

Esteso `TESTS/ACCEPTANCE/test_d17_production_newton.f90` con un gruppo NEWTON contenente `rheology_lava` e `thermal_cooling` sulle equazioni di momento x/y ed energia. Il test attiva la dipendenza reologica dalla temperatura e un raffreddamento convettivo non nullo; verifica che la derivata incrociata complex-step momento/energia sia non nulla, esegue `nonlinear_solver_2d::solve_rk_step` sul blocco unico 3x3 e controlla il residuo conservativo con l'assemblatore implicito di produzione. **Test implementato; CI da verificare.** Le regressioni legacy/native, RK, mass exchange e le due correzioni SI/IM restano aperte.


### D17.20 — Test delle due correzioni reologiche nel solver di produzione (25 ottobre 2026)

Esteso `TESTS/ACCEPTANCE/test_d17_production_newton.f90` con due chiamate a `nonlinear_solver_2d::solve_imex_stage` (percorso IMEX completo che applica il SI prima del Newton). Il caso lava usa un termine di snervamento SI deliberatamente elevato, verifica matematicamente che il contributo grezzo invertirebbe il momento e controlla che il risultato del solver non inverta il verso e non modifichi le equazioni non coinvolte. Il caso `rheology_newmodel` verifica la contrazione corretta lungo la direzione originaria dopo il Newton, confrontando anche il risultato con la soluzione implicita esatta. I test controllano gli effetti finali dei due punti di proiezione; **non strumentano i valori intermedi privati del solver** e non sostituiscono l'audit dell'ordine delle chiamate. **Esito della nuova CI ancora da verificare**. D17.20–21 ancora aperti.


### D17.20 — Gate di equivalenza legacy/native per la configurazione sorgenti (25 ottobre 2026)

Esteso `TESTS/TEST_SOURCE_REGISTRY/test_source_registry.f90`: prima di inizializzare la configurazione nativa equivalente, il test costruisce e salva separatamente la configurazione legacy a tre gruppi impliciti (temperatura ANALYTIC, pori ANALYTIC, reologia NEWTON), l'esplicito STAGE al bordo e il mass exchange STEP finale. Confronta tutti i gruppi STAGE e STEP tra le due configurazioni, per metodo/cadenza/fase, insieme ordinato degli indici delle equazioni e insieme completo dei nomi fisici canonici. Gli ID di gruppo sono esclusi intenzionalmente dal confronto perche definiti diversamente dai due input. Rimangono in essere i test negativi sul parser nativo. **Questo gate dimostra equivalenza strutturale del piano di esecuzione, non ancora equivalenza numerica delle simulazioni con input legacy/native. CI da verificare.**
