# D17 — Organizzazione unificata dei termini sorgente e dei gruppi di integrazione

**Specifica di revisione autonoma | 24 ottobre 2026 | stato: PROPOSTA DA IMPLEMENTARE**

> Questo documento riguarda **solo D17**. Sostituisce, per gli aspetti in conflitto, la precedente impostazione D17 che consentiva più solve Newton indipendenti all'interno dello stesso gruppo. Non costituisce una dichiarazione che il codice sia già allineato né un'autorizzazione a cambiare la fisica o gli schemi RK-IMEX. Le indicazioni «attuale» riportano il comportamento osservato nel repository `demichie/IMEX_SfloW2D_dev` sul branch `master` consultato il 24 ottobre 2026 e nella documentazione D17 allegata alla roadmap. Le indicazioni «richiesto» definiscono invece il contratto da implementare.

## D17.0 — Obiettivo e principio di semplificazione

D17 ha un solo obiettivo funzionale: per **ogni termine sorgente fisico attivo** decidere **(a) come integrarlo**, **(b) quando integrarlo** e **(c) a quale gruppo di integrazione appartiene**. Il raggruppamento serve a stabilire quali equazioni trattare *insieme*, in particolare quando è richiesto un solve Newton accoppiato. Nessun livello addizionale di raggruppamento o scheduling deve essere introdotto senza una necessità concreta.

Le tre proprietà sono:

1. **Metodo**: `EXPLICIT`, `NEWTON` (trattamento semi-implicito più implicito esistente), oppure `ANALYTIC` (aggiornamento con soluzione analitica disponibile).
2. **Cadenza**: `STAGE` (all'interno dell'integrazione RK-IMEX, a ciascuno stadio previsto dal tableau) oppure `STEP` (dopo il completamento dell'intero passo RK).
3. **Gruppo**: insieme di termini e relative equazioni trattati come una sola unità numerica, con metodo e cadenza uniformi.

**Regola vincolante: un gruppo = un metodo + una cadenza + una unità di soluzione.** Nei gruppi `NEWTON` questo implica **un unico sistema non lineare con tutte le equazioni target del gruppo**, mai una sequenza di solve Newton separati dentro lo stesso gruppo. Se due problemi devono essere risolti separatamente, appartengono a gruppi differenti. Un gruppo `ANALYTIC` contiene **un solo aggiornamento/solve analitico**. Un gruppo `EXPLICIT` assembla i contributi espliciti, senza un solve implicito.

**Esempio vincolante.** Velocità `NEWTON` e temperatura `ANALYTIC`, entrambe `STAGE`, **non** possono condividere il gruppo: servono due gruppi. Velocità e temperatura entrambe `NEWTON` e `STAGE` possono condividere un gruppo: allora Newton risolve **congiuntamente** il sistema per `(q_2,q_3,q_4)` (salvo layout specifico del modello). Se devono essere risolte una dopo l'altra, servono due gruppi Newton ordinati. Il raggruppamento non autorizza a trascurare derivate incrociate realmente presenti nel residuo.


**Stato di avanzamento D17.0 — IN CORSO (24 ottobre 2026).** Architettura definita nella specifica; implementazione complessiva ancora parziale. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

## D17.0a — Contratto uniforme delle leggi fisiche

**Ogni famiglia di termini sorgente è implementata nel proprio modulo fisico, con lo stesso contratto di quattro interfacce: `explicit`, `semi_implicit`, `implicit`, `analytic`.** Le interfacce esistono anche quando una particolare famiglia non supporta quel trattamento. Il catalogo dichiara esclusivamente il supporto dei tre metodi `EXPLICIT`, `NEWTON`, `ANALYTIC`; le interfacce dei metodi non disponibili devono, se chiamate erroneamente, produrre un errore diagnostico (mai un rate nullo fittizio). Il parser/validatore impedisce già prima della simulazione di selezionarle. `NEWTON` richiede **entrambe** le interfacce `semi_implicit` e `implicit`: la prima può restituire un contributo nullo se non occorre un trattamento semi-implicito. Il catalogo non deve descrivere come la legge ripartisce internamente i contributi. **Le quattro interfacce fisiche non introducono un quarto metodo numerico**.

Il requisito vale indistintamente per reologia, raffreddamento/termica, pressione dei pori, trasferimenti di massa e ogni altra legge sorgente presente o futura. `source_terms_2d` invoca solo le interfacce pertinenti alla configurazione validata e **assembla tutti i contributi**, senza eccezioni o canali laterali per una famiglia fisica. Le routine esistenti possono essere mantenute come implementazione interna delle interfacce: evitare duplicazione delle formule e wrapper inutilmente stratificati.


**Stato di avanzamento D17.0a — IN CORSO (24 ottobre 2026).** Contratto delle quattro interfacce definito; uniformazione delle leggi non completata. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

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


**Stato di avanzamento D17.1 — IN CORSO (24 ottobre 2026).** Confini di responsabilità definiti; separazione dei moduli incompleta. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

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


**Stato di avanzamento D17.2 — IN CORSO (24 ottobre 2026).** Lookup nominale e tre capabilities introdotti in forma transitoria; restano ID fisici, sottotermini e metadati legacy. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

## D17.3 — Semantica dei metodi

L’insieme dei metodi di integrazione **non deve essere considerato definitivamente chiuso** a `EXPLICIT`, `NEWTON` e `ANALYTIC`. Questi sono i tre metodi implementati in D17, ma il meccanismo di selezione e validazione deve consentire l’aggiunta successiva di altri metodi numerici, inclusi eventuali integratori basati su modelli appresi. Questa estensibilità è un requisito architetturale, non richiede di implementare ulteriori metodi in D17.

### D17.3.1 `EXPLICIT`

`source_terms_2d` assembla il vettore di rate dei termini selezionati dal gruppo, `S_G^EX(q)`. A cadenza `STAGE`, il vettore entra nella parte esplicita del tableau RK. A cadenza `STEP`, se il trattamento è supportato, il gruppo viene eseguito dopo RK una volta sull'intervallo `dt`, nel proprio ordine di esecuzione. La valutazione e l'aggiornamento degli stati al bordo, quando presenti, restano separati dal rate di cella.

### D17.3.2 `NEWTON`

`NEWTON` indica **un unico trattamento numerico di gruppo**, che preserva la distinzione già esistente tra componenti **semi-implicite** (eventuale predictor/correzione) e **implicite** (residuo e Jacobiana). `source_terms_2d` assembla `S_G^SI` e `S_G^IM` e ne limita la selezione ai termini del gruppo. `nonlinear_solver_2d` raccoglie l'unione delle equazioni target del gruppo e costruisce **un unico vettore di incognite e un unico residuo congiunto**. La Jacobiana è dimensionata all'intero gruppo e comprende gli accoppiamenti realmente presenti, valutati lungo il medesimo stato delle incognite.

Per uno stadio diagonale RK, l'intervallo implicito pertinente è `a_ii^IM * dt`, come nel percorso esistente; la rappresentazione esatta del predictor SI e del residuo deve mantenere la semantica numerica validata, senza rifattorizzazioni matematiche arbitrarie. `NEWTON` a cadenza `STEP` è un obiettivo di configurabilità **solo se** viene specificato e validato un operatore di step coerente: non è attualmente consentito dal registro; non dichiararne il supporto prima dell'implementazione e della validazione. Nel frattempo l'input `STEP+NEWTON` deve fallire con diagnostica esplicita.

### D17.3.3 `ANALYTIC`

`source_terms_2d` offre l'operatore analitico `q_after = Phi_G(q_before, delta_t)` per **un solo problema analitico**. L'operatore restituisce lo stato aggiornato per il target, non un rate o un incremento arbitrario. A `STAGE` usa l'intervallo indicato dal tableau (nella configurazione attuale `a_ii^IM*dt`); a `STEP`, solo qualora quel particolare solve abbia una semantica e un'implementazione supportate, riceve `dt` e agisce dopo RK. Fino ad allora `STEP+ANALYTIC` deve essere rifiutato all'avvio. **Non** accumulare più solve analitici indipendenti nello stesso gruppo: usare gruppi distinti.


**Stato di avanzamento D17.3 — IN CORSO (24 ottobre 2026).** Metodi già presenti in parte; interfacce uniformi e rifiuto dei metodi non supportati da completare. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

## D17.4 — Cadenza e ordine

- `STAGE`: esecuzione secondo il tableau RK-IMEX, non un semplice aggiornamento Euler ripetuto. I contributi espliciti e impliciti vengono combinati secondo i coefficienti e i pesi degli stadi. La valutazione ripetuta del rate in Newton (residuo/Jacobiana/line search) non equivale a più applicazioni temporali della sorgente.
- `STEP`: esecuzione **una sola volta** dopo l'assemblaggio RK e prima dell'avanzamento dell'orologio. Un gruppo di trasferimento di massa `STEP/EXPLICIT` conserva l’atomicità e gli aggiornamenti coerenti di stato, inventari e terreno **attraverso il medesimo esecutore dei gruppi STEP**.
- L'ordine è **tra gruppi** compatibili con la fase di esecuzione. Nessuna numerazione di solve indipendenti dentro un gruppo. L'ordine `STAGE` deve rispettare la separazione dei contributi espliciti e impliciti del tableau: un unico ordinamento numerico non deve essere usato per pretendere di cambiare l'ordine matematico del metodo RK.
- Identità di gruppo: combinazione `(cadence, group_id)` oppure ID globalmente univoco, secondo il contratto di input scelto; usare una regola unica e documentata. L'ordine non costituisce l'identità del gruppo. Conservare il namespace stage/step indipendente senza ambiguità.


**Stato di avanzamento D17.4 — IN CORSO (24 ottobre 2026).** Cadenze e ordinamento esistono; esecutore STEP comune e revisione del scheduling ancora mancanti. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

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


**Stato di avanzamento D17.5 — IN CORSO (24 ottobre 2026).** Validazione di un solo solver per gruppo introdotta; descrittori non ancora semplificati e ancora basati anche su ID. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

## D17.6 — Gruppi `ANALYTIC`

- Un solo metodo: `ANALYTIC`; una sola cadenza; **un solo solve**.
- La soluzione analitica deve essere effettivamente disponibile per tutti i termini del gruppo considerati nel loro insieme, non supposta dalla sola presenza di un singolo termine analitico.
- Non creare automaticamente un solve analitico composito concatenando più operatori individuali.
- Se velocità, temperatura e pressione dei pori sono integrate analiticamente in modo indipendente, usare gruppi differenti (anche se il vecchio input le riferiva allo stesso `IMEX_group`).
- L'operatore aggiorna soltanto le equazioni target, lasciando invariato il restante stato conservativo, a meno di una trasformazione dichiarata e validata.


**Stato di avanzamento D17.6 — IN CORSO (24 ottobre 2026).** Un solve per gruppo analitico già imposto; adeguamento al contratto minimale e regressioni complete mancanti. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

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


**Stato di avanzamento D17.7 — IN CORSO (24 ottobre 2026).** Registro, input nativo ed esecutore impongono ora un solo solver Newton per gruppo (commit 079225a0, 209a6e0, db750a5b); accoppiamento multi-target e regressione successiva da verificare. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

## D17.8 — Modulo `source_terms_2d` e API di assemblaggio

**Conservare le quattro interfacce di assemblaggio**, alimentate in modo uniforme dalle corrispondenti quattro interfacce dei moduli fisici, adattandole a un *descrittore di gruppo* già validato e rimuovendo il passaggio di un selettore di sotto-solve distinto dal gruppo quando non più necessario:

- `eval_expl_terms(..., group)` → rate esplicito delle equazioni del gruppo;
- `eval_semi_implicit_terms(..., group)` → contributo SI di un gruppo Newton;
- `eval_implicit_terms(..., group)` → contributo IM dello stesso gruppo, in REAL e COMPLEX per residuo/Jacobiana;
- `eval_analytic_terms(group, ..., interval, ...)` → stato aggiornato dall'unico operatore analitico del gruppo.

L'elenco `term_names` seleziona i contributi e il modulo assembla i vettori sulle equazioni appropriate; non deve scegliere autonomamente la cadenza o il metodo. Non introdurre un dispatcher universale con tipo di output ambiguo: rate SI/IM/EX e aggiornamento analitico hanno contratti matematici diversi. **Regola generale senza eccezioni:** ogni contributo sorgente delle equazioni, **inclusi tutti i trasferimenti di massa**, deve essere assemblato in `source_terms_2d`. In particolare erosione, deposizione, perdite della fase continua, entrainment, vaporizzazione e contributi accoppiati a quantità di moto, energia, componenti e pressione dei pori non devono essere assemblati una seconda volta nell'esecutore temporale. L'assemblatore deve produrre un risultato coerente che comprenda sia il vettore delle sorgenti delle equazioni conservative sia i rate accessori necessari per gli inventari dei materiali e per la variazione del letto. Gli eventuali termini derivati vengono inclusi automaticamente nelle equazioni pertinenti, senza essere schedulati come sorgenti indipendenti. Tenere distinta la *valutazione/assemblaggio dei rate* dall'*applicazione atomica* dello step.


**Stato di avanzamento D17.8 — DA INIZIARE (24 ottobre 2026).** Rinomina in source_terms_2d e quattro assemblatori definitivi non realizzati. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

## D17.9 — Rate espliciti e sorgenti al bordo

**Base già implementata:** routing dei rate di cella `STAGE/EXPLICIT` filtrato dal gruppo e hook al bordo separato. Nel riallineamento non duplicare le sorgenti al bordo come rate di cella e non alterare l'operatore spaziale, la curvatura o le finestre di attivazione. Una sorgente di bordo produce uno stato imposto: la sua partecipazione alla configurazione deve rispettare tale rappresentazione. Conservare i test già validi e aggiungere verifiche mirate alla struttura semplificata.


**Stato di avanzamento D17.9 — IN CORSO (24 ottobre 2026).** Routing esplicito preesistente; alcune selezioni per nome introdotte; integrazione finale nei nuovi assemblatori da realizzare. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

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


**Stato di avanzamento D17.10 — IN CORSO (24 ottobre 2026).** Assemblaggio SI/IM e Newton indicizzato già esistenti; unico solve per gruppo applicato, ma copertura congiunta, correzioni e regressioni ancora incomplete. **Aggiornamento 24 ottobre:** commit `521f752c` allinea il solver alla lista completa `group_equations` e verifica che il singolo descrittore Newton copra esattamente quelle equazioni; le valutazioni SI/IM congiunte e i confronti numerici restano da completare. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

## D17.11 — Kernel REAL/COMPLEX

**Base già realizzata nel perimetro verificato:** riutilizzo del kernel COMPLEX della correlazione lahar dal wrapper REAL. Conservare la compatibilità complex-step e prevenire la duplicazione delle formule. Dove un termine selezionabile entra in Newton, verificare la disponibilità delle sue valutazioni necessarie al residuo e alla Jacobiana. Non imporre file `.inc` dove il riuso del kernel esistente è sufficiente.


**Stato di avanzamento D17.11 — IN CORSO (24 ottobre 2026).** Kernel REAL/COMPLEX preesistenti; compatibilità del nuovo contratto e Jacobiane congiunte da verificare. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

## D17.12 — Solve analitici

**Base già implementata:** interfaccia analitica in `source_terms_2d`, con target velocità, temperatura o pressione dei pori. Semplificare il contratto rimuovendo i riferimenti a un descrittore di sotto-solve indipendente quando la scelta è già univocamente determinata dal gruppo; mantenere l'intervallo fisico e l'aggiornamento del solo target. Rifiutare gruppi che richiederebbero due operatori analitici distinti.


**Stato di avanzamento D17.12 — IN CORSO (24 ottobre 2026).** Bridge analitico preesistente; adattamento alla struttura finale e regressioni da completare. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

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


**Stato di avanzamento D17.12a — DA INIZIARE (24 ottobre 2026).** thermal_cooling_flag indipendente e matrice delle combinazioni lava/non lava non implementati. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

## D17.13 — Default e compatibilità legacy

**Requisito obbligatorio:** in assenza di `SOURCE_GROUPS` il comportamento supportato deve coincidere con l'attuale baseline legacy (stessi termini attivi, trattamento SI/IM, ordine effettivo, finestre, rate, aggiornamenti e bilanci). Il costruttore legacy deve produrre gruppi che rispettino la nuova regola **un gruppo / un solve**:

- gruppi Newton legacy con target accoppiati → **un gruppo Newton con un solo sistema congiunto**;
- gruppi Newton legacy separati → gruppi Newton distinti;
- più target analitici legacy sotto una medesima label → **gruppi analitici distinti**, nell'ordine compatibile documentato;
- rate espliciti `STAGE` → gruppi espliciti secondo la baseline;
- batch di mass-exchange `STEP` → gruppo esplicito unico, conservando atomicità e posizione temporale.

La tabella dei default deve includere la nuova attivazione termica indipendente e il comportamento in assenza del relativo flag secondo D17.12a. Deve indicare condizioni di attivazione, metodo, cadenza, gruppo, equazioni target e precedenze per ogni combinazione fisica supportata. Confrontare sia il registro normalizzato sia l'effettiva evoluzione numerica. Non riprodurre una struttura di sotto-solve solo per imitare la nomenclatura legacy.


**Stato di avanzamento D17.13 — IN CORSO (24 ottobre 2026).** Default legacy preesistenti e test metadata superati prima degli ultimi cambiamenti; equivalenza numerica end-to-end ancora da dimostrare. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

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


**Stato di avanzamento D17.14 — IN CORSO (24 ottobre 2026).** Batch mass-exchange legacy presente; moduli per processo, interfaccia comune e rimozione dell'orchestratore dedicato non realizzati. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

## D17.15 — Sequenza del timestep

Sequenza production da preservare:

1. Il driver determina il `dt` effettivo applicando i vincoli esistenti.
2. Calcolo delle medie/fattori delle finestre temporali delle sorgenti.
3. Preparazione OU condizionale.
4. Avanzamento IMEX-RK con i gruppi `STAGE`.
5. Esecuzione una sola volta dei gruppi `STEP` in ordine, con il gruppo `EXPLICIT` mass-exchange **ultimo** nella configurazione legacy, senza percorso temporale speciale.
6. Avanzamento `runtime%t`.

È necessario testare le sorgenti circolari, fissurali, radiali e al bordo, inclusi periodi inattivi, rampe, finestre parziali e wrap, nonché `dt`, output e restart nei casi pertinenti. La verifica già svolta delle cinque fasi source-window/OU/RK/step/clock rimane evidenza storica nel suo perimetro.


**Stato di avanzamento D17.15 — IN CORSO (24 ottobre 2026).** Sequenza RK/STEP legacy esistente; STEP comune definitivo e regressioni di ordine ancora da fare. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

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


**Stato di avanzamento D17.15a — DA INIZIARE (24 ottobre 2026).** Spostamento sistematico delle operazioni sullo stato da time_integration_2d a state_2d non eseguito. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

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


**Stato di avanzamento D17.16 — IN CORSO (24 ottobre 2026).** Validazione preventiva del numero unico di solver nativi aggiunta; catalogo nominale, conflitti e contratti restano da consolidare. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

## D17.17 — Separazione fisica/numerica

Conservare l'isolamento dei moduli `*_laws_2d` rispetto a RK, gruppi e infrastruttura spaziale; l'assemblaggio dei contributi nelle equazioni rimane in `source_terms_2d`. La selezione fisica deve essere guidata dall'elenco termini del **gruppo**, non da identificatori di sotto-solve non necessari. Il solver legge residuo/rate dal bridge equazionale, non importa direttamente leggi fisiche. Completare audit statico e call graph, distinguendo le API di sorgente dalle API fisiche non-sorgente dell'operatore spaziale.


**Stato di avanzamento D17.17 — IN CORSO (24 ottobre 2026).** Moduli di leggi già separati in parte; separazione finale dei quattro assemblatori e responsabilità non completata. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

## D17.18 — Adattamento delle routine esistenti

Minimizzare le modifiche alle formule: **uniformare tutti i moduli di leggi sorgente alle quattro interfacce** e adattare la selezione nei quattro assemblatori (`eval_expl_terms`, `eval_semi_implicit_terms`, `eval_implicit_terms`, `eval_analytic_terms`) alla nuova struttura semplice. **Integrare il mass-exchange in questi stessi assemblatori:** la sola interfaccia di legge attiva è ora `explicit`, che alimenta `eval_expl_terms`; le altre tre interfacce di `mass_exchange_laws_2d` rifiutano esplicitamente le chiamate. **Non conservare `eval_mass_exchange_terms` come percorso di integrazione distinto**, né `mass_exchange_step_2d` come orchestratore speciale: riutilizzare soltanto le operazioni di aggiornamento necessarie nell’esecutore ordinario dei gruppi `STEP`, senza spezzare il batch. Rimuovere eventuali chiamate duplicate, rispettare finestre e maschere, nessuna integrazione doppia fra `STAGE` e `STEP`, nessun cambiamento involontario di segno, conservazione, stato o Jacobiana.

**Riorganizzazione del modulo:** prima classificare tutte le routine pubbliche dell’attuale `equation_terms_2d` come *sorgente* o *non-sorgente*. Spostare soltanto le seconde verso il proprietario appropriato, mantenendo signature e comportamento ove possibile; rinominare quindi modulo e file in `source_terms_2d` e aggiornare import, compilazione, test e documentazione. Separare i commit/logical steps di spostamento/rinomina da quelli che modificano il comportamento dei gruppi, così da rendere le regressioni attribuibili. Non creare un ulteriore dispatcher generico.


**Stato di avanzamento D17.18 — IN CORSO (24 ottobre 2026).** Adattamenti parziali e vincoli documentati; uniformazione di tutte le leggi, rinomina e integrazione mass-exchange da fare. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

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


**Stato di avanzamento D17.18a — DA INIZIARE (24 ottobre 2026).** Procedura vincolante specificata ma prova rheology_newmodel end-to-end non implementata. **Verifica:** le precedenti CI di compilazione e del registro sono passate, ma le modifiche recenti ai gruppi Newton (commit 079225a0, 209a6e0, db750a5b) non sono ancora state sottoposte a nuovi test; non dichiarare questa sezione completata senza i gate pertinenti di D17.20–D17.21.

## D17.19 — Piano implementativo aggiornato (ordine consigliato)

| Passo | Azione | Tipo | Accettazione mirata |
|---|---|---|---|
| 1 | Congelare inventario termini e baseline dei casi legacy; registrare test storici e SHA | Conservare | Catalogo e comparatori disponibili |
| 2 | Formalizzare contratto **un gruppo/un metodo/una cadenza/un solve**; eliminare la semantica dei multi-solve Newton | **Correggere** | Specifica, test negativi per multi-solve |
| 3 | Semplificare `source_group_type` e il registro, con un unico insieme di target e nessun scheduling interno | **Correggere** | Build, test di registry e immutabilità |
| 4 | Semplificare il parser `SOURCE_GROUPS` e adapter legacy; validare subito dopo input | **Correggere** | Input validi, conflitti e messaggi negativi |
| 5 | Adattare l'assemblaggio SI/IM a **tutti i termini del gruppo** | **Correggere** | Test dei rate per gruppo, REAL/COMPLEX |
| 6 | **Riutilizzare il Newton indicizzato già presente nella dev** (`group_equations(:)`), senza riscriverlo: una sola chiamata Newton per gruppo con residuo/Jacobiana congiunti; preservare le **due correzioni numeriche** dopo SI e dopo Newton, riconoscendo i blocchi tramite indici di equazione (D17.10.1) | **Correggere** | Test sistemi 2x2, 1x1, 3x3, accoppiamento e correzioni nei due punti |
| 7 | Adattare bridge analitico all'operatore unico di gruppo | Adeguare | Un solo analytic; confronto legacy |
| 8 | Preservare routing esplicito `STAGE` e hook al bordo | Conservare/adeguare | Rate e BC senza duplicazione |
| 8a | Separare il raffreddamento dalla selezione reologica: aggiungere `thermal_cooling_flag` in `TEMPERATURE_PARAMETERS`, migrare i controlli su layout, sorgenti, coefficienti e kernel, senza cambiare la legge fisica; preservare compatibilità dei vecchi input quando il flag è assente | **Correggere** | Matrice 2x2 reologia lava/non lava e raffreddamento on/off, parser e bilanci energetici |
| 9 | Uniformare **tutte** le famiglie di leggi sorgente al contratto fisico `explicit`/`semi_implicit`/`implicit`/`analytic`, conservando i kernel esistenti; registrare per ciascuna sorgente **un solo nome, un target, tre capabilities** (`EXPLICIT`/`NEWTON`/`ANALYTIC`), con SI+IM obbligatorie per NEWTON | **Correggere** | Test per ogni famiglia: dispatch coerente, disponibilità e rifiuto di interfacce non implementate |
| 10 | Riorganizzare `src/physics/mass_exchange/` per processi, accorpare deposizione + perdita indotta, migrare i kernel esistenti e introdurre il coordinamento interno in `mass_exchange_laws_2d` senza duplicare i limitatori condivisi; integrare la sola interfaccia `explicit` in `source_terms_2d::eval_expl_terms`, eliminare il routing temporale separato e preservare il batch legacy come **ultimo gruppo `STEP/EXPLICIT`** | **Correggere** | Bilanci, test dei singoli processi, limitatore comune e dipendenze, conservazione dell’inventario, nessun doppio aggiornamento, corretto ordine STEP e build aggiornato |
| 10a | Estrarre da `time_integration_2d` i controlli/correzioni fisico-specifiche in `state_2d` con nomi neutrali (`finalize_dry_cell_state`, controlli su masse/temperatura, roundoff), riusando `state_conversion_2d` e conservando sequenza/soglie; centralizzare anche l'inizializzazione delle primitive asciutte nell'halo | **Correggere** | Audit delle chiamate e test di equivalenza sugli stati intermedi/finali, celle asciutte, errori, observer e halo di ricostruzione |
| 10b | Migrare `correct_friction_momentum` in `state_2d` come correzione della direzione della quantità di moto a due stati; conservarne le due chiamate dal solver dopo SI e dopo Newton, tramite metadati di indici | **Correggere** | Test before/after di entrambi i punti, SI attivo/inattivo, sistema accoppiato e assenza di chiamate automatiche in finalizzazione generica |
| 10c | Migrare catalogo, gruppi, input e dispatcher da ID numerici fisici a nomi canonici; accorpare le componenti delle reologie sotto un solo termine configurabile e rendere il validatore indipendente dai modelli | **Correggere** | Esempio `rheology_newmodel` D17.18a, test compatibilità e assenza di ID fisici |
| 11 | Eliminare verifiche di configurazione ridondanti dagli hot loop, non i controlli numerici | Semplificare | Audit e prestazioni senza regressioni |
| 12 | Eseguire regressioni proporzionate al rischio su legacy/nativo, finestre, RK, OU, output, restart, ifx/CI pertinenti | Verificare | Evidenze e baseline approvate |
| 13 | Chiudere D17 solo quando D17.20–D17.21 sono soddisfatti | Gate | Decisione documentata |

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
- **Stato complessivo:** **IN CORSO / NON COMPLETATO**. «Commit pubblicato», «compilazione superata», «test del registro superati» e «equivalenza numerica dimostrata» sono stati distinti.
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
