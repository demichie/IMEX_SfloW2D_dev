# IMEX-SfloW2D — Architettura e roadmap two-layer

**Stato del piano:** 15 novembre 2026 — branch `feature/two-layer`.  
**Ambito:** evoluzione dal solutore single-layer esistente a due layer fisici (L1 corrente inferiore densa, L2 nube superiore diluita). **La produzione resta single-layer: `N_LAYERS=2` è disabilitato fino al gate finale.**  
Questo documento è l'unica roadmap operativa vigente. Descrive l'architettura obiettivo, le decisioni approvate, lo stato effettivo e le verifiche ancora richieste; non è un registro dei commit. Per i dettagli tecnici di audit consultare [TWO_LAYER_CODE_AUDIT.md](TWO_LAYER_CODE_AUDIT.md).

## 1. Obiettivo e regole vincolanti

### 1.1 Semplicità e proprietà dei dati

1. **Non creare sovrastrutture se esistono già gli strumenti necessari.** Prima di aggiungere moduli, tipi, wrapper, maschere, dispatcher o interfacce, verificare e riutilizzare i proprietari e gli operatori correnti.
2. **Scegliere sempre la soluzione più semplice, diretta e definitiva.** Non introdurre generalizzazioni speculative, percorsi transitori, API ridondanti o due implementazioni dello stesso algoritmo.
3. **Una sola rappresentazione autorevole** per ogni quantità, configurazione, componente di stato e inventario. Vietate copie persistenti sincronizzate, coefficienti fisici duplicati e attivazioni globali parallele a quelle per layer.
4. **Concludere le migrazioni.** Eliminare interfacce interne e strutture sostituite, aggiornando tutti i chiamanti; mantenere soltanto le compatibilità di input espressamente concordate.
5. **Responsabilità chiare:** `physics` implementa leggi e rate fisici; `spatial` flussi e prodotti non conservativi; `numerics` geometria dei candidati, scheduling, contatti al confine numerico, integrazione Lie e rollback. Il Lie non contiene formule fisiche.
6. **Kernel cell-local:** le leggi di `physics`, in particolare `mass_exchange`, ricevono stati, configurazioni ed esposizioni locali, senza indici di griglia né accesso allo stato globale. Non ricostruiscono contatti già determinati all'esterno.
7. **Contatti minimi:** l'unico classificatore esistente distingue quattro casi sufficienti: L1+L2 presenti, entrambi assenti, solo L1, solo L2. Riutilizzarlo direttamente; non creare classificatori, cache o API equivalenti.
8. **Nessuna doppia infrastruttura two-layer:** un solo layout di componenti, un solo HP-PCCU, un solo registro Lie e un solo assemblatore mass exchange. Esattamente due layer, non un framework generico a N layer.

### 1.2 Integrità scientifica

9. Preservare equazioni, conservazione, dimensioni, metodi, splitting, regolarizzazioni, limitatori, tolleranze, ordine di valutazione e risultati single-layer. Ogni modifica fisica o numerica significativa richiede motivazione e approvazione esplicita, mai adattamenti per far passare i test.
10. Distinguere rigorosamente scambi **interni** L1↔L2, con contributi opposti e budget di coppia, dagli scambi **esterni** con atmosfera e substrato; non duplicare inventari condivisi.
11. Mantenere coerenza REAL/COMPLEX, ammissibilità, stati dry/wet, sorgenti e metodi numerici realmente supportati, rollback e riproducibilità OpenMP. Nessun fallback silenzioso fra EXPLICIT, SEMI_IMPLICIT, NEWTON/IMPLICIT e ANALYTIC.
12. `N_LAYERS=2` resta rifiutato da parser e layout di produzione finché implementazione, test accoppiati e approvazione finale non sono completi. Le prove L2 isolate non abilitano il modello.

### 1.3 Metodo di sviluppo e verifica

13. **Audit completo prima delle modifiche:** branch e revisione corrente, proprietari, dichiarazioni, `USE`, chiamanti, `.inc` REAL/COMPLEX, parser, diagnostiche, I/O/restart, TESTS, EXAMPLES e dipendenze Makefile. Non usare sostituzioni cieche.
14. Modificare gruppi di dipendenze completi, evitando patch correttive ripetute e modifiche estranee. Risolvere la causa degli errori, cercare problemi analoghi ed eliminare il codice superato.
15. Applicare `docs/FORTRAN_STYLE.md`: quattro spazi per livello lessicale, massimo 120 colonne, nessun tab o whitespace finale, direttive OpenMP e Doxygen corretti.
16. Verificare ciò che è stato effettivamente eseguito: compilazioni optimized e strict/debug; regressioni numeriche e fisiche, REAL/COMPLEX, conservazione, esempi, restart/NetCDF, OpenMP e prestazioni quando pertinenti. Un build riuscito non certifica la fisica; evidenze di commit precedenti non certificano modifiche successive.
17. Proseguire sul branch concordato `feature/two-layer`. Non creare branch né effettuare merge su `master` senza autorizzazione. Documentare soltanto lo stato verificato, distinguendo codice implementato e gate scientifici ancora aperti.

## 2. Architettura fisica e informatica

### Stato e proprietà

- `q/qp(componente,layer,x,y)` e workspace di integrazione hanno struttura rank-4. **La dimensione componente e gli indici sono comuni** a L1/L2; le variabili opzionali stocastica e Darcy occupano ciascuna un solo slot se attive in almeno un layer (OR dei flag canonici). Uno slot non attivo in un layer rimane esattamente zero in stato, trasporto, sorgenti, confini e restart.
- Le configurazioni dei processi appartengono a `L1` o `L2`; topografia `B`, geometria del terreno, substrato, `deposit_tot`, specie/proprietà solide, atmosfera, temperatura `T_ambient` e fase portante gas/liquido/lava sono **SHARED**. Non duplicare proprietà di materiale o inventari.
- Le attivazioni di reologia stocastica, trasporto conservativo stocastico e pressione interstiziale sono distinte e canoniche **per layer**. Il processo OU può essere necessario anche senza trasporto stocastico conservativo; il sorgente reologico dipende solo dal relativo flag.
- I confini laterali W/E/S/N possiedono ciascuno un array canonico `bc(component,layer)`; il workspace spaziale seleziona il layer senza duplicare i kernel. Per `N_LAYERS=1` i gruppi `*_BOUNDARY_CONDITIONS` e `*_L1_BOUNDARY_CONDITIONS` sono alternative mutuamente esclusive. Le NAMELIST `_L2_BOUNDARY_CONDITIONS` sono predisposte e riutilizzano lo stesso parsing; gli slot opzionali inattivi ricevono BC Neumann omogenee. `N_LAYERS=2` rimane disabilitato; test e integrazione accoppiata restano aperti.
- Nel single-layer, per le famiglie previste, `<FAMILY>_PARAMETERS` e `<FAMILY>_L1_PARAMETERS` sono nomi alternativi mutuamente esclusivi dello **stesso oggetto L1**; nessuna copia. Quando L2 sarà abilitato, richiedere gli input espliciti L1/L2 e rifiutare il nome breve. I nomi scalari obsoleti non sono interfacce interne da preservare.
- Il substrato possiede l'unico `T_ground`, il deposito canonico è `deposit_tot`; `T_ambient` è la sola temperatura atmosferica. Le proprietà solide sono condivise, ma il massimo packing di ammissibilità non va confuso con hindrance di deposizione o futuro packing L1.

### Geometria e operatore spaziale

- Un solo operatore HP-PCCU riutilizzato con workspace e domini wet/dry indipendenti per layer. Le superfici **derivate**, non immagazzinate, sono `eta_1=B+h1` e `eta_tot=B+h1+h2`.
- Per la correzione di pendenza forte, la metrica condivisa è `G=1/(1+B_x^2+B_y^2)`; gradienti e Hessiana provengono sempre da **B**, non da `B+h1`. La cinematica e la curvatura usano le velocità correnti del rispettivo layer.
- Nei flussi idrostatici e nella stima delle celerità entra `max(gprime_i,0)`, **senza curvatura**; la correzione di curvatura appartiene alla sorgente di pendio e, ove previsto, allo sforzo normale efficace. Preservare i termini metrici esistenti.
- Per `L=rho_m2*gprime2*h2`, le pressioni non conservative accoppiate seguono il manoscritto: `S_inter1=-G*h1*d_j(L)` e `S_inter2=-L*d_j(G*h1)`. Il fattore `G` è esterno alla derivata per L1 e interno per L2. Senza clamp la somma è `-d_j(G*h1*L)`; con clamp di buoyancy verificare e diagnosticare il difetto di bilancio, non correggerlo artificiosamente.
- Riusare l'API di quadratura/path di `nonconservative_2d`. Le velocità d'onda separate dei due layer **non** forniscono da sole il CFL accoppiato: derivare/verificare bound, iperbolicità, shear, dry/wet e clamp. Rispettare il residuo ammissibile su topografia curva indicato nel manoscritto, senza promettere well-balancing esatto non previsto.

### Sorgenti, contatti e Lie

- Un unico Lie del primo ordine con timestep globale, ordine dei gruppi, cell-outer/OpenMP, gruppo 0, candidati e rollback. Ogni gruppo avrà scope `L1`, `L2` o `INTERLAYER`. I gruppi locali chiamano gli stessi kernel monolayer con il solo stato bersaglio.
- All'ingresso **di ogni gruppo**, usare i candidati L1/L2 correnti: applicare l'attuale rimozione buoyancy L2 quando pertinente; poi classificare direttamente i quattro contatti **fuori da `physics`**. Le esposizioni terreno/atmosfera, già risolte, arrivano ai kernel locali di scambio. Nessuna cache di contatto; la scomparsa di L2 non interrompe le sorgenti L1.
- Un gruppo `INTERLAYER` valuta un solo rate `Gamma` dallo stesso snapshot per entrambi i layer, convenzione positiva L2→L1: `q1'=q1+dt*Gamma`, `q2'=q2-dt*Gamma`. Disponibilità/limiter coerenti e validazione/rollback atomici di entrambi gli stati e delle risorse condivise; conservare massa per specie, momento e variabile termica secondo il contratto fisico, distinguendo dissipazione cinetica e perdite esterne.
- Scambi con substrato e atmosfera dipendono dal contatto della **cella**: L2 solo può scambiare col terreno, L1 solo con l'atmosfera. Nessun doppio conteggio dei depositi o delle altre risorse SHARED. Cold-start L2 e comparsa/scomparsa di layer richiedono test dedicati.

**Riferimento scientifico vincolante:** *IMEX_SfloW2D_v3*, in particolare geometria delle grandi pendenze, pressione di coppia, chiusure fisiche, path-conservative e verifiche stratificate. Non introdurre nuove leggi termiche o combinazioni lava/granulare two-layer senza decisione scientifica separata.

## 3. Stato di avanzamento

La parola **completato** indica implementazione preparatoria disponibile, non certificazione dell'intero modello. **Avviato** significa codice o test parziali e gate aperti. Build e regressioni sono state confermate sullo stato precedente; il nuovo controllo degli slot COMPLEX inattivi richiede ancora CI sul commit che lo introduce.

| Ambito | Stato | Risultato e limite |
|---|---|---|
| Stato rank-4 e migrazione L1 | **Completato (strutturale)** | Stato/workspace canonici e percorsi L1 adattati; regressioni integrate ancora da completare. |
| Deposito/substrato/proprietà SHARED | **Completato (strutturale)** | `deposit_tot`, `T_ground`, specie e proprietà condivise con proprietari unici; test estesi di inventario, lava e restart aperti. |
| Configurazioni reologiche, termiche, Darcy, stochastic/OU e mass exchange L1/L2 | **Completato (preparazione)** | Ownership per layer e selezione nei percorsi locali/REAL/COMPLEX predisposte; L2 non selezionabile in produzione. |
| Classificazione dei quattro contatti | **Completato (isolato)** | Riutilizzabile senza nuove strutture; aggiornamento per gruppo con candidati accoppiati da realizzare. |
| HP-PCCU e domini indipendenti per L1/L2 | **Completato (isolato)** | Clamp idrostatico L2, test di trasporto Euler CFL-limited isolato e bilanci dedicati; nessuna pressione accoppiata. |
| Conteggio OR degli slot opzionali e conversioni REAL/COMPLEX | **Avviato** | Dimensionamento condiviso e conversioni REAL/COMPLEX verificati in CI; ulteriore controllo degli zeri COMPLEX inattivi aggiunto. Il test PCCU verifica ora l'algebra dei flussi con attivazioni L1/L2 opposte (nuova CI da eseguire); restano da verificare ricostruzione completa, sorgenti, IC/BC e restart sul layout opzionale comune. |
| Certificazione single-layer completa | **Avviato** | Build/regressioni e vari EXAMPLES già verificati; matrice scientifica, I/O/restart, OpenMP e confronto quantitativo non interamente chiusi. |
| HP-PCCU accoppiato e CFL two-layer | **Da fare** | Pressione di coppia, path, wet/dry, topografia, stratified rest e iperbolicità. |
| Lie L1/L2/INTERLAYER | **Da fare** | Scope, candidato doppio, contatti aggiornati per gruppo, transazioni condivise e rollback atomico. |
| Chiusure interlayer e scambi esterni accoppiati | **Da fare** | Settling, elutriazione, drag, entrainment, deposito/erosione, lofting e bilanci. |
| Input, output e restart two-layer | **Avviato (BC)** | BC canoniche per layer e alias L1 predisposti; parser BC L2 predisposto ma non attivabile; input L2 completo, NetCDF, diagnostiche e restart ancora da realizzare. |
| Validazione integrata e attivazione `N_LAYERS=2` | **Da fare** | Benchmark accoppiati, prestazioni, conservazione e approvazione prima dell'abilitazione. |

## 4. Roadmap operativa e condizioni di uscita

| Fase | Stato | Lavoro necessario per chiuderla |
|---|---|---|
| **A — Fondazioni single-layer** (TL00–TL03) | **Avviata, preparazione sostanzialmente realizzata** | Chiudere baseline riproducibile, matrice deposito/lava/dry, EXAMPLES, restart/NetCDF, verifiche quantitative e OpenMP; assicurare equivalenza L1 salvo delta termici approvati. |
| **B — Proprietà e fisica per layer** (TL04) | **Avviata, implementazione principale completata** | Certificare parsing/alias, tutte le configurazioni e combinazioni opzionali, budget dei processi e regressioni; nessuna attivazione L2 in produzione. |
| **C — Operatori L2 indipendenti** (TL05) | **Avviata, test isolati realizzati** | BC canoniche per layer, alias L1 e gruppi L2 predisposti; BC inattive omogenee. I test di ricostruzione HP, tracce/flussi/path/residui con slot opzionali opposti sono verificati in CI. Aggiunti test di stato iniziale asciutto su due layer e sorgenti geometriche locali (CI da verificare). I/O e restart accoppiati appartengono a TL12. |
| **D — Spaziale accoppiato** (TL06) | **Avviato** | Esteso `coupled_characteristic_audit.py` con una scansione riproducibile di stati localmente congelati: densita' uguali/differenti, contrasti forti, h1/h2, coefficiente metrico G e shear. Test distinti per autovalori complessi e per contenimento delle sole velocita' reali; le finestre di non iperbolicita' sono diagnosticate sulla griglia di shear, non spacciate per soglie analitiche (CI pendente). Non testa gradienti di densita'/composizione, closure variabili, wet/dry o curvature; `N_LAYERS=2` resta disabilitato. |
| **E — Lie accoppiato** (TL07) | **Da fare** | Scope dei gruppi, doppio candidato, gruppo 0, refresh dei contatti, buoyancy L2 e transazioni/rollback atomici, mantenendo cell-outer/OpenMP. |
| **F — Chiusure fisiche** (TL08–TL11) | **Da fare** | Settling L2→L1; elutriazione L1→L2 con cold-start; drag di coppia; scambi esterni e basali per layer esposto. Verificare risorse SHARED, disponibilità, bilanci e limiti. |
| **G — Integrazione e rilascio** (TL12–TL13) | **Da fare** | Input L2, diagnostiche, NetCDF/restart, benchmark scientifici, riproducibilità, prestazioni e regressioni single-layer; soltanto dopo approvazione abilitare `N_LAYERS=2`. |

**Ordine delle dipendenze:** verificare B/C e la baseline A mentre si prepara D; realizzare E prima degli scambi F che usano candidati e contatti aggiornati; certificare G solo dopo D–F. I gate A–C possono richiedere regressioni aggiuntive senza bloccare un lavoro preparatorio indipendente, ma **nessun gate può essere dichiarato superato senza le sue evidenze**.

## 5. Criteri di accettazione e rischi da controllare

- **Single-layer:** nessuna regressione non autorizzata in formule, risultati, ordine Lie, bagnato/asciutto, gestione del deposito e OpenMP; input brevi e `_L1` supportati ove previsti, senza duplicazioni interne.
- **Stato e slot:** per L1/L2 configurati in modo opposto, unici indici comuni, componenti inattive esattamente nulle lungo conversioni, flussi, sorgenti, IC/BC, stochastic, I/O e restart; controlli REAL/COMPLEX.
- **Spaziale:** test su terreno piano, a pendenza costante e curvo, diversi campi di velocità, solo L1, solo L2, stratified rest, clamp, identità della coppia non clampata, errori ammessi su terreno curvo, conservazione e CFL stabile.
- **Contatti/scambi:** coprire tutte le quattro combinazioni wet/dry a ogni gruppo, entrambi i confini esposti nei casi monolayer, comparsa/scomparsa, disponibilità del substrato, assenza di doppio conteggio e rollback.
- **Accoppiamento:** trasferimenti opposti per specie/momento/termica dove fisicamente conservativi, dissipazione ed eventuali difetti esterni diagnosticati; verificare modalità numeriche realmente supportate senza fallback.
- **I/O e prestazioni:** reader/writer/restart coerenti con un solo schema, test di continuazione, EXAMPLES rappresentativi, build debug/optimized, runtime e riproducibilità OpenMP.
- **Decisioni scientifiche ancora aperte:** chiusure settling/elutriazione/drag, pressioni e regolarità per dry/shear, uso Darcy in L2, packing/consolidazione, trattamento dei difetti del clamp, dettagli dei metodi impliciti e combinazioni di carrier. Non fissarle tramite convenienze implementative.

**Definizione di completamento:** l'implementazione è coerente lungo tutta la catena delle dipendenze, non restano strutture obsolete o rappresentazioni duplicate, le prove richieste sono effettivamente superate e il comportamento fisico autorizzato è preservato. Nessun merge o attivazione automatica segue dalla sola chiusura tecnica di una fase.
