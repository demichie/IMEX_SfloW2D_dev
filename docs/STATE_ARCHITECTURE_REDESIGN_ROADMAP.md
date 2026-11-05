# IMEX_SfloW2D — Roadmap di riprogettazione del sistema degli stati

**Specifica tecnica operativa** · 5 novembre 2026 · Repository `demichie/IMEX_SfloW2D_dev` · Branch `master`.

> **Ambito:** progettazione clean-slate della gestione e implementazione degli stati e del loro utilizzo in numerica e fisica. Non è una riprogettazione dell’intero codice: leggi fisiche, discretizzazioni, integratori e consumatori vengono modificati solo nella misura necessaria per usare i nuovi contratti di stato. La destinazione deve essere ricavata dai flussi di dati, dai costi computazionali e dagli invarianti numerici, non dalla suddivisione attuale dei file o dalla compatibilità con organizzazioni storiche. Il codice esistente è una fonte di formule, invarianti e casi da preservare, non il vincolo architetturale di partenza. Questo documento è una roadmap tecnica specialistica subordinata alla [roadmap generale](IMEX_SfloW2D_Development_Roadmap.md): non introduce un piano generale concorrente. Aggiornare checklist, prove e rischi qui; riportare nella roadmap generale solo stato aggregato e priorità.

> **Stato:** NUCLEO `state_new` IMPLEMENTATO; VALIDAZIONE P3 IN CORSO; MIGRAZIONE F5 AVVIATA. Il nucleo parallelo contiene tipi REAL/COMPLEX, workspace, valutazioni minime e complete della miscela, conversioni conservative/primitive, costruttori da `q`/`qp`, aggiornamento cinematico e refresh coerente dopo mutazioni termodinamiche o composizionali. L’integrazione con i consumatori di produzione e la validazione scientifica restano fasi successive.

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
- [~] **F5 — Migrazione dei sorgenti.** Aggiornare reologie, mass exchange, Darcy, termica, iniezioni e dispatcher Lie/Newton all'API definitiva, senza wrapper legacy. Iniezioni: valutazione proprietà del materiale tramite `state_new` già introdotta; allineati i rate di solidi e gas alle frazioni chiuse usate dalla miscela (bottom circular, fissural, radial). La migrazione complessiva resta in corso: mancano test di regressione sui rate in presenza di normalizzazione, overfill, e conversione degli altri sorgenti/dispatcher.

**F5 thermal analytic (5 novembre 2026):** il kernel esatto di raffreddamento riceve ora direttamente la temperatura iniziale (unico dato dello stato consumato). Il dispatcher avanza lo scalare fra modelli senza `current_state=state` e senza copia profonda degli array di composizione. Aggiornato `TEST_RHEOLOGY_CELL_LOCAL/test_lava_cooling.f90` per il nuovo contratto. Restano da eseguire i test CI e la migrazione dei tipi termici REAL/COMPLEX e dei chiamanti Lie/Newton alla nuova API; F5 non e conclusa.
- [ ] **F6 — Proprietà delle facce HP-PCCU.** Riutilizzare proprietà della traccia finale per flussi, velocità caratteristiche, pressione idrostatica e conversioni; mantenere semantica wet/dry e well-balanced.
- [ ] **F7 — Eliminazione delle API precedenti.** Eliminare routine non usate, calcoli duplicati, import, helper e documentazione obsoleti; ispezione finale delle dipendenze senza cicli.
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
| F5 | IN CORSO | Audit dipendenze 5 novembre 2026; preparazione dei contratti `state_new` e test cinematica/Richardson | Migrare reologie, processi mass exchange, iniezioni, termica, Darcy, dispatcher sorgenti e Lie/Newton; validare |
| F6 | DA FARE | — | Profilazione e percorso di faccia da consolidare |
| F7 | DA FARE | — | Rimozione API e dipendenze precedenti |
| F8 | DA FARE | — | Prove scientifiche e Actions sul commit finale |
| F9 | DA FARE | — | Benchmark e chiusura documentata |


**Regola di esecuzione:** l'ordine F0–F9 è una sequenza logica di responsabilità; F2/F3/F4 possono essere sviluppate in un unico change set consistente e validato alla sua conclusione. Non introdurre wrapper per mantenere compilabili i commit intermedi.

### Registro di migrazione F5 — primo change set (5 novembre 2026)

**Contratti consolidati in questo change set:**
- [x] Il setup valida le dimensioni di **tutti e quattro** gli array `xs`, `xg`, `alphas`, `alphag`, sia REAL sia COMPLEX.
- [x] `update_physical_state_kinematics` aggiorna anche Richardson quando cambiano `u,v`; i costruttori e `refresh_physical_state` non ne duplicano il calcolo.
- [x] `TEST_STATE_NEW` include un test dell'aggiornamento di velocita e della parita REAL/COMPLEX del Richardson.

**F5 — cambio contratti delle reologie (5 novembre 2026):** tutti i 12 moduli reologici (dispatcher e leggi) ora dichiarano direttamente `physical_state_real`/`physical_state_complex` da `state_new_types`. La promozione REAL → COMPLEX usata dalle valutazioni esplicite composte è un'operazione numerica esplicita in `state_new_types`, con buffer COMPLEX preinizializzato: non è un forwarding wrapper. Il vecchio costruttore del `state_conversion_2d` non viene più usato dalle reologie. La nuova promozione richiede aggiornamento dei chiamanti e della gestione degli scratch; compilazione e test restano NON VERIFICATI fino alla migrazione coordinata di `source_terms_2d` e dell'esecutore Lie/Newton. Non marcare F5 completata.

**Inventario delle dipendenze ancora da migrare (audit statico su master `1ca1fe6`):**
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
