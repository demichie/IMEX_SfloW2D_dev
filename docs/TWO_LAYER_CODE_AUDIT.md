# IMEX-SfloW2D — Audit tecnico di fattibilità two-layer

**Data di riferimento:** 7 novembre 2026  
**Base ispezionata:** `master` al commit `5a3778398c9ffc0452bc54c680d1583a1d5667d9` e roadmap consolidata nella PR #6.  
**Metodo:** audit statico dei sorgenti e delle interfacce connesse, dei test e dei workflow, confrontati con la specifica del manoscritto. **Non sono stati compilati né eseguiti benchmark o test dinamici in questo audit.**  
**Specifica operativa:** [Architettura e roadmap two-layer](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md).

## Esito

La proposta è **realizzabile come estensione dell'architettura attuale**, preservando il motore path-conservative e l'integratore Lie globale. Non è un'aggiunta limitata all'indice layer: richiede estensione esplicita della configurazione fisica per layer, del registro Lie, della ricostruzione accoppiata, delle leggi fisiche di L1/L2 e dei formati di output/restart. La roadmap identifica i gate che devono precedere ogni attivazione di `N_LAYERS=2`.

La fattibilità riguarda la **struttura del software**. Stabilità, iperbolicità e riproduzione dei benchmark non sono dimostrate dall'audit statico e restano condizioni scientifiche bloccanti.

## A. Risultati dell'ispezione dei sorgenti

| ID | Priorità | Riscontro verificato | Conseguenza e criterio di soluzione |
|---|---|---|---|
| **A01** | Alta | `src/core/model_layout_2d.f90` contiene `n_layers` ma `supported_layer_count` e `src/io/inpout_2d.f90` rifiutano valori diversi da 1. Gli stati sono rank-3; i workspace in `time_integration_2d`, `lie_source_executor_2d`, ricostruzione e operatori sono anch'essi monolayer. | Migrazione completa rank-4 di dati e consumer, con gate single-layer prima di abilitare L2. Non basta cambiare il controllo input. |
| **A02** | **Critica — contratto deposito** | `src/state/state_fields.f90:refresh_substrate_totals` assegna zero a `deposit_tot` quando `n_solid==0` e lo ricalcola solo dagli inventari granulari. La lava usa `lava_deposit` indipendente; `lie_source_executor_2d` lo aggiorna e chiama il refresh del totale. | `deposit_tot` canonico necessita dispatch esclusivo granulare/lava, differente modalità di accumulo, rimozione della duplicazione dei candidati e test con `n_solid=0`. Non è sufficiente sommare due campi. |
| **A03** | Alta | `lie_source_executor_2d` accumula `bed_rate_cell` da `dB_dt`, proietta il campo una volta su `B_vertex` e aggiorna la geometria prima del refresh dei totali. Lo stato possiede snapshot degli inventari. | Distinguere inventario cumulativo, rate topografico e terreno erodibile; preservare proiezione/commit e rollback atomici. |
| **A04** | **Critica — configurazione layer** | `src/core/parameters_2d.f90` espone `rheology_model` e flag di fase globali; le leggi costitutive hanno molti parametri globali. `source_configuration_2d` costruisce una sola lista di sorgenti e `source_registry_2d` rifiuta più di una reologia, sorgenti duplicate e mass-exchange in gruppi distinti. | Contratti e contesti fisici read-only per layer, input L1/L2/coppia, unicità per identità+scope; passaggio dei parametri ai kernel senza modificare globali dentro OpenMP. |
| **A05** | Alta | `src/spatial/nonconservative_2d.f90` integra già callback generiche dell'integrando su un percorso esteso con quadratura Gauss 1/2/3 e path callback opzionale. `pccu_2d` usa la stessa API. | **Nessuna modifica alla API del motore**; nuove callback fisiche e assemblaggio accoppiato. |
| **A06** | **Critica — geometria dei path** | `pccu_2d:eval_hydrostatic_path` accetta uno stato esteso `(h,Gamma,eta,G)`; la ricostruzione attuale è quella di una superficie `eta=B+h`. | Nel two-layer utilizzare `eta_1=B+h1` e `eta_tot=B+h1+h2`, una geometria condivisa delle tracce e percorsi faccia/cella accoppiati. Verificare topografia, clamp, wet/dry e identità della somma. |
| **A07** | **Critica — CFL** | `src/physics/equation_terms_2d.f90:eval_local_speeds_direction` usa `u ± sqrt(G*gprime*h)` del layer isolato; `spatial_operator_2d` prende massimi di tali velocità dalle facce. | Il bound per i modi interlayer non segue automaticamente. Derivare e verificare un bound giustificato e comportamento in perdita di iperbolicità; gate prima dei benchmark coupled. |
| **A08** | Alta | `mass_exchange_deposition_2d` implementa sedimentazione per classi, con concentrazione massima e esponente hindered prefissati nel kernel, e seppellimento del carrier mediante `coeff_porosity`. Il manoscritto two-layer introduce packing L1, gas burial, modalità `d32` e condizione di stress efficace. | Riutilizzare velocità terminali e calcolo delle specie, **ma implementare e testare** le chiusure basali specifiche del manoscritto. Nessun presupposto di equivalenza fisica a codice invariato. |
| **A09** | Alta | `source_catalog_2d`, `source_registry_2d` e `lie_source_executor_2d` descrivono le sorgenti come un'unica sequenza globale di un solo stato; il gruppo 0 legge lo snapshot pre-spatial. | Introdurre scope L1/L2/INTERLAYER, rate di trasferimento unico, controllo simultaneo e rispetto del gruppo 0. Evitare a priori Newton accoppiato generalizzato. |
| **A10** | Alta | `src/io/inpout_2d.f90` usa NetCDF `(x,y,time)` senza layer, restart binario v7 senza `n_layers` nell'header, e salva `lava_deposit` in una registrazione finale separata; rifiuta versioni non correnti. | Definire schema two-layer, output single-layer compatibile, versioning e reader v7 esplicito o rifiuto dichiarato. Test di continuazione con deposito e topografia. |
| **A11** | Alta | `src/physics/source_configuration_2d.f90` valida la mutua esclusione fra `gas_flag`, `liquid_flag`, `lava_flag`. | Il two-layer L1/L2 non deve essere interpretato come combinazione di due flag di fase continua incompatibili. Condivisione delle specie e differenza delle reologie richiedono un contratto nuovo di configurazione senza infrangere le regole monolayer. |
| **A12** | Media | `src/physics/mass_exchange/mass_exchange_entrainment_2d.f90` è già presente e coordinato da `mass_exchange_2d`. | Riuso per atmosfera→L2 con gating della copertura L1, senza duplicare formule; contabilizzare in modo distinto l'eventuale loss/lofting. |
| **A13** | Media | I test esistenti includono `TEST_NONCONSERVATIVE_PATH`, `TEST_HYDROSTATIC_PATH`, `TEST_PCCU`, `TEST_MASS_EXCHANGE`, `TEST_SOURCE_REGISTRY`, `TEST_BINARY_RESTART` e test di stato/geometry. I workflow principali `build-only`, `physics-regressions` e `solver-validation` sono a esecuzione manuale. | Riutilizzarli per il baseline e aggiungere test coupled specifici. Non attribuire certificazione a workflow non eseguiti. |

## B. Interfacce e responsabilità da rendere esplicite

### B1. Stato / accesso per layer

- Il dato globale fluido ha layout `(component,layer,x,y)`; topografia, stato di deposito ed erodibile restano condivisi.
- Le API locali monolayer mantengono vettori contigui. `q`, `qp`, stato trial, temporanei, bounds, I/O e test vengono migrati coerentemente.
- `model_layout%n_physical_variables` e `n_variables` sono distinti per nome ma **attualmente uguali per contratto**. Qualunque divergenza futura richiede test e conversioni adeguate, non un'assunzione tacita.
- L1 e L2 condividono indici di variabili; slot non validi nel layer non attivo devono essere mascherati in tutte le trasformazioni e nel sistema sorgente.
- Il dominio attivo per cella deve considerare `wet(L1) OR wet(L2)`; ogni sweep per layer seleziona il sottoinsieme valido e conserva la ricostruzione/halo necessari a entrambi. Non ignorare celle `dry L1 / wet L2`.

### B2. Sorgenti / stato fisico

- Separare chiusure di legge `physics` e commit numerico `Lie`.
- Context immutable per layer, con parametri reologici e termici differenziabili; non cambiare globali durante l'iterazione OpenMP.
- Permettere ripetizione della stessa identità sorgente in layer differenti; applicare i vincoli sull'inventario alle risorse **effettivamente condivise**.
- Assegnare al gruppo INTERLAYER uno snapshot coerente di entrambi, `Gamma` unico signed, conservazione di classi, massa, momento e energia termica; processi esterni separati, compreso il gas packing approssimato.
- Stato pre-spatial e ordine speciale del gruppo 0: testare sia iniezione in celle inizialmente asciutte sia scambi successivi.
- Gate di buoyancy per L2 indipendente dalla rimozione di L1; eseguire i controlli prima di scartare componenti del donatore in modo da non perdere inventari.

### B3. Ricostruzione / HP-PCCU

- Geometria comune e tre quote `B`, `eta_1`, `eta_tot`; tracce wet/dry per layer coerenti.
- Coppia non conservativa `-G*h1*d(K2*h2)` / `-K2*h2*d(G*h1)` nello stesso percorso esteso; clamp applicato secondo il manoscritto.
- Somma in assenza di clamp uguale a `-[G*h1*K2*h2]_L^R`; testare integrando e **operatore intero** (flusso, pressione propria, topografia, path, residuo per cella). La riduzione single-layer richiede proprietà, velocità e appropriate chiusure coincidenti.
- Separare il residuo previsto su letto curvo dall'errore numerico dello schema.
- CFL con modi barotropici/baroclinici o bound verificato, con gestione di stati non iperbolici e shear.

### B4. Deposito / evoluzione del terreno

- Il deposito lavico è un accumulo diretto di spessore; quello granulare deriva dalle classi solide e dalla porosità valida.
- Non utilizzare `n_solid==0` come segnale di deposito nullo in lava.
- Il rate geom. `dB_dt` non deve ricostruire una seconda volta l'inventario. Validazione delle proiezioni sulle quote nodali e della differenza tra deposito cumulativo ed erosione.
- Ramo `topo_change_flag=.FALSE.`: mantenere l'inventario indipendentemente dalla modifica della topografia.
- Restart v7 e output legacy devono essere controllati in fase di razionalizzazione, prima della migrazione rank-4.

### B5. Chiusure scientifiche necessarie

- Settling L2→L1: trasferimento con `u_Gamma`, thin-plug blend, amp cap, donor temperature; packing gas distinto. No deposito terreno.
- L1→terreno: basal settling con `ALPHA_PACKING_L1`, gas burial, `SETTLING_D32_FLAG_L1`, `DEP_EFF_STRESS_FLAG_L1` secondo i casi scientifici dichiarati.
- Elutriazione: Rouse, Shields, shear-power, linear legacy, co-lift, cold-start ed inventory cap; specie e gas del donatore.
- Entramento atmosferico L2, buoyant lofting e opzioni di partizione co-ignimbrite. Testare il materiale che esce dal dominio e quello che resta mobile in L1.
- Drag interfaciale Chézy: puro scambio di momento; dissipazione energetica distinta.

## C. Test e gate minimi

| Gate | Procedura | Esito richiesto |
|---|---|---|
| **G0 — baseline** | Tag del commit che ha superato build + esempi completi, prova restart; raccolta tempi/memoria e campi finali. | Evidenza riproducibile; nessuna attribuzione di test non eseguiti. |
| **G1 — deposito** | Granulare e lava, `n_solid=0`, rate di massa/energia, `deposit_tot`, `dB_dt`, erosion, topografia on/off, restart v7. | Non regressione delle leggi; nessun azzeramento lava e nessun doppio conteggio. |
| **G2 — rank-4 mono** | Test layout, conversions, source groups, geometry, spatial, output, restart, OpenMP e EXAMPLES. | Equivalenza dei risultati e footprint performance entro tolleranze dichiarate. |
| **G3 — per-layer config** | Due reologie con parametri differenti e nomi sorgente duplicati ma scope diverso; OpenMP race check; configurazione monolayer legacy. | Configurazione deterministica e indipendente dal numero di thread. |
| **G4 — HP-PCCU** | Path pair antisimmmetrico nel budget, riduzione idrostatica totale con topografia, flat/inclined/curved, dry/wet, clamp e CFL. | Identità discretizzate testate; residuo curvo interpretato; stabilità accettabile. |
| **G5 — Lie pair** | Scambio fittizio di specie/momento/termica e rollback forzato durante modifica del terreno; test gruppo 0. | Trasferimenti conservativi strettamente interni e commit atomico. |
| **G6 — chiusure** | Settling, consolidazione L1, elutriazione, drag, co-lift, lofting, pressione interstiziale e combinazioni. | Bilanci numerici e contratti fisici confrontati al manoscritto. |
| **G7 — produzione** | Benchmark integrati, long run, restart, profiling OpenMP, regressione single-layer. | Approvazione scientifica prima di abilitare `N_LAYERS=2`. |

## D. Limiti di questa verifica

L'ispezione permette di confermare le dipendenze e i vincoli architetturali, ma **non** dimostra bitwise equivalence, conservazione discreta, well-balancing, correttezza del bound CFL, scaling OpenMP o raggiungimento dei benchmark. Queste proprietà richiedono build ed esecuzioni. La roadmap evita quindi di dichiarare completati i gate.

**Valutazione:** fattibilità strutturale **positiva con prerequisiti espliciti**; nessun ostacolo che richieda riscrivere il motore `nonconservative_2d` o abbandonare il Lie globale. Le priorità architetturali più importanti sono la semantica canonica del deposito, i contesti immutabili per layer, la ricostruzione della geometria accoppiata e la CFL del sistema two-layer.
