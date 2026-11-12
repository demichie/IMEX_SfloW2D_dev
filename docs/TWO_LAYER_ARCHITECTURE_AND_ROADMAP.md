# IMEX-SfloW2D — Architettura e roadmap two-layer

**Piano di riferimento:** 12 novembre 2026. **Branch:** `feature/two-layer`. **Baseline CI verificata:** `bef39b1c` (build e regressioni). **TL04-C:** contatti e riferimento termico per-layer verificati; configurazioni non termiche e scope ancora aperti.
**Stato:** codice di produzione ancora **single-layer**; `N_LAYERS=2` resta disabilitato. La presente roadmap è il **solo piano vigente**: specifica, stato verificato, attività, gate e milestone. La cronologia dettagliata è nei commit Git, non in versioni alternative del piano.
**Documenti collegati:** [audit tecnico corrente](TWO_LAYER_CODE_AUDIT.md), [roadmap generale](IMEX_SfloW2D_Development_Roadmap.md), manoscritto *IMEX_SfloW2D_v3*.

## 1. Obiettivo, invarianti e ambito

Il modello ha due sole modalità: `N_LAYERS=1` (solutore attuale) e, **solo dopo il gate finale**, `N_LAYERS=2` con L1 corrente inferiore densa e L2 nube superiore diluita nello stesso scenario gas/solidi. Non è previsto un framework con numero arbitrario di layer, né il ripristino dell'IMEX-RK storico. Le combinazioni lava/granulare tra i due layer non sono abilitate dal semplice fatto che il single-layer supporta la lava.

1. **Una sola rappresentazione autorevole** di configurazioni, proprietà fisiche e inventari; nessun globale duplicato, copia sincronizzata, wrapper transitorio o duplicazione di formule.
2. **Compatibilità single-layer per costruzione:** stessi operatori, ordine Lie, gruppo 0, criteri wet/dry, rollback, OpenMP, convenzioni di deposito e risultati, **salvo modifiche fisiche esplicitamente approvate**. È già approvato l'uso esclusivo di `T_ambient` in sostituzione di `T_env`, anche se cambia i risultati termici.
3. **Separazione dei ruoli:** integratore Lie = algoritmi, ordine, timestep, candidati, controlli numerici e rollback; moduli `physics` = leggi, sorgenti, decisioni fisiche sui contatti; `spatial` = flussi e termini di percorso. Nessuna fisica introdotta nel Lie.
4. `N_LAYERS=2` rifiutato negli input e nel layout di produzione fino a TL13/M6. Le funzionalità parziali sono testabili tramite harness isolati, non selezionabili come modello completo.
5. Non cambiare formule, regolarizzazioni, limitatori, splitting o tolleranze implicitamente durante un refactoring. Rispettare disponibilità esplicita delle modalità EXPLICIT, SEMI_IMPLICIT, NEWTON/IMPLICIT e ANALYTIC; mai fallback silenziosi.

## 2. Audit corrente, regressioni e limiti

| Area | Implementato e riscontrato nel codice | Ancora da implementare o certificare |
|---|---|---|
| Layout | `model_layout_2d` contiene `n_layers`; `state_type%q/qp` e workspace di integrazione sono rank-4 `(componente,layer,x,y)`; i consumer di produzione selezionano esplicitamente L1. | `supported_layer_count` accetta solo 1; evoluzione L2, maschere inattive e test coupled. |
| Deposito | `deposit_tot` è l'inventario canonico dello spessore; ramo lava distinto dal granulare, senza inventario autonomo `lava_deposit`; rollback e percorsi v7 sono stati adeguati. | Gate G1 completo: casi granulari/lava, `n_solid=0`, topografia on/off, I/O e restart round-trip documentati. |
| Reologia | `rheology_config_type` contiene i coefficienti e i metadati; `rheology_l1/l2` sono due oggetti; parser legacy/L1 usa alias verso **unico L1 canonico**; kernel e dispatcher migrati a `config` esplicito e vecchi metadati globali eliminati. | Il coordinatore `source_terms_2d` seleziona ancora direttamente `rheology_l1`; il registro Lie non possiede scope L1/L2; L2 non valutato in produzione. |
| Contesto per layer | `simulation_2d%layer_physics(:)` conserva fase e reologia; il contesto referenzia ora la configurazione termica canonica senza una seconda copia persistente. La vista read-only viene passata al passo Lie. | Completare le altre configurazioni SHARED/L1/L2 e gli scope fisici solo nei gate dedicati; nessuna sincronizzazione fra copie. |
| Lie | Un solo Lie primo ordine, cell-outer/OpenMP, gruppi numerici ordinati, trial e rollback; chiamata corrente `execute_lie_source_groups(...q_trial(:,1,:,:),qp_trial(:,1,:,:)...)`. | Supportare gruppi con scope `L1`, `L2`, `INTERLAYER`, candidati entrambi i layer e commit atomico senza fisica nel driver. |
| Termica | `thermal_config_type` contiene coefficienti, metodo e `T_ground` per layer; T3 passa la configurazione in modo esplicito ai kernel REAL/COMPLEX/analitici. Il nuovo contesto L1 collega la configurazione autorevole all'integratore senza accesso diretto al globale nel modulo numerico. | Validare la nuova CI e gli esempi integrati; selezione L2 solo dopo abilitazione delle strutture accoppiate. |
| Ambiente | `T_ambient` è letta esclusivamente da `GAS_TRANSPORT_PARAMETERS`, validata finita e positiva prima degli utilizzi; `T_env` e gli ingressi INITIAL/RESTART sono rimossi. | Certificare parser, casi senza fase gas e delta termico autorizzato sul nuovo commit. |
| Sorgenti fisiche | Il registro corrente è ancora single-layer; `layer_contacts_2d` puro classifica atmosfera, terreno e contatto interlayer dai candidati locali, senza essere richiamato dal Lie. | Connessione futura a TL07/TL11 e ownership delle altre chiusure fisiche. |
| Spaziale | `nonconservative_2d` è un integratore di percorso generico già utilizzato da HP-PCCU. | Tracce accoppiate, pressione interlayer, CFL coupled, wet/dry e iperbolicità. |
| I/O | NetCDF single-layer `(x,y,time)`, reader/writer restart v7 con sola sezione L1, file legacy/probe/runout. | Schema/versione per due layer, restituzione completa e compatibilità dimostrata. |

**Evidenza CI:** su `bef39b1c` [Build optimized/debug #43306145842](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43306145842) e [Cell physics/numerical regressions #43306324027](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43306324027) risultano **success**. Questi workflow non certificano full EXAMPLES, NetCDF/restart, OpenMP o un solver two-layer. I workflow non certificano full EXAMPLES, NetCDF/restart, OpenMP o un solver two-layer.

## 3. Dati, proprietà e configurazioni autorevoli

### 3.1 Stato e substrato

Layout di produzione: `q(n_variables,n_layers,nx,ny)`; `qp(n_physical_variables,n_layers,nx,ny)`. Le dimensioni di componenti, primitive ed equazioni hanno contratti distinti. `q(:,layer,j,k)` e `qp(:,layer,j,k)` sono vettori contigui per kernel monolayer. Trial, residui, snapshot fluidi e I/O seguono lo stesso layout. Topografia, geometria, `deposit_tot`, inventari di deposito/erosione e substrato **non** hanno indice layer.

Gli indici delle equazioni sono comuni ai due layer. Le quantità significative soltanto in L1 (es. eccesso di pressione interstiziale) occupano slot inattivi in L2, con maschera coerente in conversioni, trasporto, IC/BC, restart e output. Ogni worklist attiva include celle wet L1 **oppure** wet L2: `dry L1/wet L2` non può scomparire.

`deposit_tot` è lo spessore cumulativo immobilizzato (m), non il cambiamento netto della quota: in granulare è derivato dall'inventario delle classi e dalla porosità; in lava è un accumulo diretto, anche con `n_solid=0`. `dB_dt` rappresenta la modifica geometrica; non deve duplicare l'inventario. Erosione ed erodibile restano distinti. Le leggi granulari e lava sono mutuamente esclusive, anche nel single-layer.

### 3.2 Ownership: L1, L2, INTERLAYER, SHARED

Ogni proprietà deve avere **un solo proprietario**. Le proprietà intrinseche delle specie solide e gas, i parametri dell'atmosfera, gravità, geometria e timestep restano condivisi; configurazione reologica e coefficienti di chiusure specifiche appartengono al layer; coefficienti di trasferimenti L1↔L2 appartengono alla **coppia**. I parametri non diventano per-layer per la sola presenza in una namelist.

| Famiglia esistente | Contratto target | Passo |
|---|---|---|
| `RHEOLOGY_PARAMETERS` | Unico `rheology_config_type`, istanze L1/L2, legacy alias L1 in `N_LAYERS=1`; già realizzato per le leggi. | Completare scope/consumatori al momento dell'attivazione L2. |
| `TEMPERATURE_PARAMETERS` | Unico `thermal_config_type`, istanze `thermal_l1/l2`, legacy alias verso L1; coefficienti del processo per layer. | TL04-T (specifica §4). |
| `PORE_PRESSURE_PARAMETERS` | Chiusura fisica L1; nessuna copia operativa L2 solo per simmetria. | TL04-C/TL11. |
| `SOLID_TRANSPORT_PARAMETERS` e `GAS_TRANSPORT_PARAMETERS` | Proprietà delle specie e ambiente **condivisi**; opzioni delle chiusure in L1/L2/coppia se necessarie. | TL04-C/TL08–TL11. |
| `MASS_EXCHANGE_PARAMETERS` | Distinguere processi terreno↔fluido, atmosfera↔fluido e L1↔L2; non usare una configurazione globale per processi fisici differenti. | TL07–TL11. |
| `NEW RUN_PARAMETERS`, IC e BC | Configurazione globale della simulazione e condizioni del singolo layer; evitare doppie proprietà ambiente. | TL04/TL05/TL13. |
| Sorgenti radiali/laterali/fessurali e collapse | Dati geometrici condivisi dove appropriato; destinatario e proprietà del materiale espliciti. | TL07/TL11. |
| Parametri numerici, run, diagnostiche, restart | Unici; I/O con dati per layer dove richiesto. | TL02/TL13. |

### 3.3 Input reologico già implementato

In `N_LAYERS=1` è ammessa **una sola** tra `RHEOLOGY_PARAMETERS` (oggetto `rheology` puntato a `rheology_l1`) e `RHEOLOGY_L1_PARAMETERS` (oggetto `rheology_l1`); la sintassi è quella dei componenti del tipo derivato. Niente copie sincronizzate e niente vecchi coefficienti globali. La namelist L2 viene rifiutata nel percorso di produzione attuale. Con `N_LAYERS=2`, dopo il gate finale, sono richieste L1 e L2 distinte e la legacy è vietata. Per le altre famiglie adottare tipi derivati **quando** migrano, non introdurre doppi rappresentanti in anticipo.

## 4. Contratto termico definitivo e contatti fisici

### 4.1 Namelist e temperatura atmosferica

La configurazione canonica implementata è **un solo** `thermal_config_type` con due istanze `thermal_l1`, `thermal_l2`. Il parser applica:
- `N_LAYERS=1`: **esattamente una** fra `TEMPERATURE_PARAMETERS` (`thermal%...`) e `TEMPERATURE_L1_PARAMETERS` (`thermal_l1%...`); il puntatore legacy e la L1 si riferiscono **allo stesso oggetto**, senza sincronizzazione.
- `N_LAYERS=2` (solo quando abilitato): `TEMPERATURE_L1_PARAMETERS` e `TEMPERATURE_L2_PARAMETERS`, niente legacy. Input doppi, misti, incoerenti o senza dati richiesti sono errori espliciti.
- I componenti usano la sintassi del tipo derivato, come per la reologia; fixture ed esempi attivi sono migrati in modo coerente. La temperatura del terreno si specifica con `thermal%T_ground` nel gruppo legacy, o `thermal_l1%T_ground`/`thermal_l2%T_ground` nei gruppi per-layer; il vecchio campo libero `T_GROUND` non è più accettato.
- `thermal_method`, `thermal_group` e l'abilitazione termica appartengono alla configurazione del **processo del layer**, non a globali che il parser sincronizza. Il registro Lie globale ordina comunque tutti i gruppi. TL04-T3 passa la configurazione termica read-only fino ai dispatcher fisici; la selezione `thermal_l1` avviene al confine single-layer in `time_integration_2d`, non dentro `source_terms_2d` o `thermal_laws_2d`.

**`T_env` è stata eliminata dai sorgenti di produzione** (dichiarazioni, namelist, default e consumer). **`T_ambient` è l'unica temperatura dell'atmosfera e viene letta esclusivamente dalla namelist obbligatoria `GAS_TRANSPORT_PARAMETERS`**, anche se `gas_flag=F`: sono stati rimossi gli ingressi da `INITIAL_CONDITIONS` e `RESTART_PARAMETERS` ed è richiesto un valore esplicito finito e positivo. Non introdurre alias o temperatura atmosferica per-layer. `T_ground` è un parametro di ciascun `thermal_config_type` per semplificare l'interfaccia del modello sorgente, ma rappresenta **un unico terreno fisico**: quando `N_LAYERS=2` verrà abilitato, il parser richiederà `thermal_l1%T_ground == thermal_l2%T_ground` con confronto esatto e rifiuto esplicito delle differenze. Non sincronizzare automaticamente i valori. La lettura e la verifica di coppia sono predisposte, ma attualmente inaccessibili in produzione perché `N_LAYERS=2` è ancora rifiutato.

La legge atmosferica esistente resta:
```text
Q_atm,k = -k_rad,k*(T_k**4 - T_ambient**4) - k_conv,k*(T_k - T_ambient)
k_rad,k = SBconst*emissivity_k*exp_area_fract_k
k_conv,k = atm_heat_transf_coeff_k*exp_area_fract_k
```
Nessuna nuova legge di conduzione verso il terreno, di scambio termico L1-L2 o modifica del trattamento REAL/COMPLEX/ANALYTIC. L'eventuale variazione dei risultati single-layer **derivante dalla sostituzione autorizzata di `T_env` con `T_ambient`** è ammessa e va misurata/documentata; tutte le altre differenze richiedono diagnosi. L'attuale default termico legato alla lava va migrato esplicitamente e testato, non disattivato implicitamente.

**Ownership aggiornata:** `thermal_cooling_flag`, `thermal_method`, `thermal_group`, `emissivity`, `atm_heat_transf_coeff`, `exp_area_fract`, `thermal_conductivity`, `c_p`, `enne`, `emme` e `T_ground` sono componenti di `thermal_config_type`. La duplicazione del valore `T_ground` fra due configurazioni è intenzionale e subordinata all'invariante di uguaglianza: non esiste una terza variabile globale e i kernel potranno ricevere soltanto la configurazione termica del layer. `T_ambient` resta l'unica temperatura atmosferica condivisa esterna e `SBconst` una costante universale. `thermal_conductivity`, `c_p`, `enne` ed `emme` non alimentano attualmente il raffreddamento radiativo/convettivo: conservarne l'input non autorizza una nuova legge. I coefficienti `radiative_term_coeff` e `convective_term_coeff` sono grandezze **derivate**, valutate dalle funzioni pure del modulo termico, non copie persistenti.

### 4.2 Modulo fisico dei contatti `src/physics/layer_contacts_2d.f90`

Il modulo `layer_contacts_2d` è implementato come componente **puro, cell-local, senza stato persistente**, con `evaluate_layer_contacts(thickness,contacts,valid)`: classifica `atmosphere_exposed(1:2)`, `terrain_exposed(1:2)` e `interlayer_contact` usando unicamente `dry_thickness_tolerance`. Accetta una o due altezze candidate finite e non negative; per input invalidi restituisce `valid=.FALSE.` e nessun contatto. Non introduce alcuna soglia, flag namelist o legge di trasferimento.

| Stato locale | Atmosfera | Terreno | Interfaccia L1/L2 |
|---|---|---|---|
| wet L1, wet L2 | L2 | L1 | Presente |
| wet L1, dry L2 | L1 | L1 | Assente |
| dry L1, wet L2 | L2 | L2 | Assente |
| dry L1, dry L2 | Nessuno | Nessuno | Assente |
| `N_LAYERS=1`, wet | L1 | L1 | Assente |
| `N_LAYERS=1`, dry | Nessuno | Nessuno | Assente |

La tabella è il contratto di contatto **per geometria di colonna depth-averaged**: non deduce frazioni di copertura laterale non rappresentate dallo stato. La classificazione può cambiare fra gruppi Lie quando cambiano gli spessori candidati. Il modulo non decide se applicare settling, elutriazione, attrito, conduzione o altre leggi. In particolare, l'elutriazione può creare L2 in un *cold-start*: `interlayer_contact=F` non è un divieto universale di trasferimento.

### 4.3 Confine fra Lie e fisica

**Nessuna logica dei contatti nel Lie.** Il driver rende disponibili **i due stati candidati della stessa cella al livello di assemblaggio fisico**, con dati pre-gruppo coerenti e scope numerico. Il livello `physics` calcola i contatti al momento della valutazione; il kernel monolayer riceve **soltanto lo stato del layer bersaglio** e produce **soltanto il rate dello stesso layer**; eventuali informazioni esterne (esposizione, ambiente, fondo) sono dati già preparati, mai lo stato completo dell'altro layer. Il contesto fisico `INTENT(IN)` è read-only e OpenMP-safe. I termini `INTERLAYER` sono operatori **distinti**: ricevono due stati e producono la coppia di contributi.

La preparazione dei contatti deve usare il **candidato corrente**, non una lettura da `state%q` pubblicato, che durante Lie può essere vecchio. Il gruppo 0 conserva l'eccezione pre-spatial. La nuova infrastruttura può essere aggiunta come TL04 **inerte**, ma non modificare il gating delle sorgenti single-layer esistenti fino a test di equivalenza: un `*_exposed` non sostituisce indiscriminatamente regolarizzazioni wet/dry storiche.

### 4.4 Verifica specifica della migrazione

Testare equivalenza legacy/L1, rifiuto duplicati e L2 prematuro, configurazioni termiche indipendenti e kernel REAL/COMPLEX/ANALYTIC. La suite `TEST_SOURCE_REGISTRY` comprende ora anche i quattro casi di contatto L1/L2, single-layer wet/dry, soglia canonica e input non validi. Nessun contatto è ancora applicato nel percorso Lie.

## 5. Partizione matematica e Lie

```text
dq1/dt = H1(q1,q2,B) + S_local1(q1;C1,E1) + Gamma(q1,q2) + S_external1
dq2/dt = H2(q1,q2,B) + S_local2(q2;C2,E2) - Gamma(q1,q2) + S_external2
```
`H1/H2` includono flussi, pressione, topografia e termini non conservativi **spaziali**. `S_local` è una legge monolayer; `E_k` è contesto fisico esterno già preparato; `Gamma` è il trasferimento strettamente interno fra layer. Gli scambi esterni hanno bilancio separato. Ogni legge locale riceve il solo stato del proprio layer e restituisce il solo rate di quel layer.

**Un unico integratore Lie globale**, con un solo timestep, ordinamento e rollback. Il registro usa `(source_identity,scope)` per l'unicità; ogni gruppo dichiara `L1`, `L2` o `INTERLAYER`, metodo, ordine e footprint, compresa la gestione delle risorse condivise. La stessa legge può comparire in L1 e L2. L'attuale `source_group_type` **non** dispone ancora di `scope`, e `time_integration_2d` passa al Lie soltanto `q_trial(:,1,:,:)`, `qp_trial(:,1,:,:)`: non dichiarare completata TL07.

Il Lie gestisce candidati e snapshot pre-spatial del gruppo 0, transazioni dello stato e substrato, geometria, stochastic, buone proprietà numeriche e rollback. Il livello `physics` gestisce contatti, abilitazione delle leggi, rate e trasferimenti. Con gruppi `INTERLAYER` si acquisisce uno **stesso snapshot pre-gruppo** di L1 e L2, si calcola un solo `Gamma(n_eqns)` con convenzione positiva L2→L1 e si costruiscono contemporaneamente:
```text
q1_candidate = q1 + dt*Gamma
q2_candidate = q2 - dt*Gamma
```
Un solo limiter di disponibilità condiviso; nessun limitatore indipendente che rompa l'antisimmetria. Validazione e rollback atomici di entrambi i layer, substrato, geometria e diagnostiche. La variabile termica conservativa trasferita va contabilizzata; l'energia cinetica e la dissipazione possono avere regole ulteriori, esplicite. L'approssimazione del gas di packing è un **difetto esterno diagnosticato**, non un trasferimento fittiziamente conservativo. Non imporre un Newton globale `2*n_eqns` senza evidenza di necessità.

Preservare esattamente la semantica del gruppo 0 e delle celle inizialmente asciutte, i quattro trattamenti numerici, la struttura cell-outer e gli scratch per worker; nessuna allocazione per cella o modifica di globali nell'OpenMP.

## 6. Operatore spaziale accoppiato HP-PCCU

Riutilizzare `src/spatial/nonconservative_2d.f90` **senza cambiarne l'API**: callback generica `I(W,dW/ds)` e quadratura a 1/2/3 punti. Non introdurre un secondo motore né differenze finite automatiche di derivate di percorso. Percorso e variabili estese scelti esplicitamente e coerenti su facce, interno e wet/dry.

La geometria comune è:
```text
B       = quota del substrato
eta_1   = B + h1
eta_tot = B + h1 + h2
```
Con `G=g_cell`, `K2=rho_m2*gprime2` e `L=K2*h2`, nella direzione `j` la coppia di pressione è:
```text
S_inter1 = -G*h1*d_j(L)
S_inter2 = -L*d_j(G*h1)
I1       = -G*h1*dL/ds
I2       = -L*d(G*h1)/ds
I1 + I2  = -d(G*h1*L)/ds
D1 + D2  = -[G*h1*L]_L^R    (senza clamp)
```
Valutare insieme i contributi sulla stessa geometria e con la stessa quadratura. Il fattore metrico è **fuori** dalla derivata per L1 e **dentro** per L2. Clamp di buoyancy, perdita dell'identità di somma e relative diagnosi sono espliciti. Il test della riduzione alle proprietà identiche (`K1=K2`) riguarda il **residuo idrostatico completo**, con topografia e spessori arbitrari, non soltanto la coppia. Il caso curvo può avere il residuo fisico del manoscritto `-K2*h1*h2*d_j(G)`; non promettere well-balancing non previsto.

Le velocità `u ± sqrt(G*gprime*h)` dei due layer separati **non** garantiscono un CFL valido del sistema accoppiato. Prima del rilascio derivare/verificare il bound dei modi, condizioni d'iperbolicità e comportamento in shear, clamp e wet/dry.

## 7. Chiusure fisiche da implementare

- **Settling L2→L1 (TL08):** classi, velocità terminali comuni, `u_Gamma`, entalpia del donatore, `CHI_XFER_L2`, `H_PLUG_XFER_L2`, `XFER_U2_AMP_CAP`; cap su inventario senza distruggere il momento. **Non** è deposito al terreno e non modifica `deposit_tot`.
- **Elutriazione L1→L2 (TL09):** `SHEAR_POWER`/`LINEAR_USTAR`, Rouse/Shields, disponibilità `Apack`, co-lift e densità target, cold-start di L2, limiti per classe. Presenza geometrica di interfaccia non è prerequisito assoluto per il cold-start.
- **Drag (TL10):** forza di coppia vettoriale tipo Chézy con contributi opposti ai momenti, nulla per velocità coincidenti/layer assente; separare dissipazione dell'energia cinetica da trasporto della variabile termica.
- **Scambi esterni (TL11):** entrainment atmosfera→layer esposto, deposito/erosione sul terreno per layer a contatto, consolidazione L1 con `ALPHA_PACKING_L1`, gas burial, `d32`, pressione efficace e pore pressure; lofting L2 e budget plume. Non duplicare contributi se entrambi i layer sono presenti. Preservare ramo lava single-layer (non accoppiarlo come granular+lava).
- **Pressione interstiziale:** attiva soltanto in L1; slot L2 inattivo; seeding e dissipazione da validare.
- **Coefficiente packing:** distinguere packing L1, jamming dello hindered settling e massimo packing di ammissibilità.
- **Termica:** riusare cooling atmosferico e sue varianti già presenti; eventuale conduzione col terreno o scambio termico interlayer **fuori ambito della migrazione attuale**, richiede decisione e gate fisici autonomi.

Le formule esatte, dimensioni, opzioni e limiti di validità delle chiusure TL08–TL11 vanno estratti dalle sezioni 2.5–2.7, 2.12 e 2.14 del manoscritto e verificati prima dell'implementazione; non dichiarare equivalenti leggi attuali e leggi del manoscritto per analogia dei nomi.

## 8. I/O, input e restart

Single-layer: mantenere nomi/contenuto NetCDF legacy `(x,y,time)`, file `.asc/.q_2d`, sonde e diagnostiche, salvo modifiche deliberate e documentate (`T_ambient` e deposito canonico). Il rank-4 interno non cambia da sé il formato; il restart v7 legge/scrive L1 ed inventari condivisi. Verificare byte-compatibilità reale o prevedere reader di conversione/rifiuto esplicito senza perdita silenziosa.

Two-layer: schema versionato con `n_layers`, layout, modelli e configurazioni, due stati e slot inattivi, substrato unico, topografia, variabili di continuazione, stochastic, diagnosi; NetCDF con dimensione/suffissi L1/L2, `eta_1`, `eta_tot`, Ri, densità, budget plume e `deposit_tot` una sola volta. Testare restart in assenza/presenza di scambi, deposizione, wet/dry e confronto della continuazione.

## 9. Piano operativo: task TL00–TL13

Ogni task produce codice definitivo, test e documentazione coerenti. **Implementato** non equivale a **gate scientifico chiuso**.

| Task | Stato al commit audit | Deliverable e criteri |
|---|---|---|
| **TL00 — baseline** | **Parziale** | Identificare riferimento immutabile, tag verificato, protezione branch, matrici test/tempi e output riproducibili; gate G0. |
| **TL00a — deposito canonico** | **Implementato, gate G1 aperto** | Completare test granulare/lava, `n_solid=0`, `topo_change_flag`, budget, v7/NetCDF e rollback senza doppio inventario. |
| **TL01 — rank-4 stato** | **Implementato per L1, gate G2 aperto** | Assicurare layout/contiguità e inizializzazione; test delle dimensioni, dello slot L2 inerte e dei campi asciutti. |
| **TL02 — migrazione consumer** | **Implementato per L1, audit integrato aperto** | Verificare tutte le dipendenze di state, spatial, Lie, stochastic, input/output, test, OpenMP; eliminare wrapper obsoleti. |
| **TL03 — gate single-layer** | **CI parziale superata, gate G2 aperto** | Full EXAMPLES, bilanci, wet/dry, I/O/restart, OpenMP, profiling/memoria e confronto baseline; registrare eventuale delta termico autorizzato separatamente. |
| **TL04-R — configurazione reologica** | **Implementata e CI superata** | Tipo canonico, parser legacy/L1, parametri espliciti nei kernel; completare le dipendenze dal contesto per L2 quando attivato. |
| **TL04-T — configurazione termica** | **T1/T2 e T_ground verificate; T3 compilata e regressioni superate** | Configurazione canonica L1/L2, input legacy e propagazione read-only attraverso source assembler, NEWTON REAL/COMPLEX e metodo analitico. Restano esempi integrati, confronto single-layer e scope L2. |
| **TL04-C — altre configurazioni e contatti** | **Contatti e riferimento termico per layer verificati in CI; fase complessiva parziale** | Contatti puri, descriptor L1 che referenzia il thermal canonico senza duplicarlo, percorso read-only fino all'integratore; test L1/L2 dedicati aggiornati. Restano ownership delle altre chiusure, IC/BC, slot inattivi e scope. |
| **TL05 — HP-PCCU indipendente per slot** | **Da implementare** | Trasporto/ricostruzione per L1 e L2 isolati e wet/dry combinati, `eta_1`/`eta_tot`; recupero single-layer quando assente una corrente. |
| **TL06 — HP-PCCU accoppiato** | **Da implementare** | Percorsi coppia, clamp, topografia, identità, CFL e iperbolicità; gate G4. |
| **TL07 — Lie con scope** | **Da implementare** | Registro unico con L1/L2/INTERLAYER, due candidati cell-local, rate monolayer isolati, `Gamma` di coppia, commit/rollback atomico, gruppo 0 e OpenMP; gate G5. |
| **TL08 — settling L2→L1** | **Da implementare** | Conservazione di masse/specie/momento/termica e nessun deposito terreno. |
| **TL09 — elutriazione L1→L2** | **Da implementare** | Inventario, co-lift, cold-start e compatibilità con i contatti. |
| **TL10 — drag interlayer** | **Da implementare** | Momento uguale/opposto, forza nulla ai limiti, analisi stabilità. |
| **TL11 — sorgenti esterne e basali** | **Da implementare** | Uso dei contatti da `physics` per atmosfera/terreno, consolidazione, entrainment, lofting, pore pressure; chiusure del manoscritto e gate G6. |
| **TL12 — validazione scientifica** | **Da implementare** | Benchmark coupled, resting states, budgets, wet/dry, OpenMP, CFL, restart, prestazioni, regressioni single-layer; gate G7. |
| **TL13 — integrazione e rilascio** | **Da implementare** | Contratti di input/NetCDF/restart, esempi e postprocessing, documentazione delle opzioni e limiti, decisione esplicita prima di abilitare `N_LAYERS=2`. |

**Dipendenze effettive:** G0/G1/G2 consolidano il riferimento single-layer; TL04-R è già disponibile. TL04-T/C (incluse interfacce read-only e contatti inerti) precedono la connessione ai kernel e a TL07/TL11; TL05 precede TL06; TL07 precede TL08–TL11; TL12 richiede TL06–TL11 e I/O integrati; TL13 chiude G7. Le formule e i test isolati delle chiusure possono essere preparati prima, senza collegarli al percorso di produzione.

## 10. Milestone M0–M6 e gate di integrazione

| Milestone | Copertura e risultato | Stato al commit audit | Condizione di uscita |
|---|---|---|---|
| **M0** | TL00: baseline e specifiche | **Da certificare** | G0: tag, audit e benchmark riproducibile, branch protetto verificato. |
| **M1** | TL00a: deposito unico | **Implementata, non certificata** | G1: deposizione lava/granulare, restart, topografia e output. |
| **M2** | TL01–TL03: rank-4 single-layer validato | **Parziale, non certificata** | G2: full EXAMPLES, confronto baseline, NetCDF/restart, OpenMP e numerica. |
| **M3** | TL04–TL05: infrastruttura L1/L2 isolata | **Contatti verificati; descrittore termico da certificare, resto aperto** | G3 richiede isolamento, scope, assenza race, test single-layer integrati e trasporto L2 isolato. |
| **M4** | TL06: operatore spaziale accoppiato | **Non iniziata** | G4: path, pressione completa, topografia, CFL, wet/dry, clamp e iperbolicità. |
| **M5** | TL07–TL11: Lie coppia e chiusure | **Non iniziata** | G5/G6: rate per-layer, coppia conservativa, budget esterni, rollback, contatti fisici e tutti i processi verificati. |
| **M6** | TL12–TL13: modello disponibile | **Non iniziata** | G7: benchmark scientifici, prestazioni, I/O e restart, approvazione esplicita all'attivazione. |

Le milestone sono **punti di proposta di merge**, non autorizzazioni automatiche. `master` deve rimanere utilizzabile single-layer. Per ogni PR: build optimized/debug e strict, regressioni pertinenti, campioni integrati, correttezza fisica, performance e OpenMP, evidenze dei gate, revisione e approvazione esplicita dell'utente. Nessun merge su `master` senza autorizzazione. Dopo un merge, sincronizzare consapevolmente il branch; nessuna storia riscritta o tag presunti.

## 11. Matrice di verifica minima e prove mancanti

| Gate | Prove obbligatorie | Evidenza attuale |
|---|---|---|
| **G0** | Baseline taggato, esempi riproducibili, confronti, condizioni OpenMP e prestazioni | Non certificato. |
| **G1** | `deposit_tot` granulare/lava, `n_solid=0`, rollback, erosione, topografia e restart v7 | Codice e alcuni test esistenti; matrice integrata non certificata. |
| **G2** | Rank-4 L1, IC/BC, lie, spatial, NUMERICS, full EXAMPLES, NetCDF, restart, OpenMP | CI build/regressioni su `410ad99e` **success**; resto non certificato. |
| **G3** | Configurazioni reologica/termica distinte, contatti quattro casi + mono, OpenMP e input | Modulo contatti verificato su `4e9a591a`; riferimento termico dai descrittori implementato ma non ancora certificato sul nuovo commit. Gate complessivo aperto. |
| **G4** | Idrostatica con path e topografia, coppia, residuo curvo, clamp, wet/dry, CFL | Da implementare. |
| **G5** | Sorgenti L1/L2 indipendenti, coppia `Gamma`, limiti e rollback atomico, gruppo 0 | Da implementare. |
| **G6** | Settling, elutriazione cold-start, drag, consolidazione, entrainment, lofting, termica e budget | Da implementare. |
| **G7** | Benchmark manoscritto, stabilità, restart coupled, I/O, OpenMP e single-layer confrontato | Da implementare. |

I test dei contatti non devono introdurre un gating nuovo nel single-layer, e i test dei kernel termici devono verificare REAL/COMPLEX e il trattamento analitico. L'evoluzione delle sorgenti deve usare stati candidati coerenti, non stati pubblicati obsoleti. Testare esplicitamente `wet/wet`, `wet/dry`, `dry/wet`, `dry/dry`, fronti sottili, risultati con un layer assente, rollback e transazioni di substrato.

## 12. Rischi e decisioni scientifiche da chiudere nei task

1. Catalogare tutte le proprietà termiche (`c_p`, `enne`, `emme`, `thermal_conductivity`, ecc.) prima di decidere i campi esatti del tipo canonico; proprietà condivise non vanno duplicate. Nessuna nuova conduzione è autorizzata.
2. Verificare il criterio canonico wet/dry e il significato della copertura in una colonna depth-averaged; non inferire frazioni di copertura sub-cella non presenti.
3. Definire le combinazioni ammissibili di carrier/specie e gli slot inattivi di L2, con input completi per L1/L2/coppia.
4. Derivare il percorso esteso e la CFL coupled, includendo perdita di iperbolicità, clamp e topografia curva.
5. Separare trasferimenti strettamente conservativi da gas packing e altri scambi esterni; definire bilanci di massa, momento ed energia termica.
6. Determinare le chiusure di settling, packing, entrainment, lofting e consolidazione del manoscritto prima di abilitarle.
7. Verificare il contratto del gruppo 0, Newton/analitico e rollback anche su coppie di stati.
8. Definire conversione/versioning restart e nuovi output NetCDF, con riproduzione single-layer.
9. Dimostrare performance e riproducibilità OpenMP con configurazioni read-only e nessuna fisica hard-coded nel Lie.

## 13. Regola operativa

Procedere su `feature/two-layer` con audit completo della catena dipendenze **prima** di ogni modifica, gruppi coerenti di commit, test REAL/COMPLEX quando coinvolti e registrazione di esiti sul commit effettivo. La roadmap resta un **documento consolidato dello stato e del piano vigente**: aggiornare sezioni esistenti, non accumulare note cronologiche né mantenere contemporaneamente piano vecchio e nuovo. Lo storico è nei commit; le evidenze statiche nel [rapporto di audit](TWO_LAYER_CODE_AUDIT.md).

**Riferimenti tecnici:** *IMEX_SfloW2D_v3* (sez. 2.1–2.2, 2.5–2.7, 2.12, 2.14, 3–4); `src/state/`, `src/physics/`, `src/spatial/nonconservative_2d.f90`, `src/spatial/pccu_2d.f90`, `src/numerics/lie_source_executor_2d.f90`, `docs/FORTRAN_STYLE.md`.
