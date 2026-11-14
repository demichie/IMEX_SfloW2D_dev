# IMEX-SfloW2D — Audit tecnico del branch two-layer

**Verifica statica:** 14 novembre 2026. **Branch:** `feature/two-layer`; baseline CI [`d901609`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/d9016093a4b1f5c790e083a7d137a055c890a66e) con TL05-A/B operatori e worklist per slot. **Piano vigente:** [roadmap two-layer](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md). **Scopo corrente:** applicare il principio di semplicità ai contatti: usare il classificatore puro esistente, senza filtro aggiuntivo o attivazione di L2. La sezione TL04-C0 in fondo conserva la fotografia storica precedente.

## Evidenze verificate

- Le modifiche termiche sono state preparate esclusivamente su `feature/two-layer`; nessun merge in `master` è autorizzato o effettuato in questa fase.
- [Build optimized/debug #43505177703](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43505177703) e [regressioni #43505283197](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43505283197) **success** su `d901609`, che include TL04-C3/C4 e TL05-A/B. L'utente ha rinviato la nuova esecuzione degli EXAMPLES; le prove precedenti non certificano i nuovi commit.
- L'audit ha riscontrato la propagazione reologica/termica/Darcy e la nuova ownership shared di terreno e solidi. Per TL04-C0 sono stati letti sul branch i moduli di input, il coordinatore mass exchange, le sei leggi, source registry, Lie e relativi caller, oltre a esempi/test rappresentativi. Restano da completare una ricerca esaustiva di tutti i generatori/fixture e la validazione sul checkout.
- Nessuna compilazione/esecuzione è stata lanciata da questo audit; gli esiti CI sono evidenze **esterne** limitate ai job eseguiti.

## Mappa corrente delle dipendenze

| ID | Rilevazione sul codice attuale | Lavoro aperto / vincolo |
|---|---|---|
| **A01 — Layout** | `src/core/model_layout_2d.f90:supported_layer_count` ritorna `n_layers==1`; `src/state/state_fields.f90` alloca `q/qp` rank-4. `src/numerics/time_integration_2d.f90` alloca trial rank-4 ma seleziona `(:,1,:,:)` in chiamate Lie e commit. | Test di layout e funzioni L2 isolate; non attivare `N_LAYERS=2` per eludere i gate. |
| **A02 — Config reologica** | `rheology_config_type` e `rheology_l1/l2` sono canonici; il parser accetta intenzionalmente il nome breve e `_L1` in alternativa per lo stesso oggetto. I kernel reologici REAL/COMPLEX/analitici ricevono `config` read-only. | TL04-R già verificata su `662e796`; L2 resta inattivo. Non rimuovere l'alias NAMELIST single-layer né creare copie sincronizzate. |
| **A03 — Contesto fisico** | `simulation_2d` possiede `layer_physics(:)`, `spatial_operators(:)` e `domains(:)` secondo il layout canonico. Il timestep di produzione usa L1, con geometria e inventari comuni. | L2 ancora non evoluto; non duplicare coefficienti fisici o inventari shared. |
| **A04 — Namelist termica e substrato** | `thermal_l1/l2` possiedono il cooling; `thermal => thermal_l1` mantiene intenzionalmente il nome breve in `N_LAYERS=1`. Il solo `substrate%T_ground` descrive il terreno; non esistono più copie termiche della temperatura del suolo. | Tenere indipendenti proprietà shared e closure per-layer. Validare i casi fisici dove i precedenti input avevano temperature del terreno differenti. |
| **A05 — Temperatura atmosferica** | `T_ambient` è l'unica temperatura atmosferica letta dalla namelist GAS; `T_env` non esiste più nella produzione. TL04-T1 è stata verificata in CI. | Restano verifiche integrate del delta termico autorizzato e degli esempi, in particolare lava. |
| **A06 — Termica fisica** | TL04-T3 usa configurazione read-only dal coordinamento single-layer fino al dispatcher termico, senza duplicazioni di coefficienti; cooling REAL/COMPLEX/analitico invariato nelle formule. | Validare esempi integrati lava/dry e futuro scope L2. |
| **A07 — Contatti** | `layer_contacts_2d.f90` è l'unica classificazione canonica dei quattro stati L1/L2; `TEST_SOURCE_REGISTRY` verifica già tutti i casi, le soglie e gli input invalidi. La routine ridondante `filter_contact_mass_exchange_mask` del commit `e4ef157` è stata rimossa. | TL07/TL11: usare direttamente `evaluate_layer_contacts` dentro l'assemblaggio cell-local `source_terms_2d`, senza helper o API aggiuntivi. Lie e kernel mass exchange rimangono invariati, L2 non attivo. |
| **A08 — Registro/integrazione Lie** | `source_registry_2d.f90:source_group_type` contiene id, method, source_names, equation_ids, cache; nessun `scope`. `lie_source_executor_2d.f90:execute_lie_source_groups` riceve array rank-3 di L1 ed esegue gruppi cell-outer. | Scope L1/L2/INTERLAYER, coppia candidati e trasferimento unico, rollback atomico; kernel monolayer ricevono solo q/qp e rate del layer bersaglio. |
| **A09 — HP-PCCU/path** | Workspace `spatial_operators(:)` e `domains(:)` indipendenti, CI su `d901609`. TL05-C1 aggiunge in `hyperbolic_2d` il clamp idrostatico della sola L2 (`Gamma2=max(rho_m2*gprime2,0)`), selezionato al setup dall'indice canonico di layer, senza nuovo kernel. Test gas L2 isolato integrato in `TEST_SPATIAL_OPERATOR`. | C1 richiede CI sul nuovo commit e la verifica del test gas. TL05-C rimane aperto; TL06 integra i path accoppiati `B,eta_1,eta_tot`, wet/dry e bound CFL coupled. Nessuna estensione dei flag globali della fase per attivare soluzioni improprie. |
| **A10 — Deposito e substrato** | `deposit_tot` canonico e ramo lava/granulare distinti; scratch e transazioni del Lie aggiornati. | Prove end-to-end (lava `n_solid=0`, erosione, topografia, restart v7, NetCDF) ancora non certificate dall'evidenza disponibile. |
| **A11 — Specie e processi** | Substrato e specie solide restano SHARED; TL04-C3/C4 ha migrato coefficienti e attivazioni del mass exchange in `mass_exchange_l1/l2`, mantenendo un assemblatore e i due nomi NAMELIST L1 alternativi. `gas_loss_flag` appartiene a Darcy. | EXAMPLES della migrazione differiti; budget scientifici, restart e OpenMP da completare, senza ripristinare coefficienti globali. |
| **A12 — I/O** | NetCDF single-layer `(x,y,time)`, restart v7 con accesso L1 e substrato condiviso, formati legacy. | Versione two-layer, schema del layer, continuazione e test round-trip; preservare produzione mono. |
| **A13 — Workflow** | Entrambe le Actions sono **success** su `d901609` (13 novembre); gli EXAMPLES dopo il mass exchange sono rinviati. | Nuovo clamp HP-PCCU per L2 richiede build optimized/debug e regressioni numeriche/spaziali sul nuovo commit; le Actions precedenti non certificano TL05-C1. |

## Catena di dipendenze da rispettare

1. **Configurazione:** parser `inpout_2d.f90`, `constitutive_parameters_2d.f90`, `source_input_2d.f90`, `layer_physics_2d.f90`, `simulation_2d.f90`; reset, validazione, nomi legacy, esempi, generatori, restart e diagnostiche.
2. **Kernel:** `source_terms_2d.f90`, `thermal_laws_2d.f90`, `thermal_cooling_2d.f90`, `thermal_cooling_rate.inc`, REAL/COMPLEX/analitico e Newton; configurazioni esplicite senza globali fisici condivisi mutabili.
3. **Lie:** `source_configuration_2d.f90`, `source_catalog_2d.f90`, `source_registry_2d.f90`, `lie_source_executor_2d.f90`, `time_integration_2d.f90`, workspace/scratch per worker, gruppo 0, rollback e OpenMP. Il Lie resta privo di decisioni fisiche sui contatti.
4. **Contatti:** preparazione nel livello `physics` a partire dalla coppia di stati candidati cell-local; kernel monolayer con stato/rate monolayer, coppie soltanto negli operatori INTERLAYER. Contatto ≠ abilitazione: cold-start L2 ammesso.
5. **Spaziale:** `nonconservative_2d.f90`, `pccu_2d.f90`, `hyperbolic_2d.f90`, `reconstruction_2d.f90`, `spatial_operator_2d.f90`; trace wet/dry e bound CFL.
6. **Build/tests:** `src/Makefile.am` ha dipendenze esplicite; aggiungere la dipendenza di ogni nuovo modulo ai consumer; aggiornare test unitari, test input, esempi, generatori e suite CI con ogni migrazione coerente.

## Gate e limiti della verifica

| Gate | Evidenza disponibile | Stato |
|---|---|---|
| G0 baseline/tag | Audit storico di baseline, nessuna certificazione completa verificata oggi | **Aperto** |
| G1 deposito | Implementazione e test mirati, round-trip/full matrix non dimostrati | **Aperto** |
| G2 rank-4 single-layer | Build e regressioni **success** su `d901609`; EXAMPLES dopo il mass exchange rinviati. Baseline quantitative, restart/NetCDF e OpenMP non interamente certificati. | **Parziale: aperto** |
| G3 configurazioni L1/L2 | Reologia, termica, Darcy, substrato, solidi e mass exchange per-layer verificati dalle due Actions su `32e4e33`. Contact classifier inerte; EXAMPLES e gate scientifico C5 aperti. | **Parziale: aperto** |
| G4 spazio | No coupled production | **Aperto** |
| G5 Lie coppia | No scope, no candidate pair production | **Aperto** |
| G6 chiusure | Da implementare | **Aperto** |
| G7 produzione | No full coupled benchmark | **Aperto** |

**Esito:** la preparazione strutturale (rank-4 e configurazione reologica) è reale e verificata parzialmente; **non** dimostra il two-layer. Non sono stati modificati codici fisici o risultati per questo audit. Qualunque nuovo commit richiede le sue prove. Il piano di lavoro aggiornato e le milestone si trovano **solo** nella roadmap collegata.

## TL04-C0 — chiusura dipendenze del mass exchange (13 novembre 2026)

**Baseline di audit:** `feature/two-layer` a `5e0e8ed`, che modifica soltanto la roadmap dopo la CI verde sul codice `662e796`. Questo è un audit di lettura; non costituisce implementazione C3/C4, non certifica i test aggiornati e non abilita `N_LAYERS=2`.

| Sottosistema | Implementazione riscontrata | Destinazione e obblighi per la migrazione |
|---|---|---|
| Lettura mass exchange | `src/numerics/source_input_2d.f90` dichiara `MASS_EXCHANGE_PARAMETERS` con scalari `mass_exchange_method/group`, `erosion_coeff`, `settling_flag`, `loss_rate`, `entrainment_flag`, `gas_loss_flag`, `liquid_vaporization_flag` e le opzioni lava. Il reader rifiuta `SOURCE_GROUPS`, controlla il numero delle occorrenze e valida lava, metodo e gruppo. `inpout_2d.f90` azzera in modo distribuito diversi coefficienti. | Introdurre in `mass_exchange_parameters_2d.f90` `mass_exchange_config_type` per L1/L2 e un unico METHOD/GROUP per ciascun assemblatore combinato. In single-layer leggere **una sola** tra `MASS_EXCHANGE_PARAMETERS` con `mass_exchange%...` e `MASS_EXCHANGE_L1_PARAMETERS` con `mass_exchange_l1%...`, direttamente sullo stesso oggetto canonico. Rifiutare duplicati e L2 prematuro. Nessun copia-incolla dei coefficienti fra nomi o parser di conversione. |
| Owner di processi | `constitutive_parameters_2d` possiede `erosion_coeff` e `loss_rate` (scalari `ALLOCATABLE`), `settling_flag`, `entrainment_flag`, `A_drag`, `B_drag`, `collective_settling_flag` e `gamma_steam`; `parameters_2d` possiede `liquid_vaporization_flag`. `mass_exchange_lava_deposition_2d` possiede i sette parametri di immobilizzazione. `alphastot_min` è shared. | Spostare solo parametri di processo nei sottotipi canonici erosione, deposizione, entrainment, liquid loss, vaporization e lava immobilization; eliminare i globali obsoleti insieme a tutti i consumatori. Mantener separate temperatura ambiente, materiali solidi e stato del substrato. Non confondere il packing condiviso con i coefficienti della legge di drag. |
| Darcy / gas loss | `gas_loss_flag` è letto nella NAMELIST mass exchange ma consultato in `pore_pressure_darcy_2d.f90` e inizializzato in `inpout_2d.f90`. | Il proprietario fisico è la configurazione Darcy: migrare la sua sintassi e i consumer come parte coerente del riordino input; non duplicarlo nei sottotipi del mass exchange. Preservare i criteri di attivazione gas/pore-pressure. |
| Transport e geometria | `gamma_steam` è oggi letto da `GAS_TRANSPORT_PARAMETERS` ed usato nella vaporizzazione; `liquid_vaporization_flag` comanda anche `geometry_2d:topography_zones`. `loss_rate` è controllato nel parser dei fluidi. | Allineare i consumer ed i controlli alla configurazione di processo senza cambiare le precondizioni e il timing di setup; `water_level` e le proprietà della fase restano shared. Considerare la semantica shared delle zone topografiche con possibili due closure layer-specifiche. |
| Fisica e budget | `mass_exchange_2d.f90` coordina i sei processi mediante un solo `source_mask`, limiti condivisi degli inventari carrier/solidi e aggregazione di massa, momento, energia e `dB_dt`. Le leggi sono in `mass_exchange_{erosion,deposition,carrier_liquid_loss,entrainment,vaporization,lava_deposition}_2d.f90`. | Propagare configurazioni `INTENT(IN)` lungo lo stesso coordinatore ed i kernel, senza cambiare formule, ordine di valutazione, limiti, trattamento wet/dry, process rates o scratch per-worker. `NEWTON`/`ANALYTIC` per lo scambio combinato restano non supportati. |
| Scheduling e Lie | `source_configuration_2d` legge i globali di abilitazione per iscrivere i processi nel registry e assegna lo stesso METHOD/GROUP. `simulation_2d` inizializza il descriptor L1; `time_integration_2d`, `lie_source_executor_2d`, `source_terms_2d` raggiungono l'assemblatore con la reologia ma non ancora con `mass_exchange_config_type`. | Aggiungere un riferimento non proprietario a `layer_physics_type` e la propagazione esplicita dalla simulazione fino alla valutazione fisica, senza leggere `mass_exchange_l1` dentro kernel riusabili, senza nuove sorgenti Lie separate né nuova logica fisica nel Lie. |
| Test, input, output | `TEST_MASS_EXCHANGE/test_projection_stress.f90` e `TEST_SOURCE_REGISTRY/test_source_registry.f90` impostano i vecchi globali; `TEST_SPATIAL_OPERATOR` e `TEST_WORKSPACE_LIFECYCLE` consultano ancora flag correlati. Esempi CF/TAAL/BW/DAMBREAK/ETNA_LAVA contengono variabili scalari nella NAMELIST mass exchange; alcune opzioni di vaporization/drag appartengono oggi a gruppi differenti. `inpout_2d` usa anche gli stessi flag per diagnostiche e canali di output. | Migrare in un unico gruppo coerente EXAMPLES, TESTS, generatori, fixture, verifiche di input, diagnostiche e formati di output interessati. Testare le due grafie alternative L1, collisione dei gruppi, rifiuto L2, flag non validi, compatibilità numerica dei casi esistenti e modalità esplicita. |
| Build | `src/Makefile.am` dichiara l'ordine dei moduli, dipendenze Fortran e i sorgenti mass exchange; `source_input_2d` oggi dipende perfino dal modulo della legge lava per importarne i parametri. | Collocare il nuovo modulo di parametri prima di reader, source configuration, descriptor e leggi; aggiornare tutte le dipendenze `.mod`, eliminando la dipendenza inversa input-verso-legge quando possibile. |

**Vincolo Fortran verificato localmente con gfortran:** `NAMELIST` con due grafie e alias `POINTER` verso lo stesso oggetto canonico funziona anche con componenti di sottotipi. Un derivato con componenti `ALLOCATABLE` non è leggibile nella `NAMELIST` ordinaria senza defined input/output: il compilatore rifiuta il sorgente già in compilazione. Quindi le quantità scalar `erosion_coeff`/`loss_rate` vanno rappresentate come componenti scalari non allocatable con un default definito; eventuali array allocatable di specie restano su owner condivisi adatti all'input. È un test della semantica del linguaggio, **non** una build del repository.

**Stato di TL04-C0/C3/C4 (13 novembre 2026):** audit delle dipendenze esteso ai consumer Fortran, alle fixture `.inp`/`.template`, ai generatori e ai runner standalone della geometria. La configurazione canonica L1/L2, le due grafie alternative della NAMELIST single-layer, la migrazione dei sei processi e la propagazione attraverso il Lie sono predisposte in un unico gruppo di modifiche. `GAS_LOSS_FLAG` appartiene alla configurazione Darcy; non viene duplicato nel mass exchange. **Verifica ancora pendente:** non è stato possibile compilare l'intero repository in questo ambiente; servono entrambe le Actions sul commit effettivo, le regressioni numeriche e gli EXAMPLES. La prova isolata con `gfortran` della NAMELIST nested/alias è riuscita, ma non certifica il parser o la fisica di produzione. **TL04-C5 e G3 restano aperti.**

## TL05-A — workspace spaziale per layer slot (13 novembre 2026)

**Dipendenze ispezionate:** il driver `src/IMEX_SfloW2D.f90` accede tre volte all'owner spaziale del contesto (preparazione HP e CFL). `simulation_2d` lo inizializza, lo passa due volte a `time_integration_2d` e lo finalizza; l'integratore lo usa per il residuo spaziale L1 e le diagnostiche di errore. Il backend `spatial_operator_2d` possiede già ricostruzione, flussi e velocità caratteristiche nel medesimo workspace.

**Modifica preparatoria:** `simulation_context_type%spatial_operators(:)` sostituisce lo scalare, viene allocato con `model_layout%n_layers`, inizializzato e finalizzato per elemento. `time_integration_2d` riceve l'array e continua a selezionare esclusivamente il primo elemento. Il driver utilizza `spatial_operators(1)` per la preparazione del workset HP e per il CFL. Non nasce un secondo motore spaziale e non cambiano formule, path, limiti o ordinamento numerico.

**Regressione:** `TESTS/TEST_SPATIAL_OPERATOR/test_spatial_operator.f90` confronta bit per bit i risultati CFL e HP-PCCU di due workspace distinti sullo stesso stato single-layer; `TEST_WORKSPACE_LIFECYCLE` esercita due inizializzazioni/finalizzazioni di un array di workspace per slot. La suite spaziale della workflow CI lo esegue in configurazione strict e optimized. Verificata localmente con `gfortran` la sintassi di array allocatable di tipi con procedure e passaggio assumed-shape; la build del repository sul nuovo commit resta da verificare.

**Non implementato:** trasporto indipendente L2, `eta_1`/`eta_tot`, pressione interlayer, CFL del sistema accoppiato e attivazione `N_LAYERS=2`. Restano aperti TL04-C5/G3, TL05 e tutti i gate successivi. La sezione TL04-C0 sottostante documenta una baseline precedente e non va letta come stato del codice odierno.

## TL05-B — workset spaziali indipendenti per slot (13 novembre 2026)

**Problema rilevato:** gli operatori `spatial_operators(:)` sono già separati ma il contesto possiede ancora un solo `domain_type` con maschere wet/dry, liste di celle/facce e aloni HP. L2 può avere celle bagnate/asciutte differenti da L1: condividere le worklist genererebbe una selezione errata dei flussi e dei CFL anche prima di introdurre la pressione interlayer.

**Migrazione:** `simulation_context_type%domains(:)` è il solo owner delle worklist per layer; la dimensione segue `model_layout%n_layers`. La coppia `domains(i)`/`spatial_operators(i)` conserva scratch indipendenti e viene inizializzata/finalizzata simmetricamente. Il driver, i reader/restart, lo stochastic, il passo Lie, le diagnostiche e il calcolo del CFL continuano a usare soltanto `domains(1)`. Il test-only observer `TEST_PCCU_GATE_H/build_dt_observer.py` segue lo stesso riferimento, senza cambiare algoritmo.

**Regressioni predisposte:** `TEST_SPATIAL_OPERATOR` usa due `domain_type`, due `spatial_operator_type` e uno stato temporaneo rank-4 a due slot. Verifica CFL e residuo identici per stati identici, quindi cambia soltanto l'inventario di slot 2, ricalcola la sua worklist ed esercita il kernel esistente con condizioni dry/wet differenti. Accerta che maschere e stato L1 rimangano invariati, che il CFL sia positivo/finito e che le celle di soluzione abbiano residui finiti. `TEST_WORKSPACE_LIFECYCLE` esercita l'indipendenza e la deallocazione ripetuta di entrambi gli array di oggetti. I test non dichiarano valida la differente fisica del futuro layer superiore.

**Verifica:** compilazione/esecuzione dell'intero repository e regressioni sul futuro commit non ancora effettuate. La prova isolata con `gfortran` di array allocatable di tipi derivati e procedure di inizializzazione/finalizzazione è riuscita, ma non sostituisce la build o i benchmark. Rimangono differiti gli EXAMPLES come richiesto. `N_LAYERS=2` resta esplicitamente rifiutato in produzione e G3/G4 non sono chiusi.

## TL07/TL11 — contatti diretti, senza nuove sovrastrutture (14 novembre 2026)

I quattro casi fisici (entrambi presenti, nessuno, solo L1, solo L2)
sono già rappresentati interamente in `layer_contacts_2d` e testati in
`TEST_SOURCE_REGISTRY`. `source_terms_2d` resta l'assemblatore locale:
quando TL07 fornirà entrambi i candidati di colonna, chiamerà direttamente
`evaluate_layer_contacts` e, se necessario, disabiliterà i bit dei
processi incompatibili nel `source_mask` **locale** già esistente.
Nessuna nuova routine di filtraggio o seconda rappresentazione dei
contatti. `lie_source_executor_2d` mantiene `j,k`, gestione dei
candidati, transazioni e rollback senza decidere quali leggi fisiche
attivare. Nessun kernel di `mass_exchange` deve accedere a indici di
cella o essere riscritto.

**Correzione del commit `e4ef157`:** il selettore aggiunto allora,
`filter_contact_mass_exchange_mask`, è stato rimosso con i suoi test
dedicati e la dipendenza extra in `src/Makefile.am`. Rimangono
invariati il classificatore e le verifiche dei quattro stati che
esistevano prima di quel commit. Questa correzione elimina soltanto
codice non utilizzato nel percorso produttivo: non attiva L2 e non
modifica rate o formule single-layer.

**Lavori ancora aperti:** TL05-C (trasporto L2 fisicamente coerente),
TL07 (candidati accoppiati e rollback atomico), TL11 (applicazione
diretta delle esposizioni fisiche e test dei rate/budget) e le
dipendenze di alcuni kernel da `gas_flag`/`liquid_flag` globali.
Il commit di semplificazione necessita di verifica statica;
nessuna CI nuova è implicita nella pubblicazione.

## TL05-C1 — riduzione della pressione L2 nella closure HP-PCCU condivisa (14 novembre 2026)

**Audit delle dipendenze:** `simulation_2d` possiede gia' `spatial_operators(:)` e lo inizializza con un indice per slot; `spatial_operator_2d` incapsula `hyperbolic_workspace_type`. La pressione idrostatica del backend e' calcolata da coefficienti `gamma_interface` prodotti da `eval_interface_speed_bounds`, riutilizzati nei percorsi di faccia e di cella; `hydrostatic_arithmetic_scale` usa gli stessi coefficienti per la cancellazione dei residui a riposo. `reconstruction_2d`, `pccu_2d`, conversioni REAL/COMPLEX e quadri di stato rimangono inalterati. Il manoscritto `IMEX_v3.pdf` (sez. 2.1.3 e 2.11) prescrive `max(gprime_i,0)` nei flussi idrostatici e negli autovalori, senza curvatura in pressione.

**Implementazione:** il solo operatore HP-PCCU gia' esistente riceve l'indice canonico L1/L2 in fase di inizializzazione. Per L2 `hyperbolic_2d` azzera i coefficienti `Gamma` negativi nei trace gamma usati dal CFL, nei termini di percorso ricostruiti localmente e nel limite numerico dei residui a riposo. Per L1 il flag rimane falso e tutti i calcoli restano invariati. Le velocita' caratteristiche attuali hanno gia' celerita' nulla se `gprime*h<0`; nessuna nuova formula o struttura dei flussi. `TEST_SPATIAL_OPERATOR` esegue un caso gassoso buoyant L2-only a riposo e uno dense gas-solid L2 in moto, conservando il controllo bitwise del single-layer e l'isolamento del dominio L1. L'operatore pubblico mantiene il default L1 e non viene abilitata la modalita' two-layer in produzione.

**Verifica:** controlli statici di interfacce, Fortran 120 colonne, dipendenze e copertura dei path; le due Actions devono essere eseguite e lette sul commit effettivo. Nessuna build locale del repository e nessuna prova EXAMPLE dichiarate. Rimangono aperti TL05-C (trasporto L2 completo), TL06 (pressione accoppiata e CFL due-layer) e TL07/TL11 (candidati e sorgenti secondo i contatti).
