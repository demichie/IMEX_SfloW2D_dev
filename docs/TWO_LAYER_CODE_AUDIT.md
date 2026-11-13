# IMEX-SfloW2D — Audit tecnico del branch two-layer

**Verifica statica:** 12 novembre 2026. **Branch:** `feature/two-layer`; baseline CI `bef39b1c` completata con successo; riferimento canonico termico nel contesto per-layer verificato in CI. **Piano vigente:** [roadmap two-layer](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md).

## Evidenze verificate

- Le modifiche termiche sono state preparate esclusivamente su `feature/two-layer`; nessun merge in `master` è autorizzato o effettuato in questa fase.
- [Build optimized/debug #43306145842](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43306145842) e [regressioni #43306324027](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43306324027) **success** su `bef39b1c`. La successiva migrazione delle altre configurazioni richiederà una nuova CI.
- L'audit attuale copre dichiarazioni e ownership termica, parser, dispatcher dei kernel, `source_terms_2d`, Newton REAL/COMPLEX, Lie, test diretti, Makefile e documentazione; la reologia resta fuori dall'intervento TL04-T3.
- Nessuna compilazione/esecuzione è stata lanciata da questo audit; gli esiti CI sono evidenze **esterne** limitate ai job eseguiti.

## Mappa corrente delle dipendenze

| ID | Rilevazione sul codice attuale | Lavoro aperto / vincolo |
|---|---|---|
| **A01 — Layout** | `src/core/model_layout_2d.f90:supported_layer_count` ritorna `n_layers==1`; `src/state/state_fields.f90` alloca `q/qp` rank-4. `src/numerics/time_integration_2d.f90` alloca trial rank-4 ma seleziona `(:,1,:,:)` in chiamate Lie e commit. | Test di layout e funzioni L2 isolate; non attivare `N_LAYERS=2` per eludere i gate. |
| **A02 — Config reologica** | `rheology_parameters_2d.f90` definisce `rheology_config_type` e `rheology_l1/l2`; il parser usa l'alias legacy verso L1. I kernel reologici `rheology_lava_2d`, `rheology_voellmy_2d`, `rheology_mui_2d` e il dispatcher `rheology_2d` ricevono già `config` con `INTENT(IN)`. `parameters_2d` non possiede più il vecchio `rheology_model` globale. | TL04-R: `source_terms_2d` riceve ora una configurazione esplicita lungo Lie, Newton REAL/COMPLEX, percorso analitico, semi-implicito, esplicito e mass exchange; `layer_physics_type` referenzia il dato canonico e ne espone una vista per valore. Il nuovo commit richiede build e regressioni. Nessuna nuova legge e nessuna attivazione L2. |
| **A03 — Contesto fisico** | `simulation_2d` possiede `layer_physics(:)` con fase/reologia validate; `layer_physics_type` mantiene ora un riferimento non proprietario alla configurazione termica canonica. `time_integration_2d` riceve solo la vista read-only selezionata dal contesto, senza dipendere da `thermal_l1`. | CI del descriptor verificata su `bef39b1c`; preparare scope L2 e le altre configurazioni senza duplicati sincronizzati. |
| **A04 — Namelist termica** | `thermal_config_type` contiene i coefficienti e `T_ground`; `thermal_l1/l2` sono istanze canoniche e il legacy `thermal => thermal_l1` ne è un alias. Parser single-layer verificato in CI su `b9e29e78`; validazione L1/L2 del terreno predisposta, ma produzione L2 disabilitata. | Ricertificare sul commit TL04-T3 input, esempi e comportamento fisico senza mutare namelist o soglie. |
| **A05 — Temperatura atmosferica** | `T_ambient` è l'unica temperatura atmosferica letta dalla namelist GAS; `T_env` non esiste più nella produzione. TL04-T1 è stata verificata in CI. | Restano verifiche integrate del delta termico autorizzato e degli esempi, in particolare lava. |
| **A06 — Termica fisica** | TL04-T3 usa configurazione read-only dal coordinamento single-layer fino al dispatcher termico, senza duplicazioni di coefficienti; cooling REAL/COMPLEX/analitico invariato nelle formule. | Validare esempi integrati lava/dry e futuro scope L2. |
| **A07 — Contatti** | `src/physics/layer_contacts_2d.f90` introduce classificazione pura cell-local su 1/2 spessori candidati, usando `dry_thickness_tolerance`; test integrati in `TEST_SOURCE_REGISTRY` per tutti i casi e input invalidi. Il Lie non importa il modulo. | Certificato in CI su `4e9a591a`, integrare solo nei futuri operatori TL07/TL11; non modificare il gating single-layer. |
| **A08 — Registro/integrazione Lie** | `source_registry_2d.f90:source_group_type` contiene id, method, source_names, equation_ids, cache; nessun `scope`. `lie_source_executor_2d.f90:execute_lie_source_groups` riceve array rank-3 di L1 ed esegue gruppi cell-outer. | Scope L1/L2/INTERLAYER, coppia candidati e trasferimento unico, rollback atomico; kernel monolayer ricevono solo q/qp e rate del layer bersaglio. |
| **A09 — HP-PCCU/path** | `nonconservative_2d.f90` è generico; `pccu_2d`, `hyperbolic_2d` hanno ancora fisica single-layer. | Path pressione accoppiata con `B,eta_1,eta_tot`, integrando completo, clamp, identità somma, CFL coupled e iperbolicità. |
| **A10 — Deposito e substrato** | `deposit_tot` canonico e ramo lava/granulare distinti; scratch e transazioni del Lie aggiornati. | Prove end-to-end (lava `n_solid=0`, erosione, topografia, restart v7, NetCDF) ancora non certificate dall'evidenza disponibile. |
| **A11 — Specie e processi** | `constitutive_parameters_2d.f90` contiene proprietà intrinseche e parametri globali di altre chiusure; `source_input_2d.f90` conserva metadati di pore pressure e mass exchange, **non** quelli termici ora in `thermal_config_type`. | Distinguere SHARED, L1, L2, INTERLAYER; nessuna duplicazione delle proprietà delle specie; chiusure del manoscritto da implementare con valida ownership. |
| **A12 — I/O** | NetCDF single-layer `(x,y,time)`, restart v7 con accesso L1 e substrato condiviso, formati legacy. | Versione two-layer, schema del layer, continuazione e test round-trip; preservare produzione mono. |
| **A13 — Workflow** | Le due Actions sono **success** su `bef39b1c`, inclusa la modifica del descriptor termico. | Ricertificare ogni nuovo commit con build optimized/debug e regressioni; eseguire full EXAMPLES quando disponibile. |

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
| G2 rank-4 single-layer | Build optimized/debug e regressioni fisiche/numeriche **success** su `b9e29e78`; full EXAMPLES e restart ancora da certificare | **Parziale: aperto** |
| G3 configurazioni L1/L2 | Modulo contatti e thermal descriptor superati in CI su `bef39b1c`; altri scope aperti | **Parziale: aperto** |
| G4 spazio | No coupled production | **Aperto** |
| G5 Lie coppia | No scope, no candidate pair production | **Aperto** |
| G6 chiusure | Da implementare | **Aperto** |
| G7 produzione | No full coupled benchmark | **Aperto** |

**Esito:** la preparazione strutturale (rank-4 e configurazione reologica) è reale e verificata parzialmente; **non** dimostra il two-layer. Non sono stati modificati codici fisici o risultati per questo audit. Qualunque nuovo commit richiede le sue prove. Il piano di lavoro aggiornato e le milestone si trovano **solo** nella roadmap collegata.

## TL04-C — audit delle configurazioni residue (12 novembre 2026)

**Baseline tecnica:** `feature/two-layer` a `bef39b1c`; nessuna modifica fisica o numerica in questo audit. La verifica CI del passaggio reologico TL04-R precede la migrazione Darcy L1.

| Famiglia | Owner e consumer identificati | Decisione prima della migrazione |
|---|---|---|
| Pore pressure (L1) | `constitutive_parameters_2d` conserva `hydraulic_permeability`, `dynamic_permeability_flag`, `alpha_trans`, `N_inh`, `f_inhibit_mode`; `source_input_2d` conserva METHOD/GROUP; `parameters_2d` conserva flag e `pore_pres_fract`. `inpout_2d` legge/resetta `PORE_PRESSURE_PARAMETERS` nei percorsi input; `source_configuration_2d` seleziona metodo/gruppo; `pore_pressure_darcy_2d` usa i coefficienti nel percorso analitico e nei fattori di diffusione. | Migrare come **singola configurazione canonica L1**, non creare L2 per simmetria. Valutare separatamente le quantità condivise (viscosità, specie, pressione ambiente). Mantenere `pascal_coeff` come cache derivata con ownership e lifetime espliciti; non duplicarne valori in un descriptor. Verificare tutti i percorsi di lettura/reset prima di cambiare NAMELIST. |
| Proprietà SHARED | `constitutive_parameters_2d` contiene specie solide/gas, proprietà ambientali e costanti; `inpout_2d` le legge da `SOLID_TRANSPORT_PARAMETERS`/`GAS_TRANSPORT_PARAMETERS`. | Non trasferire proprietà intrinseche nei layer; distinguere le opzioni di singole chiusure dai dati condivisi prima di definire nuovi tipi. |
| Mass exchange | `source_input_2d` legge `MASS_EXCHANGE_PARAMETERS`, che oggi aggrega `erosion_coeff`, `settling_flag`, `loss_rate`, `entrainment_flag`, gas loss, vaporization e deposizione lava; `source_configuration_2d` applica METHOD/GROUP comuni. | Non creare un unico tipo per-layer che mescoli terreno, atmosfera e scambio interlayer. Separare le famiglie solo nei rispettivi gate TL07–TL11; conservare il percorso lava single-layer. |

**Verifiche obbligatorie per la successiva migrazione pore pressure:** chiusura completa dei `USE` e delle chiamate dirette, metodi analitici e rifiuto dei metodi non supportati, generatori/fixture `EXAMPLE_2D_PP` e `EXAMPLE_DAMBREAK`, reset/restart, build optimized/debug, test `TEST_RHEOLOGY_CELL_LOCAL/test_pore_pressure_diffusion.f90`, suite numeriche e confronto dei risultati single-layer. Non modificare formule, soglie, floating-point order o il flag `N_LAYERS=2`.

**Limite dell'audit:** la ricerca GitHub indicizza il branch di default; le dichiarazioni e i frammenti dei file citati sono stati letti esplicitamente su `feature/two-layer`. Prima di effettuare sostituzioni resta necessaria una ricerca completa della dependency closure sul checkout del branch e la verifica con compilatore. Questo audit non certifica ancora tale chiusura.
