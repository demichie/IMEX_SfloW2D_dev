# IMEX-SfloW2D — Audit tecnico del branch two-layer

**Verifica statica:** 12 novembre 2026. **Repository:** `demichie/IMEX_SfloW2D_dev`; branch `feature/two-layer`; baseline CI termica confermata `b9e29e7814421cbd07a4eaa757ea065b8738e1f7`; TL04-T3 preparata, test ancora da eseguire sul nuovo commit. **Piano vigente:** [roadmap two-layer](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md). Questo documento registra riscontri correnti e lacune, non una cronologia delle vecchie architetture.

## Evidenze verificate

- Le modifiche termiche sono state preparate esclusivamente su `feature/two-layer`; nessun merge in `master` è autorizzato o effettuato in questa fase.
- [Build optimized/debug #43295920602](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43295920602) e [Cell physics/numerical regressions #43296214969](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43296214969): entrambi **completed/success** su **`b9e29e78`**; nessuna CI del nuovo TL04-T3 è ancora certificata.
- L'audit attuale copre dichiarazioni e ownership termica, parser, dispatcher dei kernel, `source_terms_2d`, Newton REAL/COMPLEX, Lie, test diretti, Makefile e documentazione; la reologia resta fuori dall'intervento TL04-T3.
- Nessuna compilazione/esecuzione è stata lanciata da questo audit; gli esiti CI sono evidenze **esterne** limitate ai job eseguiti.

## Mappa corrente delle dipendenze

| ID | Rilevazione sul codice attuale | Lavoro aperto / vincolo |
|---|---|---|
| **A01 — Layout** | `src/core/model_layout_2d.f90:supported_layer_count` ritorna `n_layers==1`; `src/state/state_fields.f90` alloca `q/qp` rank-4. `src/numerics/time_integration_2d.f90` alloca trial rank-4 ma seleziona `(:,1,:,:)` in chiamate Lie e commit. | Test di layout e funzioni L2 isolate; non attivare `N_LAYERS=2` per eludere i gate. |
| **A02 — Config reologica** | `rheology_parameters_2d.f90` definisce unico `rheology_config_type`, `rheology_l1/l2`; `inpout_2d.f90` usa `TYPE(...),POINTER :: rheology => rheology_l1` per legacy; dispatcher reologico e mass exchange ricevono config esplicita; `parameters_2d` non ha più globale `rheology_model`; `source_input_2d` non ha più `rheology_method/group`. | `source_terms_2d` utilizza ancora `rheology_l1` come selezione single-layer, da rendere layer-aware con la migrazione del chiamante; non reintrodurre bridge globali. |
| **A03 — Contesto fisico** | `simulation_2d.f90` possiede `layer_physics(:)`, inizializza L1 e congela per valore la config in `layer_physics_2d.f90`. | Contesti termici L1/L2 e chiusure; snapshot isolati, read-only e OpenMP-safe. |
| **A04 — Namelist termica** | `thermal_config_type` contiene i coefficienti e `T_ground`; `thermal_l1/l2` sono istanze canoniche e il legacy `thermal => thermal_l1` ne è un alias. Parser single-layer verificato in CI su `b9e29e78`; validazione L1/L2 del terreno predisposta, ma produzione L2 disabilitata. | Ricertificare sul commit TL04-T3 input, esempi e comportamento fisico senza mutare namelist o soglie. |
| **A05 — Temperatura atmosferica** | `T_ambient` è l'unica temperatura atmosferica letta dalla namelist GAS; `T_env` non esiste più nella produzione. TL04-T1 è stata verificata in CI. | Restano verifiche integrate del delta termico autorizzato e degli esempi, in particolare lava. |
| **A06 — Termica fisica** | `thermal_cooling_rate.inc` è invariato; coefficienti radiativo/convettivo ricavati con funzioni pure da `thermal_l1`, senza copie globali, e atmosfera da `T_ambient`; REAL/COMPLEX/analitico conservati. | Nessuna nuova legge terreno/interlayer; verificare equivalenza single-layer e percorsi di input lava prima di chiudere il gate. |
| **A07 — Contatti** | Non esiste ancora modulo condiviso `layer_contacts_2d` nel tree; il Lie non conosce scope. | Nuovo modulo `src/physics/layer_contacts_2d.f90` puro e cell-local, contatti atmosfera/terreno/coppia da spessori candidati e wet/dry canonico; niente logica fisica nel Lie; non modificare gating single-layer. |
| **A08 — Registro/integrazione Lie** | `source_registry_2d.f90:source_group_type` contiene id, method, source_names, equation_ids, cache; nessun `scope`. `lie_source_executor_2d.f90:execute_lie_source_groups` riceve array rank-3 di L1 ed esegue gruppi cell-outer. | Scope L1/L2/INTERLAYER, coppia candidati e trasferimento unico, rollback atomico; kernel monolayer ricevono solo q/qp e rate del layer bersaglio. |
| **A09 — HP-PCCU/path** | `nonconservative_2d.f90` è generico; `pccu_2d`, `hyperbolic_2d` hanno ancora fisica single-layer. | Path pressione accoppiata con `B,eta_1,eta_tot`, integrando completo, clamp, identità somma, CFL coupled e iperbolicità. |
| **A10 — Deposito e substrato** | `deposit_tot` canonico e ramo lava/granulare distinti; scratch e transazioni del Lie aggiornati. | Prove end-to-end (lava `n_solid=0`, erosione, topografia, restart v7, NetCDF) ancora non certificate dall'evidenza disponibile. |
| **A11 — Specie e processi** | `constitutive_parameters_2d.f90` contiene proprietà intrinseche e parametri globali di altre chiusure; `source_input_2d.f90` conserva metadati di pore pressure e mass exchange, **non** quelli termici ora in `thermal_config_type`. | Distinguere SHARED, L1, L2, INTERLAYER; nessuna duplicazione delle proprietà delle specie; chiusure del manoscritto da implementare con valida ownership. |
| **A12 — I/O** | NetCDF single-layer `(x,y,time)`, restart v7 con accesso L1 e substrato condiviso, formati legacy. | Versione two-layer, schema del layer, continuazione e test round-trip; preservare produzione mono. |
| **A13 — Workflow** | Due workflow success sulla baseline `b9e29e78`; T3 non ha ancora evidenza di compilazione/regressione. | Non attribuire esiti a EXAMPLES completi, restart/NetCDF, OpenMP o benchmark che non risultano eseguiti. |

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
| G3 configurazioni L1/L2 | Reologia, TL04-T1/T2 verificate in CI; migrazione `T_ground` da certificare, contatti e T3 ancora aperti | **Parziale: aperto** |
| G4 spazio | No coupled production | **Aperto** |
| G5 Lie coppia | No scope, no candidate pair production | **Aperto** |
| G6 chiusure | Da implementare | **Aperto** |
| G7 produzione | No full coupled benchmark | **Aperto** |

**Esito:** la preparazione strutturale (rank-4 e configurazione reologica) è reale e verificata parzialmente; **non** dimostra il two-layer. Non sono stati modificati codici fisici o risultati per questo audit. Qualunque nuovo commit richiede le sue prove. Il piano di lavoro aggiornato e le milestone si trovano **solo** nella roadmap collegata.
