# IMEX-SfloW2D — Audit tecnico del branch two-layer

**Verifica statica:** 12 novembre 2026. **Repository:** `demichie/IMEX_SfloW2D_dev`; branch `feature/two-layer`; HEAD ispezionato `410ad99e5c7d3071df4d4cd0cd2d8161c8ea13d1`. **Piano vigente:** [roadmap two-layer](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md). Questo documento registra **solo riscontri correnti e lacune**, non una cronologia delle vecchie architetture.

## Evidenze verificate

- Il confronto GitHub `master...feature/two-layer` mostra **98 commit avanti, zero indietro** all'HEAD ispezionato. Non è stato effettuato alcun merge in `master` nell'ambito di questo audit.
- [Build optimized/debug #43264995583](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43264995583) e [Cell physics/numerical regressions #43265344609](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43265344609): entrambi **completed/success** su **`410ad99e`**.
- L'audit ha ispezionato HEAD e commit del branch, inclusi `0e44f407` (rimozione globali residui), `410ad99e` (parser corretto), `48c6969` (alias L1), `13bd2bd` (lava/lahar), `5be747c` (lava/deposito), `80d12b0` (mu(I)), `6caffa9` (Froude/Lucas), `d4f434c` (Voellmy), `ee4d0a6` (snapshot contesto).
- Nessuna compilazione/esecuzione è stata lanciata da questo audit; gli esiti CI sono evidenze **esterne** limitate ai job eseguiti.

## Mappa corrente delle dipendenze

| ID | Rilevazione sul codice attuale | Lavoro aperto / vincolo |
|---|---|---|
| **A01 — Layout** | `src/core/model_layout_2d.f90:supported_layer_count` ritorna `n_layers==1`; `src/state/state_fields.f90` alloca `q/qp` rank-4. `src/numerics/time_integration_2d.f90` alloca trial rank-4 ma seleziona `(:,1,:,:)` in chiamate Lie e commit. | Test di layout e funzioni L2 isolate; non attivare `N_LAYERS=2` per eludere i gate. |
| **A02 — Config reologica** | `rheology_parameters_2d.f90` definisce unico `rheology_config_type`, `rheology_l1/l2`; `inpout_2d.f90` usa `TYPE(...),POINTER :: rheology => rheology_l1` per legacy; dispatcher reologico e mass exchange ricevono config esplicita; `parameters_2d` non ha più globale `rheology_model`; `source_input_2d` non ha più `rheology_method/group`. | `source_terms_2d` utilizza ancora `rheology_l1` come selezione single-layer, da rendere layer-aware con la migrazione del chiamante; non reintrodurre bridge globali. |
| **A03 — Contesto fisico** | `simulation_2d.f90` possiede `layer_physics(:)`, inizializza L1 e congela per valore la config in `layer_physics_2d.f90`. | Contesti termici L1/L2 e chiusure; snapshot isolati, read-only e OpenMP-safe. |
| **A04 — Namelist termica** | `inpout_2d.f90` definisce ancora `TEMPERATURE_PARAMETERS` con scalari globali (`thermal_cooling_flag`, emissività, conduzione, `exp_area_fract`, `c_p`, `enne`, `emme`, `T_env`, method/group). | Tipo derivato `thermal_config_type` e input legacy/L1/L2; audit proprietà e test. |
| **A05 — Temperatura atmosferica** | `T_env` è dichiarata in `constitutive_parameters_2d.f90`, default 300 K e uso in `source_terms_2d`; `T_ambient` è presente in GAS, INITIAL e RESTART, e nell'inizializzazione degli stati. | Conservare **solo `T_ambient` letta obbligatoriamente da GAS**, anche `gas_flag=F`; eliminare `T_env` e le altre sorgenti d'input. Delta fisico esplicitamente autorizzato ma da misurare. |
| **A06 — Termica fisica** | `thermal_cooling_rate.inc`: `rate=-radiative*(T**4-ambient**4)-convective*(T-ambient)`; kernel REAL/COMPLEX e analitico ricevono `ambient` esplicito. | Usare `T_ambient`; non introdurre legge di terreno o termica interlayer; verificare tutti i trattamenti e il default lava. |
| **A07 — Contatti** | Non esiste ancora modulo condiviso `layer_contacts_2d` nel tree; il Lie non conosce scope. | Nuovo modulo `src/physics/layer_contacts_2d.f90` puro e cell-local, contatti atmosfera/terreno/coppia da spessori candidati e wet/dry canonico; niente logica fisica nel Lie; non modificare gating single-layer. |
| **A08 — Registro/integrazione Lie** | `source_registry_2d.f90:source_group_type` contiene id, method, source_names, equation_ids, cache; nessun `scope`. `lie_source_executor_2d.f90:execute_lie_source_groups` riceve array rank-3 di L1 ed esegue gruppi cell-outer. | Scope L1/L2/INTERLAYER, coppia candidati e trasferimento unico, rollback atomico; kernel monolayer ricevono solo q/qp e rate del layer bersaglio. |
| **A09 — HP-PCCU/path** | `nonconservative_2d.f90` è generico; `pccu_2d`, `hyperbolic_2d` hanno ancora fisica single-layer. | Path pressione accoppiata con `B,eta_1,eta_tot`, integrando completo, clamp, identità somma, CFL coupled e iperbolicità. |
| **A10 — Deposito e substrato** | `deposit_tot` canonico e ramo lava/granulare distinti; scratch e transazioni del Lie aggiornati. | Prove end-to-end (lava `n_solid=0`, erosione, topografia, restart v7, NetCDF) ancora non certificate dall'evidenza disponibile. |
| **A11 — Specie e processi** | `constitutive_parameters_2d.f90` contiene proprietà intrinseche e parametri globali di chiusure; `source_input_2d.f90` conserva method/group globali per termica, pore pressure e mass exchange. | Distinguere SHARED, L1, L2, INTERLAYER; nessuna duplicazione delle proprietà delle specie; chiusure del manoscritto da implementare con valida ownership. |
| **A12 — I/O** | NetCDF single-layer `(x,y,time)`, restart v7 con accesso L1 e substrato condiviso, formati legacy. | Versione two-layer, schema del layer, continuazione e test round-trip; preservare produzione mono. |
| **A13 — Workflow** | Due workflow success sull'HEAD audit; i test di unità e regressione esistenti coprono una parte della fisica e della numerica. | Non attribuire esiti a EXAMPLES completi, restart/NetCDF, OpenMP o benchmark che non risultano eseguiti. |

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
| G2 rank-4 single-layer | Build optimized/debug e regressioni fisiche/numeriche **success** su `410ad99e` | **Parziale: aperto** |
| G3 configurazioni L1/L2 | Reologia migrata e CI success; termica/contatti/registri L2 mancanti | **Parziale: aperto** |
| G4 spazio | No coupled production | **Aperto** |
| G5 Lie coppia | No scope, no candidate pair production | **Aperto** |
| G6 chiusure | Da implementare | **Aperto** |
| G7 produzione | No full coupled benchmark | **Aperto** |

**Esito:** la preparazione strutturale (rank-4 e configurazione reologica) è reale e verificata parzialmente; **non** dimostra il two-layer. Non sono stati modificati codici fisici o risultati per questo audit. Qualunque nuovo commit richiede le sue prove. Il piano di lavoro aggiornato e le milestone si trovano **solo** nella roadmap collegata.
