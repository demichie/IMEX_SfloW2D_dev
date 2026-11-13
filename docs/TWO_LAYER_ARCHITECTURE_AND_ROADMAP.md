# IMEX-SfloW2D — Architettura e roadmap two-layer

**Plan updated:** 13 November 2026. **Branch:** `feature/two-layer`. **Verified CI baseline:** [`672e0a0`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/672e0a06d6465dd6f2b3e5f9e91ef7ccc1a8e36f) (build optimized/debug and cell-physics/numerical regressions successful). **TL04-C3/C4:** canonical combined mass-exchange configuration and propagation CI-verified; scientific gate C5/G3 remains open pending the deferred EXAMPLES and extended physical checks. **TL05 preparation:** spatial operator workspaces CI-verified; independent per-slot active-domain worklists and isolated numerical wet/dry tests are the next step. No two-layer execution or coupled HP-PCCU.
**Stato:** codice di produzione ancora **single-layer**; `N_LAYERS=2` resta disabilitato. La presente roadmap è il **solo piano vigente**: specifica, stato verificato, attività, gate e milestone. La cronologia dettagliata è nei commit Git, non in versioni alternative del piano.
**Documenti collegati:** [audit tecnico corrente](TWO_LAYER_CODE_AUDIT.md), [roadmap generale](IMEX_SfloW2D_Development_Roadmap.md), manoscritto *IMEX_SfloW2D_v3*.

## 1. Obiettivo, invarianti e ambito

Il modello ha due sole modalità: `N_LAYERS=1` (solutore attuale) e, **solo dopo il gate finale**, `N_LAYERS=2` con L1 corrente inferiore densa e L2 nube superiore diluita nello stesso scenario gas/solidi. Non è previsto un framework con numero arbitrario di layer, né il ripristino dell'IMEX-RK storico. Le combinazioni lava/granulare tra i due layer non sono abilitate dal semplice fatto che il single-layer supporta la lava.

1. **Una sola rappresentazione autorevole** di configurazioni, proprietà fisiche e inventari; nessun globale duplicato, copia sincronizzata, wrapper transitorio o duplicazione di formule.
2. **Equivalenza fisica e numerica single-layer, con due nomi NAMELIST ammessi per L1:** preservare operatori, formule, ordine Lie, gruppo 0, wet/dry, rollback, OpenMP, deposito e risultati, salvo modifiche fisiche esplicitamente approvate. Le famiglie di configurazione per-layer accettano in `N_LAYERS=1` sia `<FAMILY>_PARAMETERS` sia `<FAMILY>_L1_PARAMETERS`, ma mai entrambi nello stesso input; entrambi accedono direttamente allo stesso oggetto canonico L1. È già approvato l'uso esclusivo di `T_ambient` invece di `T_env`, anche se cambia i risultati termici.
3. **Separazione dei ruoli:** integratore Lie = algoritmi, ordine, timestep, candidati, controlli numerici e rollback; moduli `physics` = leggi, sorgenti, decisioni fisiche sui contatti; `spatial` = flussi e termini di percorso. Nessuna fisica introdotta nel Lie.
4. `N_LAYERS=2` rifiutato negli input e nel layout di produzione fino a TL13/M6. Le funzionalità parziali sono testabili tramite harness isolati, non selezionabili come modello completo.
5. Non cambiare formule, regolarizzazioni, limitatori, splitting o tolleranze implicitamente durante un refactoring. Rispettare disponibilità esplicita delle modalità EXPLICIT, SEMI_IMPLICIT, NEWTON/IMPLICIT e ANALYTIC; mai fallback silenziosi.
6. **Un solo oggetto autorevole per layer, con convenzione NAMELIST single-layer consolidata:** per reologia, termica, Darcy e mass exchange, `N_LAYERS=1` accetta `<FAMILY>_PARAMETERS` oppure `<FAMILY>_L1_PARAMETERS` in modo mutuamente esclusivo; il nome breve legge direttamente l'istanza L1 canonica tramite alias di input, non attraverso una copia o una traduzione. In `N_LAYERS=2` sono obbligatori `<FAMILY>_L1_PARAMETERS` e `<FAMILY>_L2_PARAMETERS` e il nome breve è vietato. Aggiornare insieme EXAMPLES/TESTS, parser, chiamanti e documenti, eliminando vecchi coefficienti scalari, wrapper di conversione e interfacce interne obsolete. Non mantenere altre compatibilità di sintassi o copie sincronizzate.

## 2. Audit corrente, regressioni e limiti

| Area | Implementato e riscontrato nel codice | Ancora da implementare o certificare |
|---|---|---|
| Layout | `model_layout_2d` contiene `n_layers`; `state_type%q/qp` e workspace di integrazione sono rank-4 `(componente,layer,x,y)`; i consumer di produzione selezionano esplicitamente L1. | `supported_layer_count` accetta solo 1; evoluzione L2, maschere inattive e test coupled. |
| Deposito | `deposit_tot` è l'inventario canonico dello spessore; ramo lava distinto dal granulare, senza inventario autonomo `lava_deposit`; rollback e percorsi v7 sono stati adeguati. | Gate G1 completo: casi granulari/lava, `n_solid=0`, topografia on/off, I/O e restart round-trip documentati. |
| Rheology | Canonical L1/L2 rheology and explicit REAL/COMPLEX/analytic source propagation are implemented and CI-verified on `23beb9c`. | L2 production scopes remain disabled; extended scientific gates are open. |
| Contesto per layer | `simulation_2d%layer_physics(:)` conserva fase e reologia; il contesto referenzia ora la configurazione termica canonica senza una seconda copia persistente. La vista read-only viene passata al passo Lie. | Completare le altre configurazioni SHARED/L1/L2 e gli scope fisici solo nei gate dedicati; nessuna sincronizzazione fra copie. |
| Lie | Un solo Lie primo ordine, cell-outer/OpenMP, gruppi numerici ordinati, trial e rollback; chiamata corrente `execute_lie_source_groups(...q_trial(:,1,:,:),qp_trial(:,1,:,:)...)`. | Supportare gruppi con scope `L1`, `L2`, `INTERLAYER`, candidati entrambi i layer e commit atomico senza fisica nel driver. |
| Thermal | Canonical L1/L2 cooling configurations and explicit kernel propagation are implemented; `T_ground` has moved out of thermal configs into the single shared `substrate%T_ground`. The accepted single-layer NAMELIST names both resolve to L1. | Extended single-layer thermal verification remains open; no new thermal law is activated. |
| Substrate | TL04-C1 owns ground temperature, porosity, erodible inventory settings and topography setup in `substrate_parameters_2d`; code, input fixtures and generators were migrated. | Extended scientific checks on terrain energy, inventory/restart and formerly distinct ground temperatures remain part of G3. |
| Ambiente | `T_ambient` è letta esclusivamente da `GAS_TRANSPORT_PARAMETERS`, validata finita e positiva prima degli utilizzi; `T_env` e gli ingressi INITIAL/RESTART sono rimossi. | Certificare parser, casi senza fase gas e delta termico autorizzato sul nuovo commit. |
| Physical sources | The registry remains single-layer; `layer_contacts_2d` is implemented and pure. Combined mass-exchange L1/L2 configurations have one canonical owner; `mass_exchange_2d.f90` receives the selected read-only config and retains shared inventory limiters. CI verified on `32e4e33`. | Complete C5 scientific verification, including deferred EXAMPLES. Activate contact dispatch and L2 source execution only in TL07/TL11. |
| Spaziale | `nonconservative_2d` resta l'integratore di percorso HP-PCCU. I workspace di ricostruzione/flussi/CFL sono per slot (CI verde su `672e0a0`); TL05-B prepara anche `domains(:)` con workset di celle, facce e ricostruzione indipendenti, mantenendo il driver L1. | Certificare le liste per slot con CI, implementare il trasporto L2 isolato con la sua closure fisica e poi `eta_1`/`eta_tot`, pressione accoppiata, CFL coupled, wet/dry e iperbolicità. |
| I/O | NetCDF single-layer `(x,y,time)`, restart v7 with L1 plus shared inventories, current `.asc/.q_2d` outputs and diagnostics. | Define one authoritative schema when I/O changes, update reader/writer/tests together and verify scientifically meaningful single-layer continuation; no superseded-format reader or backward-compatibility requirement. |

**CI evidence:** [Build (optimized and debug) #43503101813](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43503101813) and [Cell physics and numerical regressions #43503392454](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43503392454) both succeeded on [`672e0a0`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/672e0a06d6465dd6f2b3e5f9e91ef7ccc1a8e36f), after TL04-C3/C4 and the spatial-workspace preparation. One additional regression run on that commit was cancelled; a distinct successful run exists. The user deferred EXAMPLES after mass-exchange migration; earlier examples do not certify these changes. NetCDF/restart, OpenMP, conservation and two-layer scientific gates remain open. The subsequent per-slot domain worklists require new CI.

## 3. Dati, proprietà e configurazioni autorevoli

### 3.1 Stato e substrato

Layout di produzione: `q(n_variables,n_layers,nx,ny)`; `qp(n_physical_variables,n_layers,nx,ny)`. Le dimensioni di componenti, primitive ed equazioni hanno contratti distinti. `q(:,layer,j,k)` e `qp(:,layer,j,k)` sono vettori contigui per kernel monolayer. Trial, residui, snapshot fluidi e I/O seguono lo stesso layout. Topografia, geometria, `deposit_tot`, inventari di deposito/erosione e substrato **non** hanno indice layer.

Gli indici delle equazioni sono comuni ai due layer. Le quantità significative soltanto in L1 (es. eccesso di pressione interstiziale) occupano slot inattivi in L2, con maschera coerente in conversioni, trasporto, IC/BC, restart e output. Ogni worklist attiva include celle wet L1 **oppure** wet L2: `dry L1/wet L2` non può scomparire.

`deposit_tot` è lo spessore cumulativo immobilizzato (m), non il cambiamento netto della quota: in granulare è derivato dall'inventario delle classi e dalla porosità; in lava è un accumulo diretto, anche con `n_solid=0`. `dB_dt` rappresenta la modifica geometrica; non deve duplicare l'inventario. Erosione ed erodibile restano distinti. Le leggi granulari e lava sono mutuamente esclusive, anche nel single-layer.

### 3.2 Ownership: L1, L2, INTERLAYER, SHARED

**One authoritative owner per physical quantity; two mutually exclusive NAMELIST names for the single-layer L1 configuration.** Each layer-specific closure has separate canonical L1/L2 configuration; shared solid-species, carrier-phase, atmosphere and ground properties have only one owner; true L1-L2 transfer coefficients belong to the pair. The accepted short and `_L1` group names for rheology, thermal, Darcy and mass exchange must refer directly to the same canonical L1 object, not to snapshots or independently stored coefficients. Do not infer layer ownership from the NAMELIST spelling. Do not activate `N_LAYERS=2` before its scientific gates.

| Family | Canonical target | Phase |
|---|---|---|
| `RHEOLOGY_PARAMETERS` / `RHEOLOGY_L1_PARAMETERS` | Either name reads canonical `rheology_l1` when `N_LAYERS=1`; `_L1` plus `_L2` are required when L2 is enabled, with no short group in two-layer input. | Implemented; L2 scopes later |
| `TEMPERATURE_PARAMETERS` / `TEMPERATURE_L1_PARAMETERS` | Either name reads canonical `thermal_l1` when `N_LAYERS=1`; thermal L1/L2 own cooling laws and METHOD/GROUP. **`T_ground` is now owned only by the shared substrate**. `T_ambient` stays shared in GAS input. | TL04-T implemented; TL04-C1 verified |
| `PORE_PRESSURE_PARAMETERS` / `PORE_PRESSURE_L1_PARAMETERS` | Either name reads canonical Darcy `pore_pressure_l1` when `N_LAYERS=1`; L1/L2 have separate derived Pascal caches, runtime still L1-only. | Verified on `662e796`; L2 runtime remains disabled |
| `SOLID_TRANSPORT_PARAMETERS` | Shared `rho_s`, `inv_rho_s` (derived), `diam_s`, `sphericity_s`, `sp_heat_s` and `maximum_solid_packing` have one canonical owner in `core/solid_transport_parameters_2d.f90`. The input group and values are unchanged; no compatibility aliases or layer copies. | TL04-C2 implemented; CI passed on `662e796` |
| `GAS_TRANSPORT_PARAMETERS` and liquid/lava transport | Shared atmosphere, gas and carrier phase properties (including pressure and Sutherland coefficients). Only actual closure coefficients may move into law configs. | TL04-C0/C3 |
| `SUBSTRATE_PARAMETERS` | **One** `substrate_config_type` owns `T_ground`, `erodible_porosity`, `erodible_file`, initial erodible thickness, `erodible_deposit_flag`, `topography_file` and `topo_change_flag`. The allocatable `erodible_fract(:)` has the same shared module owner but remains separate from the derived type for Fortran NAMELIST support. Retired `T_erodible` and thermal ground copies are removed. | TL04-C1 implemented; CI passed on `662e796` |
| `MASS_EXCHANGE_PARAMETERS` / `MASS_EXCHANGE_L1_PARAMETERS` | In `N_LAYERS=1`, accept exactly one of these names and read the same canonical `mass_exchange_l1` object, without copies or conversion. In `N_LAYERS=2`, require `MASS_EXCHANGE_L1_PARAMETERS` and `MASS_EXCHANGE_L2_PARAMETERS` and reject the short group. `mass_exchange_config_type` owns nested erosion, deposition, entrainment, carrier loss, vaporization and lava immobilization configurations. **One METHOD/GROUP for the combined mass-exchange source of each layer.** | TL04-C3/C4 implemented; CI passed on `32e4e33`, C5 open |
| Injections, new run, IC/BC, numerical and I/O | Explicit source targets where needed; shared geometry, material data and numerical controls. `pore_pres_fract` belongs to injection, not Darcy dissipation. | TL07/TL13 |

**Ground state versus ground configuration:** `substrate_config_type` describes the one terrain, not its evolving fields. Preserve canonical `state_type%erodible`, `deposit`, `erosion`, `deposit_tot`, `erosion_tot`, geometric topography and transactional rollback/restart without layer duplicates. `coeff_porosity` derives from ground porosity; `rho_c_sub` is not input and derives from ground temperature and the carrier phase. Preserve their evaluation semantics. `maximum_solid_packing` is a shared admissibility limit, distinct from deposition hindrance `alpha_max` and future layer-specific packing closures. Audit `alphastot_min` as the current common erosion/deposition threshold; do not introduce inconsistent copies.

**One combined source:** Retain `src/physics/mass_exchange/mass_exchange_2d.f90` as the sole assembler of the active erosion, deposition, entrainment, carrier-loss, vaporization and lava-immobilization processes, with shared carrier/packing limiters and mass, momentum and energy budgets. Retain the current separate `mass_exchange_*_2d.f90` law files; add only `mass_exchange_parameters_2d.f90` for canonical nested configs. Do not create separate authoritative `erosion_l1` or multiple independent Lie integrations. Any future NEWTON or ANALYTIC treatment must solve the *combined* active processes consistently; currently only EXPLICIT mass exchange is supported and other treatments must remain explicitly rejected. `gas_loss_flag` has been moved from mass-exchange input to the canonical Darcy L1/L2 configuration without duplicate storage; the Darcy source reads its supplied configuration.

**Local contacts:** `layer_contacts_2d` already distinguishes both-wet, L1-only, L2-only and both-dry candidates. Bed erosion/deposition applies to the **terrain-exposed** layer; atmospheric entrainment to the **atmosphere-exposed** layer, each using its independent law config. L2-only may contact both boundaries; L1-only likewise. Evaluate contacts from the *current candidate* in physics, not from committed state or Lie control logic. Conservative L2-to-L1 settling (TL08) is not ground deposition; L1-to-L2 elutriation (TL09) may cold-start an absent receiver.

### 3.3 TL04-C implementation steps and acceptance gates

**Execution status (commit `32e4e33`):** TL04-C1/C2 shared ownership and TL04-C3/C4 combined exchange configurations are implemented; both CI workflows passed on `32e4e33`. The current migrated mass-exchange EXAMPLES are deferred by the user, and TL04-C5 / gate G3 remain open. This does not prevent preparatory TL05 workspace ownership, but does not authorize coupled physics or switching on L2. Preserve single-layer equations, evaluation order, Lie source order, physical defaults, wet/dry handling and rollback unless a physical change is explicitly approved. **Input contract:** retain both the short and explicit `_L1` NAMELIST names for single-layer configurations as exclusive alternatives referring to the same canonical L1 instance. Eliminate obsolete scalar readers, translation logic and duplicate storage, not the approved single-layer group names.

1. **TL04-C0 — dependency/input closure audit.** Enumerate owners, declarations, USE sites, nested includes, NAMELIST defaults/reset and readers, source registry, per-process kernels, REAL/COMPLEX/analytical callers, NetCDF/restart, tests, examples and Makefile ordering. Record actual `T_ground` versus `T_erodible` inputs/defaults and the `erodible_deposit_flag` energy branch. Define the definitive per-layer input schema, preserving approved single-layer short/`_L1` alternatives for rheology, thermal, Darcy and mass exchange; enumerate obsolete scalar syntax and conversion paths to remove, not the approved NAMELIST aliases. Audit conflicting old temperature values as a scientific decision, not a reason to support old inputs.
2. **TL04-C1 — SHARED substrate (implemented; CI verified).** `substrate_parameters_2d.f90` owns canonical `substrate_config_type`, `substrate` and the separate allocatable `erodible_fract(:)`; `SUBSTRATE_PARAMETERS` is the definitive input. `T_ground` replaces `T_erodible` and the removed thermal ground copies. Terrain temperature, porosity, erodible-species fractions, raster filename, initial thickness, deposit-to-erodible policy, `topography_file` and `topo_change_flag` have migrated; general run controls remain outside the substrate. Derived `rho_c_sub`/`coeff_porosity` are non-input, and ground setup does not depend exclusively on L1 erosion. All affected example/test inputs and generators were migrated without old readers. The latest CI and user-reported examples passed; detailed thermal-energy comparisons, lava defaults, erodible-deposit semantics, substrate inventory/restart and other G3 checks still require explicit evidence.
3. **TL04-C2 — SHARED solid and phase properties (implemented; CI verified).** `MAXIMUM_SOLID_PACKING` remains in `SOLID_TRANSPORT_PARAMETERS`; `rho_s`, derived `inv_rho_s`, `diam_s`, `sphericity_s`, `sp_heat_s` and `maximum_solid_packing` now have exactly one owner in `core/solid_transport_parameters_2d.f90`. Darcy, mass-exchange, rheology, state-mixture and input consumers import this owner explicitly. Both CI workflows passed on `662e796`; extended scientific consistency with `alpha_trans` remains part of G3. Leave atmospheric/carrier data shared; distinguish process-specific deposition drag parameters `A_drag`, `B_drag`, `collective_settling_flag`, hindrance parameters and vaporization heat fraction from material constants. Preserve the current common `alphastot_min` semantics.
4. **TL04-C3 — per-layer composite mass-exchange configuration.** Introduce `mass_exchange_parameters_2d.f90` with one `mass_exchange_config_type` per layer (`mass_exchange_l1/l2`), owning nested erosion/deposition/entrainment/liquid-loss/vaporization/lava-deposition configs and a **single METHOD/GROUP per layer**. In `N_LAYERS=1`, support **either** `MASS_EXCHANGE_PARAMETERS` (short single-layer name with `mass_exchange%...` fields) **or** `MASS_EXCHANGE_L1_PARAMETERS` (explicit `mass_exchange_l1%...` fields), never both; each reads the **same canonical L1 instance**, e.g. via a NAMELIST pointer alias, without separate coefficients, copy-back or conversion. When `N_LAYERS=2` is enabled, require both `MASS_EXCHANGE_L1_PARAMETERS` and `MASS_EXCHANGE_L2_PARAMETERS` and reject `MASS_EXCHANGE_PARAMETERS`; reject any L2 group while two-layer production is disabled. Retire the old flat scalar entries in `MASS_EXCHANGE_PARAMETERS`, not the approved group name; migrate all EXAMPLES/TESTS inputs, generators and documentation to derived-type components using the two supported single-layer spellings. Check Fortran NAMELIST allocatable-component behavior before committing to the type layout. Preserve current rheology-dependent deposition selection until a separate model-selection decision and validation.
5. **TL04-C4 — reuse the existing combined assembler.** Thread config read-only through `simulation_2d/layer_physics_2d`, source configuration, source terms, Lie/Newton interfaces and the existing `mass_exchange_2d.f90`; each distinct physical law stays in its current file. Keep common inventory limiters and rate aggregation inside this one assembler; no additional coordinator and no separately integrated Lie source per child law. Preserve unsupported NEWTON/ANALYTIC behavior until genuinely coupled algorithms for **all selected processes** exist.
6. **TL04-C5 — verification of single-layer preparation.** Check acceptance of either short or `_L1` single-layer group (but never both), rejection of `_L2` in single-layer and rejection of the short group in two-layer input; reject retired scalar syntax and conversion wrappers; verify both spellings share the exact L1 configuration, plus updated EXAMPLES/TESTS fixtures and generators, SHARED ownership, default temperatures, solid packing and sediment packing separation, source order, wet/dry, mass/momentum/energy and substrate inventory accounting; compile strict/debug and optimized; run `TEST_MASS_EXCHANGE`, `TEST_SOURCE_REGISTRY`, pertinent REAL/COMPLEX suites, `EXAMPLE_TAAL`, `EXAMPLE_TAAL_BOTTOM`, `EXAMPLE_CF`, pore-pressure, lava and liquid examples, plus restart/NetCDF and OpenMP reproducibility where affected. Record evidence against the actual commit; no CI from older commits validates these edits.
7. **TL07/TL11 — deferred layer activation.** Add L1/L2/INTERLAYER scheduling and paired candidate-state/atomic rollback infrastructure first; only then connect `layer_contacts_2d` inside physics to select the applicable layer-specific external exchange. Exercise wet/wet, wet/dry, dry/wet and dry/dry, including L2-only bed erosion/deposition and L1-only atmospheric entrainment; no double counting of common terrain/atmosphere resources. Interlayer transfers remain separate pair operators under TL08–TL10, with coherent global budgets.

### 3.4 Rheology input: current state and required cleanup

**Current and target single-layer contract:** `RHEOLOGY_PARAMETERS` and `RHEOLOGY_L1_PARAMETERS` are intentional, mutually exclusive names for the **same canonical `rheology_l1` object**; `RHEOLOGY_L2_PARAMETERS` remains rejected while two-layer production is disabled. When `N_LAYERS=2` is activated, require explicit `_L1` and `_L2` groups and reject the short name. Preserve the analogous short/`_L1` contract for thermal and pore pressure, and apply it to mass exchange. Keep the NAMELIST pointer aliases where needed for direct shared-object reading, but never add duplicate configuration ownership, sync code or obsolete scalar entries. Preserve the verified physical formulas and numerical kernels.

Il passaggio della configurazione reologica dal contesto di simulazione fino a `source_terms_2d` è stato implementato in TL04-R, senza selezione diretta di `rheology_l1` nei consumer riutilizzabili. La propagazione copre sorgenti esplicite, semi-implicite, implicite REAL/COMPLEX, soluzione analitica e mass exchange reologico, con aggiornamento coerente di chiamanti, test e dipendenze. `rheology_2d` e le singole leggi reologiche (come `rheology_lava_2d`) usano già `TYPE(rheology_config_type), INTENT(IN) :: config`: mantenerle inalterate salvo dipendenze residue dimostrate dall'audit. La configurazione canonica mantiene un solo owner; `layer_physics_type` possiede un riferimento non proprietario e restituisce una vista per valore, senza snapshot persistenti da sincronizzare. Occorre completare le verifiche CI e integrali del nuovo gruppo di modifiche. Il percorso `N_LAYERS=1` deve conservare identici ordine numerico, formule e risultati; la selezione L2 resta inattiva fino ai gate dedicati.

## 4. Contratto termico definitivo e contatti fisici

### 4.1 Namelist e temperatura atmosferica

La configurazione canonica implementata è **un solo** `thermal_config_type` con due istanze `thermal_l1`, `thermal_l2`. Il parser applica:
- **Definitive single-layer contract:** `N_LAYERS=1` accepts either `TEMPERATURE_PARAMETERS` (`thermal%...`) or `TEMPERATURE_L1_PARAMETERS` (`thermal_l1%...`) through a NAMELIST alias to the exact same canonical L1 object; providing both is an error. Preserve both group names; do not introduce separate thermal coefficient storage.
- When `N_LAYERS=2` is eventually enabled: require `TEMPERATURE_L1_PARAMETERS` and `TEMPERATURE_L2_PARAMETERS`; reject `TEMPERATURE_PARAMETERS`, duplicate groups, missing required groups and invalid L2 parameters. Obsolete scalar syntax remains unsupported.
- Ground temperature is now stored only in `substrate%T_ground`, read from `SUBSTRATE_PARAMETERS` after TL04-C1. The old `T_erodible` input and `thermal_l1/l2%T_ground` components were removed, and affected example/test inputs and generators were migrated without old ground-temperature readers.
- `thermal_method`, `thermal_group` e l'abilitazione termica appartengono alla configurazione del **processo del layer**, non a globali che il parser sincronizza. Il registro Lie globale ordina comunque tutti i gruppi. TL04-T3 passa la configurazione termica read-only fino ai dispatcher fisici; la selezione `thermal_l1` avviene al confine single-layer in `time_integration_2d`, non dentro `source_terms_2d` o `thermal_laws_2d`.

**`T_env` was eliminated from production.** `T_ambient` is the only shared atmospheric temperature, read from mandatory `GAS_TRANSPORT_PARAMETERS` and validated positive/finite. **Current state after TL04-C1:** `substrate%T_ground` is the single terrain temperature shared across layers; the removed thermal copies are no longer compared in the input parser. The affected examples/tests use the canonical substrate input. Retain a scientific audit of cases where historically independent erodible and thermal ground temperatures could have differed; the single-temperature ownership is implemented, but CI success alone is not proof of full energy-equivalence across all scenarios.

La legge atmosferica esistente resta:
```text
Q_atm,k = -k_rad,k*(T_k**4 - T_ambient**4) - k_conv,k*(T_k - T_ambient)
k_rad,k = SBconst*emissivity_k*exp_area_fract_k
k_conv,k = atm_heat_transf_coeff_k*exp_area_fract_k
```
Nessuna nuova legge di conduzione verso il terreno, di scambio termico L1-L2 o modifica del trattamento REAL/COMPLEX/ANALYTIC. L'eventuale variazione dei risultati single-layer **derivante dalla sostituzione autorizzata di `T_env` con `T_ambient`** è ammessa e va misurata/documentata; tutte le altre differenze richiedono diagnosi. L'attuale default termico legato alla lava va migrato esplicitamente e testato, non disattivato implicitamente.

**Current thermal ownership:** cooling activation, method/group, emissivity, atmospheric transfer, exposed area and thermal closure coefficients remain in `thermal_config_type`; `T_ground` has moved to the shared substrate in TL04-C1. `T_ambient` remains SHARED and `SBconst` universal. `thermal_conductivity`, `c_p`, `enne` and `emme` do not currently affect atmospheric radiation/convection: do not introduce new laws. Radiative/convective coefficients remain computed, not persistent duplicates.

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

Testare le due NAMELIST single-layer ammesse (breve e `_L1`) sullo stesso oggetto canonico, il rifiuto di gruppi duplicati e di L2 prematuro, i vecchi scalari non supportati, le configurazioni termiche indipendenti e i kernel REAL/COMPLEX/ANALYTIC. La suite `TEST_SOURCE_REGISTRY` comprende ora anche i quattro casi di contatto L1/L2, single-layer wet/dry, soglia canonica e input non validi. Nessun contatto è ancora applicato nel percorso Lie.

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

**TL05-A — workspace HP-PCCU per slot (CI verificata su `672e0a0`):** `simulation_context_type` possiede `spatial_operators(:)` dimensionato sul layout canonico, con inizializzazione/finalizzazione per elemento. `time_integration_2d` riceve l'array ma usa soltanto `spatial_operators(1)`; anche il CFL produttivo resta L1. `TEST_SPATIAL_OPERATOR` confronta bit per bit due operatori indipendenti sullo stesso stato.

**TL05-B — workset attivi per slot (implementazione in verifica):** il contesto possiede `domains(:)` dimensionato sul layout canonico, anziché un solo `domain`. Ciascun dominio mantiene separatamente maschere wet/dry, liste di celle e facce, aloni e workset HP, aggiornati da `domain_type%check_solve`. Il driver, il Lie, lo stochastic e I/O/restarter selezionano esplicitamente `domains(1)` per preservare il single-layer. In `TEST_SPATIAL_OPERATOR`, una coppia di stati conservativi rank-4 e due domini/operatori indipendenti verifica l'identità per input uguali e l'assenza di interferenze quando cambia il wet/dry dello slot 2. `TEST_WORKSPACE_LIFECYCLE` verifica allocazione/finalizzazione ripetuta delle coppie. Il test usa ancora il kernel monolayer sulla stessa closure solo come prova di isolamento numerico, **non** come validazione della fisica L2. La futura selezione per-layer delle sorgenti e della fase, i path `eta_1`/`eta_tot` e il trasporto L2 fisicamente corretto restano aperti. TL05 e i gate scientifici non sono conclusi.

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

Single-layer: preserve physical and numerical behavior and the established dual NAMELIST naming convention for per-layer physics (`<FAMILY>_PARAMETERS` or `<FAMILY>_L1_PARAMETERS`, exclusive and bound to the same L1 configuration). Update input components and example/test fixtures when their owners migrate, without obsolete flat scalar readers or conversion wrappers. NetCDF output, `.asc/.q_2d`, probes and diagnostics retain current behavior unless an explicitly scoped I/O change is needed; the internal rank-4 layout does not itself mandate an output-format change. The existing restart v7 is the current implementation, not a permanent backward-compatibility commitment; any future schema migration defines one authoritative supported format and updates writers/readers/tests together, without supporting superseded formats unless explicitly authorized.

Two-layer: schema versionato con `n_layers`, layout, modelli e configurazioni, due stati e slot inattivi, substrato unico, topografia, variabili di continuazione, stochastic, diagnosi; NetCDF con dimensione/suffissi L1/L2, `eta_1`, `eta_tot`, Ri, densità, budget plume e `deposit_tot` una sola volta. Testare restart in assenza/presenza di scambi, deposizione, wet/dry e confronto della continuazione.

## 9. Piano operativo: task TL00–TL13

Ogni task produce codice definitivo, test e documentazione coerenti. **Implementato** non equivale a **gate scientifico chiuso**.

| Task | Stato al commit audit | Deliverable e criteri |
|---|---|---|
| **TL00 — baseline** | **Parziale** | Identificare riferimento immutabile, tag verificato, protezione branch, matrici test/tempi e output riproducibili; gate G0. |
| **TL00a — deposito canonico** | **Implementato, gate G1 aperto** | Completare test granulare/lava, `n_solid=0`, `topo_change_flag`, budget, v7/NetCDF e rollback senza doppio inventario. |
| **TL01 — rank-4 stato** | **Implementato per L1, gate G2 aperto** | Assicurare layout/contiguità e inizializzazione; test delle dimensioni, dello slot L2 inerte e dei campi asciutti. |
| **TL02 — migrazione consumer** | **Implementato per L1, audit integrato aperto** | Verificare tutte le dipendenze di state, spatial, Lie, stochastic, input/output, test, OpenMP; eliminare wrapper obsoleti. |
| **TL03 — gate single-layer** | **CI attuale superata, EXAMPLES riportati OK; gate G2 aperto** | Documentare matrice EXAMPLES e confronto baseline; completare bilanci, wet/dry, I/O/restart, OpenMP, profiling/memoria e delta termico autorizzato. |
| **TL04-R — rheology configuration** | **Implemented; CI verified on `662e796`** | Explicit canonical L1 configuration propagated through source assembly, Lie, REAL/COMPLEX and analytic callers; both single-layer NAMELIST names intentionally retained. L2 production selection and complete scientific gate remain open. |
| **TL04-T — thermal configuration** | **Implemented; CI verified on `662e796`** | Canonical L1/L2 cooling configs and explicit propagation; TL04-C1 relocated shared ground temperature. Both single-layer NAMELIST names intentionally retained; extended thermal scientific gate open. |
| **TL04-C — remaining configurations and contacts** | **C1–C4 implemented; both CI workflows green on `32e4e33`; C5/G3 open** | Preserve one combined mass-exchange assembler and the approved single-layer NAMELIST aliases. Complete deferred EXAMPLES, extended physics, restart and OpenMP verification before certifying C5. |
| **TL05 — HP-PCCU indipendente per slot** | **A: operatori per slot CI verificati su `672e0a0`; B: domain worklists per slot da validare** | Ogni slot ha workspace HP-PCCU e workset di celle/facce separati; il calcolo produttivo resta L1. Testare l'isolamento wet/dry numerico e poi implementare l'effettivo trasporto L2, `eta_1`/`eta_tot` e la riduzione single-layer. |
| **TL06 — HP-PCCU accoppiato** | **Da implementare** | Percorsi coppia, clamp, topografia, identità, CFL e iperbolicità; gate G4. |
| **TL07 — Lie con scope** | **Da implementare** | Registro unico con L1/L2/INTERLAYER, due candidati cell-local, rate monolayer isolati, `Gamma` di coppia, commit/rollback atomico, gruppo 0 e OpenMP; gate G5. |
| **TL08 — settling L2→L1** | **Da implementare** | Conservazione di masse/specie/momento/termica e nessun deposito terreno. |
| **TL09 — elutriazione L1→L2** | **Da implementare** | Inventario, co-lift, cold-start e compatibilità con i contatti. |
| **TL10 — drag interlayer** | **Da implementare** | Momento uguale/opposto, forza nulla ai limiti, analisi stabilità. |
| **TL11 — sorgenti esterne e basali** | **Da implementare** | Uso dei contatti da `physics` per atmosfera/terreno, consolidazione, entrainment, lofting, pore pressure; chiusure del manoscritto e gate G6. |
| **TL12 — validazione scientifica** | **Da implementare** | Benchmark coupled, resting states, budgets, wet/dry, OpenMP, CFL, restart, prestazioni, regressioni single-layer; gate G7. |
| **TL13 — integrazione e rilascio** | **Da implementare** | Contratti di input/NetCDF/restart, esempi e postprocessing, documentazione delle opzioni e limiti, decisione esplicita prima di abilitare `N_LAYERS=2`. |

**Dependencies:** Canonical rheology, thermal, Darcy, inert contacts, shared substrate, solid-species/packing ownership and per-layer combined mass-exchange configuration are implemented; CI on `32e4e33` passed. Finish TL04-C5 scientific certification when the deferred EXAMPLES and extended physical checks are available; in parallel prepare slot-local TL05 workspaces without enabling coupled physics. Complete isolated TL05 transport before TL06; TL07 provides paired candidates and rollback before TL11 contact dispatch and TL08–TL10 transfers; TL12 validates TL06–TL11; TL13 governs final `N_LAYERS=2` activation. Earlier green CI does not validate later commits.

## 10. Milestone M0–M6 e gate di integrazione

| Milestone | Copertura e risultato | Stato al commit audit | Condizione di uscita |
|---|---|---|---|
| **M0** | TL00: baseline e specifiche | **Da certificare** | G0: tag, audit e benchmark riproducibile, branch protetto verificato. |
| **M1** | TL00a: deposito unico | **Implementata, non certificata** | G1: deposizione lava/granulare, restart, topografia e output. |
| **M2** | TL01–TL03: rank-4 single-layer validato | **Parziale, non certificata** | G2: full EXAMPLES, confronto baseline, NetCDF/restart, OpenMP e numerica. |
| **M3** | TL04–TL05: isolated L1/L2 infrastructure | **Shared substrate, solids, L1/L2 exchange configs e workspace HP-PCCU L1 CI-verified on `672e0a0`; workset di dominio per slot in verifica; L2 fisico isolato e C5 aperti** | G3: single authoritative SHARED properties, two canonical per-layer exchange configs, one mass-exchange assembler, tests, OpenMP, inactive L2 and isolated L2 transport. |
| **M4** | TL06: operatore spaziale accoppiato | **Non iniziata** | G4: path, pressione completa, topografia, CFL, wet/dry, clamp e iperbolicità. |
| **M5** | TL07–TL11: Lie coppia e chiusure | **Non iniziata** | G5/G6: rate per-layer, coppia conservativa, budget esterni, rollback, contatti fisici e tutti i processi verificati. |
| **M6** | TL12–TL13: modello disponibile | **Non iniziata** | G7: benchmark scientifici, prestazioni, I/O e restart, approvazione esplicita all'attivazione. |

Le milestone sono **punti di proposta di merge**, non autorizzazioni automatiche. `master` deve rimanere utilizzabile single-layer. Per ogni PR: build optimized/debug e strict, regressioni pertinenti, campioni integrati, correttezza fisica, performance e OpenMP, evidenze dei gate, revisione e approvazione esplicita dell'utente. Nessun merge su `master` senza autorizzazione. Dopo un merge, sincronizzare consapevolmente il branch; nessuna storia riscritta o tag presunti.

## 11. Matrice di verifica minima e prove mancanti

| Gate | Prove obbligatorie | Evidenza attuale |
|---|---|---|
| **G0** | Baseline taggato, esempi riproducibili, confronti, condizioni OpenMP e prestazioni | Non certificato. |
| **G1** | `deposit_tot` granulare/lava, `n_solid=0`, rollback, erosione, topografia e restart v7 | Codice e alcuni test esistenti; matrice integrata non certificata. |
| **G2** | Rank-4 L1, IC/BC, Lie, spatial, NUMERICS, full EXAMPLES, NetCDF, restart, OpenMP | Build e regressioni su `672e0a0` **success**; gli EXAMPLES dopo il mass exchange sono rinviati. Confronto baseline, NetCDF/restart, OpenMP e matrice scientifica non interamente certificati. |
| **G3** | Rheology/thermal/Darcy configs, one substrate, shared solid packing, per-layer combined mass exchange, four contact cases and no L2 activation | C1–C4 code verified by both CI workflows on `32e4e33`; inert contact classification tested. Migrated EXAMPLES, detailed budgets, I/O, OpenMP and C5 scientific certification remain pending. Physical contact dispatch remains deferred. |
| **G4** | Idrostatica con path e topografia, coppia, residuo curvo, clamp, wet/dry, CFL | Da implementare. |
| **G5** | Sorgenti L1/L2 indipendenti, coppia `Gamma`, limiti e rollback atomico, gruppo 0 | Da implementare. |
| **G6** | Settling, elutriazione cold-start, drag, consolidazione, entrainment, lofting, termica e budget | Da implementare. |
| **G7** | Benchmark manoscritto, stabilità, restart coupled, I/O, OpenMP e single-layer confrontato | Da implementare. |

Contact tests must never change single-layer gating. Verify candidate-state `wet/wet`, `wet/dry`, `dry/wet`, `dry/dry` and both exposed boundaries when only one layer is wet, plus shared-substrate transactions and rollback. Thermal and exchange kernels require REAL/COMPLEX checks where supported; do not silently enable unsupported methods.

## 12. Rischi e decisioni scientifiche da chiudere nei task

1. Catalogare tutte le proprietà termiche (`c_p`, `enne`, `emme`, `thermal_conductivity`, ecc.) prima di decidere i campi esatti del tipo canonico; proprietà condivise non vanno duplicate. Nessuna nuova conduzione è autorizzata.
2. Verificare il criterio canonico wet/dry e il significato della copertura in una colonna depth-averaged; non inferire frazioni di copertura sub-cella non presenti.
3. Definire le combinazioni ammissibili di carrier/specie e gli slot inattivi di L2, con input completi per L1/L2/coppia.
4. Derivare il percorso esteso e la CFL coupled, includendo perdita di iperbolicità, clamp e topografia curva.
5. Separare trasferimenti strettamente conservativi da gas packing e altri scambi esterni; definire bilanci di massa, momento ed energia termica.
6. Determinare le chiusure di settling, packing, entrainment, lofting e consolidazione del manoscritto prima di abilitarle.
7. Verificare il contratto del gruppo 0, Newton/analitico e rollback anche su coppie di stati.
8. Definire formato autorevole di restart/NetCDF e testare continuazione e risultati single-layer dopo la migrazione; nessun reader di formati superati senza richiesta esplicita.
9. Dimostrare performance e riproducibilità OpenMP con configurazioni read-only e nessuna fisica hard-coded nel Lie.

**TL04-C scientific/input risks:** Shared `T_ground` is now implemented in the authoritative substrate and repository inputs were migrated; explicitly audit historical cases where ground and erodible temperatures differed, since one temperature can change physical assumptions even when current CI and EXAMPLES pass. Retain `erodible_deposit_flag` energy semantics and historical lava defaults. Do not equate `maximum_solid_packing`, deposition hindrance `alpha_max` and future L1 packing coefficients. Keep current `alphastot_min` shared between erosion/deposition until an authorized change. Do not claim implicit or analytical combined mass exchange until coupled rate/residual, limits, supported treatment and conservation have all been demonstrated.

## 13. Regola operativa

Procedere su `feature/two-layer` con audit completo della catena dipendenze **prima** di ogni modifica, gruppi coerenti di commit, test REAL/COMPLEX quando coinvolti e registrazione di esiti sul commit effettivo. **Regola permanente: una sola configurazione autorevole per layer e nessuna API interna obsoleta; la coppia di nomi NAMELIST `<FAMILY>_PARAMETERS`/`<FAMILY>_L1_PARAMETERS` è un contratto supportato solo per `N_LAYERS=1`, con lettura diretta dello stesso oggetto, non una copia sincronizzata.** Eliminare lettori scalari superati, wrapper e ogni altra sintassi transitoria; mantenere coerenza di EXAMPLES, TESTS e parser. La conservazione dei risultati single-layer non richiede la conservazione di vecchi componenti o strutture di input. La roadmap resta un **documento consolidato dello stato e del piano vigente**: aggiornare sezioni esistenti, non accumulare note cronologiche né mantenere contemporaneamente piano vecchio e nuovo. Lo storico è nei commit; le evidenze statiche nel [rapporto di audit](TWO_LAYER_CODE_AUDIT.md).

**Riferimenti tecnici:** *IMEX_SfloW2D_v3* (sez. 2.1–2.2, 2.5–2.7, 2.12, 2.14, 3–4); `src/state/`, `src/physics/`, `src/spatial/nonconservative_2d.f90`, `src/spatial/pccu_2d.f90`, `src/numerics/lie_source_executor_2d.f90`, `docs/FORTRAN_STYLE.md`.
