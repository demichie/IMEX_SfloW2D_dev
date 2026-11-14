# IMEX-SfloW2D — Architettura e roadmap two-layer

**Plan updated:** 14 November 2026. **Branch:** `feature/two-layer`. **Verified CI baseline:** [`9427826`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/942782656d43f4e2e3fe2bd6659a95e4ab83ba82) (build optimized/debug and cell-physics/numerical regressions successful). **TL04-C3/C4:** canonical mass-exchange configuration remains verified; scientific gate C5/G3 remains open. **TL05-A/B/C1:** independent domains/operators and L2 buoyancy-pressure clamp are verified; TL05-C2 adds isolated one-step conservation checks, with new CI pending. No coupled two-layer execution.
**Stato:** codice di produzione ancora **single-layer**; `N_LAYERS=2` resta disabilitato. La presente roadmap è il **solo piano vigente**: specifica, stato verificato, attività, gate e milestone. La cronologia dettagliata è nei commit Git, non in versioni alternative del piano.
**Documenti collegati:** [audit tecnico corrente](TWO_LAYER_CODE_AUDIT.md), [roadmap generale](IMEX_SfloW2D_Development_Roadmap.md), manoscritto *IMEX_SfloW2D_v3*.

## 1. Obiettivo, invarianti e ambito

Il modello ha due sole modalità: `N_LAYERS=1` (solutore attuale) e, **solo dopo il gate finale**, `N_LAYERS=2` con L1 corrente inferiore densa e L2 nube superiore diluita nello stesso scenario gas/solidi. Non è previsto un framework con numero arbitrario di layer, né il ripristino dell'IMEX-RK storico. Le combinazioni lava/granulare tra i due layer non sono abilitate dal semplice fatto che il single-layer supporta la lava.

1. **Semplicità e una sola rappresentazione autorevole** di configurazioni, proprietà fisiche e inventari. Riutilizzare direttamente ciò che esiste: non introdurre sovrastrutture, wrapper, livelli di dispatch o stati intermedi se i moduli attuali sono sufficienti. Nessun globale duplicato, copia sincronizzata o duplicazione di formule.
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
| Physical sources | The registry remains single-layer; the existing pure `layer_contacts_2d` already describes all four possible L1/L2 contact states. Canonical L1/L2 mass-exchange configurations and the combined rate/limiter backend are CI verified on `32e4e33`. | At TL07/TL11, determine contacts outside `physics` at each Lie group entry from current candidates, after existing buoyancy removal; pass only the resulting local exposures to physics. Do not add a contact-mask API or duplicate physical kernels. Complete source scopes, candidate-state integration and phase-selection audit; C5/EXAMPLES verification remains open. |
| Spaziale | Un solo motore HP-PCCU; `spatial_operators(:)` e `domains(:)` indipendenti (CI su `d901609`). Clamp idrostatico di L2 TL05-C1 verificato dalle due Actions su `9427826`. TL05-C2 estende direttamente il test esistente a un passo Euler CFL-limited L2 gas/solidi, con controlli di ammissibilità e budget, CI da verificare. | Il driver evolve ancora solo L1. TL05-C resta aperto; seguono interfacce `eta_1`/`eta_tot`, pressione accoppiata, CFL two-layer e wet/dry in TL06. |
| I/O | NetCDF single-layer `(x,y,time)`, restart v7 with L1 plus shared inventories, current `.asc/.q_2d` outputs and diagnostics. | Define one authoritative schema when I/O changes, update reader/writer/tests together and verify scientifically meaningful single-layer continuation; no superseded-format reader or backward-compatibility requirement. |

**CI evidence:** [Build (optimized and debug) #43584710511](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43584710511) and [Cell physics and numerical regressions #43584822116](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43584822116) both succeeded on [`9427826`](https://github.com/demichie/IMEX_SfloW2D_dev/commit/942782656d43f4e2e3fe2bd6659a95e4ab83ba82), covering TL05-C1 and earlier preparation. The user deferred extended EXAMPLES after the mass-exchange migration; older examples do not certify C1/C2. NetCDF/restart, complete OpenMP, coupled conservation and two-layer scientific gates remain open. TL05-C2 requires new CI on its own commit.

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
7. **TL07/TL11 — deferred layer activation.** Add L1/L2/INTERLAYER scheduling and paired candidate-state/atomic rollback infrastructure first; then call the existing contact classifier from the numerical cell/group loop, outside `physics`, and pass the resulting local exposures to layer-specific external exchange. Exercise wet/wet, wet/dry, dry/wet and dry/dry, including L2-only bed erosion/deposition and L1-only atmospheric entrainment; no double counting of common terrain/atmosphere resources. Interlayer transfers remain separate pair operators under TL08–TL10, with coherent global budgets.

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

### 4.3.1 Principio di implementazione minima dei contatti (14 novembre 2026)

La sola classificazione necessaria e' quella gia' implementata in
`layer_contacts_2d`. Per gli spessori candidati correnti `(h1,h2)`
esistono **esattamente quattro stati**: entrambi presenti
(L1 a contatto con il terreno, L2 con l'atmosfera), nessuno presente
(nessun contatto), solo L1 presente (L1 con entrambi), solo L2
presente (L2 con entrambi). La tabella in sezione 4.2 e' il riferimento
autorevole; non introdurre altre classificazioni o rappresentazioni.

Quando TL07 rendera' disponibili i candidati L1/L2 della stessa cella,
il Lie chiamera' **direttamente** `evaluate_layer_contacts` fuori da
`physics`, dopo il controllo di buoyancy e all'inizio di ogni gruppo,
per decidere quali processi siano ammissibili per il layer bersaglio. Erosione e deposizione sul terreno
richiedono `terrain_exposed(i)`; entrainment atmosferico richiede
`atmosphere_exposed(i)`. Gli altri processi conservano le rispettive
condizioni fisiche; la vaporizzazione legata al deposito dipende
dalla deposizione effettiva. Contatto necessario non significa
attivazione automatica: valgono ancora configurazione, fase, inventari,
limitatori e condizioni di wet/dry.

`lie_source_executor_2d` conserva indici di griglia, scelta del layer,
scheduling, stati candidati, classificazione dei contatti e rollback.
In `physics` si lavora solo con dati **cell-local**, senza indici `j,k`,
identita' del layer o calcolo di `exposed`: la fisica riceve i contatti
gia' determinati all'esterno. Poiche' `mass_exchange_2d`
accetta gia' `source_mask`, il suo filtro locale sara' applicato
**direttamente nell'assemblaggio esistente**, senza nuova routine
pubblica, wrapper, ulteriore tipo o secondo coordinatore.
`mass_exchange_2d` e le sei leggi fisiche restano invariati.

Per L2 restano da rimuovere le dipendenze dalla fase globale
(`gas_flag`/`liquid_flag`) dove la fase e' specifica del layer.
Il registro resta single-layer e i controlli per L2 non si attivano
prima dei candidati accoppiati e rollback atomico di TL07/TL11.
Test richiesti in quella fase: quattro stati, cambi di esposizione tra
gruppi Lie, budget del substrato/atmosfera, dry/cold-start e regressioni L1.

**Correzione della preparazione precedente:** la routine
`filter_contact_mass_exchange_mask` aggiunta in `e4ef157` e' stata
rimossa, insieme ai relativi test e alla dipendenza Makefile non piu'
necessaria. `layer_contacts_2d` e i suoi test originari sono sufficienti
per questa fase; nessuna modifica alle leggi o al percorso Lie produttivo.

### 4.4 Verifica specifica della migrazione

Testare le due NAMELIST single-layer ammesse (breve e `_L1`) sullo stesso oggetto canonico, il rifiuto di gruppi duplicati e di L2 prematuro, i vecchi scalari non supportati, le configurazioni termiche indipendenti e i kernel REAL/COMPLEX/ANALYTIC. La suite `TEST_SOURCE_REGISTRY` comprende ora anche i quattro casi di contatto L1/L2, single-layer wet/dry, soglia canonica e input non validi. Nessun contatto è ancora applicato nel percorso Lie.

## 5. Partizione matematica e Lie

```text
dq1/dt = H1(q1,q2,B) + S_local1(q1;C1,E1) + Gamma(q1,q2) + S_external1
dq2/dt = H2(q1,q2,B) + S_local2(q2;C2,E2) - Gamma(q1,q2) + S_external2
```
`H1/H2` includono flussi, pressione, topografia e termini non conservativi **spaziali**. `S_local` è una legge monolayer; `E_k` è contesto fisico esterno già preparato; `Gamma` è il trasferimento strettamente interno fra layer. Gli scambi esterni hanno bilancio separato. Ogni legge locale riceve il solo stato del proprio layer e restituisce il solo rate di quel layer.

**Un unico integratore Lie globale**, con un solo timestep, ordinamento e rollback. Il registro usa `(source_identity,scope)` per l'unicità; ogni gruppo dichiara `L1`, `L2` o `INTERLAYER`, metodo, ordine e footprint, compresa la gestione delle risorse condivise. La stessa legge può comparire in L1 e L2. L'attuale `source_group_type` **non** dispone ancora di `scope`, e `time_integration_2d` passa al Lie soltanto `q_trial(:,1,:,:)`, `qp_trial(:,1,:,:)`: non dichiarare completata TL07.

Il Lie gestisce candidati e snapshot pre-spatial del gruppo 0, transazioni dello stato e substrato, geometria, stochastic, buone proprietà numeriche e rollback. Il livello numerico determina i contatti; `physics` riceve l'esposizione locale gia' selezionata e valuta leggi, rate e trasferimenti. Con gruppi `INTERLAYER` si acquisisce uno **stesso snapshot pre-gruppo** di L1 e L2, si calcola un solo `Gamma(n_eqns)` con convenzione positiva L2→L1 e si costruiscono contemporaneamente:
```text
q1_candidate = q1 + dt*Gamma
q2_candidate = q2 - dt*Gamma
```
Un solo limiter di disponibilità condiviso; nessun limitatore indipendente che rompa l'antisimmetria. Validazione e rollback atomici di entrambi i layer, substrato, geometria e diagnostiche. La variabile termica conservativa trasferita va contabilizzata; l'energia cinetica e la dissipazione possono avere regole ulteriori, esplicite. L'approssimazione del gas di packing è un **difetto esterno diagnosticato**, non un trasferimento fittiziamente conservativo. Non imporre un Newton globale `2*n_eqns` senza evidenza di necessità.

**Pre-TL07 architecture preparation (14 November 2026, CI pending).**
The existing contact classifier is moved unchanged from `physics/` to
`numerics/`, with no second implementation or change to its four-case contract.
The one-time source-registry configuration builder is also owned by `numerics/`
and receives the selected canonical rheology, thermal, pore-pressure and
mass-exchange configurations explicitly, rather than reading `_l1` inside
its module. The production caller still selects L1; L2 sources and the
per-group contact refresh remain deferred. The existing physical kernels,
phase flags and single-layer numerical path are not changed by this step.

**TL07 — contact refresh at every Lie group entry (approved contract).**
Keep the existing cell-outer / group-inner traversal. At the start of **each**
Lie group, use the current post-advection or previous-group accepted L1/L2
candidate states. First apply the already implemented L2 buoyancy-removal
criterion (`gprime2 <= 0`): clear L2's conservative and primitive local
candidate, not only its thickness. This is the existing physical sink, not
an additional limiter or a newly invented conservation transfer. Unlike the
single-layer early exit, removal of L2 must not terminate the cell's Lie
sequence while L1 may still require source updates. Next invoke the **existing**
`evaluate_layer_contacts` on the two current candidate thicknesses,
**outside** `physics`, then pass the resulting terrain/atmosphere
exposure and interface-contact information to the relevant cell-local physics.
No previous-contact cache, change-detection flag, duplicate classifier,
extra wrapper, or additional runtime contact representation is needed.

Contacts are fixed during evaluation of one group and refreshed before the
next; accepted deposition, entrainment, cold-start, loss of a layer and the
post-advection state therefore affect the subsequent group automatically.
Apply final buoyancy removal after the last group as in the current Lie.
Preserve group-0 pre-spatial input semantics, strict candidate validation,
shared-substrate budgets and atomic rollback. A source treatment requiring
contact changes *inside* one implicit group needs a separate documented
numerical treatment; do not silently reinterpret group semantics. This is
a TL07/TL11 design decision, **not implemented**; `N_LAYERS=2` remains disabled.

Preservare esattamente la semantica del gruppo 0 e delle celle inizialmente asciutte, i quattro trattamenti numerici, la struttura cell-outer e gli scratch per worker; nessuna allocazione per cella o modifica di globali nell'OpenMP.

## 6. Operatore spaziale accoppiato HP-PCCU

Riutilizzare `src/spatial/nonconservative_2d.f90` **senza cambiarne l'API**: callback generica `I(W,dW/ds)` e quadratura a 1/2/3 punti. Non introdurre un secondo motore né differenze finite automatiche di derivate di percorso. Percorso e variabili estese scelti esplicitamente e coerenti su facce, interno e wet/dry.

La geometria comune è:
```text
B       = quota del substrato
eta_1   = B + h1
eta_tot = B + h1 + h2
```
### Approved two-layer geometric contract — IMEX_v3.pdf

**Decision (13 November 2026, before TL05-C):** adopt the terrain-based geometry of
*IMEX_SfloW2D_v3: A two-layer shallow-water model for pyroclastic density currents with
pore-pressure transport, polydisperse settling and shear-driven elutriation* as the
scientific reference for large-slope correction, curvature and interlayer pressure.
The model uses a global Cartesian frame and vertically measured layer thicknesses.
This records the chosen equations; it does **not** mean they are already implemented
or scientifically validated in the current repository.

- **SHARED reference geometry:** the ground elevation `B(x,y)`, `grad(B)` and
  `B_xx, B_xy, B_yy` have one authoritative owner, with no L1/L2 copies.
  The large-slope metric is `G = g_cell = 1/(1 + B_x**2 + B_y**2)` for **both**
  layers. Do **not** replace `B` by `B+h1` when forming metric or terrain curvature.
- **Layer-local kinematics:** when the large-slope correction is enabled, each
  layer's bed-tangent vertical velocity is `w_i = u_i*B_x + v_i*B_y`;
  use its own `(u_i,v_i)` in velocity magnitudes of the applicable closures.
  Its curvature acceleration is
  `kappa_curv_i = B_xx*u_i**2 + 2*B_xy*u_i*v_i + B_yy*v_i**2`.
  Thus terrain derivatives are SHARED but `w_i` and `kappa_curv_i` are
  computed per layer from the **current** layer velocities; no stored copies.
- **Curvature placement:** distinguish the nonnegative reduced gravity
  `max(gprime_i,0)`, which enters hydrostatic pressure fluxes and the
  single-layer wave-speed estimates, from
  `gtilde_i = max(gprime_i,0) + kappa_curv_i`, used by the topographic
  slope source and (for dense L1) the effective normal stress of rheology.
  Do **not** add curvature to pressure fluxes or CFL celerities.
  Preserve the existing explicit metric-gradient source contributions.
- **Moving interface, not a replacement terrain:** `eta_1=B+h1` and
  `eta_tot=B+h1+h2` are derived from the two-layer state, not stored
  alternative geometries. The upper flow feels the gradient of `B+h1`
  through the hydrostatic pressure/coupling terms (Eqs. 6–7, 14–17),
  **not** through a separate slope metric or a Hessian of `B+h1`.
  The lower flow feels the spatial gradient of the overlying load.
- **Nonconservative pair:** use the two products and metric placements given
  directly below (Sect. 2.7.1 of the manuscript): `G` lies **outside**
  the L1 derivative and **inside** the differentiated L2 expression.
  When `L = rho_m2*gprime2*h2` is nonnegative and no clamp is active,
  their sum is `-d_j(G*h1*L)`. The ash-cloud-side pressure load is clamped
  to zero for buoyant L2; the lower-side product retains the raw load
  gradient as specified in the manuscript. Track any resulting pair-budget
  defect explicitly; do not silently force antisymmetry in the clamped case.
- **Well-balancing qualification:** exact stratified rest on a flat bed,
  and the paired product identity without clamping, are required checks.
  On curved terrain the manuscript accepts a residual proportional to
  `-rho_m2*gprime2*h2*h1*d_j(G)`, at the single-layer curved-bed
  balance floor. Do not claim exact curved-bed equilibrium or change
  the metric placements solely to remove this residual.
- **Activation and scope:** these decisions constrain TL05-C (isolated L2
  transport) and TL06 (coupled HP-PCCU). They do not enable `N_LAYERS=2`,
  introduce interlayer curvature forces, or authorize changes to existing
  single-layer flags, regularizations or numerical evaluation order.
  Any layer-specific reinterpretation of `SLOPE_CORRECTION_FLAG` or
  `CURVATURE_TERM_FLAG` requires a separate explicit decision.

**Scientific reference:** IMEX_v3.pdf, Sect. 2.1.2 (Eqs. 6–7), Sect. 2.3
(Eq. 34), Sect. 2.7.1 (metric placements and buoyancy clamp), Sect. 2.11
(Eq. 73), Sect. 3.0.2 (path-conservative products), and Sect. 4.0.3
(stratified lake-at-rest verification).

**Required TL05/TL06 tests:** flat bed (`G=1`, zero curvature), constant-slope
terrain, curved terrain with distinct `(u1,v1)` and `(u2,v2)`, a nonuniform
`h1` with `B` fixed, L1-only/L2-only dry reductions, stratified rest,
the unclamped pair-sum and clamped buoyancy-budget cases, coupled
reconstruction/wet-dry and CFL/hyperbolicity diagnostics. Preserve
single-layer bitwise regressions whenever the new features are inactive.

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

**TL05-B — workset attivi per slot (CI verificata su `d901609`):** il contesto possiede `domains(:)` dimensionato sul layout canonico, anziché un solo `domain`. Ciascun dominio mantiene separatamente maschere wet/dry, liste di celle e facce, aloni e workset HP, aggiornati da `domain_type%check_solve`. Il driver, il Lie, lo stochastic e I/O/restarter selezionano esplicitamente `domains(1)` per preservare il single-layer. In `TEST_SPATIAL_OPERATOR`, una coppia di stati conservativi rank-4 e due domini/operatori indipendenti verifica l'identità per input uguali e l'assenza di interferenze quando cambia il wet/dry dello slot 2. `TEST_WORKSPACE_LIFECYCLE` verifica allocazione/finalizzazione ripetuta delle coppie. Il test usa ancora il kernel monolayer sulla stessa closure solo come prova di isolamento numerico, **non** come validazione della fisica L2. La futura selezione per-layer delle sorgenti e della fase, i path `eta_1`/`eta_tot` e il trasporto L2 fisicamente corretto restano aperti. TL05 e i gate scientifici non sono conclusi.

**TL05-C1 — L2 isolated gas hydrostatic pressure (CI verified on `9427826`, 14 November 2026).**
The existing `spatial_operator_type%initialize(layer)` accepts the material
slot (1 by default; 2 for the upper dilute layer). It passes a single
nonnegative-pressure choice to its owned `hyperbolic_workspace_type`.
For L2, the face CFL calculation stores
`Gamma2 = rho_m2*max(gprime2,0)` instead of a negative `Gamma2`;
the same value is used in face hydrostatic paths, in-cell paths and
stationary-pressure roundoff scaling. For L1 the existing sequence,
formula and signed `Gamma1` are unchanged. Negative `gprime2` already
produced zero celerity in the existing speed bounds, so no new
characteristic-speed calculation or second physics kernel is introduced.
The canonical mixture and shared terrain derivatives are reused.

`simulation_2d` passes its existing loop index at spatial-workspace
initialization. `TEST_SPATIAL_OPERATOR` keeps its positive-gravity
bitwise L1/slot-2 comparison and additionally exercises a gas L2-only
state: hot buoyant gas at rest (zero pressure/transport residual)
and dense moving gas/solids (positive reduced gravity, finite spatial
term and CFL). The single-layer liquid fixture is preserved for the
reference comparison, and no producer enables `N_LAYERS=2`.
This remains a **partial TL05-C implementation**, not certification of
coupled pressure, L2 sources, paired positivity or two-layer CFL.
Both Actions are verified on `9427826`; extended EXAMPLES and full
scientific gates are still open.

**TL05-C2 — isolated L2 transport invariants (CI pending, 14 November 2026).**
Extend `TEST_SPATIAL_OPERATOR` in place, without adding production
modules, numerical kernels or a second integrator. Starting from a
fully wet dense gas/solid L2 and dry L1, use the existing per-slot
worklists and HP-PCCU CFL to take one test-only forward-Euler
conservative transport step. An inactive/constant boundary collar
eliminates external advective fluxes; verify positive energy,
admissible component masses and conservation of total mixture mass,
solid mass and thermal-energy variable. The moving dense case must
produce nonzero transport on the reference 2D grid; L1 and the
single-layer reference arrays must remain unchanged. The existing
buoyant-gas zero-force and identical-layer bitwise tests are retained.
These isolated tests **do not** activate L2 in production or prove
the coupled two-layer CFL, source transfers, shared substrate budgets,
or curved-bed stratified well-balancing.

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
| **TL05 — HP-PCCU indipendente per slot** | **A/B/C1 CI verificata su `9427826`; C2 test Euler in verifica, C aperto** | Uno stesso HP-PCCU serve entrambi gli slot. L2 usa `max(gprime2,0)` e i test gas isolati includono un aggiornamento conservativo CFL-limited con bilanci di massa, solidi ed energia. CI C2 e completamento fisico L2 ancora aperti, prima di TL06. |
| **TL06 — HP-PCCU accoppiato** | **Da implementare** | Percorsi coppia, clamp, topografia, identità, CFL e iperbolicità; gate G4. |
| **TL07 — Lie con scope** | **Da implementare; contratto contatti tra gruppi definito** | Registro unico L1/L2/INTERLAYER, candidati cell-local e rollback atomico; loop esterno celle, interno gruppi. All'ingresso di ogni gruppo: rimozione buoyancy L2 gia' esistente, poi classificazione dei quattro contatti sui candidati correnti **fuori da `physics`**; esposizioni passate alle leggi locali. Non interrompere i gruppi L1 se L2 scompare; preservare gruppo 0, `Gamma` di coppia e OpenMP; gate G5. |
| **TL08 — settling L2→L1** | **Da implementare** | Conservazione di masse/specie/momento/termica e nessun deposito terreno. |
| **TL09 — elutriazione L1→L2** | **Da implementare** | Inventario, co-lift, cold-start e compatibilità con i contatti. |
| **TL10 — drag interlayer** | **Da implementare** | Momento uguale/opposto, forza nulla ai limiti, analisi stabilità. |
| **TL11 — sorgenti esterne e basali** | **Da implementare; classificatore dei contatti esistente** | Usare i contatti aggiornati dal Lie all'ingresso di ogni gruppo e passarli a `source_terms_2d`, filtrando i bit del `source_mask` nell'assemblaggio gia' presente, senza una nuova API. Mantenere i kernel di mass exchange invariati, completare le dipendenze di fase per layer, sorgenti esterne e bilanci SHARED; gate G6, sezione 4.3.1. |
| **TL12 — validazione scientifica** | **Da implementare** | Benchmark coupled, resting states, budgets, wet/dry, OpenMP, CFL, restart, prestazioni, regressioni single-layer; gate G7. |
| **TL13 — integrazione e rilascio** | **Da implementare** | Contratti di input/NetCDF/restart, esempi e postprocessing, documentazione delle opzioni e limiti, decisione esplicita prima di abilitare `N_LAYERS=2`. |

**Dependencies:** Canonical rheology, thermal, Darcy, inert contacts, shared substrate, solid-species/packing ownership and per-layer combined mass-exchange configuration are implemented; CI on `32e4e33` passed. Finish TL04-C5 scientific certification when the deferred EXAMPLES and extended physical checks are available; in parallel prepare slot-local TL05 workspaces without enabling coupled physics. Complete isolated TL05 transport before TL06; TL07 provides paired candidates and rollback before TL11 contact dispatch and TL08–TL10 transfers; TL12 validates TL06–TL11; TL13 governs final `N_LAYERS=2` activation. Earlier green CI does not validate later commits.

## 10. Milestone M0–M6 e gate di integrazione

| Milestone | Copertura e risultato | Stato al commit audit | Condizione di uscita |
|---|---|---|---|
| **M0** | TL00: baseline e specifiche | **Da certificare** | G0: tag, audit e benchmark riproducibile, branch protetto verificato. |
| **M1** | TL00a: deposito unico | **Implementata, non certificata** | G1: deposizione lava/granulare, restart, topografia e output. |
| **M2** | TL01–TL03: rank-4 single-layer validato | **Parziale, non certificata** | G2: full EXAMPLES, confronto baseline, NetCDF/restart, OpenMP e numerica. |
| **M3** | TL04–TL05: isolated L1/L2 infrastructure | **Shared substrate, solids, L1/L2 exchange configs and spatial operators/domains per slot CI-verified on `d901609`; isolated physical L2 transport and C5 remain open** | G3: single authoritative SHARED properties, two canonical per-layer exchange configs, one mass-exchange assembler, tests, OpenMP, inactive L2 and isolated L2 transport. |
| **M4** | TL06: operatore spaziale accoppiato | **Non iniziata** | G4: path, pressione completa, topografia, CFL, wet/dry, clamp e iperbolicità. |
| **M5** | TL07–TL11: Lie coppia e chiusure | **Non iniziata** | G5/G6: rate per-layer, coppia conservativa, budget esterni, rollback, contatti fisici e tutti i processi verificati. |
| **M6** | TL12–TL13: modello disponibile | **Non iniziata** | G7: benchmark scientifici, prestazioni, I/O e restart, approvazione esplicita all'attivazione. |

Le milestone sono **punti di proposta di merge**, non autorizzazioni automatiche. `master` deve rimanere utilizzabile single-layer. Per ogni PR: build optimized/debug e strict, regressioni pertinenti, campioni integrati, correttezza fisica, performance e OpenMP, evidenze dei gate, revisione e approvazione esplicita dell'utente. Nessun merge su `master` senza autorizzazione. Dopo un merge, sincronizzare consapevolmente il branch; nessuna storia riscritta o tag presunti.

## 11. Matrice di verifica minima e prove mancanti

| Gate | Prove obbligatorie | Evidenza attuale |
|---|---|---|
| **G0** | Baseline taggato, esempi riproducibili, confronti, condizioni OpenMP e prestazioni | Non certificato. |
| **G1** | `deposit_tot` granulare/lava, `n_solid=0`, rollback, erosione, topografia e restart v7 | Codice e alcuni test esistenti; matrice integrata non certificata. |
| **G2** | Rank-4 L1, IC/BC, Lie, spatial, NUMERICS, full EXAMPLES, NetCDF, restart, OpenMP | Build e regressioni su `d901609` **success**; gli EXAMPLES dopo il mass exchange sono rinviati. Confronto baseline, NetCDF/restart, OpenMP e matrice scientifica non interamente certificati. |
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
