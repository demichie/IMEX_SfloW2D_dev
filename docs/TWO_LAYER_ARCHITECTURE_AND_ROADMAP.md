# IMEX-SfloW2D: architettura single-layer / two-layer e roadmap

**Data:** 7 novembre 2026  
**Stato:** PROPOSTA, nessun codice implementato  
**Repository:** demichie/IMEX_SfloW2D_dev, branch master  
**Documento specialistico:** subordinato a [IMEX_SfloW2D_Development_Roadmap.md](IMEX_SfloW2D_Development_Roadmap.md), che rimane il tracker operativo generale.

## 1. Obiettivi e decisioni concordate

Estendere il modello depth-averaged esistente alle sole due configurazioni `N_LAYERS=1` e `N_LAYERS=2`, senza un framework generico a N layer. L1 è il flusso inferiore denso, L2 il flusso superiore diluito; le equazioni e le chiusure two-layer devono seguire il manoscritto *IMEX_SfloW2D_v3* fornito dall'utente, adattando la sua precedente integrazione IMEX-RK all'attuale Lie di primo ordine.

Decisioni:
1. Introdurre **subito, come primo intervento futuro**, la dimensione layer come **secondo indice**: `q(n_variables,n_layers,nx,ny)` e `qp(n_physical_variables,n_layers,nx,ny)`. Durante questa fase `n_layers=1` rimane l'unica configurazione ammessa.
2. Conservare indici di variabili/equazioni identici nei due layer. Il vettore di L2 può contenere slot inattivi (es. pressione interstiziale, fisicamente propria di L1), da mascherare esplicitamente.
3. Unico integratore **Lie globale**, timestep, rollback e pubblicazione dello stato; ogni gruppo opera su L1, L2 o sulla coppia.
4. Kernel fisici dei singoli layer invariati e **cell-local**; nessun indice spaziale j,k dentro physics. Stato fisico locale e workspace Newton restano monolayer; per un gruppo interlayer si preparano i due stati locali necessari.
5. Ogni scambio interno conservativo L1/L2 produce **un solo vettore di trasferimento firmato** `Gamma(n_eqns)`: positivo per L2 -> L1, applicazione `+Gamma` a L1 e `-Gamma` a L2, con un'unica limitazione di inventario e verifica di ammissibilità per entrambi.
6. Il **settling L2 -> L1 non forma deposito**: tutto il materiale resta nei layer mobili. La **deposizione L1 -> terreno** è un processo distinto, aggiorna inventari del substrato e può modificare topografia.
7. L'**entrainment atmosferico** è già in `src/physics/mass_exchange/mass_exchange_entrainment_2d.f90`: riusarlo come scambio esterno verso L2, senza nuova directory o duplicazioni.
8. Nuova directory `src/physics/interlayer/` per settling, elutriazione, drag e assemblaggio fisico locale accoppiato; il gradiente idrostatico non conservativo interlayer va invece in `src/spatial/`.
9. Conservare il single-layer con formule, ordine Lie e risultati invariati; lavoro precedente di rifattorizzazione/ottimizzazione **stateNEW in pausa**, separato da questo piano.

Non obiettivi: N>2, nuova architettura IMEX-RK, Newton due-layer generale preventivo, migrazione MPI/GPU o integrazione ML in questa fase, file wrapper senza fisica, cambiamenti gratuiti di formule.

## 2. Situazione attuale verificata nel repository

| Componente | Stato verificato | Impatto previsto |
|---|---|---|
| `src/core/model_layout_2d.f90` | Esiste già `n_layers`, ma `supported_layer_count` ammette soltanto 1. | Conservare il vincolo fino alla certificazione two-layer. |
| `src/state/state_fields.f90` | `state_type%q` e `qp` a rango 3; inventari e diagnostiche separati. | Passare a rango 4 solo per i campi fluidi pertinenti. |
| `src/numerics/time_integration_2d.f90` | Lie di primo ordine, `spatial_rate`, `qp_before`, `q_trial` e `qp_trial` a rango 3; rollback di substrato e geometria. | Adeguare i workspace, preservare transazione globale. |
| `src/numerics/lie_source_executor_2d.f90` | Gruppi ordinati, scratch REAL/COMPLEX per worker, trasferimenti e controlli locali. | Aggiungere in una seconda fase il supporto a gruppi di scope L1/L2/INTERLAYER. |
| `src/numerics/source_registry_2d.f90` | Gruppi con metodo, ID, footprint; gruppo 0 con semantica speciale. | Estendere metadata senza alterare il gruppo 0. |
| `src/physics/source_terms_2d.f90` | Valutazione locale, assemblaggio e limitazione dei rate. | Riutilizzare, distinguendo scambi interlayer dagli esterni. |
| `src/physics/mass_exchange/` | Sono presenti entrainment, deposizione ed erosione. | Nessuna duplicazione. |
| `src/spatial/` | HP-PCCU, ricostruzione e termini non conservativi esistenti. | Aggiungere coppia di pressione idrostatica interlayer. |
| `src/Makefile.am` | Elenco esplicito dei sorgenti. | Aggiornare solo quando nasceranno nuovi moduli. |

**Nota dimensionale:** `n_variables` e `n_physical_variables` sono metadata distinti. Non assumere che q e qp abbiano necessariamente lo stesso primo indice massimo. Audit preliminare obbligatorio su tutti gli accessi: allocazioni, snapshot, halo, BC, ricostruzione, input, NetCDF, restart, diagnostiche, tests.

## 3. Equazioni e partizione degli operatori

Per le due variabili conservative locali `q1`, `q2`:

```text
d q1/dt = H1(q1,q2,B) + S_local1(q1) + Gamma(q1,q2) + S_external1
d q2/dt = H2(q1,q2,B) + S_local2(q2) - Gamma(q1,q2) + S_external2
```

`H`: flussi spaziali, gravità/topografia e pressione accoppiata path-conservative. `S_local`: reologia, termica, pore pressure e processi propri del layer. `Gamma`: solo trasferimenti L1/L2 realmente interni e conservativi. `S_external`: atmosfera, terreno, iniezioni e altri ingressi/uscite. La partizione è una scelta di architettura del codice; non sostituisce letteralmente le etichette S1/S2/S3 del manoscritto.

### 3.1 Idrostatica interlayer — spatial

La sezione 2.7 del manoscritto usa la coppia, per una direzione orizzontale j:

```text
L1: -gcell*h1 * d_j(rho_m2*gprime2*h2)
L2: -(rho_m2*gprime2*h2) * d_j(gcell*h1)
```

Con trattamento specifico della buoyancy/clamp nel secondo membro. Sono prodotti non conservativi spaziali: discretizzazione coerente con flussi, path, velocità alle facce, fattori metrici e topografia. Implementare in `spatial/` e **non** in un gruppo Lie locale. Verificare lago stratificato fermo su letto piano, residuo noto su letto curvo, bilanci e velocità caratteristiche; non rivendicare bilancio esatto dove il clamp o il residuo curvo lo impediscono.

### 3.2 Scambi interlayer — physics/interlayer

- **Settling L2 -> L1:** per specie solida `s`, un `Gamma_s>0` assegnato come guadagno a L1 e perdita a L2; modificare coerentemente massa totale, momento con un'unica velocità di interfaccia, energia termica con proprietà/temperatura del donatore e specie trasportate. Nessun aggiornamento di substrato o topografia.
- **Elutriazione L1 -> L2:** `Gamma_s<0`; velocità e temperatura del donatore L1; filtro e legge di shear secondo manoscritto; gas co-lift se attivo.
- **Drag:** solo momento, `D=f*rho_m2*|u2-u1|*(u2-u1)` verso L1 e `-D` verso L2, con gate di spessore; valutare esplicito prima di considerare solver accoppiato implicito.
- **Energia termica:** il rate della variabile termica conservativa è opposto nei due layer per i trasferimenti interni. Non confondere con conservazione dell'energia cinetica e dissipazione, che richiederebbero una chiusura separata.

**Eccezione fisica esplicita:** nel manoscritto il gas interstiziale aggiunto a L1 per la formazione del materiale a packing da settling L2 può non essere interamente sottratto a L2. Non va nascosto in un `Gamma` nominalmente antisimmmetrico: mantenere un correttivo identificabile e diagnosticarne il deficit, oppure scegliere e validare scientificamente una chiusura conservativa alternativa. La deposizione L1 -> terreno e il gas sepolto sono separati. Evitare di chiamare "deposito" lo scambio L2 -> L1.

### 3.3 Sorgenti esterne

- **Atmosfera -> L2:** riuso `mass_exchange_entrainment_2d`; aria ambiente ferma apporta massa ed entalpia, non momento diretto. Quando L2 copre L1 non applicare anche l'entrainment atmosferico diretto a L1. Con `N_LAYERS=1` mantenere comportamento esistente.
- **L1 <-> substrato:** deposizione, erosione e inventari/tessitura del terreno rimangono nei moduli `mass_exchange` attuali, con rollback di substrato e geometria.
- **Lofting di L2:** distinguere materiale esportato in plume, solidi eventualmente reindirizzati verso L1 nelle opzioni del manoscritto, e budget co-ignimbrite; non rimuovere automaticamente anche L1.
- **Reologie:** stesse leggi in `rheology_laws/` con parametri e scelta indipendenti per layer, senza directory reologiche duplicate.
- **Pore pressure:** applicabile a L1; stabilire esplicitamente lo slot inattivo L2 e la semantica di ricostruzione/output.

## 4. Layout, prestazioni e proprietà locali

Layout canonico:

```fortran
q(n_variables, n_layers, nx, ny)
qp(n_physical_variables, n_layers, nx, ny)
```

In Fortran `q(:,layer,j,k)` è contiguo; `q(:,:,j,k)` rappresenta un blocco contiguo della cella. Questo preserva le interfacce fisiche locali con un vettore 1D. Non introdurre la dimensione layer su topografia, geometria e substrato, comuni ai due flussi.

Gli array globali `q_trial`, `qp_trial`, `spatial_rate` e gli snapshot di fluido pertinenti diventeranno coerenti con la rappresentazione rank-4. Lo stato locale del Lie potrà essere `q_cell(n_eqns,n_layers)`; il workspace per worker manterrà uno stato fisico monolayer per gruppo L1/L2 e ne userà due solo quando serve una sorgente interlayer. Newton resta monolayer per le reologie indipendenti: nessun `2*n_eqns` imposto a tutte le celle.

Performance gate: rilevare temporanei Fortran, contiguità, conversioni ripetute, allocazioni e scaling OpenMP. Nessuna generalizzazione n-layer o propagazione dell'indice layer nei kernel non accoppiati.

## 5. Lie globale: gruppi e transazione

Un solo ordine di fasi Lie e un unico `dt` per entrambi i layer. Il registro identifica per ogni gruppo **metodo numerico** e **scope** (`L1`, `L2`, `INTERLAYER`). Le interfacce fisiche EXPLICIT/SEMI_IMPLICIT/IMPLICIT/ANALYTIC restano disponibili secondo il contratto corrente; per un metodo non implementato si produce errore esplicito.

- Gruppo L1: evolve solo q(:,1,j,k).
- Gruppo L2: evolve solo q(:,2,j,k).
- Gruppo INTERLAYER: legge **lo stesso stato pre-gruppo** dei due layer, calcola un unico `Gamma(n_eqns)`, limita rispetto a tutti gli inventari coinvolti, costruisce simultaneamente `q1* = q1 + dt*Gamma` e `q2* = q2 - dt*Gamma`.
- Processi concorrenti che consumano la stessa specie nello stesso gruppo: limiter con budget condiviso, non due limitatori indipendenti.
- Conferma del timestep solo dopo validazione di entrambi i layer e di substrato/geometria/stocastico; fallimento di uno => rollback completo.
- Le sorgenti L1 e L2 realmente indipendenti possono condividere una fase Lie, perché i loro operatori commutano; la **prima implementazione** può usare una sequenza globale semplice e ottimizzare la fusione in seguito.
- Audire rigorosamente il **gruppo 0** attuale, che combina sorgenti esplicite e operatore spaziale, prima di decidere collocazione dei nuovi gruppi.
- Stato asciutto gestito per layer: wet/wet, wet/dry, dry/wet, dry/dry; rimozione per buoyancy e controlli di ammissibilità non devono distruggere l'altro layer.

## 6. Struttura delle directory prevista

```text
src/
  physics/
    mass_exchange/                 # esistente: entrainment, deposito, erosione ...
    rheology_laws/                 # esistente: scelta per-layer
    thermal_laws/                  # esistente
    pore_pressure/                 # esistente: fisica L1
    interlayer/                    # nuova, solo fisica locale nuova
      interlayer_settling_2d.f90
      interlayer_elutriation_2d.f90
      interlayer_drag_2d.f90
      interlayer_exchange_2d.f90  # solo se ha vero lavoro di assemblaggio/limiting
  spatial/
    interlayer_pressure_2d.f90    # nuovo se utile separarlo da moduli esistenti
  numerics/
    source_catalog_2d.f90
    source_registry_2d.f90
    lie_source_executor_2d.f90
    time_integration_2d.f90
  state/
    state_fields.f90
    state_types.f90
```

Niente nuovi moduli di solo forwarding o duplicazione di formule delle leggi esistenti; mantenere separati calcolo fisico e applicazione numerica.

## 7. Roadmap operativa dettagliata

Tutti i blocchi sono **DA INIZIARE**. Le attività TL00–TL03 costituiscono un gate obbligatorio prima di abilitare due layer.

| ID | Blocco e interventi | Criteri di accettazione |
|---|---|---|
| **TL00** | Baseline: censire riferimenti rank-3, allocazioni, snapshot, condizioni al contorno, output, NetCDF, restart, diagnostiche, test; fissare risultati e tempi baseline. | Mappa dipendenze verificata; esempi 1-layer e build di riferimento riproducibili. |
| **TL01** | Stato rank-4 con solo layer 1: `model_layout` e `state_fields`; conversioni, allocazione, init dry, owner dello stato. | Compilazione e tests locali; formula invariata, `N_LAYERS=2` ancora respinto. |
| **TL02** | Migrazione dei consumer: Lie, spatial rate, flux/reconstruction, sorgenti, driver, BC/IC, I/O/diagnostiche. Compatibilità restart gestita esplicitamente. | EXAMPLES, NetCDF, restart round-trip, bounds e OpenMP superati. |
| **TL03** | Certificazione single-layer con confronti numerici, bitwise ove realistico, mass/energy budgets, runtime e memoria, varianti stocastiche e rollback. | Nessuna regressione inspiegata. **Gate bloccante** prima del two-layer. |
| **TL04** | Attivazione controllata stato L1+L2 senza accoppiamento: variabili comuni/inattive, parametri e reologie indipendenti, IC/BC, output. | Round-trip di entrambi gli stati, dry L1/wet L2 e viceversa, accessi verificati. |
| **TL05** | Spatial per layer: applicare HP-PCCU a entrambi, CFL globale, worklist e ricostruzioni per layer, senza cross-pressure. | Test di singolo layer isolato in ciascuno slot, positività e stabilità. |
| **TL06** | Coppia idrostatica interlayer path-conservative: pressione, metriche, clamp e contributi alle interfacce; audit del bound CFL accoppiato. | Lago stratificato a riposo su piano, caratterizzazione residuo su curva, momentum budget e dam-break. |
| **TL07** | Registry Lie con scope L1/L2/INTERLAYER e trasferimento singolo `Gamma`; limiter comune, stati pre-gruppo, rollback di coppia. | Test manufactured exchange conservativo, donor positivity, full rollback, OpenMP. |
| **TL08** | Settling L2 -> L1 con trasferimento specie, massa, momento e termica; riuso delle formule di settling, gestione thin/dry e gas associato. | Conservazione interna delle componenti strettamente scambiate; nessun deposito al suolo; conto esplicito dell'eventuale gas approssimato. |
| **TL09** | Elutriazione L1 -> L2: legge shear e filtri, gas co-lift, stesso rate su due lati. | Bilanci classi/energia/momento e controlli disponibilità del donatore. |
| **TL10** | Drag interfaciale opposto nei due momenti; valutazione esplicita della rigidezza e possibile soluzione implicita solo se giustificata. | Momento totale invariato, drag nullo quando appropriato, stabilità della rilassazione. |
| **TL11** | Riutilizzo entrainment atmosfera -> L2, deposizione/erosione L1 <-> terreno, pore pressure L1, reologie distinte, lofting e budget plume. | Nessun doppio entrainment; bilanci ambiente/terreno coerenti; output diagnostici chiari. |
| **TL12** | Integrazione scientifica: benchmark del manoscritto, colonne chiuse/aperte, rest states, topografia curva, sorgenti simultanee, wet/dry, multi-thread. | Matrice benchmark documentata, differenze scientifiche spiegate, regressioni single-layer ancora verdi. |
| **TL13** | Documentazione finale e configurazione: input, output/restart versionati, esempi, limitazioni, aggiornamento roadmap generale e rilascio di `N_LAYERS=2`. | Review e autorizzazione esplicita prima di rendere la modalità two-layer pubblicamente disponibile. |

**Dipendenze:** TL00 → TL01 → TL02 → TL03 (gate) → TL04 → TL05 → TL06; TL07 deve precedere i termini di scambio TL08–TL10; TL11 integra sorgenti locali/esterne; TL12 richiede entrambi i percorsi spatial e Lie; TL13 conclude. È possibile sviluppare singoli test di fisica in parallelo, ma non saltare i gate di integrazione.

## 8. Verifiche scientifiche e di regressione

1. In `N_LAYERS=1`: equazioni, reologie, ordine e metodi Lie, densità/buoyancy, dry state, output e crash/restart devono riprodurre il baseline. Eseguire le famiglie EXAMPLES effettivamente presenti (inclusi ETNA/ETNA_LAVA e TAAL/TALL secondo nome verificato), fino al tempo finale; documentare le differenze agli ultimi bit se non evitabili.
2. In `N_LAYERS=2`: L2 assente deve ridurre a L1-only; verificare anche L1 assente/L2 presente secondo le closure consentite.
3. Lago stratificato su piano, letto curvo, fronte, dam-break con partner vuoto, CFL rispetto allo spettro accoppiato e well-balancing.
4. Colonna di solo settling interlayer: `Delta(q1+q2)=0` per le componenti interne e assenza assoluta di modifica a terreno; limiti di inventario per ogni specie.
5. Solo deposito L1 -> terreno: inventario terra cresce e L1 diminuisce secondo porosità/convenzioni; L2 inalterato.
6. Elutriazione: massa e momento in coppia, entalpia del donatore, gas co-lift con budget esplicito delle eccezioni.
7. Drag: effetto nullo per velocità uguali o layer assente, momento totale invariato.
8. Entrainment L2: massa ed entalpia ambientali aggiunte, momento diretto nullo, nessun doppio entrainment L1 coperto.
9. Rifiuto di qualsiasi candidato: nessun commit parziale L1/L2/substrato/geometria/stocastico.
10. Fortran memory/bounds checking, assenza di temporanei imprevisti, profiling OpenMP e riproducibilità.

## 9. Questioni aperte da risolvere prima dei blocchi corrispondenti

- Layout effettivo di pressure-of-pore e altri scalari validi in un solo layer, footprint delle equazioni inattive, primitive supplementari e diagnostiche.
- Eventuali differenze deliberatamente non conservative del manoscritto (gas di packing), senza violare il contratto antisimmmetrico dei trasferimenti strettamente interni.
- Semantica e ordine del gruppo 0 Lie quando le due reologie e gli scambi hanno metodologie differenti.
- Criterio CFL/hyperbolicity del sistema accoppiato, che potrebbe non essere garantito dalla sola velocità d'onda per layer.
- Posizione/spessore del layer superiore sopra L1 parzialmente asciutto, topografia condivisa e boundary conditions.
- Input per la scelta delle reologie e dei gruppi indipendenti e retrocompatibilità input single-layer.
- Schema NetCDF two-layer e versioning dei restart senza cambiare silenziosamente il formato storico.
- Differenza tra stato sorgente conservativo, trasferimento di entalpia termica ed eventuale dissipazione di energia cinetica.
- Necessità di metodo implicito interlayer se lo schema esplicito risulta troppo rigido.
- Scelta finale delle opzioni scientifiche del manoscritto (lofting, co-lift, deposizione/porosità) prima di considerare completa l'implementazione.

## 10. Regole di esecuzione

- **La creazione di questo file non autorizza modifiche al codice.**
- Uno o pochi sub-blocchi verificabili per commit; non unire migrazione rank-4 e modifica della fisica.
- Non dichiarare completato TL03 senza confronto dei risultati e test realmente eseguiti.
- Le equazioni del manoscritto sono riferimento scientifico, ma i suoi schemi IMEX-RK storici non sostituiscono l'architettura Lie corrente.
- Dopo ciascun blocco completato, registrare SHA, files cambiati, test/esiti, differenze e attività ancora aperte; riportare il progresso sintetico nel tracker generale.
- Il lavoro `stateNEW` rimane in pausa: non riavviare refactoring non richiesto.

## 11. Riferimenti

- Manoscritto fornito dall'utente: *IMEX_SfloW2D_v3: A two-layer shallow-water model for pyroclastic density currents with pore-pressure transport, polydisperse settling and shear-driven elutriation*, sezioni 2.1, 2.5–2.7, 3–4.
- [Roadmap generale](IMEX_SfloW2D_Development_Roadmap.md).
- Moduli correnti: `src/core/model_layout_2d.f90`, `src/state/state_fields.f90`, `src/numerics/time_integration_2d.f90`, `src/numerics/lie_source_executor_2d.f90`, `src/numerics/source_registry_2d.f90`, `src/physics/source_terms_2d.f90`, `src/physics/mass_exchange/`, `src/spatial/`.

**Prossima attività, solo dopo revisione e approvazione del piano:** TL00, audit/baseline. Nessuna modifica sorgente eseguita.
