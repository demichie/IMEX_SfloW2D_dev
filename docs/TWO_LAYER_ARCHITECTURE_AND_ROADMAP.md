# IMEX-SfloW2D — Architettura single-layer / two-layer e roadmap

**Riferimento:** 7 novembre 2026  
**Ambito:** architettura di produzione e piano di implementazione per `N_LAYERS=1` e `N_LAYERS=2`.  
**Stato implementativo:** il codice in produzione è single-layer; la modalità two-layer non è ancora implementata. I blocchi TL00–TL13 sono aperti.  
**Tracciamento:** [Roadmap generale](IMEX_SfloW2D_Development_Roadmap.md) per lo stato complessivo del progetto; questo documento è la specifica tecnica del percorso two-layer.

## 1. Obiettivo e configurazioni del modello

IMEX-SfloW2D supporta attualmente un modello depth-averaged a un layer. L'architettura di riferimento prevede due sole configurazioni:

- **`N_LAYERS=1`**: modello esistente, con fisica, discretizzazione HP-PCCU, Lie splitting e risultati di riferimento da preservare.
- **`N_LAYERS=2`**: **L1**, corrente inferiore densa (undercurrent), e **L2**, corrente superiore diluita (ash cloud), secondo il manoscritto *IMEX_SfloW2D_v3*. La sua formulazione fisica viene integrata nell'attuale architettura numerica Lie di primo ordine, senza reintrodurre lo schema IMEX-RK storico.

Non è previsto un framework a un numero arbitrario di layer, né una riscrittura generalizzata delle leggi costitutive. Il progetto mantiene la separazione tra stato, fisica locale, operatori spaziali e integrazione temporale. Il workstream indipendente di ulteriori ottimizzazioni dello stato resta sospeso durante questa implementazione.

## 2. Stato corrente del codice

| Componente | Architettura corrente |
|---|---|
| `src/core/model_layout_2d.f90` | Metadati comprendenti `n_layers`; `supported_layer_count` ammette ancora soltanto un layer. |
| `src/state/state_fields.f90` | Stato del flusso `q`, `qp` attualmente a rango 3; inventari del substrato, diagnostiche e snapshot posseduti dallo stato. |
| `src/numerics/time_integration_2d.f90` | Integrazione Lie di primo ordine; unico ciclo temporale e rollback di stato, substrato e geometria. |
| `src/numerics/source_catalog_2d.f90`, `source_registry_2d.f90`, `lie_source_executor_2d.f90` | Catalogo delle sorgenti e relativi metodi; gruppi Lie ordinati; scratch locali per worker, controlli e transazioni. |
| `src/physics/source_terms_2d.f90` e `src/physics/mass_exchange/` | Kernel fisici cell-local, sorgenti conservative, deposizione, erosione, entrainment e coordinamento degli inventari. |
| `src/spatial/nonconservative_2d.f90` | Motore **già generale** per integrali di percorso: callback dell'integrando, stato esteso e tangente al percorso; Gauss a 1, 2 o 3 punti; percorso personalizzabile. |
| `src/spatial/pccu_2d.f90` | Idrostatica single-layer e assemblaggio central-upwind/path-conservative; usa il motore generico. |
| `src/spatial/hyperbolic_2d.f90`, `reconstruction_2d.f90`, `spatial_operator_2d.f90` | Tracce ricostruite, percorsi alle interfacce e interni alle celle, velocità caratteristiche, CFL e assemblaggio del residuo spaziale. |
| `src/Makefile.am` | Elenco esplicito dei moduli compilati. |

La distinzione tra `n_variables`, `n_physical_variables` e `n_equations` resta esplicita: le rispettive dimensioni non devono essere considerate intercambiabili.

## 3. Stato globale, stato fisico locale e substrato

### 3.1 Layout dei campi fluidi

Il layout globale è:

```fortran
q (n_variables,          n_layers, nx, ny)
qp(n_physical_variables, n_layers, nx, ny)
```

La dimensione di componente precede quella di layer. In Fortran, `q(:,layer,j,k)` e `qp(:,layer,j,k)` sono vettori contigui, utilizzabili dagli attuali kernel locali. L'accesso con entrambe le dimensioni di componente e layer permette di gestire la coppia di stati di una cella senza rendere i kernel monolayer dipendenti dalla griglia.

I workspace globali di evoluzione (`q_trial`, `qp_trial`, `spatial_rate`, cache/snapshot fluidi pertinenti) seguono lo stesso principio. Topografia, geometria e substrato sono **unici e condivisi**: non acquistano una dimensione layer.

La migrazione al rango 4 viene eseguita inizialmente con `n_layers=1`; `N_LAYERS=2` resta rifiutato fino alla validazione della nuova modalità. La migrazione comprende tutte le dipendenze: allocazioni, inizializzazione, conversioni, condizioni iniziali e al contorno, halo, ricostruzione, operatori, diagnostiche, NetCDF, restart e test.

### 3.2 Variabili dei due layer

Le posizioni delle variabili e delle equazioni rimangono identiche fra L1 e L2. Una grandezza definita fisicamente soltanto in L1, come l'eccesso di pressione interstiziale, occupa uno slot inattivo in L2: maschera, inizializzazione asciutta, conversioni, trasporto, output e metodi numerici devono rispettare questa inattività.

Le chiusure e i parametri reologici possono differire fra L1 e L2; si riutilizzano i moduli fisici esistenti, senza copie per-layer. Gli stati fisici, le valutazioni di una singola cella e gli scratch Newton restano monolayer. Solo i processi accoppiati richiedono contemporaneamente due stati fisici locali.

### 3.3 Deposito unico indipendente dal tipo di flusso

**`deposit_tot(:,:)` è lo spessore canonico del deposito, espresso in metri, per tutte le configurazioni fisiche supportate.** Le configurazioni di flusso granulare e di lava sono incompatibili nella stessa simulazione: **non** esiste una somma contemporanea di deposito granulare e lavico.

- **Flusso granulare:** `deposit(:,:,:)` conserva l'inventario volumetrico delle classi solide; lo spessore totale del deposito è derivato mediante la porosità del substrato, `deposit_tot = SUM(deposit,DIM=3)/(1-erodible_porosity)`, secondo la convenzione esistente.
- **Flusso di lava:** `deposit_tot(:,:)` rappresenta direttamente lo spessore della lava immobilizzata, senza conversione mediante la porosità granulare.
- **`lava_deposit(:,:)`**, presente nel codice corrente, è ridondante rispetto alla rappresentazione canonica dello spessore e viene eliminato solo dopo un audit completo di dipendenze e conservazione: Lie, snapshot/rollback, NetCDF, restart, diagnostiche e topografia.
- **`erosion` / `erosion_tot`:** mantengono la semantica distinta dell'erosione del substrato. L'eventuale evoluzione della quota topografica non coincide automaticamente con lo spessore cumulativo del deposito: evitare doppio conteggio.

Restano **distinte le due leggi fisiche**: deposizione granulare per sedimentazione e immobilizzazione della lava per rate prescritto o dipendente da temperatura/viscosità. Entrambe producono gli opportuni trasferimenti di massa ed energia e l'aggiornamento geometrico mediante il coordinamento esistente. La variabile comune è lo **spessore risultante**, non la legge di trasferimento. `deposit_tot` ha un proprietario e una convenzione coerenti in tutta la simulazione, inclusi rollback e output.

La razionalizzazione del deposito è un blocco autonomo di compatibilità single-layer, senza cambiare le formule fisiche.

## 4. Partizione matematica e integrazione temporale

Per i due vettori conservativi mobili `q1` e `q2`:

```text
dq1/dt = H1(q1,q2,B) + S_local1(q1) + Gamma(q1,q2) + S_external1
dq2/dt = H2(q1,q2,B) + S_local2(q2) - Gamma(q1,q2) + S_external2
```

- **`H1,H2`:** flussi inerziali, pressione idrostatica, topografia e prodotti spaziali non conservativi, gestiti da HP-PCCU.
- **`S_local`:** reologia, termica, pressione interstiziale e altre sorgenti che agiscono su un solo layer.
- **`Gamma`:** trasferimenti **strettamente interni e conservativi** fra i due layer.
- **`S_external`:** scambi con atmosfera, substrato, iniezioni e rimozione/esportazione dal dominio.

Il codice usa **un solo Lie splitting globale**, con sequenza ordinata di gruppi, timestep condiviso, scelta CFL globale, validazione e rollback atomico. Ogni gruppo ha un metodo e uno scope `L1`, `L2` o `INTERLAYER`. Il gruppo 0 e la combinazione fra operatori spaziali e sorgenti mantengono la loro semantica attuale fino alla verifica delle dipendenze.

Le quattro interfacce di trattamento fisico (EXPLICIT, SEMI_IMPLICIT, IMPLICIT, ANALYTIC) mantengono il contratto di disponibilità: le modalità non implementate producono errore esplicito, non fallback impliciti o no-op. I gruppi monolayer indipendenti possono essere eseguiti separatamente; un'eventuale fusione delle fasi commutanti è un'ottimizzazione successiva.

## 5. Termini non conservativi e discretizzazione HP-PCCU

### 5.1 Motore numerico esistente

Il motore canonico è **`src/spatial/nonconservative_2d.f90`**, attraverso:

```fortran
CALL eval_path_contribution( &
    direction, stateL, stateR, integrand_provider, &
    path_contribution, n_quad=3, path_provider=... )
```

La sua callback riceve lo **stato esteso** `W(s)` e la **tangente** `dW/ds` lungo il percorso, e restituisce l'integrando vettoriale. Il modulo non deve sapere se `W` contiene conservative, primitive, tracce ricostruite, geometria o variabili di entrambi i layer.

La forma `f(W) * d(g(W))/dx` è un caso particolare già rappresentabile:

```text
I(W,dW/ds) = f(W) * d(g(W))/ds
D_LR       = integral_0^1 I(W(s),dW/ds) ds
```

Lo stesso vale per somme `SUM_m f_m(W) * d(g_m(W))/ds`. **Non è prevista una nuova API fattorizzata né un secondo motore di quadratura:** le callback fisiche restituiscono direttamente il prodotto completo. Le derivate lungo il percorso sono definite analiticamente dalla callback e dalla scelta dello stato/percorso, senza differenze finite automatiche.

Il percorso è una scelta numerica esplicita. Interpolare variabili conservative, primitive oppure grandezze estese non è in generale equivalente. Per ogni termine occorre utilizzare la stessa convenzione fra tracce di faccia, percorsi interni e trattamento dry/wet.

### 5.2 Coppia idrostatica interlayer

Definiti `G = gcell`, `K2 = rho_m2 * gprime2` e `L = K2*h2`, i due termini sorgente non conservativi lungo la direzione orizzontale `j` sono:

```text
L1:  S_inter1 = -G*h1 * d_j(L)
L2:  S_inter2 = -L    * d_j(G*h1)
```

Lungo il percorso `W(s)`:

```text
I1 = -G*h1 * dL/ds
I2 = -L    * d(G*h1)/ds
```

Si calcolano **insieme, sulle stesse tracce e con la stessa quadratura**, disponendo nei rispettivi componenti del momento di L1 e L2. In assenza di clamp:

```text
I1 + I2 = -d(G*h1*L)/ds
D1 + D2 = -[G*h1*L]_L^R
```

La seconda identità costituisce un vincolo di bilancio per il calcolo dei percorsi, con l'accuratezza della quadratura prescelta. Il clamp di galleggiamento sul membro di L2 rompe deliberatamente questa identità; la sua applicazione e i relativi difetti di bilancio devono essere espliciti e diagnosticabili, secondo il manoscritto.

Il fattore metrico ha una collocazione **non intercambiabile**: all'esterno della derivata in L1, all'interno in L2. La topografia deve essere inclusa nella verifica dell'intero operatore idrostatico, non solo dei due contributi di accoppiamento.

Con le medesime proprietà idrostatiche, la somma delle equazioni continue deve ridursi al termine di pressione single-layer di spessore `h=h1+h2`, **compreso il contributo topografico**. La riduzione non richiede `h1=h2`, e non vale generalmente per densità diverse. L'equivalenza delle equazioni continue non comporta automaticamente uguaglianza bitwise dei residui discreti.

Secondo la formulazione scelta nel manoscritto, l'ash cloud su letto curvo può presentare un residuo di riposo associato a `-K2*h1*h2*d_j(G)`. La verifica distingue fondo piano, pendenza costante e topografia curva e non promette un well-balancing esatto laddove la specifica scientifica indica un residuo.

### 5.3 Integrazione nello schema esistente

`pccu_2d.f90` e `hyperbolic_2d.f90` già utilizzano `eval_path_contribution` per la pressione single-layer; il percorso di produzione single-layer rimane invariato durante la migrazione.

Il two-layer richiede:
1. ricostruzione coerente delle tracce L1 e L2 **alla medesima interfaccia**;
2. uno stato esteso accoppiato e una callback di integrando per la coppia di pressione interlayer, senza indici `j,k` nella fisica;
3. valutazione dei percorsi di faccia e interni alle celle e assemblaggio dei contributi nei due momenti;
4. trattamento esplicito delle quattro combinazioni wet/dry, del clamp, delle proprietà termodinamiche e del percorso sulla topografia;
5. velocità caratteristiche e vincolo CFL coerenti con il sistema accoppiato, non limitati a un semplice massimo acritico di due CFL monolayer.

I nuovi contributi idrostatici appartengono all'**operatore spaziale**, non a un gruppo sorgente Lie locale. Il nucleo di quadratura `nonconservative_2d.f90` non richiede modifiche.

## 6. Sorgenti fisiche accoppiate ed esterne

### 6.1 Trasferimenti strettamente interni

Un trasferimento conservativo fra L1 e L2 produce **un solo vettore firmato** `Gamma(n_eqns)`, positivo da L2 a L1. Il coordinatore calcola i due candidati dallo **stesso stato pre-gruppo**:

```text
q1_candidate = q1 + dt*Gamma
q2_candidate = q2 - dt*Gamma
```

Il limite di inventario è condiviso: non si applicano limitatori indipendenti che rompano l'antisimmmetria. Accettazione o rifiuto sono atomici per entrambi i layer e per gli eventuali inventari condivisi. La conservazione riguarda massa, specie, momento e **variabile termica conservativa** effettivamente trasferita; non implica automaticamente conservazione dell'energia cinetica o una legge per la dissipazione.

- **Settling L2 → L1:** trasferimento delle classi solide, con velocità di interfaccia e temperatura/entalpia del donatore, limitato dal materiale disponibile; entrambi i layer restano mobili. **Non** aggiorna `deposit_tot`, il substrato o la topografia.
- **Elutriazione L1 → L2:** trasferimento da shear, con selezione granulometrica, velocità ed entalpia del donatore L1 e gas eventualmente sollevato; usa la stessa contabilità di coppia.
- **Drag interfaciale:** `D = f*rho_m2*|u2-u1|*(u2-u1)`, applicato come `+D` al momento di L1 e `-D` a quello di L2; nullo se manca uno dei due strati e quando le velocità coincidono. La scelta iniziale è esplicita, con trattamento implicito solo se richiesto dalla rigidezza.

**Eccezione scientifica:** l'approssimazione di co-trasporto del gas interstiziale a packing descritta nel manoscritto può produrre una piccola componente non esattamente antisimmmetrica. Essa va esplicitata come correzione esterna o difetto diagnostico distinto; non va nascosta in `Gamma` e non deve rendere false le prove di conservazione dei termini strettamente interni.

### 6.2 Scambi con l'ambiente e il terreno

- **Entramento atmosferico:** riuso di `src/physics/mass_exchange/mass_exchange_entrainment_2d.f90` per l'atmosfera → L2. Aria ambiente ferma apporta massa ed entalpia, senza momento orizzontale diretto. Con L2 sopra L1, L1 non riceve un secondo entrainment atmosferico; il caso monolayer conserva il comportamento esistente.
- **Deposizione ed erosione L1 ↔ substrato:** riuso delle leggi e degli inventari granulometrici esistenti; `deposit_tot` è il riferimento unico per lo spessore del deposito. Nessuna modifica al terreno dal solo settling L2 → L1.
- **Lava → substrato:** distinta legge di immobilizzazione, ma lo stesso `deposit_tot` canonico; nessuna somma fisica con depositi granulari.
- **Pore pressure:** chiusura propria di L1 e slot non attivo L2.
- **Reologia e termica:** moduli esistenti riutilizzati per ciascun layer con parametri distinti.
- **Lofting L2:** rimozione/esportazione e possibile reindirizzamento dei solidi secondo le opzioni del manoscritto, senza eliminare impropriamente L1; bilanci di massa e plume identificabili.

## 7. Proprietà del Lie globale e parallelismo

Il registro associa a ciascun gruppo identità fisica, metodo numerico, footprint delle equazioni e scope (`L1`, `L2`, `INTERLAYER`). I gruppi interlayer acquisiscono lo stesso snapshot dei due stati prima del calcolo; il rifiuto di un singolo candidato comporta rollback dell'intero timestep.

L'implementazione conserva:
- kernel in `src/physics/` **cell-local**, senza `j,k` o accessi a campi globali;
- scratch fisici e Newton per worker, riutilizzati e non allocati per ogni cella;
- assenza di un solver Newton globale `2*n_eqns` introdotto senza necessità;
- un solo ordine Lie, un solo `dt`, un solo commit dello stato;
- transazioni coerenti di fluido, substrato, geometria, diagnostiche e variabili stocastiche;
- OpenMP sui loop più ampi possibili, con controlli di contiguità e temporanei Fortran.

Per ragioni di robustezza e riproducibilità, la modalità `N_LAYERS=1` è un percorso esplicitamente verificato in tutte le fasi di sviluppo; nessuna generalizzazione inutile della fisica locale.

## 8. Organizzazione dei moduli

```text
src/
  state/
    state_fields.f90              # fluidi, substrato comune, snapshot
    state_types.f90               # stati fisici locali
    state_conversions.f90         # conversioni conservative/primitive
  physics/
    mass_exchange/               # leggi esistenti: deposito, erosione, entrainment...
    rheology_laws/               # leggi esistenti, configurabili per layer
    thermal_laws/                # leggi esistenti
    interlayer/                  # nuova fisica locale accoppiata
      interlayer_settling_2d.f90
      interlayer_elutriation_2d.f90
      interlayer_drag_2d.f90
      interlayer_exchange_2d.f90 # solo se contiene assemblaggio sostanziale
  spatial/
    nonconservative_2d.f90       # motore generico esistente, invariato
    pccu_2d.f90                  # percorsi e algebra PCCU
    hyperbolic_2d.f90            # ricostruzioni, facce e residui
    interlayer_pressure_2d.f90   # callback fisica, se utile come modulo separato
  numerics/
    source_catalog_2d.f90
    source_registry_2d.f90
    lie_source_executor_2d.f90
    time_integration_2d.f90
```

Non introdurre wrapper di solo inoltro, duplicazioni delle leggi fisiche, un nuovo framework per prodotti non conservativi o estrazioni di formule elementari senza beneficio reale. I moduli nuovi sono inseriti nell'elenco esplicito di `src/Makefile.am` solo quando vengono effettivamente implementati.

## 9. Roadmap operativa

Tutti i blocchi seguenti sono **DA INIZIARE**. Ogni blocco è un intervento separato e verificabile; i primi gate proteggono il single-layer funzionante.

| ID | Obiettivo tecnico | Criteri di accettazione |
|---|---|---|
| **TL00 — baseline e audit** | Identificare commit single-layer validato, fissare tag annotato; definire branch di sviluppo separato; inventario delle dipendenze rank-3 e dei consumatori di `deposit_tot`/`lava_deposit`; misurare risultati, output, tempi e memoria degli esempi. | Riferimento immutabile, mappa dipendenze e test riproducibili. La protezione di `master` e lo stato del branch sono verificati, non presunti. |
| **TL00a — deposito canonico** | Rendere `deposit_tot` lo spessore unico per granulari **oppure** lava; mantenere le leggi fisiche separate; audit e rimozione dell'inventario `lava_deposit` solo se davvero ridondante; allineare rollback, topografia, I/O e restart. | Deposito/erosione e geometria invariati fisicamente; EXAMPLE_ETNA_LAVA e casi granulari confrontati; nessun doppio conteggio o conversione di porosità nel ramo lava. |
| **TL01 — layout rank-4, un layer** | Introdurre il secondo indice layer nei campi fluidi globali e nei metadati, lasciando ammesso soltanto `N_LAYERS=1`; preservare i kernel cell-local. | Build e test di stato/convertitori, nessuna modifica intenzionale delle formule. |
| **TL02 — migrazione consumer globali** | Adeguare Lie, workspace, spatial, ricostruzione, BC/IC, driver, output, NetCDF, restart, diagnosi, snapshot e test. | EXAMPLES, restart round-trip, bounds, OpenMP e formati documentati. |
| **TL03 — gate single-layer** | Confrontare evoluzioni, campi, budget, performance, memoria, rami wet/dry, buoyancy, rollback e risultati finali con il baseline taggato. | Nessuna regressione inspiegata; bitwise quando realisticamente ottenibile. Nessuna abilitazione two-layer prima della chiusura del gate. |
| **TL04 — stato two-layer inerte** | Slot L1/L2, variabili inattive, configurazioni/reologie distinte, dry state, IC/BC e output per layer, senza accoppiamento. | Round-trip e validità per wet/wet, wet/dry, dry/wet, dry/dry. |
| **TL05 — HP-PCCU per layer** | Ricostruzione e flussi di ciascun layer, dominio/halo coerenti, timestep globale preliminare. | Ciascuno slot riproduce il comportamento del singolo flusso in assenza di accoppiamento. |
| **TL06 — pressione interlayer** | Utilizzare **l'attuale `nonconservative_2d` invariato** con callback della coppia idrostatica; assemblare percorsi di faccia e interni, topografia, metriche, clamp e CFL accoppiato. | Identità dei due path, riduzione single-layer con topografia, equilibrio su fondo piano/inclinato, residuo curvo caratterizzato, casi wet/dry, test di stabilità. |
| **TL07 — Lie con scope e trasferimento di coppia** | Registro L1/L2/INTERLAYER, `Gamma` unico, limiter condiviso, due candidati simultanei, validazione e rollback atomici, rispetto del gruppo 0. | Test manufactured di conservazione, positività dei donatori, rollback e OpenMP. |
| **TL08 — settling L2 → L1** | Riutilizzare chiusure di settling; trasferire specie, massa, momento ed energia con donatore L2; considerare il gas packing e i limiti thin/dry. | Conservazione dei trasferimenti interni, eccezioni diagnosticate, **nessun deposito al terreno**. |
| **TL09 — elutriazione L1 → L2** | Implementare legge shear, filtro granulometrico, co-lift e rate di coppia. | Bilanci conservativi, limite d'inventario e trasporto termico/momento coerenti. |
| **TL10 — drag interfaciale** | Forza quadratica uguale e opposta sui momenti; verificare rigidità esplicita. | Momento totale conservato, forza nulla nei casi limite, stabilità. |
| **TL11 — sorgenti esterne e chiusure** | Entrainment atmosfera→L2, L1↔terreno, pressione interstiziale L1, reologie differenziate, lofting L2, budget plume, corretto uso del deposito canonico. | Nessun doppio entrainment o modifica del deposito da soli scambi interlayer; budget e output diagnostici coerenti. |
| **TL12 — validazione scientifica** | Benchmark del manoscritto, resting states, topografia curva, coppie di processi simultanei, wet/dry, CFL, restart, OpenMP e confronto con single-layer. | Matrice test e differenze scientifiche documentate, baseline 1-layer ancora superato. |
| **TL13 — finalizzazione** | Documentazione di configurazione, output/restart e limiti scientifici; aggiornamento tracker generale; verifica integrata e abilitazione `N_LAYERS=2`. | Review finale, test superati e decisione esplicita di abilitazione. |

**Dipendenze principali:** TL00 → TL00a → TL01 → TL02 → TL03 (gate) → TL04 → TL05 → TL06; TL07 precede TL08–TL10; TL11 completa l'interazione con ambiente e substrato; TL12 e TL13 concludono. Eventuali test fisici isolati possono essere preparati prima senza saltare i gate di integrazione.

## 10. Matrice di validazione

### Single-layer
- Build, Actions autorizzate, test di unità, conversioni e gruppi Lie.
- Esecuzione completa degli esempi disponibili, inclusi scenari granulari/gas e lava; verifica dei campi finali, massa, energia, volume, runout, topografia e deposito.
- Verifica esplicita del ramo lava con `deposit_tot` e della porosità solo dove applicabile.
- NetCDF, restart, rollback, OpenMP e confronto prestazionale con baseline immutabile.

### Two-layer: operatori spaziali
- Un layer vuoto: recupero della formulazione pertinente di layer isolato.
- Identità di percorso `D1+D2 = -[G*h1*K2*h2]_L^R` senza clamp.
- `K1=K2`: somma dei termini idrostatici **inclusi i contributi topografici** equivalente al single-layer `h=h1+h2`, con spessori arbitrari.
- Riposo su topografia piana e a pendenza costante; residuo di well-balancing misurato su fondo curvo.
- Clamp per buoyancy, asciutto/bagnato, dam-break e bound CFL del sistema accoppiato.

### Two-layer: scambi e transazioni
- Settling L2→L1 con invarianti di massa/specie/momento/termica e `deposit_tot` invariato.
- Elutriazione L1→L2 con filtro e gas co-lift; difetto gas packing separatamente contabilizzato dove previsto.
- Drag con momento totale invariato e limiti corretti a velocità uguali o layer assente.
- Entrainment atmosferico aggiunto una sola volta a L2; nessun momento diretto dell'aria ferma.
- Deposizione L1→substrato e lava→substrato su configurazioni **mutuamente esclusive**, attraverso il medesimo `deposit_tot`.
- Fallimento di una cella/gruppo: nessun commit parziale dei due layer, del substrato o della geometria.
- Verifiche di temporanei, contiguità, scalabilità OpenMP e bilanci anche in scenari combinati.

## 11. Vincoli scientifici e questioni da chiudere nei rispettivi blocchi

1. Definire esattamente le variabili L2 inattive e le loro semantiche per input, trasporto e restart.
2. Specificare il percorso esteso e la rappresentazione delle grandezze termodinamiche nella quadratura accoppiata, senza modificare il motore `nonconservative_2d`.
3. Verificare velocità caratteristiche e iperbolicità del sistema two-layer, inclusa l'influenza delle interfacce e del clamp.
4. Distinguere identità dell'operatore continuo, precisione di quadratura e well-balancing del residuo discretizzato su topografia curva.
5. Trattare esplicitamente le eventuali chiusure non conservative del gas di packing.
6. Definire condizioni al contorno e topografia condivisa nelle combinazioni wet/dry e negli stati con layer superiore presente senza L1.
7. Stabilire versioni e retrocompatibilità di restart/NetCDF e la semantica di `deposit_tot` dopo la razionalizzazione del deposito.
8. Verificare l'ordine dei gruppi Lie (in particolare gruppo 0) e il limite di stabilità per eventuali termini interlayer rigidi.
9. Confermare le opzioni scientifiche specifiche del manoscritto per lofting, co-lift, porosità e deposizione prima della validazione finale.

## 12. Protocollo operativo

- Ogni blocco ha un commit o una PR circoscritta, con test/esiti documentati e riferimento al baseline.
- `master` è il riferimento single-layer: lo sviluppo two-layer avviene su un branch dedicato; la protezione effettiva del branch e il tag del baseline vengono verificati in TL00.
- Nessun cambiamento di formula fisica è introdotto soltanto per adattare il layout. Gli aggiustamenti scientifici richiedono giustificazione e validazione autonoma.
- Le ottimizzazioni non pertinenti al two-layer e il precedente workstream di ulteriore refactoring dello stato restano fuori ambito.
- La compilazione da sola non chiude un gate scientifico; i benchmark e il confronto degli output hanno valore separato.
- Il progresso di TL00–TL13 viene riportato sinteticamente nella roadmap generale, mantenendo questo file come specifica tecnica consolidata.

## 13. Riferimenti

- Manoscritto scientifico: *IMEX_SfloW2D_v3: A two-layer shallow-water model for pyroclastic density currents with pore-pressure transport, polydisperse settling and shear-driven elutriation*, in particolare sezioni 2.1–2.2, 2.5–2.7, 3–4.
- [Roadmap generale del progetto](IMEX_SfloW2D_Development_Roadmap.md).
- Sorgenti esistenti: `src/state/`, `src/physics/mass_exchange/`, `src/spatial/nonconservative_2d.f90`, `src/spatial/pccu_2d.f90`, `src/spatial/hyperbolic_2d.f90`, `src/numerics/lie_source_executor_2d.f90`.
