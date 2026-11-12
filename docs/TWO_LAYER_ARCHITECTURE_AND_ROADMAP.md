# IMEX-SfloW2D — Architettura single-layer / two-layer e roadmap

**Riferimento:** 7 novembre 2026  
**Ambito:** architettura di produzione e piano di implementazione per `N_LAYERS=1` e `N_LAYERS=2`.  
**Stato implementativo (11 novembre 2026):** TL04 avviato con una prima configurazione fase/reologia per-layer, ancora senza coefficienti fisici locali; `N_LAYERS=2` non è ancora abilitato. Sul branch di sviluppo `state_type%q/qp` e i workspace fluidi di Lie sono migrati a rank-4 con uso esplicito del layer 1 nei consumer principali, senza wrapper; il solver rimane funzionalmente single-layer: la build CI optimized/debug del commit `155a0b6` è riuscita, mentre le regressioni funzionali rimangono da eseguire. Il deposito unico lavico è stato verificato dall'utente sul desktop prima della migrazione rank-4; compilazione optimized/debug e regressioni di fisica e numerica sono superate su GitHub. Gli esempi sono stati verificati dall'utente sul desktop dopo la migrazione; restano aperti test documentati di restart, NetCDF, OpenMP e confronto baseline. I gate di accettazione TL01–TL03 rimangono aperti.  
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
| `src/state/state_fields.f90` | Stato del flusso `q`, `qp` migrato a rango 4 sul branch, con `N_LAYERS=1`; inventari del substrato, diagnostiche e snapshot posseduti dallo stato. |
| `src/numerics/time_integration_2d.f90` | Integrazione Lie di primo ordine; trial e residuo spaziale rank-4, esecuzione del layer 1 attraverso sezioni cell-local, unico ciclo temporale e rollback di stato, substrato e geometria. |
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
- **`lava_deposit(:,:)`** non è più un inventario autonomo nell'implementazione TL00a sul branch; gli usi sono migrati a \`deposit_tot\`. La verifica di conservazione, rollback, NetCDF, restart, diagnostiche e topografia costituisce un gate ancora aperto. I file NetCDF creati dal nuovo schema espongono `deposit_tot`; l'append a un NetCDF precedente privo del campo è rifiutato esplicitamente, senza perdita silenziosa della storia.
- **`erosion` / `erosion_tot`:** mantengono la semantica distinta dell'erosione del substrato. L'eventuale evoluzione della quota topografica non coincide automaticamente con lo spessore cumulativo del deposito: evitare doppio conteggio.

Restano **distinte le due leggi fisiche**: deposizione granulare per sedimentazione e immobilizzazione della lava per rate prescritto o dipendente da temperatura/viscosità. Entrambe producono gli opportuni trasferimenti di massa ed energia e l'aggiornamento geometrico mediante il coordinamento esistente. La variabile comune è lo **spessore risultante**, non la legge di trasferimento. `deposit_tot` ha un proprietario e una convenzione coerenti in tutta la simulazione, inclusi rollback e output.

La razionalizzazione del deposito è un blocco autonomo di compatibilità single-layer, senza cambiare le formule fisiche.

### 3.4 Contratto di aggiornamento e derivazione di \`deposit_tot\`

\`deposit_tot\` è lo **spessore cumulativo positivo di materiale immobilizzato** rispetto allo stato iniziale, non il cambiamento netto di quota topografica e non l'inventario erodibile. È un campo del substrato privo dell'indice layer. L'erosione ha il suo inventario e la propria semantica; l'eventuale incorporazione del deposito nell'erodibile deve seguire \`erodible_deposit_flag\` e le convenzioni attuali.

Nell'implementazione TL00a \`state_fields:refresh_substrate_totals\` mantiene \`deposit_tot\` come inventario autorevole quando \`lava_flag=T\`, anche con \`n_solid=0\`; nel ramo granulare lo deriva da \`deposit\`. Il Lie incrementa direttamente il candidato \`deposit_tot\`. Il record finale del restart v7 conserva la posizione e il significato binario storico, ma per lava viene letto/scritto mediante \`deposit_tot\`. Questa compatibilità deve essere convalidata con un round-trip eseguito. Per ottenere il comportamento canonico sono quindi indispensabili **due rami fisici esclusivi e verificabili**:

\`\`\`text
granular carrier: deposit_tot = SUM(deposit_volumes)/(1-erodible_porosity)
lava carrier:     deposit_tot = accumulated immobilized lava thickness
\`\`\`

Nel ramo lava, \`n_solid=0\` è valido e non implica deposito nullo; non leggere né validare \`erodible_porosity\` per costruire il deposito lavico. Non sommare \`lava_deposit\` a \`deposit_tot\` e non aggiornare entrambi come inventari indipendenti. Il rate \`dB_dt\` è il **rate della geometria**: aggiornamento dell'inventario e proiezione sulle quote nodali \`B_vertex\` devono avvenire una volta sola, evitando conversioni/cumulazioni doppie. Il dato canonico deve essere coerente anche se \`topo_change_flag=.FALSE.\` e il materiale si deposita.

La sorgente lava attuale fornisce \`lava_deposit_rate\` e \`dB_dt\`; l'integratore usa \`lava_deposit\` come candidato cell-local, controlla il suo valore e lo include nel rollback. La migrazione sostituisce tali usi con il candidato \`deposit_tot\` **solo dove il ramo è lava**; nel granulare le classi solide continuano a costituire l'inventario autorevole e \`deposit_tot\` è derivato. Valutare se sostituire il campo \`lava_deposit_rate\` del contratto \`mass_exchange_rates_type\` con una denominazione generica del trasferimento geometrico, oppure mantenerlo come rate identificato dalla legge lavica, senza trasformarlo in secondo inventario.

Per rollback si salvano tutti gli inventari autorevoli modificabili e la geometria; per restart lo spessore canonico deve essere ripristinato senza perdita, preservando la lettura esplicita del formato storico v7 o documentandone il rifiuto. I test confrontano separatamente casi granulari e lava, anche con \`n_solid=0\`, prima/dopo restart e con/senza \`topo_change_flag\`.

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

### 4.1 Configurazione per layer: stato del codice e contratto

Oggi \`parameters_2d\` e \`constitutive_parameters_2d\` contengono numerosi parametri globali, inclusi \`rheology_model\`, flag di fase, termica e proprietà costitutive. \`source_configuration_2d\` costruisce un catalogo unico di nomi attivi; \`source_registry_2d\` rifiuta una seconda sorgente reologica, la ripetizione dello stesso nome e gruppi di scambio massico in più gruppi; \`source_input_2d\` usa namelist e gruppo singoli. **Aggiungere il campo scope non basta a configurare due reologie indipendenti.**

Il modello two-layer usa una configurazione **immutabile per layer** risolta in inizializzazione, con: carrier e specie condivise compatibili, scelta reologica, coefficienti della sua legge, modalità e ordine sorgenti, termica, settling L1/L2, eventuali variabili fisiche inattive, opzioni di interfaccia e validazione di combinazioni consentite. La scelta di una reologia per-layer può riutilizzare gli stessi ID del catalogo, ma la regola di unicità va applicata a \`(source_identity, scope)\`; ogni gruppo può avere footprint e parametri diversi, con un unico ordinamento Lie globale.

Le interfacce cell-local ricevono la configurazione del proprio layer come dato **read-only** oppure usano un contesto esplicitamente posseduto dalla simulazione. È vietato selezionare L1/L2 cambiando temporaneamente parametri globali in un loop OpenMP, perché provocherebbe race e cambierebbe la semantica del single-layer. Non si duplicano moduli di reologia; si rendono passabili i parametri effettivamente diversi. La migrazione dei parametri viene pianificata per dipendenze reali, senza convertire indiscriminatamente tutte le costanti del modello.

Per `N_LAYERS=1` sono ammessi `RHEOLOGY_PARAMETERS` oppure `RHEOLOGY_L1_PARAMETERS`, ma non entrambi. La namelist `RHEOLOGY_PARAMETERS` legge un oggetto derivato `rheology` che viene copiato in `rheology_l1` prima della validazione e della costruzione dei contesti fisici. Per `N_LAYERS=2` sono ammesse solo le namelist distinte L1/L2. \`N_LAYERS=2\` dispone di input per L1, L2 e coppia (settling, elutriazione, drag, pressure coupling, lofting), con controllo di completezza dei dati, proprietà per specie e combinazioni fisiche consentite. La modalità two-layer riguarda una corrente densa e una nube diluita nello **stesso scenario gas/solidi**: il fatto che il codice supporti anche lava monolayer non autorizza combinazioni lava+granulare tra layer.

### Contratto input: tipi derivati nelle NAMELIST (decisione 11 novembre 2026)

**Scelta vincolante per l'implementazione TL04:** definire una sola volta tutti i coefficienti e i metadati reologici nel tipo derivato `rheology_config_type`; usare direttamente una variabile di tale tipo in ciascuna namelist, senza ripetere l'elenco dei componenti:

```fortran
TYPE(rheology_config_type) :: rheology, rheology_l1, rheology_l2
NAMELIST /rheology_parameters/ rheology
NAMELIST /rheology_l1_parameters/ rheology_l1
NAMELIST /rheology_l2_parameters/ rheology_l2
```

Esempio del **nuovo formato da implementare**, non ancora supportato dal parser:

```fortran
&RHEOLOGY_L1_PARAMETERS
  RHEOLOGY_L1%RHEOLOGY_MODEL = 'voellmy',
  RHEOLOGY_L1%MU = 0.3,
  RHEOLOGY_L1%XI = 300.0,
  RHEOLOGY_L1%RHEOLOGY_METHOD = 'NEWTON',
  RHEOLOGY_L1%RHEOLOGY_GROUP = 1,
/
&RHEOLOGY_L2_PARAMETERS
  RHEOLOGY_L2%RHEOLOGY_MODEL = 'froude',
  RHEOLOGY_L2%MU_0 = 0.1,
  RHEOLOGY_L2%MU_INF = 0.5,
  RHEOLOGY_L2%FR_0 = 1.0,
  RHEOLOGY_L2%RHEOLOGY_METHOD = 'ANALYTIC',
  RHEOLOGY_L2%RHEOLOGY_GROUP = 2,
/
```

I coefficienti sono puramente illustrativi. Con `N_LAYERS=1` il file può contenere **esattamente una** tra `RHEOLOGY_PARAMETERS` e `RHEOLOGY_L1_PARAMETERS`. Nel primo caso si legge `rheology` e si esegue `rheology_l1 = rheology`; nel secondo si legge direttamente `rheology_l1`. In entrambi i casi, dopo questa normalizzazione **l'unica configurazione usata dal resto del codice è `rheology_l1`**. Con `N_LAYERS=2` sono richieste `RHEOLOGY_L1_PARAMETERS` e `RHEOLOGY_L2_PARAMETERS`; `RHEOLOGY_PARAMETERS` è vietata. L'input single-layer legacy deve adottare anch'esso la sintassi del tipo derivato (`RHEOLOGY%MU`, `RHEOLOGY%XI`, ecc.): non è mantenuta la lettura degli scalari nudi, né una terza definizione dei campi. Le namelist con formato errato, duplicato o mescolato sono rifiutate con messaggi espliciti. La conversione degli esempi e dei test va fatta assieme al parser. La configurazione letta è indipendente per layer; viene validata una sola volta attraverso routine comuni e poi convertita in un contesto fisico per-layer, senza mutare i globali dei coefficienti nei loop OpenMP. Per le altre famiglie di parametri adottare tipi derivati nelle namelist quando si effettua la rispettiva migrazione, non prima.

**Stato:** decisione aggiornata il 12 novembre 2026; lettura basata su tipi derivati non ancora implementata. La conservazione da verificare è quella del comportamento numerico single-layer, non della sintassi dei file storici.

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

### 5.4 Specifica dei percorsi e della CFL accoppiata

La routine corrente \`pccu_2d:eval_hydrostatic_path\` usa lo stato esteso single-layer \`(h,Gamma,eta,G)\` con \`eta=B+h\`; le sue tracce sono preparate da \`reconstruction_2d\` e i percorsi sono assemblati da \`hyperbolic_2d\`. Nel two-layer si definiscono esplicitamente **tre superfici geometriche**:

\`\`\`text
B          = quota del substrato
eta_1      = B + h1            (interfaccia L1/L2)
eta_tot    = B + h1 + h2       (superficie superiore)
\`\`\`

La ricostruzione accoppiata usa la stessa geometria alle facce, mantiene \`h1,h2>=0\` e fornisce a ciascun layer le proprie tracce e le superfici fisiche corrette. Non è corretto utilizzare \`eta=B+h2\` come superficie del secondo layer quando L1 è presente. Se L1 è assente, l'interfaccia coincide con B. Occorre definire la regolarizzazione delle proprietà di L2 agli estremi asciutti, analogamente alla regolarizzazione di \`Gamma\` nelle facce single-layer.

La callback dell'accoppiamento riceve le grandezze estese e le tangenti, già fornite dal chiamante; il percorso (primitive, conservative o esteso) è specificato insieme alla coppia di interfacce. Il prodotto completo include gli effetti di \`dG/ds\`; usare la stessa quadratura per L1 e L2, senza sovrascrivere i contributi idrostatici propri o la topografia. L'identità esatta di somma è un test del path integrato, mentre la riduzione del **residuo idrostatico completo** al single-layer richiede verificare anche il posizionamento di \`G\` nei flussi e nelle sorgenti.

Le velocità caratteristiche correnti (\`equation_terms_2d:eval_local_speeds_direction\`) usano \`u±sqrt(G*gprime*h)\` **per una sola corrente**. Non costituiscono da sole un bound dimostrato per il sistema accoppiato. Prima della modalità di produzione, identificare il Jacobiano/sistema dei modi, le ipotesi di iperbolicità, un bound documentato e un fallback in caso di shear o perdita di iperbolicità. CFL, positività e dry/wet devono essere verificati con scenari in cui entrambe le correnti sono presenti; non assumere che \`min(dt_L1,dt_L2)\` sia sempre sufficiente.

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

### 6.3 Specifica delle chiusure fisiche del manoscritto

Il riuso dei moduli attuali identifica le **leggi di base disponibili**, non implica che tutte le chiusure previste per L1 siano già implementate:

- **Sedimentazione L2→L1:** velocità terminali disponibili in \`mass_exchange_deposition_2d\`; configurazione per-layer di hindered settling, limiti per classe, trasferimento di momento a \`u_Gamma\` coerente con interfaccia e opzioni \`CHI_XFER_L2\`, \`H_PLUG_XFER_L2\`, \`XFER_U2_AMP_CAP\`, salvo esclusione scientificamente documentata. Il limite alla velocità del donatore deve preservare il momento trasferendo l'eventuale differenza a L1.
- **Consolidazione L1→terreno:** il manoscritto prevede \`ALPHA_PACKING_L1\`, seppellimento di gas interstiziale, opzioni per settling con \`d32\` comune alle classi e inibizione per stress efficace/pressione interstiziale. Il codice esistente usa una legge di deposito con coefficiente di porosità e fattori di hindered settling fissati nel kernel; il comportamento del manoscritto richiede quindi un'estensione **fisica esplicita** di tali chiusure (senza duplicarne inutilmente la velocità terminale). L'inventario del substrato resta distinto dallo scambio L2→L1.
- **Elutriazione:** modalità \`SHEAR_POWER\` con stress d'interfaccia e di superficie, filtro Rouse e Shields, factor di availability \`Apack\`, opzione \`LINEAR_USTAR\`, cold-start quando L2 è sottile, limite di disponibilità per specie e co-lift di gas determinato dalla densità target.
- **Drag:** gate di spessore e formula vettoriale Chézy completa; l'eventuale rimozione/dissipazione dell'energia cinetica resta distinta dal trasferimento della variabile termica.
- **Lofting:** soglia di buoyancy per L2, esportazione di gas/specie/momento/energia e diagnostica per classe; opzioni di trasferimento di solidi verso L1 (tutti o ripartizione in funzione della velocità di salita del plume). In presenza di solo L2 buoyant non cancellare L1.
- **Parametro packing:** distinguere chiaramente \`ALPHA_PACKING_L1\`, la concentrazione di jamming dello hindered settling e il packing massimo di ammissibilità; sono coefficienti fisici differenti, non sinonimi.
- **Pressione interstiziale:** il seeding di pressione nei nuovi sedimenti mobili L1 e la sua successiva dissipazione richiedono termini specifici; non attivare la stessa equazione in L2.

Le formule, le unità, i limiti di validità e i default devono essere mappati puntualmente alle sezioni 2.5–2.7, 2.12 e 2.14 del manoscritto in un **catalogo di chiusure** prima di implementare TL08–TL11. Le opzioni non implementate non vengono attivate silenziosamente sotto il nome del modello completo.

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

### 8.1 Formati di input, NetCDF, restart e compatibilità

L'I/O corrente è prevalentemente single-layer e risiede in \`src/io/inpout_2d.f90\`:
- il NetCDF usa \`dimids(3)=(x,y,time)\`, variabili come \`h\`, \`eta\`, \`u\`, \`v\`, \`T\` e \`deposit__XX\` senza dimensione layer;
- il restart binario **v7** registra \`q\`, topografia nodale e inventari, con una registrazione finale separata di \`lava_deposit\`, e rifiuta versioni diverse; nell'header dimensionale non compare \`n_layers\`;
- \`N_LAYERS=2\` è attualmente respinto sia dall'inizializzazione del layout sia dall'I/O;
- sono presenti anche formati legacy \`.asc\` / \`.q_2d\` e diagnostiche/probe/runout, che richiedono uno schema di compatibilità definito.

**Contratto single-layer:** per \`N_LAYERS=1\` preservare il contenuto e i nomi dei campi NetCDF e dei file legacy, salvo una modifica deliberata e testata del deposito canonico. La forma interna rank-4 non deve automaticamente produrre una dimensione layer negli output monolayer. Per restart v7, scegliere in modo esplicito se mantenere un percorso writer/reader byte-compatibile o introdurre un reader di conversione v7→nuovo schema. Non affermare retrocompatibilità se non è supportata da test reali.

**Contratto two-layer:** identificare nell'header \`n_layers\`, layout, tipo di modello, flag/chiusure e versione; memorizzare i due stati fluidi e gli eventuali slot inattivi, substrato condiviso, topografia nodale, numeriche di continuazione, stato stocastico e diagnostiche. Per NetCDF definire la dimensione layer o suffissi L1/L2, le superfici \`eta_1\` e \`eta_tot\`, la densità, il Ri, il budget co-ignimbrite e \`deposit_tot\` **una sola volta** come campo del substrato. Documentare il trattamento di plotter, postprocessing e output storici.

I test I/O sono effettuati sia senza scambi sia con scambi e modifica del terreno; un restart alla stessa istantanea deve riprodurre una continuazione identica entro il contratto di riproducibilità dichiarato.


**Avanzamento TL01/TL02 (migrazione diretta):** `state_type%q`, `state_type%qp`, i trial di Lie e il residuo spaziale utilizzano ora il layout `(component,layer,x,y)`; i chiamanti single-layer in driver, Lie, inizializzazione, stochastic e I/O accedono esplicitamente a `layer=1`. Il reader/writer del restart v7 continua a trasferire la sola sezione `(:,1,:,:)`, per mantenere il formato single-layer. Non sono introdotti wrapper né copie compatibili rank-3; il contenitore preparatorio `layered_fluid_storage_2d` viene rimosso. Restano da eseguire compilazione, test e audit degli altri consumer e dei workspace globali: TL01/TL02 e G2 non sono dichiarati completati.

### 8.2 Registro di avanzamento, verifiche e commit (11 novembre 2026)

Questa tabella separa le modifiche presenti nel branch dalle verifiche effettivamente superate. Lo stato non equivale a certificazione.

| Blocco | Stato | Evidenza / attività residua |
|---|---|---|
| TL00 — baseline e audit | **Parziale** | Audit statico in `TWO_LAYER_CODE_AUDIT.md`; base `master` identificata come `5a377839`. Tag immutabile, baseline riproducibile e protezione branch ancora da certificare. |
| TL00a — deposito canonico | **Implementato; verificato dall'utente sul desktop** | La funzionalità del deposito unico per lava è stata confermata dall'utente. Restano non documentati i risultati specifici della matrice completa G1 (restart v7, NetCDF, topografia on/off, casi granulari). |
| TL01 — stato fluido rank-4 | **In lavorazione** | `state%q/qp` a rank-4 e indice layer esplicito; build optimized/debug **superata** su `155a0b6` (Action #42990883485). Test di contiguità, layout e regressioni ancora da verificare. |
| TL02 — migrazione consumer | **In lavorazione** | Driver, Lie, inizializzazione, stochastic, input/output e test toccati; primo build fallito per `SIZE(state%qp(:,1,:,:),4)` in `inpout_2d.f90`; corretto in `155a0b6`, con build optimized/debug **superata** (Action #42990883485). Audit consumer, NetCDF, restart e OpenMP ancora da completare. |
| TL03 — regressione single-layer | **In lavorazione; regressioni CI e esempi desktop superati** | [Action #43024081071](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43024081071) **success** su `a7227c6`; l'utente conferma che gli esempi sul desktop sono OK dopo la migrazione. Restano da documentare esatti casi e condizioni, restart binario, NetCDF, OpenMP e confronto col baseline; G2/G3 aperti. |
| TL04 — contesti fisici e stato inerte | **In lavorazione** | Input deciso: un solo tipo derivato per tutte le namelist; `RHEOLOGY_PARAMETERS` (oggetto `rheology`) ammessa solo in single-layer e copiata in `rheology_l1`, oppure `RHEOLOGY_L1_PARAMETERS` direttamente; per two-layer L1/L2 obbligatorie. Parser/input ancora da migrare. Primo blocco: tipo `layer_physics_type` con fase, selezione reologia validata e identificativo canonico, istanze possedute dal `simulation_context_type` e test di configurazione indipendente. I coefficienti fisici, scope del registry, variabili inattive, IC/BC/output e attivazione dei due layer restano da implementare; `N_LAYERS=2` continua a essere rifiutato. Build e test del nuovo blocco non ancora verificati. |
| TL05–TL13 — accoppiamento e validazione | **Non iniziati** | Restano subordinati ai gate previsti e al completamento di TL04. |

**Quarto tentativo CI (11 novembre 2026):** [Cell physics and numerical regressions #43024081071](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43024081071) sul commit `a7227c6`: **success**. La correzione del ciclo begin/commit/rollback del test è stata verificata dalla CI; si chiude il sottogate delle regressioni di fisica e numerica. Non equivale alla chiusura di TL03, che richiede ancora esempi integrati e I/O/restart.

**Terzo tentativo CI (11 novembre 2026):** [Cell physics and numerical regressions #43020328048](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43020328048) al commit `1f32de0` raggiunge la chiamata `check_lie_exchange_transaction` ma fallisce nel guard `geometry_2d:apply_local_bed_rate`: il test chiamava direttamente l'esecutore Lie senza avviare la transazione geometrica richiesta. Corretto il test con `begin_geometry_transaction(sparse=.TRUE.)` prima dell'esecutore, `commit_geometry_transaction` in caso di successo e `rollback_geometry_transaction` in caso di fallimento. Non è stata modificata la logica di produzione. **Nuovo rerun necessario.**

**Secondo tentativo CI (11 novembre 2026):** la [regressione #43019726403](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43019726403) sul commit `9f1de1e` supera la compilazione di `TEST_MASS_EXCHANGE`, ma il controllo runtime Fortran fallisce in `allocate_grid` (linea 491): `state%deposit_tot` era allocato sia nella precedente `ALLOCATE` (linea 489) sia nella successiva. Il duplicato è stato rimosso mantenendo un solo `ALLOCATE` e la corrispondente singola `DEALLOCATE` in `release_grid`. **La verifica CI della correzione è ancora aperta.**

**Risultato CI successivo al build:** [Cell physics and numerical regressions #43018138981](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43018138981), commit `78bd1bd`: **failed** durante la compilazione di `TEST_MASS_EXCHANGE/test_projection_stress.f90`; `state%deposit_tot` è elencato due volte nella `DEALLOCATE` di `release_grid`. La correzione elimina soltanto il duplicato, senza cambiare fisica o layout. Il rerun del test è necessario per avanzare G2/TL03.

**Decisione sul parser TL04 (aggiornata 12 novembre 2026):** un solo tipo derivato per i coefficienti. In single-layer il nome legacy `RHEOLOGY_PARAMETERS` è ammesso, ma il suo oggetto `rheology` è copiato in `rheology_l1` subito dopo la lettura; in alternativa si accetta la namelist L1. Le vecchie assegnazioni dirette a scalari non sono conservate. In two-layer solo L1 e L2. Validazione unificata dopo la normalizzazione; modifiche al parser e agli input non ancora implementate.

**Avvio TL04 (11 novembre 2026):** introdotta una configurazione per-layer separata, con campi privati e lettura tramite query, validazione della fase e della reologia attraverso il catalogo esistente, senza modifica dei globali durante OpenMP. Il contesto di simulazione possiede la configurazione del layer singolo e la inizializza prima dei workspace. Non si considera ancora completata la parte di *parametri fisici* per layer: vanno trasferiti esplicitamente i coefficienti e aggiornati i kernel prima di abilitare due layer, senza cambiare globali tra chiamate. Aggiunto il test `TEST_LAYER_PHYSICS` al gruppo `state` della Action di regressione, da eseguire. Il codice è stato committato ma **non** compilato localmente, perché il repository remoto non è accessibile dalla macchina di compilazione di questa sessione.

**Registro commit rilevanti (ordine cronologico):**

| Commit | Attività | Esito |
|---|---|---|
| `208db5c` | Estensione test del deposito canonico e roadmap TL00a | Committato; test non certificato. |
| `b110abf` | Registrazione verifica desktop dell'utente sul deposito unico lavico | Committato. |
| `f6d75e0`, `dbb0af9` | Preparazione contenitore rank-4 indipendente e allineamento precisione | Preparazione, poi sostituita dalla migrazione diretta. |
| `5653681`, `c448432`, `6c4ea8e`, `1452259` | Migrazione diretta dello stato e consumer, correzioni e test del layout | Committati, prima verifica di compilazione ancora incompleta. |
| `36dba58`, `813622f`, `c4171db` | Migrazione dei workspace Lie, correzione allocatable e allineamento roadmap | Committati, build da verificare. |
| `155a0b6` | Correzione del rango intrinseco `SIZE` nell'output NetCDF | **Build optimized/debug superata**, [Action #42990883485](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/42990883485). |
| `9f1de1e` | Eliminazione doppia `DEALLOCATE` nel test `TEST_MASS_EXCHANGE` | Compilazione test superata, ma regressione fallita a runtime nella `ALLOCATE` duplicata. |
| `1f32de0` | Eliminazione doppia `ALLOCATE(state%deposit_tot)` in `allocate_grid` | Correzione confermata: regression CI raggiunge il test Lie, fallendo sul contratto geometrico. |
| `a7227c6` | Test `check_lie_exchange_transaction`: gestione esplicita del ciclo della transazione geometrica | **Regressione superata:** [Action #43024081071](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/43024081071). |
| `c3a89bb` | Registro della Action di regressione riuscita | Documentazione aggiornata. |
| `fb965ba` | TL04: tipo per-layer con contesto di simulazione, validazione di reologia/fase e test | Implementato; build e test da verificare. |
| **Commit successivo a `fb965ba`** | Integrazione di `TEST_LAYER_PHYSICS` nella CI e pulizia documentazione TL04 | Test da eseguire con la prossima Action. |

**Dettaglio build:** il compilatore ha segnalato 12 istruzioni `ALLOCATE` in `src/io/inpout_2d.f90`, alle linee 7686–7696 e 7824, che usavano `SIZE(state%qp(:,1,:,:),3/4)`. Poiché la sezione è rank-3, il parametro `DIM=4` è invalido. Gli accessi sono ora `SIZE(state%qp,3)` e `SIZE(state%qp,4)`, riferiti all'array originale rank-4, senza wrapper e senza copie. La build optimized/debug è stata eseguita su GitHub e risulta completata con successo sul commit `155a0b6` il giorno 11 novembre 2026 ([Action #42990883485](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/42990883485)); i test di regressione e le simulazioni non risultano ancora verificati in questo registro.

## 9. Roadmap operativa

I blocchi TL00–TL13 descrivono il percorso completo. Su `feature/two-layer` sono presenti le modifiche al catalogo lava, l'input `LAVA_TRANSPORT_PARAMETERS` e l'implementazione di TL00a; la certificazione dei test e i gate di milestone restano aperti. Ogni blocco è un intervento separato e verificabile; i primi gate proteggono il single-layer funzionante.

| ID | Obiettivo tecnico | Criteri di accettazione |
|---|---|---|
| **TL00 — baseline e audit** | Identificare commit single-layer validato, fissare tag annotato; definire branch di sviluppo separato e verificarne la protezione; inventario di rank-3, globali di fisica e configurazione, variabili di deposito, NetCDF/legacy/restart v7, geometria, domain worklist, sorgenti e test. Catalogare formule/opzioni del manoscritto, test baseline e costi. | Riferimento immutabile, mappa dipendenze e test riproducibili. La protezione di `master` e lo stato del branch sono verificati, non presunti. |
| **TL00a — deposito canonico** | Rendere `deposit_tot` lo spessore unico per granulari **oppure** lava, con derivazione granulare e incremento lavico esplicitamente distinti; eliminare `lava_deposit` come inventario autonomo dopo migrazione degli usi; allineare rate fisici, rollback, topografia, I/O e restart v7. | Deposito/erosione e geometria invariati fisicamente; EXAMPLE_ETNA_LAVA e casi granulari confrontati anche con `n_solid=0`, topografia on/off e restart v7; nessun azzeramento del deposito lavico, doppio conteggio o uso della porosità granulare nel ramo lava. |

**Avanzamento TL00a:** estesa la regressione statica `TEST_WORKSPACE_LIFECYCLE` per verificare snapshot/rollback del deposito lavico autorevole con `n_solid=0`, refresh dopo rollback e azzeramento del totale in granulare con `n_solid=0`. **Verifica dell'utente (desktop, 11 novembre 2026):** il deposito unico, utilizzato anche per la lava, risulta funzionante nelle prove effettuate localmente. I comandi, gli esiti dei singoli test e le evidenze di restart v7, NetCDF ed esempi completi non sono stati forniti; non sono quindi certificati in questa roadmap. G1 resta aperto fino alla raccolta di tali evidenze. La migrazione rank-4 TL01/TL02 non è iniziata.
| **TL01 — layout rank-4, un layer** | Introdurre il secondo indice layer nei campi fluidi globali e nei metadati, lasciando ammesso soltanto `N_LAYERS=1`; preservare i kernel cell-local. | Build e test di stato/convertitori, nessuna modifica intenzionale delle formule. |
| **TL02 — migrazione consumer globali** | Adeguare Lie, workspace, spatial, ricostruzione, BC/IC, driver, output, NetCDF, restart, diagnosi, snapshot e test. | EXAMPLES, restart round-trip, bounds, OpenMP e formati documentati. |
| **TL03 — gate single-layer** | Confrontare evoluzioni, campi, budget, performance, memoria, rami wet/dry, buoyancy, rollback e risultati finali con il baseline taggato. | Nessuna regressione inspiegata; bitwise quando realisticamente ottenibile. Nessuna abilitazione two-layer prima della chiusura del gate. |
| **TL04 — stato two-layer inerte** | Slot L1/L2, variabili inattive, **contesti fisici read-only per layer** (nessun cambio di globali dentro OpenMP), proprietà specie comuni, input/validazione delle due reologie, registry con univocità per `(identità,scope)`, dry state, IC/BC e output, senza accoppiamento. | Round-trip e validità per wet/wet, wet/dry, dry/wet, dry/dry. |
| **TL05 — HP-PCCU per layer** | Ricostruzione e flussi di ciascun layer, dominio/halo coerenti, timestep globale preliminare. | Ciascuno slot riproduce il comportamento del singolo flusso in assenza di accoppiamento; `eta_1=B+h1`, `eta_tot=B+h1+h2` e trattamenti dry/vuoto verificati. |
| **TL06 — pressione interlayer** | Utilizzare **l'attuale `nonconservative_2d` invariato** con callback della coppia idrostatica; assemblare percorsi di faccia e interni, topografia, metriche, clamp e CFL accoppiato. | Identità dei due path, riduzione del **residuo idrostatico completo** a single-layer con topografia, equilibrio su fondo piano/inclinato, residuo curvo caratterizzato, clamp, casi wet/dry; bound CFL accoppiato giustificato e perdita di iperbolicità gestita. |
| **TL07 — Lie con scope e trasferimento di coppia** | Registro L1/L2/INTERLAYER, riuso configurabile degli stessi termini in layer diversi, `Gamma` unico, limiter condiviso, due candidati simultanei, validazione e rollback atomici, rispetto del gruppo 0, ordine e footprint per-layer. | Test manufactured di conservazione, positività dei donatori, rollback e OpenMP. |
| **TL08 — settling L2 → L1** | Riutilizzare velocità terminale e chiusure di settling; trasferire specie, massa, momento ed energia con velocità d'interfaccia `u_Gamma`, blend per thin-donor e cap di amplificazione conservativo; contabilizzare separatamente l'approssimazione gas packing. | Conservazione dei trasferimenti interni, eccezioni diagnosticate, **nessun deposito al terreno**. |
| **TL09 — elutriazione L1 → L2** | Implementare SHEAR_POWER / LINEAR_USTAR, stress da slip e superficie, Rouse/Shields, cold-start con L2 sottile, availability factor, co-lift verso densità target e limiter di coppia. | Bilanci conservativi, limite d'inventario e trasporto termico/momento coerenti. |
| **TL10 — drag interfaciale** | Forza quadratica uguale e opposta sui momenti; verificare rigidità esplicita. | Momento totale conservato, forza nulla nei casi limite, stabilità. |
| **TL11 — sorgenti esterne e chiusure** | Entrainment atmosfera→L2, consolidazione/deposizione L1→terreno con packing/gas burial/d32/stress efficace secondo opzioni, eventuale erosione, pressione interstiziale L1, reologie differenziate, lofting L2 e budget plume, deposito canonico. | Nessun doppio entrainment o modifica del deposito da soli scambi interlayer; budget e output diagnostici coerenti. |
| **TL12 — validazione scientifica** | Benchmark del manoscritto, resting states, topografia curva, coppie di processi simultanei, wet/dry, CFL, restart, OpenMP e confronto con single-layer. | Matrice test e differenze scientifiche documentate, baseline 1-layer ancora superato. |
| **TL13 — finalizzazione** | Documentazione delle chiusure, loro default, configurazione per layer, NetCDF/restart versionato e retrocompatibilità v7, limiti scientifici; aggiornamento tracker generale; verifica integrata e abilitazione `N_LAYERS=2`. | Review finale, test superati e decisione esplicita di abilitazione. |

**Dipendenze principali:** TL00 → TL00a → TL01 → TL02 → TL03 (gate) → TL04 → TL05 → TL06; TL07 precede TL08–TL10; TL11 completa l'interazione con ambiente e substrato; TL12 e TL13 concludono. Eventuali test fisici isolati possono essere preparati prima senza saltare i gate di integrazione.

## 10. Milestone di integrazione in `master`

Le milestone sono **punti di rilascio incrementale** lungo la roadmap tecnica TL00–TL13. Il completamento di un singolo blocco TL non autorizza automaticamente il merge: ciascuna milestone raggruppa blocchi che insieme producono uno stato del codice **coerente, funzionante, verificato e mantenibile in `master`**.

**Principio fondamentale:** `master` deve continuare a compilare, eseguire gli esempi single-layer e mantenere la propria interfaccia pubblica, a prescindere dal fatto che il two-layer sia ancora incompleto. Fino alla milestone di rilascio **M6**, `N_LAYERS=2` rimane disabilitato nelle build di produzione e viene eventualmente attivato **solo negli harness diagnostici espliciti** delle funzionalità già implementate. Non esporre una modalità two-layer che restituisce risultati fisicamente incompleti senza avvertimento.

### 10.1 Milestone, confini e risultati

| Milestone | Blocchi coperti | Risultato integrabile in `master` | Condizione specifica di merge |
|---|---|---|---|
| **M0 — Baseline e specifiche** | TL00 (documentazione e audit) | Baseline single-layer identificato, tag e documentazione architetturale verificati; branch e regole di integrazione definiti. | Il riferimento baseline è realmente riproducibile; l'aggiornamento documentale non descrive come implementate funzioni ancora assenti. **Nessuna modifica fisica.** |
| **M1 — Deposito canonico single-layer** | TL00a | Un solo `deposit_tot` per lo spessore del deposito granulare **oppure** lavico, con leggi indipendenti e rollback/I/O coerenti. | Esempi granulari e lava, test `n_solid=0`, erosione, topografia on/off e restart; nessuna regressione fisica o ambiguità sul formato storico. |
| **M2 — Layout layer-aware certificato, ancora single-layer** | TL01–TL03 | Stato e consumer fluidi con indice layer; la configurazione pubblica continua a supportare soltanto `N_LAYERS=1`. | Gate TL03 chiuso: build, regressioni, full EXAMPLES, NetCDF, restart, OpenMP, memoria/performance; confronto con baseline M1. `N_LAYERS=2` ancora rifiutato. |
| **M3 — Infrastruttura two-layer isolata** | TL04–TL05 | Stato L1/L2, configurazioni fisiche per layer, ricostruzioni e trasporti indipendenti verificabili tramite test dedicati; nessuna simulazione coupled pubblica. | Test a entrambi gli slot, dry/wet in ogni combinazione, assenza di race su parametri, proprietà di isolamento e compatibilità single-layer; l'accoppiamento fisico resta disattivato. |
| **M4 — Operatore spaziale accoppiato verificato** | TL06 | Idrostatica interlayer HP-PCCU, geometria delle due superfici, percorsi non conservativi e criterio CFL accoppiato, tutti testati separatamente. | Identità di coppia e riduzione single-layer dell'operatore completo, well-balancing nei casi previsti, clamp, wet/dry, bound CFL e iperbolicità documentati; feature ancora non selezionabile in produzione. |
| **M5 — Accoppiamento locale e chiusure fisiche** | TL07–TL11 | Gruppi Lie per layer e coppia, limiter/rollback atomici, settling, elutriazione, drag, entrainment, deposito basale e lofting, con diagnostiche e configurazioni definite. | Test di conservazione e casi limite per ogni chiusura, budget esterni espliciti, restart/output coupled, nessun doppio scambio o perdita di inventario, regressioni single-layer; codice two-layer ancora protetto dal gate di produzione. |
| **M6 — Rilascio del modello two-layer** | TL12–TL13 | Modalità `N_LAYERS=2` abilitata all'utente con input, output, restart, documentazione delle chiusure e benchmark scientifici. | Benchmark del manoscritto e matrice completa di validazione superati, limiti scientifici dichiarati, performance e stabilità verificate, approvazione esplicita al rilascio. |

**Regola sui blocchi non ancora completamente attivati:** è ammissibile integrare un componente interno non pubblico soltanto se il codice di produzione non lo raggiunge per `N_LAYERS=1`, i test dedicati ne verificano l'implementazione e non esistono fallback silenziosi o configurazioni apparentemente supportate. Non utilizzare `N_LAYERS=2` come test di produzione prima di M6.

### 10.2 Gate comuni obbligatori per ogni merge

Ogni milestone deve soddisfare **tutti** i requisiti comuni, oltre a quelli specifici indicati nella tabella:

1. **Ambito chiuso:** tutti i blocchi TL richiesti dalla milestone sono implementati; nessun cambiamento non correlato viene incluso per comodità.
2. **Build e test realmente eseguiti:** compilazione pulita con configurazioni di riferimento; test di unità e integrazione interessati; suite regressiva single-layer completa in base al piano di accettazione, con risultati registrati. Una Action non avviata non vale come test superato. Gli esempi di lunga durata possono richiedere esecuzione manuale documentata.
3. **Fisica e numerica:** nessun cambiamento intenzionale non documentato di equazioni, schemi HP-PCCU o ordine Lie; bilanci, positività, wet/dry e rollback verificati nei processi interessati. Nessuna nuova modalità resa disponibile se non supera i propri test scientifici.
4. **Compatibilità:** configurazioni `N_LAYERS=1`, input esistenti, output attesi, diagnostiche e restart gestiti secondo il contratto della milestone; ogni incompatibilità deliberata è versionata, documentata e testata.
5. **Prestazioni e concorrenza:** assenza di regressioni non spiegate in memoria, temporanei e OpenMP; ogni nuovo percorso usa parametri read-only e scratch appropriati.
6. **Tracciabilità:** commit di partenza, file modificati, risultati/test, differenze dal baseline, limiti conosciuti e stato dei blocchi registrati nel tracker generale; roadmap specialistica aggiornata come **stato corrente consolidato**, non come cronologia di patch.
7. **Revisione prima del merge:** Pull Request da `feature/two-layer` a `master` esaminata, conflitti risolti, gate verificati e **approvazione esplicita**. Il passaggio non è automatico al semplice completamento dei test.

Se un gate fallisce, la milestone **non** è raggiunta: il lavoro rimane sul branch di sviluppo. Se dopo la fusione emerge un problema, si interviene con una correzione isolata o con rollback del commit di integrazione; non si lascia `master` in uno stato non funzionante.

### 10.3 Ciclo di sviluppo e sincronizzazione dei branch

1. Il lavoro in corso procede su **`feature/two-layer`**; i blocchi TL possono essere committati e verificati singolarmente sul branch prima della milestone.
2. Al raggiungimento di M0–M6, preparare una **PR circoscritta** verso `master`, contenente tutti e soltanto i cambiamenti non ancora integrati. I test e le evidenze di accettazione sono allegati alla PR.
3. Dopo il merge, considerare `master` come nuovo riferimento stabile; registrare il commit di integrazione e, quando utile, un tag di milestone (`two-layer-M1`, `two-layer-M2`, ...). I tag si creano **solo su commit realmente validati e integrati**, non preventivamente.
4. Continuare su `feature/two-layer` **sincronizzandolo con `master` dopo ogni merge**, senza riscrivere la storia pubblicata, cancellare commit o reintrodurre modifiche già integrate. Se la strategia di merge scelta produce una storia divergente (per esempio squash merge), allineare consapevolmente i rami prima di riprendere lo sviluppo: non presumere che un successivo PR diff sia automaticamente pulito.
5. Mantenere `master` protetto da push diretti e force-push, secondo le regole effettivamente disponibili e verificate nel repository. Le milestone non richiedono nuovi branch permanenti: eventuali branch temporanei sono strumenti di revisione, non workstream paralleli.
6. I lavori che non appartengono al modello two-layer, inclusi possibili nuovi metodi del secondo ordine per le sorgenti, **non** entrano in questa roadmap né nelle sue milestone.

### 10.4 Stato delle milestone

| Milestone | Stato | Merge in `master` | Evidenza |
|---|---|---|---|
| M0 | Da certificare | Non eseguito | Roadmap e audit preparati; baseline e gate ancora da verificare. |
| M1 | In verifica | Non eseguito | Codice TL00a presente sul branch; manca la certificazione dei test e del gate M0. |
| M2 | Non iniziata | Non eseguito | TL01–TL03 aperti. |
| M3 | Non iniziata | Non eseguito | TL04–TL05 aperti. |
| M4 | Non iniziata | Non eseguito | TL06 aperto. |
| M5 | Non iniziata | Non eseguito | TL07–TL11 aperti. |
| M6 | Non iniziata | Non eseguito | TL12–TL13 aperti. |

La tabella indica lo stato verificato della pianificazione, non un presunto avanzamento del codice. Aggiornare milestone, esiti e SHA **soltanto dopo** aver eseguito e registrato le verifiche corrispondenti.


## 11. Matrice di validazione

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

## 12. Vincoli scientifici e questioni da chiudere nei rispettivi blocchi

1. Definire esattamente le variabili L2 inattive e le loro semantiche per input, trasporto e restart.
2. Specificare il percorso esteso e la rappresentazione delle grandezze termodinamiche nella quadratura accoppiata, senza modificare il motore `nonconservative_2d`.
3. Verificare velocità caratteristiche e iperbolicità del sistema two-layer, inclusa l'influenza delle interfacce e del clamp.
4. Distinguere identità dell'operatore continuo, precisione di quadratura e well-balancing del residuo discretizzato su topografia curva.
5. Trattare esplicitamente le eventuali chiusure non conservative del gas di packing.
6. Definire condizioni al contorno e topografia condivisa nelle combinazioni wet/dry e negli stati con layer superiore presente senza L1.
7. Stabilire versioni e retrocompatibilità di restart/NetCDF e la semantica di `deposit_tot` dopo la razionalizzazione del deposito.
8. Verificare l'ordine dei gruppi Lie (in particolare gruppo 0) e il limite di stabilità per eventuali termini interlayer rigidi.
9. Confermare le opzioni scientifiche specifiche del manoscritto per lofting, co-lift, porosità e deposizione prima della validazione finale.

## 13. Protocollo operativo

- Ogni blocco è verificato sul branch di sviluppo; le PR verso `master` avvengono **alle milestone M0–M6** con gate comuni e specifici superati, test/esiti documentati e riferimento al baseline.
- `master` è il riferimento single-layer: lo sviluppo two-layer avviene su un branch dedicato; la protezione effettiva del branch e il tag del baseline vengono verificati in TL00.
- Nessun cambiamento di formula fisica è introdotto soltanto per adattare il layout. Gli aggiustamenti scientifici richiedono giustificazione e validazione autonoma.
- Le ottimizzazioni non pertinenti al two-layer e il precedente workstream di ulteriore refactoring dello stato restano fuori ambito.
- La compilazione da sola non chiude un gate scientifico; i benchmark e il confronto degli output hanno valore separato.
- Il progresso di TL00–TL13 e delle milestone M0–M6 viene riportato sinteticamente nella roadmap generale, mantenendo questo file come specifica tecnica consolidata.
- Le evidenze e i rischi emersi dall'ispezione statica del sorgente sono catalogati nel [rapporto di audit del codice e della roadmap](TWO_LAYER_CODE_AUDIT.md). L'audit non equivale a una compilazione o all'esecuzione dei benchmark.

## 14. Riferimenti

- Manoscritto scientifico: *IMEX_SfloW2D_v3: A two-layer shallow-water model for pyroclastic density currents with pore-pressure transport, polydisperse settling and shear-driven elutriation*, in particolare sezioni 2.1–2.2, 2.5–2.7, 3–4.
- [Roadmap generale del progetto](IMEX_SfloW2D_Development_Roadmap.md).
- Sorgenti esistenti: `src/state/`, `src/physics/mass_exchange/`, `src/spatial/nonconservative_2d.f90`, `src/spatial/pccu_2d.f90`, `src/spatial/hyperbolic_2d.f90`, `src/numerics/lie_source_executor_2d.f90`.
