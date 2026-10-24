# D17 — Full source-roadmap section (preserved verbatim)

**Source:** [Integrated development roadmap](./IMEX_SfloW2D_Development_Roadmap_WITH_LEGACY_LABELS_2026-10-15%20(1).md), section D17.  
**Preservation:** The following D17 section is reproduced **verbatim, without abridgment or technical edits**, from the branch document as of 2026-10-24.  
**Execution status and commits:** [D17_SOURCE_REFACTOR_HANDOFF.md](./D17_SOURCE_REFACTOR_HANDOFF.md).

> **IMPORTANT:** This is the full *historical* D17 text retained to prevent loss of details. Several points (notably numeric source IDs, catalog metadata and repeated independent Newton solves) were superseded by later agreed design choices. See the binding architecture section in the handoff. Do **not** silently mix the two versions: a full reconciled D17 specification must explicitly revise the historical passages, while retaining all unaffected details, examples, tests, defaults, acceptance gates and implementation steps.

---

### D17 — Organizzazione unificata dei termini sorgente, gruppi e separazione tra fisica e numerica

Stato ufficiale: IN CORSO — i metadati del registro, il bridge analitico,
l’interfaccia esplicita D17.9 e il routing filtrato dei rate semi-impliciti e
impliciti D17.10, incluso lo sweep ordinato dei solve Newton indipendenti nello
stesso gruppo, sono implementati e verificati localmente. Restano da completare
i requisiti successivi; D17 nel complesso non è completo. D17.9 è stato
committato in `3fc18813`; il lavoro D17.10 è ancora locale.

Baseline della fase A: `1310511e` (23 ottobre 2026). Specifica di compatibilità:
[D17_SOURCE_CONTRACT_2026_10_23.md](acceptance/D17_SOURCE_CONTRACT_2026_10_23.md);
inventario con 30 contributi fisici e 4 assemblaggi derivati:
[D17_SOURCE_TERMS_MANIFEST_2026_10_23.json](acceptance/D17_SOURCE_TERMS_MANIFEST_2026_10_23.json).
La fase A è una baseline documentale storica. Il catalogo e il registro
definiscono i 34 ID Fortran, i descrittori dei gruppi/solve e i registri
validati posseduti dalla simulazione, con adattatore dei gruppi IMEX legacy.
Rapporto:
[D17_B_SOURCE_REGISTRY_2026_10_23.md](acceptance/D17_B_SOURCE_REGISTRY_2026_10_23.md),
[evidenza machine-readable](acceptance/results_d17_b_source_registry_2026_10_23.json).
Il dispatch dei solve analitici è documentato in
[D17_C_ANALYTIC_BRIDGE.md](acceptance/D17_C_ANALYTIC_BRIDGE.md); il routing
esplicito stage, incluso l’hook distinto al bordo, è documentato in
[D17_9_EXPLICIT_STAGE_ROUTING.md](acceptance/D17_9_EXPLICIT_STAGE_ROUTING.md).
Il routing semi-implicito/implicito per gruppo e solver è documentato in
[D17_10_GROUP_FILTERED_IMPLICIT_RATES.md](acceptance/D17_10_GROUP_FILTERED_IMPLICIT_RATES.md).

Per il registro passano, in strict e optimized, 864 configurazioni legacy,
14.857 asserzioni e 28 casi invalidi rifiutati. I test D17-C e D17.9 verificano
localmente i percorsi analitici ed espliciti e le regressioni IMEX/spaziali
indicate nei rispettivi report. I risultati di parser e integratore riportati
nel rapporto B sono evidenza storica di quella fase, non una validazione
completa del D17 attuale. La normalizzazione conserva i blocchi Newton
congiunti e ordina i gruppi analitici legacy velocità/temperatura/pressione.
Validazione e copie complete del registro restano fuori dai loop di cella e
Newton.

Nessuna nuova accettazione Intel ifx/Linux/CI è attestata. Lo stato per ciascun
requisito numerato è:

| Requisito | Stato al 23 ottobre 2026 | Evidenza o lavoro ancora necessario |
|---|---|---|
| D17.1 | **Parziale** | Bridge analitico e dispatcher esplicito collegati; routing comune semi-implicito/implicito e mass-exchange da completare. |
| D17.2–D17.5 | **Implementati localmente** | Catalogo, ID, descrittori e registri `stage`/`step`; coperti dai test del registro. |
| D17.6 | **Implementato per i gruppi legacy supportati; builder nativo locale** | Validazione del solve unico e normalizzazione legacy; il builder nativo controlla copertura dei termini attivi e coerenza dei target (D17.16). |
| D17.7 | **Implementato localmente per i gruppi stage supportati** | Sweep ordinato dei solve indipendenti; blocco Newton legacy congiunto preservato come solve unico. L’equivalenza end-to-end resta aperta nei requisiti pertinenti. |
| D17.8 | **Parziale** | Le quattro interfacce ricevono metadati di gruppo/solve nei percorsi supportati; restano i percorsi numerici e di cadenza non ancora integrati. |
| D17.9 | **Implementato e committato** | Rate espliciti e hook al bordo selezionati dal gruppo; test strict/optimized e build GNU superati (`3fc18813`). |
| D17.10 | **Implementato localmente** | Rate semi-impliciti/impliciti filtrati per gruppo e solver, percorsi impliciti REAL/COMPLEX e sweep runtime ordinato dei solve Newton multipli; report e regressioni locali in D17.10. |
| D17.11 | **Implementato localmente per i kernel REAL/COMPLEX duplicati individuati** | La correlazione di viscosità lahar condivide ora il kernel COMPLEX con il wrapper REAL; valore e derivata sono testati nei due rami. Per il perimetro verificato non serve un file `.inc`; report D17.11. |
| D17.12 | **Implementato localmente** | Bridge analitico validato; vedere il report D17-C. |
| D17.13–D17.16 | **Parziali** | D17.13 resta parziale: equivalenza per tutte le attivazioni non completata. D17.14 instrada il batch mass-exchange atomico dal registro `step` una sola volta dopo RK. D17.15 verifica l’ordine delle cinque fasi production source-window/OU/RK/step-source/clock e un aggiornamento OU attivo con seed fissato; il test nativo end-to-end verifica una sorgente circolare attiva a metà intervallo e il bilancio integrato nel passo production. Test equation-term separati coprono impulso parziale, rampa lineare, wrap, inattività, finestra fissurale e rampa cosinusoidale radiale. Selezione di `dt`, output end-to-end, trasporto stocastico e tutte le combinazioni di sorgente restano aperti. D17.16 implementa localmente `SOURCE_GROUPS`, la distinzione esplicito/default legacy, i controlli di conflitto e la copertura completa dei termini attivi; il flusso completo `read_param` è coperto per una configurazione nativa rappresentativa, mentre l’equivalenza per tutte le combinazioni fisiche resta da validare. Evidenze D17.13–D17.16. |
| D17.17–D17.18 | **Parziali** | Dispatch analitico, esplicito e rate Newton filtrati sono dietro `equation_terms_2d`; il batch post-RK `step` è separato da `time_integration_2d` in `mass_exchange_step_2d`, mantenendo validazione del registro, atomicità e posizione temporale. Il test statico copre tutti i sorgenti `numerics/` e `spatial/`, proibisce import diretti dei moduli fisici e controlla l’allowlist delle API equazionali usate dallo spazio; i moduli `*_laws_2d` non dipendono dall’orchestrazione numerica/spaziale. Restano audit del call graph/runtime, integrazione term-selective e matrice D17.20. |
| D17.19 | **In corso** | Stato puntuale dei 24 elementi nella sequenza implementativa sotto. |
| D17.20 | **Parziale** | Le regressioni locali strict/optimized verificano il registro/parser, i solve analitici e Newton, i percorsi REAL/COMPLEX, gli hook spaziali, il batch post-RK e ora una finestra source attiva nel timestep production con budget di massa. Restano l’equivalenza su tutte le combinazioni fisiche, copertura end-to-end delle altre finestre/sorgenti, selezione di `dt`/output del driver, approvazione delle baseline e gate cross-platform/CI. |
| D17.21 | **Non soddisfatto** | Criteri finali non raggiunti; D17 non può essere dichiarato completato. |

#### D17.0 — Logica generale e obiettivo

Il D17 introduce una struttura comune per tutti i termini sorgente delle equazioni e separa:

1. leggi fisiche;
2. calcolo dei rate;
3. assemblaggio dei contributi nelle equazioni;
4. integrazione RK-IMEX;
5. solve impliciti;
6. solve analitici.

La struttura generale deve essere:

```text
time_integration_2d
        |
        v
nonlinear_solver_2d
        |
        v
equation_terms_2d
        |
        v
moduli fisici e leggi costitutive
```

Le leggi fisiche non devono conoscere RK, IMEX, gruppi o il clock globale. I kernel analitici ricevono un intervallo esplicito, senza conoscere come sia stato ottenuto. I rate ordinari restano quantità per unità di tempo; i proposal limitati da intervallo/inventario mantengono un contratto esplicito distinto, senza eliminare i cap esistenti.

Le routine numeriche devono invece gestire:

- stage;
- time-step;
- cadenze;
- gruppi;
- ordine dei solve;
- integrazione dei rate.

`equation_terms_2d.f90` costituisce il livello di collegamento tra i due ambiti.

---

#### D17.1 — Contratto tra fisica e numerica

I moduli fisici devono:

- definire le leggi costitutive;
- calcolare rate;
- calcolare eventuali soluzioni analitiche;
- fornire kernel reali e complessi dove richiesto dalle capability del percorso;
- non conoscere la struttura RK-IMEX.

`equation_terms_2d.f90` deve:

- ricevere un gruppo;
- selezionare i termini fisici;
- chiamare le leggi appropriate;
- assemblare i contributi nelle equazioni;
- restituire rate o risultati analitici.

Non deve reimplementare le leggi fisiche né eseguire un integratore temporale generale. Delega i solve analitici ai kernel fisici e converte il risultato nel blocco conservativo dichiarato; il coordinatore numerico compone i gruppi e ricava gli incrementi/rate effettivi.

`time_integration_2d.f90` deve gestire il ciclo RK-IMEX e l’applicazione dei gruppi.

`nonlinear_solver_2d.f90` deve occuparsi esclusivamente dei solve numerici impliciti e non deve chiamare direttamente moduli come:

```fortran
thermal_laws_2d
pore_pressure_laws_2d
rheology_2d
```

---

#### D17.2 — Termini fisici elementari

Ogni componente fisica deve essere identificata da un `term_id`.

Esempi:

```fortran
SOURCE_TERM_RHEOLOGY_VOELLMY_COULOMB
SOURCE_TERM_RHEOLOGY_VOELLMY_TURBULENT_DRAG
SOURCE_TERM_THERMAL_RADIATIVE_COOLING
SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION
SOURCE_TERM_PORE_PRESSURE_DISSIPATION
SOURCE_TERM_MASS_EXCHANGE_EROSION
SOURCE_TERM_MASS_EXCHANGE_DEPOSITION
SOURCE_TERM_MASS_EXCHANGE_ENTRAINMENT
SOURCE_TERM_BOTTOM_CIRCULAR_INJECTION
SOURCE_TERM_BOTTOM_FISSURAL_INJECTION
SOURCE_TERM_RADIAL_INJECTION
SOURCE_TERM_BOUNDARY_INJECTION
```

Ogni termine deve avere:

- equazioni interessate;
- modulo fisico responsabile;
- interfacce disponibili;
- metodo;
- cadenza;
- gruppo.

##### Implementazione

Devono essere definiti tutti i `term_id` in un modulo comune. La selezione dei termini non deve dipendere dalla posizione delle chiamate nel codice.

---

#### D17.3 — Metodi e cadenze

Ogni termine deve essere organizzato mediante:

```text
method  = 'explicit' | 'newton' | 'analytic'
cadence = 'stage'   | 'step'
group   = intero positivo
```

Il metodo determina il trattamento:

| Metodo | Interfacce |
|---|---|
| `explicit` | esplicita |
| `newton` | semi-implicita e implicita |
| `analytic` | analitica |

Non deve essere introdotto un parametro indipendente `treatment`.

La cadenza determina il momento della valutazione:

- `stage`: dentro il ciclo RK-IMEX;
- `step`: una volta per time-step, fuori dal ciclo RK.

##### Implementazione

Devono essere definiti enum o costanti per:

- metodo;
- cadenza;
- `term_id`;
- variabili target;
- solve impliciti;
- solve analitici.

---

#### D17.4 — Appartenenza ai gruppi

Ogni contributo fisico attivo deve essere registrato una sola volta, senza duplicarlo tra `stage` e `step`. Gli assemblaggi derivati (calore/composizione/pressione associati a una sorgente) seguono il padre e non costituiscono nuove applicazioni indipendenti.

Sono distinti:

```text
(stage, group=1)
(step,  group=1)
```

All’interno di un gruppo deve essere presente un solo metodo. L’identità è `(cadence,id)`; un ID `stage` non può indicare contemporaneamente un gruppo esplicito e uno analitico. ID e ordine di esecuzione sono distinti: conservare fase temporale e ordine mediante metadati espliciti.

##### Implementazione

Devono essere gestiti separatamente:

```fortran
stage_groups(:)
step_groups(:)
```

L’inizializzazione deve verificare:

- unicità dell’appartenenza;
- uniformità del metodo;
- validità del gruppo;
- coerenza tra termine e cadenza.

---

#### D17.5 — Struttura del gruppo

La struttura minima raccomandata è:

```fortran
TYPE :: source_group_type
   INTEGER :: id
   INTEGER :: cadence
   INTEGER :: method
   INTEGER :: execution_order
   INTEGER :: legacy_id = 0

   INTEGER, ALLOCATABLE :: term_ids(:)
   INTEGER, ALLOCATABLE :: equation_ids(:)

   INTEGER, ALLOCATABLE :: implicit_solver_ids(:)
   INTEGER, ALLOCATABLE :: analytic_solver_ids(:)
END TYPE source_group_type
```

Devono essere presenti due registri distinti:

```fortran
TYPE(source_group_type), ALLOCATABLE :: stage_groups(:)
TYPE(source_group_type), ALLOCATABLE :: step_groups(:)
```

Deve inoltre essere disponibile una funzione del tipo:

```fortran
LOGICAL FUNCTION group_contains(group, term_id)
```

---

I descrittori dei solve devono associare ID, blocco target, termini e ordine. Il catalogo dei termini deve distinguere rate di cella, stato al bordo e assemblaggio derivato, con capability REAL/COMPLEX e bundle ammessi. Questi dati non introducono un nuovo selettore temporale `treatment`.

#### D17.6 — Gruppi `analytic`

Un gruppo `analytic` può contenere un solo solve analitico.

Il solve aggiorna un solo blocco target, mantenendo congelati gli altri: temperatura, pressione dei pori oppure velocità vettoriale `(u,v)`. Le due componenti della velocità appartengono allo stesso solve; non sono due solve scalari indipendenti.

Valido:

```text
gruppo analytic:
    solve analitico della temperatura
    contributo radiativo
    contributo convettivo
```

Non valido:

```text
gruppo analytic:
    solve analitico della temperatura
    solve analitico della pressione dei pori
```

Nel nuovo formato nativo i due solve devono essere assegnati a gruppi distinti. Gli input legacy già validi che condividono un gruppo analitico restano accettati: l’adattatore li espande in gruppi interni distinti mantenendo l’ordine attuale velocità → temperatura → pressione, con ricostruzione dello stato tra i solve. I nuovi ID interni non ridefiniscono l’ordine né riscrivono i campi legacy nei backup.

##### Implementazione

Per ogni gruppo analitico deve essere verificato:

```text
SIZE(analytic_solver_ids) == 1
```

Deve inoltre essere definita la variabile target.

---

#### D17.7 — Gruppi `newton`

Un gruppo Newton può contenere più solve impliciti.

Esempio:

```text
gruppo newton:
    solve implicito della velocità
    solve implicito della temperatura
    solve implicito della pressione dei pori
```

I solve realmente indipendenti possono essere eseguiti secondo un ordine definito. Un blocco accoppiato resta però un unico solve Newton con tutte le equazioni e derivate incrociate: le configurazioni legacy con velocità/energia/pressione nello stesso gruppo non devono essere trasformate in sweep scalari. Metadati espliciti distinguono blocchi accoppiati e solve indipendenti; sovrapposizioni ambigue dei target sono rifiutate.

##### Implementazione

Ogni gruppo Newton deve contenere uno o più:

```fortran
implicit_solver_ids(:)
```

Ogni solve deve avere:

- identificatore;
- variabile target;
- termini associati;
- ordine di esecuzione;
- interfaccia numerica.

---

#### D17.8 — Interfacce comuni

Le interfacce generiche devono essere definite in `equation_terms_2d.f90`:

```fortran
eval_expl_terms
eval_semi_implicit_terms
eval_implicit_terms
eval_analytic_terms
```

Ogni routine di rate/solve deve:

1. ricevere il gruppo e, quando necessario, lo specifico solver_id;
2. inizializzare il risultato;
3. selezionare i termini;
4. chiamare le leggi fisiche;
5. assemblare il risultato;
6. restituire rate o risultato analitico.

Firme indicative:

```fortran
SUBROUTINE eval_expl_terms(group, q, qp, time, rate)
```

```fortran
SUBROUTINE eval_semi_implicit_terms(group, q, qp, time, rate)
```

```fortran
SUBROUTINE eval_implicit_terms(group, q, qp, time, rate)
```

```fortran
SUBROUTINE eval_analytic_terms(group, solver_id, q_before, context, interval, q_after)
```

`q_after` è lo stato conservativo completo: copiato da `q_before`, modifica solo il blocco target. Il wrapper converte temperatura/pressione/velocità primitive in Q4/Qp/Q2,Q3. Il coordinatore numerico passa `interval=dt*a_ii` negli stage e ricava il contributo IM effettivo `(q_stage-q_known)/(dt*a_ii)` dopo le correzioni; il percorso con diagonale zero non divide per tale intervallo. `analytic` conserva l’exact-flow a coefficienti congelati, non diventa un solve backward Euler.

Le firme sono indicative: prima del codice va definito un contesto read-only con geometria, Hessiana, gravità, stato, valore stocastico congelato, coperture/arco delle sorgenti, snapshot delle finestre e tempo attuale. L’interfaccia al bordo `eval_source_bdry` resta distinta dai vettori di rate: restituisce uno stato primitivo alla ricostruzione e quindi al flusso spaziale.

---

#### D17.9 — Interfaccia esplicita

L’interfaccia esplicita deve assemblare tutti i contributi espliciti richiesti.

Per esempio:

\[
R_{\mathrm{explicit}}
=
R_{\mathrm{Coulomb}}
+
R_{\mathrm{turbulent\ drag}}.
\]

Le formule fisiche non devono essere duplicate.

##### Implementazione

Modificare `eval_expl_terms` affinché:

- inizializzi il vettore `rate`;
- selezioni i `term_id`;
- chiami i kernel fisici;
- sommi i contributi;
- distingua i gruppi `stage` dai gruppi `step`.

---

#### D17.10 — Interfacce semi-implicite e implicite

Un gruppo Newton può contenere più contributi semi-impliciti e più solve impliciti.

Esempio:

```text
semi-implicito:
    Coulomb friction

implicito:
    turbulent drag
    thermal relaxation
```

##### Implementazione

Modificare:

```fortran
eval_semi_implicit_terms
eval_implicit_terms
```

per selezionare il gruppo, il termine e lo specifico `implicit_solver_id`.

Devono essere supportate versioni reali e complesse quando richieste dal solutore.

I valutatori selezionano i termini presenti sia nel gruppo sia nel descrittore
del solver; il coordinatore esegue in ordine i solve Newton registrati nel
gruppo. Evidenza e perimetro verificato: [D17.10 — group-filtered implicit
rates](acceptance/D17_10_GROUP_FILTERED_IMPLICIT_RATES.md).

---

#### D17.11 — Riutilizzo dei kernel REAL/COMPLEX

Le formule comuni non devono essere duplicate tra versioni reali e complesse. Riutilizzare kernel complessi con wrapper reali o interfacce generiche già presenti; usare `.inc` matematici ristretti quando utili, senza imporre una nuova relocation delle leggi D16 o una complessificazione dei percorsi oggi reali.

Esempi:

```text
rheology_voellmy_coulomb_rate.inc
rheology_voellmy_turbulent_drag_rate.inc
thermal_radiative_cooling_rate.inc
```

I file `.inc` devono contenere soltanto il kernel matematico comune.

I wrapper devono gestire:

- tipo reale o complesso;
- dichiarazioni;
- argomenti;
- selezione del contributo;
- risultato.

---

#### D17.12 — Interfacce analitiche

Le routine analitiche devono rimanere nei moduli fisici, ma devono essere chiamate attraverso `eval_analytic_terms`.

Il flusso deve essere:

```text
time_integration_2d
        ↓
equation_terms_2d
        ↓
legge fisica analitica
```

##### Implementazione

- aggiungere `eval_analytic_terms`;
- rimuovere le chiamate dirette dai moduli numerici;
- assicurare un solo solve analitico per gruppo;
- utilizzare gruppi distinti per variabili analitiche diverse.

---

#### D17.13 — Configurazione di default compatibile con il codice attuale

In assenza di configurazioni esplicite nel file di input, il codice deve costruire automaticamente una configurazione equivalente al comportamento attuale.

La configurazione di default deve essere determinata dall’inventario delle chiamate effettive presenti nel codice, non soltanto dal nome dei moduli.

**Stato locale:** i default IMEX compilati sono verificati attraverso
`initialize_configured_sources` nei casi senza sorgenti attive e con i solve
analitici lava/energia/pressione abilitati; il test controlla ID, metodo,
ordine, termini e target. Il test grouped-stage esegue inoltre la sequenza
default energia analitica → pressione analitica → quantità di moto Newton;
oracoli exact-flow, residui Newton e rate effettivi verificano gli aggiornamenti
dei target. Resta aperta l’equivalenza per tutte le attivazioni e i percorsi
temporali. Vedere
[D17.13 — default configuration](acceptance/D17_13_DEFAULT_CONFIGURATION.md).

##### D17.13.1 — Sorgenti esplicite `stage`

Il gruppo esplicito `stage` proposto ha **ID 4**, per non collidere con gli ID
impliciti legacy 1–3. Contiene curvatura, iniezione dal fondo circolare,
iniezioni fissurali e iniezione radiale conservativa di cella. Calore,
composizione e pressione dell’iniezione sono assemblaggi dipendenti del padre.

L’iniezione laterale al bordo è catalogata nello stesso percorso esplicito ma
mantiene la rappresentazione **stato primitivo → ricostruzione → flusso** tramite
`eval_source_bdry`: non viene aggiunta come ulteriore rate di cella.
Pressione idrostatica e pendenza ordinaria restano nel path spaziale.

Preparare le finestre mediate una volta per step prima del RK; valutare i rate
EX sullo stato dello stage richiesto, mantenendo ordine, `need_explicit_stage`
e argomento temporale attuali. Lo stato al bordo conserva il suo gate temporale
puntuale. I contributi inattivi sono omessi o nulli senza cambiare le equazioni.

---

##### D17.13.2 — Termini Newton `stage`

Il default della velocità è `stage / group=3 / method=newton`. Il catalogo
elenca per ogni reologia i contributi SI e impliciti effettivi, inclusi i
bundle aggregati (per esempio drag viscoso/turbolento lahar). La tabella
machine-readable è la specifica, non una dicitura generica «termini attuali».

Mantenere SI sullo stato pre-Newton, correzione dissipativa **prima** di Newton,
solve del blocco implicito, correzione successiva e rate effettivo di stage.
Il fix `9b2a4989` non deve essere riaperto dal refactoring. Residuo e Jacobiano
usano lo stesso filtro gruppo/solve e conservano tolleranze e semantica COMPLEX.
La partizione della quantità di moto resta implicita anche senza reologia.

Le configurazioni legacy possono mettere energia/pressione nei gruppi Newton,
anche insieme alla velocità: preservare il solve congiunto e le derivate
incrociate, oltre all’ordine tra gruppi distinti.

---

##### D17.13.3 — Solve analitici

Il default attuale è **energia / gruppo 1 / analytic**, poi
**pressione dei pori / gruppo 2 / analytic**, poi velocità Newton nel gruppo 3.
I primi due sono attivi soltanto quando le rispettive equazioni lo richiedono.
La pressione usa lo stato aggiornato dopo la temperatura.

Sono disponibili anche i solve analitici della velocità delle reologie, quando
selezionati tramite input. Ogni solve analitico consuma il bundle fisico previsto
una sola volta: non aggiungere di nuovo gli stessi rate SI/impliciti.
Le selezioni parziali di bundle non implementate vanno rifiutate.

Nel formato nativo ogni gruppo ha un solo solve e un blocco target. L’adattatore
legacy normalizza i gruppi analitici condivisi senza rifiutare configurazioni
oggi supportate e testate; la specifica allegata descrive ordine e nuovi ID.

---

##### D17.13.4 — Termini espliciti `step`

I termini attualmente applicati fuori dal ciclo RK devono rimanere inizialmente `step`.

Esempio:

```text
step group 1:
    SOURCE_TERM_MASS_EXCHANGE_EROSION
    SOURCE_TERM_MASS_EXCHANGE_DEPOSITION
    SOURCE_TERM_MASS_EXCHANGE_ENTRAINMENT
    SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS
    SOURCE_TERM_MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS
    SOURCE_TERM_MASS_EXCHANGE_CARRIER_DARCY_LOSS
    SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION
```

Configurazione equivalente:

```text
method  = explicit
cadence = step
group   = 1
```

L’elenco completo è nel manifest allegato: erosione, deposizione, perdita
liquida configurata, perdita associata alla deposizione, perdita Darcy,
entrainment e vaporizzazione. La prima migrazione conserva un unico batch
atomico: le perdite del carrier sono sommate prima del limitatore condiviso;
deposizione, vaporizzazione, inventari e topografia restano collegati. I cap
dipendenti da `dt`, il cutoff `alphastot_min`, la maschera delle sorgenti e le
ricostruzioni non possono essere sostituiti da applicazioni indipendenti per ID.

##### D17.13.5 — Inventario obbligatorio dei default

Prima dell’implementazione deve essere costruita una tabella definitiva contenente, per ogni termine:

| Campo | Contenuto |
|---|---|
| `term_id` | identificatore fisico |
| modulo | file responsabile |
| routine attuale | chiamata esistente |
| metodo default | `explicit`, `newton`, `analytic` |
| cadenza default | `stage` o `step` |
| gruppo default | numero del gruppo |
| solve target | variabile aggiornata |
| ordine | posizione nell’elaborazione |
| equazioni | componenti modificate |
| rappresentazione/capability | rate di cella, stato al bordo, assemblaggio derivato; REAL/COMPLEX disponibili |
| stato | stato pre/post solve usato, variabili congelate e refresh richiesti |
| batch | dipendenze, cap condivisi, accumulo e applicazione atomica |
| preparazione/applicazione | snapshot temporale e fase effettiva, distinti dalla sola cadenza |
| side effect | inventari, geometria e controlli applicati dalla numerica |

Questa tabella deve essere considerata la specifica di compatibilità con il codice attuale.

##### D17.13.6 — Criterio di equivalenza

La configurazione di default deve riprodurre:

1. ordine delle chiamate;
2. momento di valutazione;
3. stato utilizzato per il rate;
4. momento di applicazione;
5. aggiornamento degli inventari;
6. aggiornamento della topografia;
7. controlli di positività;
8. ricostruzione delle variabili fisiche.

Il solo fatto di assegnare lo stesso metodo non è sufficiente a garantire la compatibilità.

---

#### D17.14 — Separazione tra `stage` e `step`

I gruppi `stage` vengono valutati nel ciclo RK-IMEX.

I gruppi `step` vengono valutati una sola volta dopo il completamento del RK.

##### Implementazione

Modificare `time_integration_2d.f90` affinché:

- elabori i gruppi `stage` durante il ciclo RK;
- elabori i gruppi `step` in una fase separata;
- non ricalcoli i gruppi `step` negli stage;
- applichi i risultati nell’ordine definito dalla configurazione.

##### Stato ed evidenza locale

Il driver passa il registro validato al fractional step post-RK. La routine
seleziona i gruppi `step` nell'ordine convalidato e applica una sola volta il
batch mass-exchange indivisibile; un registro senza termini `step` non modifica
stato, inventari o topografia. Il batch non viene invocato dal ciclo degli
stage. I gruppi post-RK non mass-exchange sono rifiutati esplicitamente finché
non esiste un kernel che li implementi.

Regressione e limiti: [D17.14 — post-RK step routing](acceptance/D17_14_STEP_ROUTING.md).

---

#### D17.15 — Ordine del time-step

L’ordine compatibile è quello inventariato nella specifica allegata:

1. scelta di `dt`, finestre delle sorgenti mediate sullo step e preparazione OU;
2. per ogni stage: stato noto da contributi precedenti, controlli e conversione;
3. sweep implicita ordinata (default temperatura, pressione, velocità), incluse SI/correzioni;
4. contributo IM effettivo dopo le correzioni;
5. valutazioni EX richieste sullo stato risolto e operatore spaziale con hook al bordo;
6. completamento RK, controlli di positività/asciutto e ricostruzione già presenti;
7. proposal comuni di mass exchange sullo stato post-RK e maschera di copertura;
8. applicazione unica a stato/inventari, controlli/ricostruzione, aggiornamento del letto;
9. avanzamento del tempo, diagnostiche e output.

Non anticipare indiscriminatamente tutti i gruppi EX ai solve dello stage e non
posticipare i controlli già presenti fino al termine dello step. La numerica
coordina; il main conserva inizialmente l’ordine esterno attuale, senza imporre
uno spostamento simultaneo del driver per introdurre i registri.

##### Stato ed evidenza locale

L’ordine del driver è centralizzato in `simulation%advance_timestep` e
conserva: medie delle finestre, preparazione OU condizionale, avanzamento RK
(con solve impliciti ordinati e EX/spaziale sullo stato risolto), assembly
finale e proiezioni, batch mass-exchange con ricostruzione e aggiornamento del
letto, quindi avanzamento dell’orologio. La scelta/limitazione di `dt` resta
nel driver prima della routine; diagnostiche/output restano dopo.

Le regressioni IMEX verificano traiettorie e assembly stage per N_RK=2/3/4,
gruppi misti analitici/Newton e hook spaziali; le regressioni mass-exchange
verificano l’aggiornamento aggregato, inventari, maschere, topografia e il
routing post-RK. Entrambe le suite passano localmente in profili strict e
optimized.

**D17.15 resta parziale:** `TEST_IMEX_STAGES` verifica con un observer le
cinque fasi della routine production e l’avanzamento finale dell’orologio.
La fixture abilita un aggiornamento OU con seed fissato e ne verifica
l’effetto sul campo stocastico dopo la fase dedicata; inoltre applica una
sorgente circolare attiva per metà del timestep e verifica il bilancio integrato
di massa e la positività. La suite `TEST_CURVATURE_SOURCE` verifica a livello
equation-term le finestre per sorgenti circolari, fissurali e radiali. Restano
la selezione/limitazione di `dt`, l’output del driver, il trasporto stocastico e
l’equivalenza di tutte le sorgenti/configurazioni. Evidenza e limiti:
[D17.15 — timestep ordering](acceptance/D17_15_TIMESTEP_ORDER.md).

---

#### D17.16 — Validazione dell’input e dei gruppi

Durante l’inizializzazione devono essere controllati:

- metodo valido;
- cadenza valida;
- gruppo positivo;
- contributo fisico registrato una sola volta globalmente, senza duplicazione stage/step;
- metodo unico nel gruppo;
- interfacce richieste disponibili;
- un solo solve analitico nei gruppi `analytic`;
- uno o più solve impliciti nei gruppi `newton`;
- variabili target definite;
- ordine dei solve non ambiguo;
- termine fisico esistente.

Gli errori devono essere rilevati prima dell’avvio della simulazione. Conservare
la compatibilità dei campi `IMEX_*` attuali e distinguere input esplicitamente
fornito da default inizializzati. Rifiutare la combinazione ambigua di nuova
configurazione nativa e configurazione legacy esplicita, capability mancanti,
bundle parziali non supportati e split di limitatori condivisi non validati.
La namelist nativa opzionale si chiama `SOURCE_GROUPS`; cadenze, metodi,
termini e target sono indicati con nomi simbolici.

##### Stato ed evidenza locale

La validazione esistente controlla i campi legacy `IMEX_*` letti dalla
namelist prima dell’inizializzazione dei workspace: ID gruppo nell’intervallo
1–3 e metodi ANALYTIC/NEWTON supportati. `source_registry_2d` valida inoltre
identità e ordine dei gruppi, cadenza/fase/metodo, appartenenza dei termini,
duplicazioni globali, target/equazioni, coverage e ordine dei solver, bundle
atomici e descrittori mancanti o non validi; l’inizializzazione fallita libera
i dati parziali.

Il parser rileva namelist duplicate e input malformati e rifiuta la
combinazione con campi legacy espliciti; i sentinels permettono di distinguere
questi ultimi dai default. Il builder nativo richiede ogni termine fisico
attivo esattamente una volta, verifica i target impliciti e valida gli indici
equazione prima di costruire partizione e registro.

`TEST_SOURCE_REGISTRY/run_test.sh` passa in strict e optimized: 864
configurazioni legacy, 14.893 asserzioni e 28 configurazioni invalide
rifiutate per profilo; include casi parser validi, assenti, malformati,
duplicati e in conflitto. `TEST_MASS_EXCHANGE/run_test.sh` compila i moduli
di produzione, incluso `inpout_2d`, e passa in entrambi i profili. Report:
[D17.16 — input and source-schedule validation](acceptance/D17_16_VALIDATION.md).

`TEST_IMEX_STAGES/run_test.sh` passa in strict e optimized con una fixture
nativa che attraversa `init_param`/`read_param` e costruisce il registry per
un blocco velocity e un termine mass-exchange attivo.

**D17.16 resta parziale:** la fixture end-to-end copre una configurazione
rappresentativa; equivalenza per tutte le combinazioni fisiche e accettazione
cross-platform/CI restano aperte.

---

#### D17.17 — Riorganizzazione dei moduli

I moduli fisici devono essere organizzati, per esempio, in:

```text
rheology_2d.f90
thermal_laws_2d.f90
pore_pressure_laws_2d.f90
```

Essi non devono contenere:

- ciclo RK;
- gestione degli stage;
- gestione delle cadenze;
- selezione dei gruppi;
- solutore Newton;
- integrazione generale del time-step.

Le routine di assemblaggio devono rimanere in `equation_terms_2d.f90`.

##### Implementazione

- spostare le formule fisiche nei moduli `*_laws_2d`;
- spostare l’assemblaggio in `equation_terms_2d`;
- eliminare le formule duplicate;
- verificare che i moduli fisici non dipendano dalla numerica temporale.

**Evidenza locale D17.17:** il batch atomico mass-exchange post-RK è stato
estratto da `time_integration_2d` in `mass_exchange_step_2d`; l’integratore
conserva il tableau RK e la sequenza degli stage, mentre la routine condivisa
del timestep ne mantiene l’invocazione dopo RK e prima del clock. Rate,
maschere, aggiornamento degli inventari e topografia conservano l’ordine
precedente. Un controllo statico dei confini tra solver/integratore, bridge
equazionale e leggi fisiche è aggiunto ai test ACCEPTANCE. L’audit degli altri
percorsi e del call graph completo resta aperto:
[D17.17 — post-RK step-source separation](acceptance/D17_17_STEP_SOURCE_SEPARATION.md).

---

#### D17.18 — Adattamento delle routine esistenti

Devono essere riorganizzate:

```fortran
eval_expl_terms
eval_semi_implicit_terms
eval_implicit_terms
eval_mass_exchange_terms
```

La riorganizzazione deve:

- utilizzare i `term_id`;
- rispettare gruppi e cadenze;
- distinguere rate da aggiornamenti temporali;
- mantenere inizialmente il comportamento fisico attuale;
- eliminare chiamate duplicate;
- preservare l’ordine dei default.

Il trasferimento di massa deve essere integrato progressivamente nella struttura generale, mantenendo inizialmente la cadenza attuale fuori dal RK.

Il batch `step` post-RK è implementato localmente in `mass_exchange_step_2d`
ed è invocato da `simulation%advance_timestep`; la valutazione dei contributi
fisici resta in `equation_terms_2d`. La selezione term-selective resta
intenzionalmente differita finché i kernel e i limitatori atomici non sono
validati singolarmente.

---

#### D17.19 — Sequenza implementativa

Stato: IN CORSO. Le voci completate indicano l’implementazione locale e
l’evidenza disponibile, non il completamento dell’accettazione D17.

1. Inventariare tutte le chiamate attuali ai termini sorgente (fase A documentata nella specifica/manifest allegati). **COMPLETATO.**
2. Definire i `term_id`. **COMPLETATO.**
3. Definire metodi, cadenze e identificatori dei solve. **COMPLETATO nei metadati.**
4. Costruire la tabella completa dei default. **COMPLETATO come baseline e adattatore legacy; la costruzione dei default è verificata nei casi inattivo e lava/energia/pressione attivi, mentre l’equivalenza runtime complessiva resta da verificare (D17.13).**
5. Definire `source_group_type`. **COMPLETATO.**
6. Creare i registri `stage_groups` e `step_groups`. **COMPLETATO.**
7. Definire le firme delle quattro interfacce. **COMPLETATO per i percorsi supportati:** semi-implicita e implicita accettano anche gruppo e descrittore solver, mantenendo i call site legacy.
8. Implementare i gruppi di default. **COMPLETATO nei metadati; il bridge predefinito è testato nei casi inattivo e attivo riportati in D17.13; dispatch ed equivalenza numerica completa con il comportamento legacy restano parziali.**
9. Implementare il parsing dell’input. **IMPLEMENTATO localmente:** namelist opzionale `SOURCE_GROUPS`, rilevamento duplicati, errori sintattici e conflitti con campi legacy espliciti; parser e fixture rappresentativa end-to-end `init_param`/`read_param` testati in strict/optimized (D17.16).
10. Implementare i controlli di consistenza. **PARZIALE:** validazione legacy e nativa; 864 configurazioni legacy, 14.893 asserzioni e 28 configurazioni invalide per profilo (D17.16). Il builder richiede ogni termine fisico attivo esattamente una volta e verifica gli obiettivi impliciti; equivalenza per tutte le combinazioni fisiche e accettazione cross-platform/CI restano aperte.
11. Modificare `eval_expl_terms`. **COMPLETATO localmente per il routing esplicito `stage` e l’hook al bordo (D17.9).**
12. Modificare `eval_semi_implicit_terms`. **IMPLEMENTATO localmente:** seleziona i termini semi-impliciti del gruppo e del solver Newton indicati.
13. Modificare `eval_implicit_terms`. **IMPLEMENTATO localmente:** seleziona rate reologici, raffreddamento e dissipazione della pressione dei pori dal gruppo/solver; routing REAL/COMPLEX mantenuto.
14. Aggiungere `eval_analytic_terms`. **COMPLETATO localmente (D17.12).**
15. Modificare `time_integration_2d.f90`. **PARZIALE:** dispatch dei gruppi `stage` completato localmente (D17.9); la routine condivisa production `simulation%advance_timestep` e il test dell’ordine delle fasi sono implementati localmente (D17.15). Il batch post-RK `step` è delegato a `mass_exchange_step_2d`; la sorgente circolare a finestra parziale è verificata anche nel timestep production, mentre le medie equation-term coprono impulsi, rampe, wrap, inattività e le sorgenti fissurali/radiali. Restano da verificare selezione di `dt`, output driver e tutte le combinazioni fisiche.
16. Modificare `nonlinear_solver_2d.f90`. **PARZIALE:** descrittori gruppo/solver propagati al residuo, alla Jacobiana complessa, al line search e al rate semi-implicito; sweep ordinato dei solver di un gruppo implementato localmente (D17.10). Restano i percorsi numerici e di cadenza ancora non integrati.
17. Rimuovere le chiamate fisiche dirette dai moduli numerici. **PARZIALE, audit statico completato per i confini dei moduli:** nessun sorgente in `numerics/` o `spatial/` importa direttamente moduli `*_laws_2d`; l’uso spaziale di `equation_terms_2d` è limitato a un’allowlist di API per flussi, coefficienti idrostatici, velocità caratteristiche, hook al bordo e limiter di flusso. Le leggi non importano moduli di orchestrazione numerica/spaziale. Restano audit del call graph e verifica runtime/term-selective.
18. Riutilizzare kernel REAL/COMPLEX comuni; introdurre `.inc` solo dove utili. **IMPLEMENTATO localmente per la duplicazione individuata:** la correlazione di viscosità lahar ha un solo kernel COMPLEX; non è stato introdotto un `.inc` non necessario. Evidenza: [D17.11 — REAL/COMPLEX kernel reuse](acceptance/D17_11_REAL_COMPLEX_KERNEL_REUSE.md).
19. Creare wrapper reali e complessi. **IMPLEMENTATO localmente per la viscosità lahar:** il wrapper REAL conserva l’interfaccia generica e riusa il kernel COMPLEX; test di valore e derivata nei due rami della correlazione. La matrice D17.20 resta aperta.
20. Organizzare i solve impliciti multipli. **IMPLEMENTATO localmente:** il coordinatore esegue i solver Newton nell’ordine registrato e limita ciascun solve ai propri termini/equazioni; il test runtime controlla i residui termico e di pressione dei pori nello stato sequenziale. Il blocco legacy accoppiato resta un unico solve.
21. Organizzare i solve analitici singoli per gruppo. **COMPLETATO localmente per i gruppi supportati e verificato dal report D17-C.**
22. Integrare il trasferimento di massa. **IMPLEMENTATO localmente per il batch post-RK:** il registro `step` valida il batch mass-exchange atomico; valutazione dei rate, maschere, aggiornamento di stato/inventari e topografia sono eseguiti da `mass_exchange_step_2d` dopo RK e prima dell’avanzamento del clock. Il dispatch term-selective e la matrice completa restano aperti.
23. Verificare l’ordine delle operazioni. **PARZIALE:** ordine dei gruppi analitici legacy, routing esplicito e sequenza production source-window/OU/RK/step-source/clock verificati; aggiornamento OU attivo, media temporale equation-term di sorgenti circolari/fissurali/radiali e sorgente circolare a finestra parziale nel timestep production verificati; selezione di `dt` e ordine driver completo restano da validare.
24. Eseguire i test di regressione. **PARZIALE:** suite locali mirate B, C, D17.9, D17.10, D17.15 e D17.17 superate, incluse reologia REAL/COMPLEX, registro, mass-exchange e IMEX strict/optimized; il controllo architetturale statico passa. Le regressioni source-window verificano impulso, rampa lineare/cosinusoidale, attraversamento del periodo, intervallo inattivo, contributo fissurale e bilancio production. Matrice D17.20 completa e gate CI/piattaforme non completati.

---

#### D17.20 — Test obbligatori

Stato ufficiale: **PARZIALE**. Le verifiche locali aggiornate coprono molte
condizioni strutturali e numeriche elencate sotto; i report B, C, D17.9, D17.10,
D17.15 e D17.17 attestano solo i percorsi e le asserzioni effettivamente
eseguiti. Restano aperte l’equivalenza dei default per tutte le combinazioni
fisiche, le finestre temporali di tutte le sorgenti nel percorso production,
la selezione di `dt` e l’output del driver, nonché l’accettazione delle
piattaforme/CI e l’approvazione scientifica delle baseline.

**Evidenza locale eseguita il 2026-10-23:** `TEST_SOURCE_REGISTRY` passa in
strict/optimized con 864 combinazioni legacy, 14.893 asserzioni e 28
configurazioni invalide; `TEST_RHEOLOGY_CELL_LOCAL` passa con 255/255 fixture
per profilo, inclusi rate REAL/COMPLEX e solve Newton ordinati; `TEST_IMEX_STAGES`
passa in strict/optimized con 1/4 thread, tableau, gruppi nativi e misti,
ordine production e budget di massa della sorgente temporizzata; `TEST_MASS_EXCHANGE`
passa in strict/optimized con 1/4 thread; `TEST_CURVATURE_SOURCE` e
`TEST_STATE_CONVERSION` passano; i controlli architetturali D17 e i contratti
Doxygen passano. Questi risultati non equivalgono all’accettazione end-to-end
completa della matrice.

Devono essere verificati:

1. assenza di configurazioni nell’input;
2. sorgenti esplicite `stage`;
3. sorgenti esplicite `step`;
4. gruppo Newton con un solve implicito;
5. gruppo Newton con più solve impliciti;
6. gruppo analitico con un solo solve;
7. rifiuto di due solve analitici nello stesso gruppo;
8. gruppi `stage` e `step` con lo stesso numero;
9. termine assegnato a due gruppi con la stessa cadenza;
10. gruppo con metodi diversi;
11. metodo Newton incompleto;
12. metodo analitico senza solve analitico;
13. verifica reale/complessa;
14. sorgenti circolari, fissurali e al bordo;
15. trasferimento di massa fuori dal RK;
16. conservazione della massa;
17. positività;
18. ordine dei solve;
19. equivalenza con i default attuali;
20. indipendenza dei moduli numerici dai moduli fisici;
21. conservazione dei gruppi Newton accoppiati e del residuo congiunto;
22. compatibilità dei gruppi analitici legacy condivisi e ordine dei refresh;
23. assenza di collisione tra gruppo EX e gruppi IM nello stesso namespace;
24. rifiuto della doppia registrazione stage/step e dell’input nativo/legacy ambiguo;
25. hook al bordo senza doppio rate di cella, e sorgente di curvatura;
27. preservazione della correzione SI pre-Newton e dei rate effettivi corretti;
28. batch mass-exchange, limitatori condivisi, inventari e topografia;
29. rifiuto di selezioni parziali non supportate e confronto default/legacy normalizzati.

---

#### D17.21 — Criteri di completamento

Stato ufficiale: **DA FARE** — nessun criterio complessivo è dichiarato
soddisfatto finché restano aperti i requisiti e i test indicati sopra.

Il D17 sarà completato quando:

- ogni termine fisico sarà identificato da un `term_id`;
- la tabella dei default sarà completa;
- ogni contributo fisico sarà registrato una volta senza doppia integrazione tra cadenze;
- ogni gruppo avrà un solo metodo;
- ogni gruppo `analytic` avrà un solo solve analitico;
- ogni gruppo `newton` potrà contenere più solve impliciti;
- i solve avranno variabili target e ordine definiti;
- tutte le chiamate ai termini passeranno da `equation_terms_2d.f90`;
- `nonlinear_solver_2d.f90` non conoscerà i modelli fisici;
- `time_integration_2d.f90` non conterrà formule costitutive;
- le leggi fisiche non conosceranno RK, IMEX o gruppi;
- i rate reali e complessi saranno ottenuti senza duplicazioni;
- i gruppi `stage` e `step` saranno indipendenti;
- default e input legacy supportati riprodurranno il comportamento attuale, inclusi gruppi condivisi e accoppiati;
- l’ordine temporale sarà documentato e compatibile con correzioni, finestre e limitatori condivisi;
- gli errori saranno rilevati all’inizializzazione;
- i test di regressione saranno superati.

