# IMEX-SfloW2D — Architettura e roadmap two-layer

**Stato:** 16 novembre 2026 · branch `feature/two-layer` · `N_LAYERS=2` disabilitato in produzione.  
**Obiettivo:** completare il modello con L1 inferiore e L2 superiore, conservando il comportamento single-layer. Questa e' la roadmap operativa; non e' un registro dei commit. L'audit dettagliato e' in [TWO_LAYER_CODE_AUDIT.md](TWO_LAYER_CODE_AUDIT.md).

## 1. Principi e architettura vincolante

- **Una sola rappresentazione autorevole:** `q/qp(component,layer,x,y)` e workspace rank-4; indici di componente comuni ai due layer; slot opzionali presenti se richiesti da almeno un layer (OR dei flag) e **esattamente nulli** nel layer che non li usa. Nessuna copia persistente sincronizzata, coefficiente globale duplicato o architettura transitoria.
- **Proprieta' chiare:** configurazioni e attivazioni reologiche, termiche, pore/Darcy, stochastic/OU e mass exchange appartengono al layer; topografia `B`, geometria, atmosfera, substrato, `deposit_tot`, `T_ground`, `T_ambient` e proprieta' delle specie solide sono condivisi. Distinguere packing di ammissibilita', hindrance e packing fisico.
- **Due layer, strumenti esistenti:** riutilizzare HP-PCCU, ricostruzione, quadratura dei path, registro Lie, classificazione dei contatti e assemblatore mass exchange; workspace distinti solo dove serve uno stato per layer. Non introdurre framework a N layer, dispatcher, cache o wrapper equivalenti a oggetti gia' presenti.
- **Responsabilita':** `physics` contiene leggi e rate cell-local, senza indici di griglia; `spatial` contiene flussi e termini non conservativi; `numerics` gestisce workset, contatti, candidati, splitting Lie e controllo finale. Nessuna formula fisica nel Lie.
- **Contatti:** riutilizzare direttamente i quattro casi locali (entrambi, nessuno, solo L1, solo L2). Prima di ogni gruppo usare i candidati correnti, determinare all'esterno di `physics` le esposizioni a terreno/atmosfera e passarle alle leggi locali. Anche L2 solo puo' toccare il terreno e L1 solo l'atmosfera. Non duplicare contatti o inventari.
- **Compatibilita':** mantenere i nomi di input single-layer concordati (`<FAMILY>_PARAMETERS` oppure `<FAMILY>_L1_PARAMETERS`, mutuamente esclusivi sullo stesso oggetto L1); per two-layer richiedere configurazioni esplicite L1/L2. BC canoniche `bc(component,layer)` su W/E/S/N; BC omogenee per slot inattivi. Nessuna interfaccia interna obsoleta conservata per compatibilita'.

### Contratto fisico e numerico

- Riferimento: *IMEX_SfloW2D_v3*. Preservare equazioni, dimensioni, conservazione, path, splitting, limitatori, regolarizzazioni e ordine numerico senza approvazione per modifiche scientifiche.
- Superfici **derivate**: `eta1=B+h1` ed `eta2=B+h1+h2`. La metrica condivisa `G=1/(1+B_x^2+B_y^2)` e la curvatura derivano da `B`, non da `B+h1`; velocita' e sorgenti cinematiche restano specifiche del layer.
- Pressioni interlayer, con `L=rho_m2*gprime2*h2`: `S_inter1=-G*h1*d_j(L)`, `S_inter2=-L*d_j(G*h1)`. Evitare il doppio conteggio del contributo gia' incluso nel path idrostatico L2. Il clamp L2 della pressione buoyant resta quello autorizzato; non compensare artificialmente eventuali difetti di bilancio.
- Un solo Lie del primo ordine con timestep comune, gruppo 0, ordine dei gruppi, cell-outer/OpenMP e arresto fatale su errore, senza rollback o retry. I controlli complessivi si eseguono a fine integratore; i controlli intermedi restano solo se richiesti dal limitatore, dalle conversioni o dal solver. Scope `L1`, `L2`, `INTERLAYER`; ogni scambio interno usa un unico rate sullo stesso snapshot, trasferimenti opposti e disponibilita' condivise e controllo conclusivo dei candidati e delle risorse SHARED.

#### Semantica vincolante del gruppo 0 (spaziale + sorgenti esplicite)

Il **gruppo 0 comprende l'operatore spaziale e le sorgenti esplicite assegnate al gruppo 0**.
L'operatore spaziale e' eseguito per primo. I rate delle sorgenti del gruppo 0 vengono
**valutati sullo stato iniziale** `q^n/qp^n`, mentre il loro contributo viene applicato
**al candidato dopo lo stadio spaziale** `q^*`. Questi due stati non sono intercambiabili.

Con `L` definito come contributo spaziale con il segno gia' adottato dall'integratore,
la sequenza concettuale e':

```text
q*       = q^n + dt * L(q^n)
R0       = S0(q^n, qp^n, ... )
R0_limit = Limit(R0, q*, dt, inventari e vincoli fisici)
q^(0)    = q* + dt * R0_limit
```

**Il limitatore interviene prima dell'aggiornamento**: considera congiuntamente i
rate proposti, la disponibilita' nel candidato destinatario `q^*`, il passo
temporale e gli inventari pertinenti. Non e' corretto limitare i rate utilizzando
solo `q^n`, ne' valutarli su `q^*` al posto di `q^n`. Per i processi accoppiati,
la limitazione deve mantenere coerenti rate e trasferimenti associati, senza
alterare i bilanci di massa, specie, momento ed energia previsti dalle leggi fisiche.
Restano applicabili i vincoli specifici sui prelievi e sulle disponibilita'
preesistenti gia' stabiliti dai kernel, incluso il divieto di consumare nello
stesso passo materiale appena iniettato, dove previsto: non sostituirli con
un nuovo limitatore generico.

Per `N_LAYERS=2`, la stessa sequenza vale per entrambi i candidati L1/L2 sul
**medesimo timestep**. Le sorgenti locali usano i rispettivi stati iniziali
e le configurazioni di layer; i trasferimenti interlayer usano un unico rate
calcolato sul medesimo snapshot, con contributi opposti e limitazione
**congiunta** rispetto ai candidati destinatari di entrambi i layer e alle
risorse condivise. I contatti per i gruppi successivi si ricavano invece dai
candidati correnti, dopo l'accettazione del gruppo precedente. Gli stati locali sono aggiornati in sequenza; un errore arresta il programma senza ripristini o tentativi di recupero.

I gruppi `1,2,...` seguono in ordine, valutando le proprie sorgenti sullo
stato candidato aggiornato dai gruppi precedenti e usando il metodo configurato.
La separazione tra `spatial` e `lie_source_executor` e' organizzativa, non
cambia l'appartenenza matematica dello spaziale al gruppo 0.
**La migrazione all'esecutore unico L1/L2 non deve cambiare questa semantica,
il limitatore esistente, la sequenza di valutazione o l'ordine floating-point
senza un'esplicita decisione numerica e le relative regressioni.**
- Verificare ogni modifica con audit delle dipendenze (chiamanti, `.inc` REAL/COMPLEX, Makefile, input, I/O, test), stile `docs/FORTRAN_STYLE.md`, build optimized/strict e regressioni pertinenti. Non correggere la fisica per far passare i test. Nessun merge in `master` senza autorizzazione.

## 2. Stato effettivo e lavoro di implementazione

**Priorita' attuale: completare il percorso two-layer end-to-end.** Le verifiche strutturali e di integrita' restano obbligatorie durante lo sviluppo; l'analisi scientifica estesa e' una fase successiva. "Predisposto" non significa abilitato o validato in produzione.

| Fase | Stato | Completamento richiesto |
|---|---|---|
| **TL00–TL04 — Stato, configurazioni e fisica per layer** | Predisposti | Layout canonico e ownership L1/L2, migrazioni REAL/COMPLEX e SHARED presenti. Chiudere i percorsi non ancora esercitati: slot opzionali, alias/parser, sorgenti, IC/BC e regressioni single-layer. |
| **TL05 — OperatorI L2 isolati** | Predisposti e testati | Ricostruzione HP, BC, pressioni e trasporto L2 isolati; completare copertura integrata di stati dry, sorgenti e variabili opzionali senza creare un secondo operatore. |
| **TL06 — Spaziale accoppiato** | **Implementazione parziale verificata** | `evaluate_coupled`, workset OR, bound preliminari con `h1+h2`, CFL comune, path e residui interlayer sono presenti e testati isolatamente (quattro stati wet/dry, stratified rest, geometrie Q1). Il trial Euler verifica l'ammissibilita' dei due candidati; un caso separato a supporto compatto verifica bilanci globali di massa, energia e specie senza flussi di contorno. **Collegare l'operatore alla transazione Lie completa** (TL07), quindi verificare i bilanci negli stadi con sorgenti e i casi dry/BC. |
| **TL07 — Lie accoppiato** | **Avviato** | Il candidato spaziale per cella/layer usa un unico kernel privato dell'integratore; il test D20 ora verifica lo stop fatale dopo un errore tardivo, senza ripristinare lo stato. La transazione Lie richiama ora `evaluate_coupled` e verifica entrambi i candidati sul workset OR, con tracce CFL L1/L2 gia' preparate; il percorso L1 resta invariato. Il percorso L2 si interrompe esplicitamente prima delle sorgenti e del commit, fino a TL07-B. La pubblicazione finale dei campi fluidi comprende tutti gli slot attivi in un'unica transazione, ma resta bloccata per L2 fino all'esecutore scoped completo. Il coordinatore mass exchange riceve ora direttamente i due flag locali di contatto dai candidati correnti: terreno per erosione, deposizione esterna, perdita di carrier, vaporizzazione e immobilizzazione lava; atmosfera per entrainment. Per la produzione single-layer le autorizzazioni rimangono equivalenti al comportamento preesistente. Una regressione copre le quattro configurazioni locali verificando i rate di perdita di carrier ed entrainment. Restano da integrare e validare i gruppi scoped nel percorso end-to-end, gli scambi `INTERLAYER`, la buoyancy L2 e i controlli conclusivi; i descrittori Lie ora distinguono L1/L2 e accettano la stessa identita' di gruppo nei due layer; il configuratore registra entrambi i layer usando lo stesso percorso e `simulation_2d` inizializza i descrittori fisici L1/L2 quando richiesti, l'esecutore ha ora candidati rank-4, usa un unico ciclo cell-outer per gruppi L1/L2 e conserva lo stesso limitatore del gruppo 0; la chiamata two-layer dall'integratore resta bloccata fino al passaggio dei contatti alle leggi di scambio, alla fisica interlayer e alla validazione end-to-end; nessun fallback silenzioso tra metodi. |
| **TL08–TL11 — Scambi e chiusure** | Da integrare | Settling L2→L1, elutriazione L1→L2 e cold-start, drag interlayer, entrainment e scambi esterni con terreno/atmosfera. Condividere disponibilita' del substrato e depositi; verificare specie, momento, energia/termica e limiti di ammissibilita'. |
| **TL12 — I/O e integrazione** | Parziale (input/BC) | Input two-layer completo, output/diagnostiche, NetCDF, restart e continuazione; casi representative EXAMPLES; preservare input e risultati single-layer. |
| **TL13a — Gate funzionale** | Non superato | Percorso end-to-end coerente, nessuna struttura obsoleta, controlli essenziali del §3 superati in build strict/optimized, OpenMP e regressioni. Solo dopo, valutare l'abilitazione **controllata** di `N_LAYERS=2` per test integrati; non equivale a validazione scientifica o rilascio. |
| **TL13b — Validazione scientifica e rilascio** | Successiva al gate funzionale | Benchmark fisici, convergenza spazio/tempo, regimi critici, iperbolicita' e shear, prestazioni, limiti delle chiusure e verifica di soluzioni fisicamente affidabili. Eventuali modifiche alle equazioni richiedono decisione scientifica separata. Abilitazione generale solo dopo approvazione. |

**Sequenza:** chiudere TL06 e TL07, poi TL08–TL11; completare in parallelo le lacune TL00–TL05 e TL12 senza bloccare lavori indipendenti. Successivamente TL13a e TL13b. Non dichiarare chiusa una fase sulla base dei soli test di un kernel isolato.

### Contratto degli scambi `INTERLAYER` (non ancora eseguibili)

Lo scope `INTERLAYER` e' gia' definito, ma **nessuna sorgente fisica
interlayer e' registrata nel catalogo**. Il registro deve rifiutare
con diagnostica esplicita un gruppo con tale scope, anche se vuoto
o contenente una legge locale L1/L2. Non introdurre gruppi segnaposto
o dispatcher senza leggi fisiche.

Per attivare uno scambio occorrono:

1. Legge **cell-local** definita scientificamente per settling L2->L1,
   elutriazione L1->L2 o drag, compresi parametri, unita' e comportamento
   quando un layer e' inizialmente asciutto.
2. **Un solo rate per trasferimento**, con contributi opposti alle
   grandezze conservate nei due layer; budget di momento ed energia,
   compresa la dissipazione, esplicitamente definiti.
3. **Limitazione congiunta** prima dell'applicazione in funzione
   dei due candidati destinatari, del timestep e degli inventari
   condivisi. Per il gruppo 0 i rate si valutano su `q^n`, ma
   si limitano e applicano al candidato post-spaziale `q^*`.
4. Un'unica iterazione cell-outer nell'ordine Lie, con aggiornamento
   dei due candidati locali, nessun rollback e stop fatale su errore.
5. Regressioni su conservazione di massa/specie, momento, energia,
   wet/dry nei quattro casi, cold-start, risorse saturate,
   REAL/COMPLEX e OpenMP, senza alterare i risultati single-layer.

Il solo scope nel registro **non** equivale a implementare uno
scambio interlayer, e non autorizza `N_LAYERS=2`.

## 3. Verifiche obbligatorie durante l'implementazione

1. **Single-layer invariato:** test quantitativi, EXAMPLES rappresentativi, input, ordine Lie, sorgenti, deposizione, restart/NetCDF e OpenMP per le parti toccate; nessuna variazione fisica o di floating-point non autorizzata.
2. **Ammissibilita' e wet/dry:** tutti e quattro i casi locali, ingressi Dirichlet su celle asciutte, comparsa/scomparsa dei layer, nessuna profondita' o composizione non valida; slot inattivi nulli in REAL/COMPLEX, flussi, sorgenti, IC/BC e I/O.
3. **Bilanci e limitazione del gruppo 0:** massa totale e per specie, trasporto, momento e variabili termiche secondo il contratto fisico; verificare separatamente il rate valutato su `q^n`, la limitazione rispetto a `q^*` e agli inventari, e il candidato `q^(0)` dopo l'applicazione. Per gli scambi, controllare rate unico, trasferimenti interlayer opposti, disponibilita' condivise, assenza di doppio conteggio di `deposit_tot` e substrato, dissipazione e scambi esterni distinti. Verificare stop fatale e risorse condivise.
4. **Spaziale e tempo:** stabilita' e finitezza nei casi di sviluppo ammessi, CFL condiviso, topografia piana/inclinata/curva, stratified rest e limiti wet/dry; rispettare il residuo ammesso su terreno curvo. Il bound `h1+h2` e' **provvisorio**, non una garanzia universale di iperbolicita'.
5. **Evidenze:** build optimized e strict/debug, regressioni fisiche, test REAL/COMPLEX, OpenMP 1/4 thread e casi integrati pertinenti. Registrare SHA e test effettivamente superati; CI precedente non valida un nuovo commit.

## 4. Validazione scientifica differita e decisioni aperte

La perdita di iperbolicita' e' un limite potenziale del **sistema continuo**, non un difetto risolvibile in generale riducendo il CFL. Le scansioni preliminari a coefficienti congelati evidenziano regimi con autovalori complessi e forniscono controlli utili, ma **non certificano** le equazioni con composizione variabile, shear, clamp, geometria curva o interfacce dry/wet.

**Durante lo sviluppo:** mantenere i bound preliminari solo per i test controllati; rilevare stati non finiti/non ammissibili e non attribuire validita' fisica a un semplice run riuscito. Non introdurre nuovi solver, viscosita' artificiale, attrito correttivo, soglie arbitrarie o regolarizzazioni per mascherare i casi non iperbolici.

**In TL13b:** studiare autovalori e iperbolicita' nei regimi fisici rilevanti, robustezza e convergenza al raffinamento, limiti di shear, contrasti di densita', wet/dry e pressione/clamp; confrontare con benchmark e manoscritto. Solo allora decidere, con approvazione scientifica, se servano limitazioni di applicabilita', scambi dissipativi motivati fisicamente o una formulazione differente.

Restano da definire dove non gia' specificati dal modello: dettagli di settling, elutriazione, drag e budget energetici; eventuale Darcy in L2, packing/consolidazione, carrier combinati e contratti dei metodi impliciti. Nessuna di queste scelte va introdotta solo per convenienza implementativa.

**Definizione di completamento:** architettura semplice e definitiva, catena delle dipendenze aggiornata, test pertinenti effettivamente superati e separazione esplicita tra **software funzionalmente completo** e **modello scientificamente validato**. `N_LAYERS=2` resta disabilitato fino al gate e alla relativa autorizzazione.
