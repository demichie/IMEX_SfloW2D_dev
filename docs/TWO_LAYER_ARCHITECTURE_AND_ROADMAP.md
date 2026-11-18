# IMEX-SfloW2D — Architettura e roadmap two-layer

**Stato:** 17 novembre 2026 · branch `feature/two-layer` · `N_LAYERS=2` disabilitato in produzione.
**Obiettivo:** completare L1 inferiore e L2 superiore preservando comportamento, input e risultati single-layer.
Questa e' la specifica operativa delle scelte e del lavoro residuo, non una cronologia dei commit.
L'evidenza tecnica delle fasi precedenti e' in [TWO_LAYER_CODE_AUDIT.md](TWO_LAYER_CODE_AUDIT.md).

## 1. Architettura e invarianti

- **Stato unico:** `q/qp(component,layer,x,y)` rank-4, indici conservativi comuni ai layer; layout opzionale
  definito dall'OR delle attivazioni L1/L2, con componenti inattive esattamente nulle. Nessuna copia
  sincronizzata, coefficiente globale duplicato o struttura transitoria.
- **Proprieta' fisiche:** reologia, termica, Darcy/pore, stochastic/OU e scambi esterni configurati
  per layer. Topografia `B`, geometria, atmosfera, substrato, `deposit_tot`, `T_ground`,
  `T_ambient` e proprieta' delle specie solide sono condivisi. I parametri interlayer
  appartengono allo scambio L1/L2, non sono duplicati nelle due configurazioni.
- **Responsabilita':** `physics` contiene leggi e rate **cell-local**, senza indici di griglia;
  `spatial` contiene flussi e termini non conservativi; `numerics` gestisce candidati,
  workset, contatti, ordine Lie e controlli. Riutilizzare gli operatori, le strutture dati
  e i kernel matematici esistenti, senza framework paralleli o formule duplicate.
- **Contatti:** `layer_contacts_2d` classifica entrambi wet, entrambi dry, solo L1, solo L2.
  L2 isolato puo' toccare il terreno e L1 isolato puo' essere esposto all'atmosfera.
  I contatti derivano dai candidati correnti prima dei gruppi successivi e sono passati
  alle leggi interessate; non creare duplicati o flag persistenti sincronizzati.
  L'assenza di contatto L1/L2 non vieta il **cold-start** da un trasferimento interlayer.
- **Compatibilita':** mantenere le alternative input single-layer `<FAMILY>_PARAMETERS` e
  `<FAMILY>_L1_PARAMETERS` sullo stesso oggetto L1, mutuamente esclusive; per two-layer
  richiedere la configurazione esplicita dei layer. BC `bc(component,layer)` sui quattro lati,
  con slot inattivi omogenei. Non conservare interfacce interne obsolete.

### Spaziale accoppiato

Le superfici sono derivate: `eta1=B+h1`, `eta2=B+h1+h2`; metrica
`G=1/(1+B_x^2+B_y^2)` e curvatura derivano da `B`, non da `B+h1`.
Per `L=rho_m2*gprime2*h2`, i contributi di pressione interlayer sono
`S_inter1=-G*h1*d_j(L)` e `S_inter2=-L*d_j(G*h1)`; non duplicare il termine
gia' compreso nel path idrostatico di L2. Conservare il clamp buoyant L2
previsto dal modello, senza introdurre compensazioni artificiali.
La pressione interlayer e' un termine **spaziale**, non una sorgente Lie.

### Contratto dell'integratore Lie

Un solo timestep, un solo esecutore **cell-outer/OpenMP**, gruppi ordinati con
`scope=L1|L2|INTERLAYER` e metodo `EXPLICIT|NEWTON|ANALYTIC` secondo le capacita'
del catalogo. La directory fisica non impone ne' metodo ne' gruppo: sorgenti dello
stesso gruppo condividono il metodo; sorgenti in gruppi diversi possono usarne
di distinti. Il drag non deve condividere necessariamente il metodo del mass exchange.

**Gruppo 0:** lo spaziale opera per primo; le sorgenti assegnate al gruppo 0
sono **valutate su `q^n/qp^n`**, limitate rispetto ai candidati post-spaziali
`q^*` e agli inventari pertinenti, quindi applicate a `q^*`:

```text
q*       = q^n + dt * L(q^n)
R0       = S0(q^n, qp^n, ...)
R0_limit = Limit(R0, q*, dt, inventories)
q^(0)    = q* + dt * R0_limit
```

I gruppi successivi leggono i candidati aggiornati dai gruppi precedenti.
Ogni scambio interlayer legge uno stesso snapshot per i due layer e applica
contributi accoppiati ai due candidati. Il limitatore opera **prima** di ogni
aggiornamento; niente clipping ex post, rollback o retry. Un errore invalida
il passo e provoca stop fatale. I controlli complessivi sono a fine integratore,
con verifiche intermedie solo dove necessarie a limiter, conversioni o solver.
Preservare l'ordine di valutazione single-layer e le regole esistenti sui
prelievi dal substrato e dal materiale preesistente.

## 2. Scambi fisici interlayer

### Organizzazione dei moduli e dei gruppi

`src/physics/mass_exchange/` rimane responsabile delle sorgenti locali di
un layer verso terreno e atmosfera. Gli scambi interni L1/L2 appartengono a
`src/physics/interlayer_mass_exchange/` per settling L2→L1 ed elutriazione L1→L2.
Il drag e' fisicamente distinto e appartiene a
`src/physics/interlayer_momentum_exchange/`, con il coordinatore
`interlayer_momentum_exchange_2d.f90` e la legge
`interlayer_drag_chezy_2d.f90`. Il NAMELIST canonico
`INTERLAYER_PARAMETERS` include i due oggetti indipendenti
`interlayer_mass_exchange` e `interlayer_momentum_exchange`, con metodi e
gruppi Lie separati; `INTERLAYER_EXCHANGE_PARAMETERS` rimane un alias storico
per il solo mass exchange e le due sintassi sono mutuamente esclusive.
Il backup produce soltanto la forma canonica. Aggiungere nuove leggi soltanto quando
servono, senza interfacce transitorie. Condividere i kernel fisici gia' esistenti se equivalenti,
senza assimilare la deposizione nel substrato al trasferimento interlayer.

**Settling ed elutriazione** producono rate fisici indipendenti e sono
coordinati da `interlayer_mass_exchange_2d.f90` per i limiti di inventario e
packing, analogamente al coordinatore `mass_exchange_2d.f90` esistente.
Il coordinatore ha una responsabilita' effettiva, non introduce un nuovo
integratore o una seconda sorgente autorevole. Si prevede un **gruppo
esplicito comune per gli scambi di massa**, cosi' i rate condividono
lo snapshot e sono limitati congiuntamente; verificare al momento
dell'integrazione l'interazione con gli altri gruppi del registro.

**Drag** e' un processo distinto, senza scambio di massa. La legge
quadratica Chézy dell'equazione (63) e' integrata analiticamente
in un sottopasso a masse e coefficienti costanti, conservando il momento
totale e lasciando invariati massa, composizione ed energia termica
di ogni layer, come approvato. Il coordinatore possiede l'unica
configurazione pair-wide e la legge e' un kernel puro; la soglia di
contatto e' `h1,h2 > 1e-3 m`. Il drag usa un gruppo `ANALYTIC` separato
con id positivo: il gruppo zero non e' previsto per questa integrazione
nonlineare esatta. L'input, il registro e l'esecutore Lie condividono
ora il coordinatore fisico, e le regressioni verificano le configurazioni
separate. Resta da validare l'intero timestep two-layer prima di
attivare la fisica in produzione.

### Solidi, carrier e composizione della corrente elutriata

Le leggi del riferimento `IMEX_v3.pdf` determinano i flussi solidi
per specie di settling (L2→L1) ed elutriazione (L1→L2). Si conserva la
massa di ogni specie **con un solo rate per processo** e contributi
opposti nei due layer; momento ed energia termica sono trasferiti
secondo le velocita' e temperature fisiche dei rispettivi processi.

**Elutriazione:** il parametro di input
`elutriation_solid_volume_fraction`, con `0 < alpha_lift <= 1`, indica
la frazione volumetrica solida **della corrente trasferita**, non
la concentrazione finale di L2. Dato un volume solido proposto
`dV_s`, il volume di carrier richiesto e'
`dV_c_req = (1-alpha_lift)/alpha_lift * dV_s`; `alpha_lift=1`
corrisponde a nessun carrier associato. Non fissare ancora un
default numerico per `alpha_lift` senza giustificazione.

Il secondo input `elutriation_carrier_mode` seleziona:

- `LIMIT_SOLIDS`: in assenza di carrier sufficiente, ridurre insieme il
  flusso solido e quello di carrier per rispettare la composizione prescritta.
- `ALLOW_CONCENTRATED`: mantenere il flusso solido per quanto ammissibile,
  trasferendo il carrier disponibile; la **corrente** puo' essere piu'
  concentrata del valore prescritto, ma gli stati finali devono sempre
  rispettare packing, inventari e ammissibilita'. Non creare carrier.

I modi non cambiano la legge fisica che propone il rate solido e non
introducono un vincolo di composizione finale su L2. La diluizione
atmosferica resta un processo esterno distinto.

Le due leggi disponibili nel percorso cell-local comprendono il
settling derivato dal kernel condiviso della deposizione e l'elutriazione
`LINEAR_USTAR` o `SHEAR_POWER` del riferimento. I coefficienti Shields
e di lunghezza di sollevamento della seconda chiusura richiedono input
esplicito e calibrazione; non dedurre valori fisici dal solo successo dei test.
Le sorgenti termiche usano il calore specifico e la temperatura del donatore;
il momento del settling usa la velocita' interfaciale comune, quello
dell'elutriazione la velocita' del donatore. Il vettore accoppiato cell-local
e' collegato allo scope Lie `INTERLAYER` nell'esecutore cell-outer,
mentre il percorso di produzione L2 rimane bloccato fino al gate funzionale.

**Settling:** trasferire i solidi e soltanto il carrier L2→L1
eventualmente necessario a non superare il packing di L1, tenendo
conto **anche dell'elutriazione simultanea e del suo co-lift**.
Non aggiungere gas interstiziale senza il corrispondente prelievo
da L2 (approssimazione non conservativa descritta nel riferimento).
Il packing riguarda gli stati risultanti; il carrier associato
all'elutriazione riguarda invece la composizione della **corrente**.
Distinguere sempre packing di ammissibilita', packing della chiusura
fisica e fattore di hindered settling.

### Coordinatore e limitazione congiunta

Come in `mass_exchange`, le singole leggi restituiscono rate,
poi **un solo coordinatore cell-local** valuta prima dell'applicazione
la competizione per gli inventari di solidi e carrier e limita i
trasferimenti congiuntamente. I due layer sono simultaneamente
donatori e destinatari: non e' sufficiente applicare due limiter
monolayer indipendenti.

- Calcolare i flussi **lordi** di solidi per specie nelle due direzioni,
  la richiesta di carrier del co-lift e il carrier eventualmente
  necessario al settling. Per il gruppo 0 rispettare la distinzione
  rate su `q^n` / disponibilita' sul candidato `q^*`.
- Imporre disponibilita' dei donatori, non negativita', packing e
  ammissibilita' di **entrambi** i candidati, usando lo stesso
  timestep; considerare le disponibilita' condivise senza sottrarre
  piu' volte lo stesso materiale. Non compensare automaticamente
  apporti e prelievi lordi prima della verifica delle disponibilita'.
  La possibile riutilizzabilita' simultanea del materiale trasferito
  va decisa esplicitamente, senza dipendenza nascosta dall'ordine.
- Consentire fattori di riduzione **distinti** per settling,
  elutriazione e, in `ALLOW_CONCENTRATED`, carrier, ma calcolarli
  nello stesso problema di ammissibilita'. In `LIMIT_SOLIDS`
  il rate di co-lift mantiene la proporzione prescritta.
- Applicare solo i trasferimenti accettati, con debiti e crediti
  uguali per le **masse** di specie/carrier. Non ridurre il problema
  a un unico flusso *netto* di carrier quando i due flussi lordi
  trasportano momento o energia termica diversi. Convertire
  volumi e masse con le densita' coerenti con lo stato termodinamico.
- Coprire i quattro stati wet/dry: settling puo' formare L1 da L2;
  elutriazione puo' formare L2 da L1 senza contatto interlayer
  preesistente; drag richiede entrambi i layer presenti.
  Nessuna correzione arbitraria dello stato dopo il trasferimento.

La scelta tra `LIMIT_SOLIDS` e `ALLOW_CONCENTRATED` e' un'opzione
fisica di input; non autorizza a violare vincoli conservativi o di
ammissibilita'. Il contratto delle soglie wet/dry e delle chiusure
termodinamiche deve essere verificato nelle implementazioni.

## 3. Stato e sequenza di lavoro

**Priorita': completare il percorso two-layer end-to-end.** La
predisposizione architetturale e i test isolati non equivalgono ad
abilitazione in produzione.

| Fase | Stato | Lavoro residuo / gate |
|---|---|---|
| **TL00–TL04 — Stato e fisica per layer** | Base predisposta | Verificare i percorsi integrati ancora scoperti: slot opzionali, IC/BC, input, output e regressioni L1/L2. |
| **TL05 — Operatori L2** | Test isolati presenti | Integrare copertura dry, BC, sorgenti e variabili opzionali senza duplicare operatori. |
| **TL06 — Spaziale accoppiato** | Kernel e test parziali | Collegare completamente i path, i trial e i controlli al passo Lie two-layer; verificare bilanci, CFL comune, dry/BC e terreno curvo. |
| **TL07 — Lie accoppiato** | Percorso spaziale-Lie-commit collegato, da verificare | Il trial L1/L2 entra nell'esecutore Lie condiviso e pubblica i due layer sul workset comune; il gate di layout resta attivo. Gruppo 0 da snapshot pre-spaziale, inventari post-spaziali e transazione interlayer sono implementati. Restano regressioni integrate del timestep, optional slots, OpenMP e gate L2. La fisica stochastic L2 e' esplicitamente rifiutata finche' manca il workspace OU L2. |
| **TL08–TL11 — Chiusure e scambi** | Kernel e assemblatore cell-local predisposti | Rate di settling/elutriazione, co-lift in due modalita', limitatore congiunto, trasformazione dei volumi con densita' di carrier distinte e sorgenti accoppiate di massa/specie/momento/energia termica con regressioni isolate. Restano validazione della transazione Lie effettiva, variabili opzionali per layer, stati termodinamici finali, drag e test two-layer end-to-end. |
| **TL12 — I/O e integrazione** | Parziale | Input L1/L2, NetCDF, restart/round-trip, diagnostiche e EXAMPLES; preservare input e output single-layer. |
| **TL13a — Gate funzionale** | Non superato | Percorso end-to-end, assenza di architetture obsolete, build strict/optimized, regressioni, OpenMP e invarianti essenziali. Solo dopo, eventuale abilitazione **controllata** per test. |
| **TL13b — Validazione scientifica** | Successiva | Benchmark fisici, convergenza, robustezza, iperbolicita', shear, wet/dry, pressioni e limiti delle chiusure. Rilascio solo dopo approvazione scientifica. |

**TL07 - collegamento del timestep (17 novembre 2026):** il vecchio arresto
post-spaziale in `time_integration_2d%advance` e' stato rimosso: lo stesso
percorso sorgenti e commit viene ora attraversato da L1 e L2. Il layout
continua a rifiutare `N_LAYERS=2` nei run produttivi. Per evitare l'uso
incrociato del rumore L1, l'eventuale reologia/trasporto stochastic L2
non e' ancora abilitata. Questo collegamento richiede test di passo
accoppiato con invarianti fisici prima di essere considerato concluso.

**TL07 - regressione integrata del timestep:** `TEST_SPATIAL_OPERATOR`
esercita ora direttamente `time_integration_workspace_type%advance`
su un supporto compatto a due layer e sull'unico workset comune.
Per ciascuna delle quattro configurazioni wet/dry confronta lo stesso
passo con coefficiente Chézy nullo e positivo, controllando il bilancio
locale del momento tra i layer e l'invarianza di massa, solidi ed
energia termica per effetto del solo drag. I profili strict/optimized
usano lo script esistente, che confronta ora anche i payload del passo
accoppiato tra OpenMP 1 e 4 thread. La regressione include ora anche
settling ed elutriazione sul timestep effettivo: il confronto con le
sorgenti disattivate verifica il bilancio di tutte le equazioni
conservative e il cold-start di L1 e L2, senza duplicare l'integratore.
La prova non abilita `N_LAYERS=2` e resta da verificare in CI.
Restano aperti i test completi del gruppo 0 interlayer, sorgenti di
terreno, variabili opzionali e restart.

**Sequenza:** completare TL06/TL07 e la fisica TL08–TL11,
coltivando in parallelo i test ancora mancanti TL00–TL05 e l'I/O TL12.
Poi TL13a e TL13b. Non dichiarare una fase conclusa sulla base
del solo test di un kernel.

## 4. Verifica e questioni scientifiche aperte

- **Regressione:** nessuna variazione non autorizzata dei risultati
  single-layer, delle formule, del floating-point, dei metodi o
  degli input esistenti. Test su quattro configurazioni wet/dry,
  cold-start L1/L2, limiti packing, solidi per classe, carrier nei
  due modi, inventari esauriti, substrato/deposito e assenza di
  creazione o distruzione di massa.
- **Budget:** conservazione interlayer di massa totale, specie e
  carrier; momento e trasferimenti termici coerenti con ogni legge;
  contabilizzare la dissipazione del drag, senza conversioni
  termiche inventate. Test di rate, limiter e candidati distinti,
  inclusa la semantica `q^n`/`q^*` del gruppo 0.
- **Integrazione:** build optimized e strict/debug, consistenza
  REAL/COMPLEX dove richiesta, OpenMP 1/4 thread, CFL, geometrie
  piano/inclinato/curvo, stratified rest, BC e casi EXAMPLES.
  Audit completo di chiamanti, `.inc`, Makefile, parsing,
  diagnostiche, I/O e restart per ogni nuova interfaccia.
- **Validazione:** il bound preliminare `h1+h2` non certifica
  l'iperbolicita' del sistema continuo; il CFL non risolve
  eventuali autovalori complessi. In TL13b verificare regimi con
  shear, contrasti di densita', composizioni variabili e wet/dry.
  Non introdurre soglie, viscosita' o forze correttive per
  mascherare problemi fisici.

**Da definire scientificamente prima dell'implementazione:**
default numerico di `elutriation_solid_volume_fraction`,
possibile utilizzo di apporti simultanei nei budget lordi,
chiusura dei volumi di carrier comprimibile, parametri effettivi
di settling/elutriazione, soglie e dissipazione del drag,
ordine definitivo dei gruppi e vincoli di eventuali metodi
impliciti. Resta da valutare la copertura Darcy L2,
packing/consolidazione e carrier combinati. Fare riferimento
alle formulazioni di `IMEX_v3.pdf` distinguendo equazioni
fisiche e approssimazioni non conservative dell'implementazione
precedente; nessuna modifica scientifica implicita.

**Definizione di completamento:** architettura definitiva senza
duplicazioni, dipendenze coerenti, test realmente superati sul
commit interessato e distinzione tra software funzionalmente
integrato e modello scientificamente validato. `N_LAYERS=2`
resta disabilitato fino al gate e alla relativa autorizzazione.
