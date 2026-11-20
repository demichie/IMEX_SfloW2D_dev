# IMEX-SfloW2D — Architettura e roadmap two-layer

**Stato:** 20 novembre 2026 · branch `feature/two-layer` · `N_LAYERS=2` rifiutato per default; opt-in `IMEX_TWO_LAYER_VALIDATION=1` soltanto per test TL13a.
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
`interlayer_drag_chezy_2d.f90`. Tutte le cinque unita' operative
(`interlayer_settling`, `interlayer_elutriation`, coordinatore
`interlayer_mass_exchange`, legge `interlayer_drag_chezy` e coordinatore
`interlayer_momentum_exchange`) espongono le quattro interfacce fisiche
`EXPLICIT`, `SEMI_IMPLICIT`, `IMPLICIT`, `ANALYTIC`: soltanto le
formulazioni effettivamente implementate vengono eseguite. Le altre
sollevano un errore fatale esplicito, come nei moduli storici di
`physics/mass_exchange`. `TEST_SOURCE_REGISTRY` verifica tutti i quindici
rifiuti in entrambe le build strict e optimized. Il NAMELIST canonico
`INTERLAYER_PARAMETERS` include i due oggetti indipendenti
`interlayer_mass_exchange` e `interlayer_momentum_exchange`, con metodi e
gruppi Lie separati. `INTERLAYER_PARAMETERS` e' l'unica NAMELIST
pair-wide: il parser e il backup usano esclusivamente questa forma.
Non esiste un secondo input per il solo scambio di massa. Aggiungere nuove leggi soltanto quando
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
| **TL00–TL04 — Stato e fisica per layer** | Slot opzionali e OU L1/L2 validati nei blocchi 1-3 | L1/L2 condividono il layout opzionale ma mantengono configurazioni, campi OU e flussi casuali indipendenti. Restano IC/BC e I/O/restart completi nei blocchi successivi. |
| **TL05 — Operatori L2** | Percorso spaziale e BC L1/L2 validati nel blocco 3 | Ricostruzione HP e BC per layer, selezione di velocita' L1/L2 e componenti inattive riusano gli operatori canonici; regressioni estese a tutti i lati e alle quattro occupazioni wet/dry. |
| **TL06 — Spaziale accoppiato** | Regresioni accoppiate validate nel blocco 3 | Il passo Lie usa lo stesso workset L1/L2, le velocita' di faccia comuni sono propagate sulle sole facce attive. Regresioni per BC indipendenti sui quattro lati, quattro stati wet/dry, budget e CFL condiviso su topografie Q1; la validazione degli autovalori accoppiati con shear resta nel gate scientifico TL13b. |
| **TL07 — Lie accoppiato** | Integrato; blocco 1 verificato | Trial L1/L2, gruppo zero, transazione interlayer, drag e scambi esterni di erosione/deposizione sono esercitati sul timestep produttivo nei quattro stati wet/dry. Il blocco 2 aggiunge trasporto Darcy L1/OU L2 e cooling L2 nel passo accoppiato. Il gate funzionale generale resta chiuso. |
| **TL08–TL11 — Chiusure e scambi** | Blocco 4 verificato nelle Actions strict/optimized | Settling/elutriazione e co-lift con due modalita' sono limitati congiuntamente per inventari lordi e packing finale. I tre fattori (settling, elutriazione, carrier in ALLOW_CONCENTRATED) possono differire. I test esercitano densita' e temperature diverse, specie multiple, trasferimenti opzionali e sequenza massa/drag. Il gate scientifico e gli end-to-end restano TL13. |
| **TL12 — I/O e integrazione** | Blocco 5 implementato; regressione spaziale/I/O verificata su `7f8189a76` | Cold-start da due file .q_2d, snapshot e raster L2, NetCDF q/qp L1+L2, restart v8 con OU/RNG L2, budget per layer e totale e regressione di round-trip. Il gate produttivo e gli EXAMPLES two-layer rimangono TL13a. |
| **TL13a — Gate funzionale** | Smoke Etna ed esteso (restart, conservazione, OpenMP two-layer) verdi su `847f7906`; chiusura CI complessiva pendente | Ottenere sul commit finale workflow optimized e strict/debug interamente verdi, oltre ai test two-layer. Gate produttivo non autorizzato. |
| **TL13b — Validazione scientifica** | Successiva | Benchmark fisici, convergenza, robustezza, iperbolicita', shear, wet/dry, pressioni e limiti delle chiusure. Rilascio solo dopo approvazione scientifica. |

**TL07 - collegamento del timestep (17 novembre 2026):** il vecchio arresto
post-spaziale in `time_integration_2d%advance` e' stato rimosso: lo stesso
percorso sorgenti e commit viene ora attraversato da L1 e L2. Il layout
continua a rifiutare `N_LAYERS=2` nei run produttivi. Il blocco 2 sostituisce il precedente rifiuto del campo stochastic L2
con un workspace e un flusso RNG indipendenti per layer. La funzione
produttiva two-layer resta comunque bloccata dal gate di layout fino
all'accettazione funzionale e scientifica.

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
La regressione copre ora anche le due schedulazioni dello scambio
interlayer, nel gruppo 0 e in un gruppo successivo. Per ogni scenario
wet/dry confronta il candidato accettato dal timestep con la
transazione fisica cell-local, applicata a q* e valutata su q^n
soltanto per il gruppo 0. Un caso wet/wet combina poi mass exchange
del gruppo 0 e drag analitico nel gruppo successivo, verificandone
l'ordine e la conservazione. I payload di questi passi rientrano nel
confronto bitwise OpenMP 1/4 thread. La prova resta da verificare in CI.

**TL07 - copertura aggiuntiva del gruppo zero:** la regressione integrata
esercita anche sorgenti di curvatura indipendenti per L1 e L2 nello
stesso timestep. Per ogni cella confronta i contributi espliciti
valutati sullo stato pre-spaziale con il candidato post-spaziale,
senza introdurre un secondo integratore. Il caso e' limitato alle
fixture con curvatura attiva ed e' incluso nel payload OpenMP
strict/optimized. La verifica CI sul commit che aggiunge il caso
resta necessaria; non sostituisce i test integrati delle sorgenti
di erosione/deposizione dal substrato.

**TL07 - sorgenti esterne e scambi interlayer nello stesso timestep:**
`TEST_SPATIAL_OPERATOR` comprende ora anche i gruppi espliciti
`mass_exchange_deposition` e `mass_exchange_erosion` per ciascun layer,
seguiti dallo scambio interlayer e dal gruppo drag. Gli stessi quattro
scenari wet/dry controllano i depositi nel substrato condiviso, i
prelievi erodibili, l'ammissibilita' dei candidati e il bilancio della
massa solida complessiva (L1+L2+substrato), confrontando ciascun caso
con un passo senza scambio esterno. Sono usati il registro, il workset,
i contatti e l'esecutore Lie di produzione, senza un integratore
alternativo. I payload sono confrontati bitwise tra 1 e 4 thread
OpenMP. I test richiedono ancora esecuzione CI sul commit corrente;
il gate di produzione L2 resta chiuso.

Il contratto attuale rimane di errore fatale senza retry: q/qp vengono
pubblicati solo dopo la riuscita del Lie globale, ma il substrato puo'
essere aggiornato per cella durante il sweep. Non dichiarare un rollback
globale recuperabile del substrato senza una migrazione coerente della
transazione. La copertura delle sorgenti di terreno e' ora inclusa nella regressione
integrata del timestep, ma necessita della verifica CI finale. Variabili
opzionali, BC e restart restano nel perimetro dei blocchi successivi.
Il gate produttivo `N_LAYERS=2` rimane chiuso.

**Blocco 2 - variabili opzionali e OU per layer (18 novembre 2026):**
`simulation_context_type` possiede un array di workspace OU, uno per
layer, senza copie dei campi e senza riutilizzare il rumore di L1 per L2.
L1 mantiene la successione RNG globale, inclusa quella ripristinata
dal restart binario storico; L2 possiede una successione indipendente
e ripristina il generatore globale dopo ogni prelievo. Ciascun workspace
usa il layer assegnato per leggere massa, velocita' e componenti conservate. Il
timestep prepara i campi dopo lo spaziale per ciascun layer attivo;
l'esecutore Lie legge il corrispondente `effective_Z` per ogni gruppo.
La legge termica e i coefficienti Darcy continuano a usare le
configurazioni e le cache di layer gia' canoniche. Il parser accetta
`STOCHASTIC_L2_PARAMETERS` per due layer e calcola `std_min` e la scala
di velocita' con lo stesso kernel di validazione L1.

`TEST_STOCHASTIC_EVOLUTION` verifica evoluzione, trasporto, seed e
indipendenza del campo L2; `TEST_SPATIAL_OPERATOR` esercita il
timestep integrato con Darcy passivo L1, OU L2 e una sorgente termica
esplicita L2, verificando che non modifichi L1. Il payload
`optional_lie.bin` confronta anche i campi OU L2 con rumore seminato e
lo stato dei due layer bitwise tra 1 e 4 thread OpenMP. La compilazione e le
regressioni sul commit del blocco 2 devono ancora essere osservate.
Il restart two-layer, il NetCDF e l'attivazione produttiva rimangono
nei blocchi successivi; nessuna ipotesi fisica e' stata modificata.

**Blocco 3 - operatori spaziali, wet/dry e BC L1/L2 (18 novembre 2026):**
`spatial_operator_2d%evaluate_coupled` propaga ora verso L2 gli stessi
limiti di velocita' totale soltanto sulle facce del workset OR realmente
attive; coefficienti di pressione e ricostruzioni restano propri di
ciascun layer. La vista `qp(idx_h,:,:)` del layer inferiore viene
accettata dalla preparazione geometrica L2 con shape assunto, senza
copia temporanea del dominio e senza cambiare `B`, `G` o le
formule idrostatiche e interlayer.

`TEST_SPATIAL_OPERATOR` verifica ora in modo integrato i quattro lati
W/E/S/N con dati Dirichlet distinti per L1/L2 e quattro selettori di
traccia HP, anche a dominio asciutto. Verifica inoltre tutte le
occupazioni wet/wet, solo L1, solo L2, dry/dry usando il workset OR,
il CFL comune e la conservazione della massa con il solver spaziale
di produzione, su topografie piane e Q1 curve esistenti. Il payload
`coupled_bc_wetdry.bin` viene confrontato bitwise tra OpenMP 1 e 4
thread in strict e optimized. Le regressioni del gruppo 0 per la
curvatura e dell'equilibrio stratificato rimangono attive.
Questo blocco non modifica la chiusura degli autovalori accoppiati;
la verifica scientifica di shear/iperbolicita' resta TL13b. Il gate
produttivo `N_LAYERS=2` rimane chiuso fino all'accettazione prevista.

**Blocco 4 - scambi conservativi L1/L2 (18 novembre 2026):**
Il coordinatore fisico esistente `interlayer_mass_exchange_2d`
continua a usare i rate prodotti dalle leggi cell-local `interlayer_settling`
ed `interlayer_elutriation`, senza spostare formule o introdurre un
secondo integratore. Il limitatore impone la disponibilita' di ogni
donatore sulle masse e sui volumi lordi pre-trasferimento, preserva
le proporzioni tra le classi di ogni processo e verifica il packing
finale in entrambi i layer. Quando il packing limita un processo,
non riduce inutilmente anche l'altro: ricerca scale separate di
settling ed elutriazione a partire da un candidato congiuntamente
ammissibile. In `ALLOW_CONCENTRATED` verifica anche se l'intero
trasferimento solido diventa ammissibile riducendo soltanto il
carrier ascendente; `LIMIT_SOLIDS` mantiene invece la composizione
della corrente imposta dall'input.

I trasferimenti conservativi usano le densita' e capacita'
termiche dei donatori e la conversione di volume ai destinatari.
`TEST_SOURCE_REGISTRY` aggiunge casi vincolati con fattori
solidi diversi, carrier limitato indipendentemente, densita'
donatore/destinatario diseguali, due classi solide e temperature
differenti. `test_interlayer_lie_transaction` esercita inoltre
la sequenza scambio di massa e drag Chézy sui candidati aggiornati,
e controlla che i componenti conservativi opzionali Darcy/OU
restino nei rispettivi layer. Il controllo finale azzera i rate
provvisori se una transazione non e' ammissibile.

La modifica riguarda il limite dei rate interlayer effettivamente
vincolati; le formule di settling, elutriazione, drag e il
comportamento single-layer non cambiano. Compilazione optimized,
strict/debug, regressioni e OpenMP sono ancora da verificare sul
commit finale del blocco. Il gate produttivo `N_LAYERS=2` resta
chiuso, insieme alle verifiche scientifiche TL13b.

**Blocco 5 - schema I/O definitivo e senza retrocompatibilita':**
Il codice usa ora un unico checkpoint binario non versionato per uno o
due layer; il suo header contiene la configurazione di griglia, le
componenti e il numero di layer. Il reader valida il layout e ripristina
`q(component,layer,x,y)`, topografia e substrato condivisi, RNG/OU e
storie diagnostiche. Non esistono percorsi di lettura o scrittura per
formati precedenti: i vecchi checkpoint devono essere rigenerati.

NetCDF usa lo stesso schema per `N_LAYERS=1` e `N_LAYERS=2`:
`q_layers(component,layer,x,y,time)` e `qp_layers` sono gli unici
campi di stato di layer; `b`, `h_total`, `eta_top`,
`mass_total`, inventari del substrato e diagnostiche aggregate sono
campi condivisi. Non si conservano in parallelo i precedenti campi
singolo-layer. Il conservativo ASCII senza suffisso appartiene al
single-layer; nel two-layer ogni stato e' esplicito `_L1_` o `_L2_`.
I cold-start two-layer richiedono due `.q_2d` coerenti. L'output ESRI
usa un solo DEM condiviso.

Le NAMELIST **senza** suffisso rimangono intenzionalmente parte della
configurazione single-layer; con due layer si usano `_L1` e `_L2`.
Non sono alias di retrocompatibilita' e non devono essere rimossi.

Le regressioni integrate controllano il round-trip del checkpoint e
delle variabili NetCDF, i test di restart e i lettori C/Python sono
migrati al formato canonico e lo script DAMBREAK legge `qp_layers`.
Le azioni CI e la validazione scientifica rimangono da eseguire per il
commit di questa migrazione. Il gate di produzione `N_LAYERS=2`
rimane chiuso fino a TL13a/TL13b.

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


**TL13a — rimozione di ripristino substrato non produttivo (18 novembre 2026):**
L'API `substrate_snapshot_type`, `snapshot_substrate` e
`restore_substrate` era usata soltanto da due fixture di test; il
timestep produttivo usa un errore fatale senza retry e non dispone
di un rollback recuperabile globale. La pulizia elimina tale
interfaccia e allinea i test al contratto attuale, mantenendo la
transazione della geometria e i controlli aggregati lava/granulare.
Le Actions sul commit di pulizia devono essere verificate
separatamente. Il gate `N_LAYERS=2` resta disabilitato.


**TL13a — integrazione funzionale e continuazione (18 novembre 2026):**
`TEST_SPATIAL_OPERATOR` controlla ora anche l'append NetCDF dopo
il restart: due record temporali e stati L1/L2 conservati in entrambi.
Il workflow `solver-validation.yml` include il test strutturale two-layer
nelle build strict/debug e optimized, insieme ai test produttivi di
restart e budget; il confronto OpenMP del driver single-layer
(1/2/4 thread) e' abilitato per default.
Questo non costituisce ancora un end-to-end produttivo a due layer:
`supported_layer_count` e il parser mantengono intenzionalmente il
gate chiuso. L'esecuzione delle Actions sulla revisione finale e
l'approvazione di un test controllato two-layer sono ancora necessarie.

## Etna dense-to-dilute example (preparation, not production activation)

`EXAMPLES/EXAMPLE_ETNA_TWO_LAYER` prepares the unchanged Etna dense avalanche
as the lower layer and an **exactly dry upper layer**. The prepared L1
conservative snapshot is copied byte-for-byte from a t=0 run of the current
single-layer solver with the original `EXAMPLE_ETNA` input and ellipsoidal
preprocessing. The L2 snapshot shares its grid coordinates and has zero
conservative components in every cell. No separate particle cloud is
initially imposed: L2 can form only through configured L1→L2 elutriation.

The example uses the one canonical Etna DEM archive, supplies explicit
per-layer boundary and closure NAMELISTs and the canonical pair-wide
`INTERLAYER_PARAMETERS`, and keeps exploratory elutriation coefficients
separate from scientifically validated Etna calibration. Its lightweight
snapshot-generation tests are part of `make check`. The solver still rejects
`N_LAYERS=2`; this example is **not** evidence of two-layer production
readiness. Future acceptance must verify cold-start physics, conservative
budgets, restart, netCDF and OpenMP through the actual driver.


**TL13a — opt-in del driver e smoke test (20 novembre 2026):** il layout
ammette ora entrambi i layer, ma il parser rifiuta `N_LAYERS=2` senza
`IMEX_TWO_LAYER_VALIDATION=1`. Il driver prepara il workset OR,
le due ricostruzioni CFL e il bound spaziale accoppiato attraverso il kernel
canonico. `EXAMPLES/EXAMPLE_ETNA_TWO_LAYER/run_validation.py` prepara gli
stati e verifica il gate di default, quindi tenta un breve run integrato
con output conservativo, ESRI, NetCDF e checkpoint. La Action
`solver-validation.yml` esegue questo test solo con input manuale
`two_layer_smoke=true`. **Questa nota descrive la fase iniziale di attivazione:**
le esecuzioni successive sono documentate nell'audit di chiusura TL13a.
TL13a/TL13b restano aperte. Il bound preliminare
non certifica l'iperbolicita' con shear e non autorizza uso produttivo.


**TL13a — audit finale del workset (20 novembre 2026):** il driver e
l'integratore utilizzavano esclusivamente `simulation%domains(1)`, mentre
`simulation_2d` allocava anche un `domain_type` L2 mai usato. Il contesto
possiede ora un unico `simulation%domain` (workset OR L1/L2), riusando i
kernel e gli operatori spaziali per-layer esistenti. Non sono cambiati
lo stato conservativo, le condizioni al contorno, le leggi sorgente o
il calcolo CFL. La chiusura TL13a richiede ancora regressioni integrate
strict/optimized completate su quest'ultimo commit. TL13b resta il gate
scientifico separato.
