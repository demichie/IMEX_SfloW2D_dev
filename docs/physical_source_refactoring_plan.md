# Refactoring delle sorgenti fisiche: architettura di riferimento

## Obiettivo
Quattro famiglie: rheology, mass exchange, mass injection, thermal.
Tre livelli necessari e sufficienti: assemblaggio delle sorgenti
(`source_terms_2d`), un coordinatore per famiglia, un modulo per
ogni modello fisico. Nessun wrapper o livello di dispatch intermedio.

## Regole vincolanti
- Ogni coordinatore e ogni modello espone quattro trattamenti subito dopo
  `CONTAINS` e sempre nell'ordine: explicit, semi_implicit, implicit, analytic.
- Il coordinatore seleziona direttamente i modelli, somma i tassi calcolati
  sullo stesso stato, applica vincoli condivisi realmente necessari oppure
  compone in ordine gli operatori analitici. La composizione non coincide
  automaticamente con la soluzione esatta della somma delle sorgenti.
- Ogni modulo specifico contiene la propria fisica. Gli helper privati
  devono eseguire calcoli reali, non limitarsi a inoltrare chiamate.
- **Tutte** le interfacce esplicite, sia dei coordinatori sia dei
  modelli specifici, ricevono sempre `physical_state_type`.
  La firma esplicita non accetta mai `complex_physical_state_type`.
- Se un esplicito combina la parte semi-implicita e quella implicita,
  converte lo stato reale in `complex_physical_state_type` **al proprio
  interno**, poi richiama direttamente l'interfaccia implicita complessa.
  La conversione non va spostata nel coordinatore e non richiede wrapper.
  Se l'esplicito valuta direttamente la sorgente, usa lo stato reale.
- Le interfacce semi-implicite e analitiche ricevono lo stato reale;
  quelle implicite ricevono lo stato complesso.
- Le conversioni effettuate dentro gli espliciti compositi sono
  intenzionali. Devono usare stati temporanei locali e thread-safe;
  non duplicare conversioni ulteriori a monte.
- Non passare separatamente temperatura, velocità, densità e altre
  grandezze intrinseche già contenute o ricavabili dallo stato.
  Mantenere invece gli argomenti realmente esterni: tempo, geometria,
  condizioni ambientali, substrato e parametri della legge.
  Valutare il costo dei ricalcoli per non introdurre regressioni.
- Gli stati reali e complessi contengono le grandezze cinematiche derivate
  `w`, `mod_vel_hor`, `mod_vel2`, `mod_vel` e
  `curvature_acceleration`. Le derivate topografiche non sono campi dello
  stato: sono argomenti obbligatori di `update_physical_state_kinematics`.
  Ogni modifica di `u` o `v` richiede di ricalcolare queste grandezze con
  la geometria locale appropriata. Le valutazioni complex-step
  impiegano somme di quadrati complessi, mai `ABS` delle velocità.
  La geometria deve essere passata esplicitamente per la cella interessata;
  i valori zero iniziali non sono sostituti di una topografia non piana.
- Assemblare lo stato una volta quando possibile, senza conversioni
  reali/complesse superflue. Gli stati temporanei devono essere locali
  e thread-safe. Nessun nuovo tipo specifico di famiglia.
- Le capacità dei quattro trattamenti non vanno falsamente uniformate:
  un flusso esatto restituisce un risultato aggiornato, mentre una
  sorgente istantanea restituisce un tasso.
- Mass exchange conserva la responsabilità comune dei vincoli di
  riserva e packing, non distribuibili ingenuamente ai modelli.

## Stato delle attività
- [x] Moduli fisici separati e quattro ingressi ordinati.
- [x] Rimossi livelli privati di dispatch ridondanti in thermal e rheology.
- [x] Rimossa la routine di forwarding esplicita del mass exchange.
- [ ] Eliminare gli argomenti derivabili dallo stato laddove possibile.
- [x] Aggiungere i campi cinematici ai due stati e aggiornamento esplicito.
- [x] I costruttori degli stati accettano la geometria locale opzionale e
      le valutazioni STAGE reali/complesse la trasmettono immediatamente.
- [x] Eliminata la copia completa dello stato nel coordinatore reologico.
- [x] Eliminato il ricalcolo duplicato della cinematica dal coordinatore:
      usa i campi dello stato e riceve la geometria gia' applicata.
- [x] Rimosse le derivate topografiche dai due tipi di stato; il refresh
      cinematico riceve cinque argomenti geometrici obbligatori.
- [ ] Audit degli stati costruiti senza geometria: i costruttori mantengono
      parametri opzionali per i test esistenti e assumono zero se assenti.
      Nel propagatore reologico analitico la geometria STAGE viene passata
      dal chiamante; verificare gli altri percorsi.
- [ ] Completare l'audit degli aggiornamenti diretti degli stati, inclusi i test.
- [x] Confermare il contratto: espliciti sempre reali; conversione locale negli espliciti compositi.
- [ ] Verificare tutte le firme e le conversioni degli espliciti compositi.
- [ ] Semplificare selezione/configurazione dei processi senza duplicazione.
- [ ] Eliminare altri helper che siano solo inoltri.
- [ ] Riutilizzare gli stati costruiti in source assembly.
- [ ] Verificare Build only e Physics regressions (D17/D18) al termine
      dell'intero blocco, poi Solver validation.

## Decisione approvata: futuro refactoring pore pressure (DA ESEGUIRE)
- Creare `src/physics/pore_pressure_laws/` con coordinatore
  `pore_pressure_laws_2d.f90` e modello
  `pore_pressure_dissipation_2d.f90`, ognuno con quattro ingressi
  ordinati (explicit, semi_implicit, implicit, analytic), stato reale
  salvo implicit complesso. La famiglia comprende esclusivamente
  dissipazione/diffusione della pressione interstiziale.
- La fisica dovuta a scambio di massa e iniezione NON risiede nella
  famiglia pore pressure e NON va calcolata nei coordinatori
  `mass_exchange_2d.f90` e `mass_injection_2d.f90`.
  Ogni modulo di processo specifico calcola tutti i rate delle equazioni
  conservative che il processo modifica, incluso il rate di pressione
  interstiziale; ogni coordinatore seleziona, limita e somma.
- Spostare `eval_pore_pressure_gas_loss` nel modello Darcy di mass
  exchange; `eval_pore_pressure_inlet_source` nei modelli specifici
  di iniezione; distribuire
  `eval_pore_pressure_mass_exchange_source` nei modelli specifici
  di mass exchange. Le proprietà costitutive veramente condivise
  possono restare in un helper fisico separato, non nel coordinatore.
- Particolare attenzione allo scambio di massa: oggi il contributo
  di pressione è determinato dal rate netto dopo i limitatori.
  Verificare matematicamente linearità e ripartizione per processo;
  estendere `mass_exchange_rates_type` ai rate di pressione se necessario,
  scalando tali rate assieme al processo. Preservare il totale ottenuto
  dopo vincoli carrier/packing e impedire doppi conteggi.
- Aggiornare Makefile.am, import, test e chiamanti; preservare le
  verifiche complex-step e propagazione analitica. Nessun nuovo tipo
  di stato, wrapper o formula di coordinamento.
- Questa sezione registra una decisione progettuale, NON una modifica
  già implementata. Procedere prima con il refactoring in corso e
  convalidare le modifiche con Build e Physics Regressions.

## Baseline di validazione
- L'utente conferma Build e Physics Regressions entrambe riuscite
  dopo il blocco di pulizia delle sorgenti fino al commit `279f5b7`.
  Solver validation e confronto prestazionale ancora da eseguire.

## Chiusura strutturale Rheology (2026-10-31)
- Verificati 11 modelli: quattro trattamenti nell'ordine prescritto e tipi di stato coerenti (reale per explicit, semi-implicit, analytic; complesso per implicit).
- Verificate conversioni locali negli espliciti compositi; nessun marcatore di continuazione isolato e nessuna riga oltre 132 caratteri nei file dei modelli.
- Rimosse dal coordinatore variabili cinematiche temporanee inutilizzate; i modelli leggono la cinematica dallo stato.
- **Refactoring strutturale Rheology concluso**. La validazione Build e Physics Regressions positiva riferita dall'utente precede l'ultima pulizia del coordinatore: ripetere entrambe sul commit finale prima di dichiarare la validazione completa. Benchmark prestazionale separato.
- Prossima famiglia da chiudere: Mass exchange. La futura revisione Pore pressure rimane registrata sopra ed e' distinta.

## Audit Mass exchange (2026-10-31)
- Verificati il coordinatore e i sette modelli specifici: ciascuno espone
  le quattro interfacce ordinate (explicit, semi_implicit, implicit,
  analytic), con stato reale/complesso conforme; i trattamenti non
  disponibili sono esplicitamente rifiutati.
- Uniformate le interfacce di erosione e di entrainment rimuovendo
  argomenti gia' disponibili nello stato o non utilizzati.
- Il coordinatore mantiene correttamente limitatori condivisi di
  carrier e packing e assembla i processi dopo la loro applicazione.
- Nel coordinatore restano Bx e By nell'interfaccia esplicita,
  ormai inutilizzati internamente: eliminare in fase di aggiornamento
  congiunto di source_terms_2d.f90 e dei chiamanti/test, per non
  modificare adesso il contratto pubblico esterno.
- La migrazione dei contributi alla pressione di poro nei singoli
  modelli fisici resta differita alla ristrutturazione pore pressure
  gia' descritta; deve preservare i rate dopo i limitatori.
- Dopo queste modifiche eseguire Build e Physics Regressions; non
  attribuire all'ultimo commit i test riusciti sul commit precedente.
