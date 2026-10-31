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
- [x] Confermare il contratto: espliciti sempre reali; conversione locale negli espliciti compositi.
- [ ] Verificare tutte le firme e le conversioni degli espliciti compositi.
- [ ] Semplificare selezione/configurazione dei processi senza duplicazione.
- [ ] Eliminare altri helper che siano solo inoltri.
- [ ] Riutilizzare gli stati costruiti in source assembly.
- [ ] Verificare Build only e Physics regressions (D17/D18) al termine
      dell'intero blocco, poi Solver validation.
