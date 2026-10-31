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
- Le interfacce esplicite **dei modelli** che sommano le componenti
  semi-implicita e implicita ricevono `complex_physical_state_type`.
  Quando l'esplicito costruisce direttamente la sorgente può ricevere
  `physical_state_type`. Semi-implicito e analitico: stato reale;
  implicito: stato complesso.
- In un coordinatore che supporta modelli espliciti eterogenei il punto
  di conversione va definito una volta sola, non ripetuto all'interno
  di ogni modello; non introdurre un framework di adattatori.
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
- [ ] Consolidare la conversione complessa per gli espliciti compositi.
- [ ] Semplificare selezione/configurazione dei processi senza duplicazione.
- [ ] Eliminare altri helper che siano solo inoltri.
- [ ] Riutilizzare gli stati costruiti in source assembly.
- [ ] Verificare Build only e Physics regressions (D17/D18) al termine
      dell'intero blocco, poi Solver validation.
